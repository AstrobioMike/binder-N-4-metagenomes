library(shiny)
library(shinythemes)
library(KEGGREST)
library(tidyverse)
library(rclipboard)

# loading data and functions
source("data/data-and-functions.R")

# Use a fluid Bootstrap layout
ui <- navbarPage("Guerrero Negro mat metagenomes - Nitrogen focus", inverse = TRUE, theme = shinytheme("spacelab"), position = "fixed-top",
                 tags$style(type="text/css", "body {padding-top: 70px;}"), selected = "KO coverage plots",

        tabPanel("KO coverage plots",

            fluidPage(

                titlePanel("KO coverage plots"),

                HTML("<p><b>Nitrogen focus is based on KEGG Orthologs (KOs) from <a href='https://www.genome.jp/kegg-bin/show_pathway?map00910' target='_blank'>KEGG's Nitrogen metabolism pathway</a>.</b><br>
                     Coverage per million is normalized to within this functional grouping.</p>"),

                hr(),

                # Generate a row with a sidebar
                sidebarLayout(

                    # Define the sidebar with one input
                    sidebarPanel(
                        selectInput(inputId = "KO_ID", label = HTML("<div style='color: black'>Select KO ID to view:</div>"),
                                    choices=N_KO_norm_cov_tab$label),

                                    helpText(HTML("<p style='color: black'>This list contains the KEGG N-metabolism KOs recovered in our metagenomes.<br><br>
                                                  Choosing one will:</p>
                                                  <ul>
                                                        <li style='color: black'>plot that KOs normalized coverage to the right</li>
                                                        <li style='color: black'>plot those coverages based on the phyla of the genes contributing to that KO</li>
                                                        <li style='color: black'>adjust the table below to show the gene-level data for the selected KO</li>
                                                  </ul>
                                                  <br>
                                                  <p style='color: black'>The underlying data can be seen and downloaded from the <b><i>KO-grouped data</i></b> and <b><i>Gene-level data</i></b> pages accessible from the menu bar above.<br>
                                                  <br>
                                                  Individual gene sequences can also be accessed from the <b><i>Gene-level data</i></b> page.<br><br>

                                                  A static image of all N-related KO coverages is at the bottom of this page.<br><br>

                                                  These only include KO terms detected in our data. See the <b><i>All N-metabolism KO info</i></b> page for all KO terms included in KEGG's Nitrogen metabolism pathway.<br>
                                                  </p>")),
                        width = 5,

                    ),

                    # Create a spot for the barplot

                    column(width = 6,
                        mainPanel(
                            plotOutput("plot", height = "500px")
                        ), offset = 1
                    )

                ),

                hr(),

                plotOutput("taxa_plot", height = "500px"),

                HTML("<p style='text-align: right; padding-right: 200px'><b>*</b>Other are all those contributing less than 1% of all depths.</p>"),


                hr(),

                HTML("<p><b>This table depicts information for the contributing genes of the above-selected KO.</b><br>
                     The full underlying data is downloadable from the <b><i>Gene-level data</i></b> page.</p>"),

                DT::dataTableOutput("taxa_table"),

                # all bars
                hr(),

                HTML("<p><b>Static plot of normalized coverages of all detected N-related KO terms.</b></p>"),

                img(src = "all-N-plots.png")

            )
        ),

        tabPanel("KO-grouped data",

            titlePanel("KO-grouped data"),

            HTML("<b>This table shows combined coverage and stats information for each N-related KO recovered in our dataset.</b><br>"),
            HTML("The table can be filtered by searching on the top right."),

            hr(),

            DT::dataTableOutput("KO_grouped_table"),

            downloadButton("download_grouped_KO_data", " Download full data table")

        ),

        tabPanel("Gene-level data",

            titlePanel("Gene-level data"),

            HTML("<b>This page provides access to gene sequences and the coverage and taxonomy information for the 1,305 unique N-related genes we recovered.</b><br>"),
            HTML("The table can be filtered by searching on the top right."),

            hr(),

            rclipboardSetup(),

            textInput("gene_ID", label = h3("Get gene sequence"), value = "74"),

            h4("AA sequence:"),
            fluidRow(column(5, verbatimTextOutput("AA_seq")), column(6, uiOutput("copy_AA_seq")), column(7, HTML("<p><b>Open <a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome' target='_blank'>NCBI BLASTP</a></b></p>"))),

            h4("DNA sequence:"),
            fluidRow(column(5, verbatimTextOutput("DNA_seq")), column(6, uiOutput("copy_DNA_seq")), column(7, HTML("<p><b>Open <a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome' target='_blank'>NCBI BLASTN</a></b></p>"))),

            hr(),

            DT::dataTableOutput("gene_level_table"),

            downloadButton("download_gene_level_data", " Download full data table")

        ),

        tabPanel("All N-metabolism KO info",

            titlePanel("All N-metabolism KO info"),

            HTML("<b>This page lists all KOs from <a href='https://www.genome.jp/kegg-bin/show_pathway?map00910'>KEGG's Nitrogen metabolism pathway</a>.</b><br>"),
            HTML("The final column value is 'true' if that particular KO was recovered in our data."),

            hr(),

            DT::dataTableOutput("all_ko_table"),

            downloadButton("download_N_KO_info", " Download full data table")

        )

    )

# Define server logic required to draw a histogram
server <- function(input, output) {

    # Fill in the spot we created for a plot
    output$plot <- renderPlot({

        plot_KO(input$KO_ID)

    })

    output$KO_grouped_table <- DT::renderDataTable(DT::datatable({
        N_KO_norm_cov_tab[, 1:12]
    }, rownames = FALSE, options = list(pageLength = 20)))

    output$download_grouped_KO_data <- downloadHandler(
        filename = "N-KO-normalized-coverage.tsv",
        content = function(file) {
            write.table(N_KO_norm_cov_tab[, 1:12], file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
    )

    output$all_ko_table <- DT::renderDataTable(DT::datatable({
        N_KO_info_tab
    }, rownames = FALSE, options = list(pageLength = 20)))

    output$AA_seq <- renderText({

        if (! input$gene_ID %in% N_gene_level_tab$gene_ID ) {
            paste("Requested gene ID not present :(", sep = "\n")
        } else {
            get_seq(input$gene_ID)
        }
    })

    output$DNA_seq <- renderText({

        if (! input$gene_ID %in% N_gene_level_tab$gene_ID ) {
            paste("...", sep = "\n")
        } else {
            get_seq(input$gene_ID, "DNA")
        }
    })

    output$gene_level_table <- DT::renderDataTable(DT::datatable({
        N_gene_level_tab[, 1:16] %>% mutate(across(2:5, round, 3))
    }, rownames = FALSE, options = list(pageLength = 10)))

    output$download_gene_level_data <- downloadHandler(
        filename = "N-gene-normalized-master-tab.tsv",
        content = function(file) {
            write.table(N_gene_level_tab, file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
    )


    # Add clipboard buttons
    output$copy_AA_seq <- renderUI({
        rclipButton("clipbtn", " Copy", icon = icon("clipboard"), get_seq(input$gene_ID))
    })

    output$copy_DNA_seq <- renderUI({
        rclipButton("clipbtn", " Copy", icon = icon("clipboard"), get_seq(input$gene_ID, "DNA"))
    })


    output$taxa_plot <- renderPlot({
        plot_KO_tax(input$KO_ID, tab = N_gene_level_tab, info_tab = N_KO_norm_cov_tab)
    })

    output$taxa_table <- DT::renderDataTable(DT::datatable({
        subset_tax_tab(input$KO_ID, tab = N_gene_level_tab, info_tab = N_KO_norm_cov_tab) %>% mutate(across(2:5, round, 3))
    }, rownames = FALSE, options = list(pageLength = 10)))

    output$download_N_KO_info <- downloadHandler(
        filename = "All_N_KO_info.tsv",
        content = function(file) {
            write.table(N_KO_info_tab, file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
    )

}

# Run the application
shinyApp(ui = ui, server = server)
