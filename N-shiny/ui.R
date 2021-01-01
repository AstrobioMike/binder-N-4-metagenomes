library(shiny)
library(shinythemes)
library(shinydashboard)

# loading data and functions
source("data-and-functions.R")


# Use a fluid Bootstrap layout
ui <- navbarPage("Guerrero Negro mat metagenomes - N focus", inverse = TRUE, theme = shinytheme("spacelab"), position = "fixed-top",
                 tags$style(type="text/css", "body {padding-top: 70px;}"), selected = "Plots",

        tabPanel("Plots",

            fluidPage(

                titlePanel("Plots"),

                # Generate a row with a sidebar
                sidebarLayout(

                    # Define the sidebar with one input
                    sidebarPanel(
                        selectInput("KO_ID", "Select KO ID to plot:",
                                    choices=N_KO_norm_cov_tab$label),
                                    hr(),
                                    helpText("KEGG N-metabolism KOs recovered in our metagenomes.
                                             Bla bla bla.
                                             More bla."),
                        width = 5
                    ),

                    # Create a spot for the barplot

                    column(width = 6,
                        mainPanel(
                            plotOutput("plot")
                        ), offset = 1
                    )

                ),

                # all bars
                hr(),

                img(src = "all-N-plots.png")

            )
        ),

        tabPanel("Data",

            titlePanel("Data"), hr(),

            DT::dataTableOutput("table")

        ),

        tabPanel("All N-metabolism KOs",

            titlePanel("All N-metabolism KOs"), hr(),

            DT::dataTableOutput("all_ko_table")

        )

    )