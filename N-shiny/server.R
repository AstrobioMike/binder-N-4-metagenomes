library(shiny)
library(shinythemes)
library(shinydashboard)

# loading data and functions
source("data-and-functions.R")

# Define server logic required to draw a histogram
server <- function(input, output) {

    # Fill in the spot we created for a plot
    output$plot <- renderPlot({

        # Render a barplot
        plot_KO(input$KO_ID)

    })

    output$table <- DT::renderDataTable(DT::datatable({
        N_KO_norm_cov_tab[, 1:12]
    }, options = list(pageLength = 20)))

    output$all_ko_table <- DT::renderDataTable(DT::datatable({
        N_KO_info_tab
    }, options = list(pageLength = 20)))

}
