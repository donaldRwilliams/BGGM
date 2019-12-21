#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#' @import shiny
#'
#'
#' @examples
shiny_bggm.explore <- function(object,...){
  x <- object

  input_dots <- list(...)
  if(is.null( input_dots$node_groups )){
    cell_width <- c("50%", "50%", "50%")

  } else{
    cell_width <- c("50%", "50%", "65%")
  }
  ui <- fluidPage(
    sliderInput(
      "cut_off", "Bayes Factor Threshold:",
      min = 1, max = 100, value = 3, step = 1
    ),
    mainPanel("",
              fluidRow(
                splitLayout(cellWidths = cell_width,
                            plotOutput("alt_plot"),
                            plotOutput("null_plot"),
                            plotOutput("incon_plot"))
              )
    )
    # plotOutput("alt_plot")
  )


  server <- shinyServer(


    function(input,output){

      inputs_dots <- reactive({
        input_dots <- list(...)

      })

      sel <- reactive({
        sel <- select(x, BF_cut = input$cut_off)
      })


      incon <- reactive({
        incon <- ifelse(sel()$Adj_10 + sel()$Adj_01 == 0, 1, 0)

      })

      plts <- reactive({
        plts <- plot(sel(), ...)
      })


      # output$value <- renderText({ input$caption })

      output$alt_plot <- renderPlot({
        if(is.null(inputs_dots()$node_groups)){
          plts()$plt + ggtitle("Conditional Dependence") +
            theme(legend.title = element_blank())
        } else {
          plts()$plt + ggtitle("Conditional Dependence") +
            theme(legend.position  = "none")

        }

      })

      output$null_plot <- renderPlot({
        if(is.null(inputs_dots()$node_groups)){
          plts()$plt_null +
            ggtitle("Conditional Independence") +
            theme(legend.title = element_blank())

        } else {
          plts()$plt_null + ggtitle("Conditional Independence") +
            theme(legend.position  = "none")

        }


      })

      output$incon_plot <- renderPlot({
        plot_adjacency(incon(),...) +  ggtitle("Inconclusive") +
          theme(legend.title = element_blank())
      })


    })

  shiny::shinyApp(ui = ui, server = server)

}


#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
shiny_bggm <- function(object,...){
  UseMethod("shiny_bggm", object)
}

