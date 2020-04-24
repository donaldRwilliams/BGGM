#' Shiny App for \code{explore} Objects
#'
#' @param object object of class \code{explore}
#' @param ... currently ignored
#'
#' @export
#' @import shiny
#'
#'
#' @examples
#' \donttest{
#' fit <- explore(tas[,1:20],
#'                prior_sd = 0.20,
#'                 iter = 5000)
#'
#' sel <- select(fit)
#'
#' shiny_bggm(fit)
#' }
shiny_bggm.explore <- function(object,...){
  x <- object
  input_dots <- list(...)



  choices <- c("two.sided","greater","less")
  names(choices) <- c("\u03C1 \u2260 0","\u03C1 > 0","\u03C1 < 0")
  choices



  ui <-  fluidPage(
    titlePanel("BGGM"),
    selectInput("hyp","Hypothesis",choices=choices, selected = "greater"),
    fluidRow(column(10, sliderInput("cut_off",
                                    label = "Bayes Factor Threshold:",
                                    min = 3, max = 100, value = 3))),
    plotOutput("plots"),
    downloadButton('Download','Download Plot')
  )






  server <- shinyServer(
    function(input,output,session){
      inputs_dots <- reactive({
        input_dots <- list(...)
      })
      sel <- reactive({
        sel <- select(x, BF_cut = input$cut_off, alternative = input$hyp)
      })
      incon <- reactive({

        if(input$hyp != "two.sided"){
          incon <- ifelse(sel()$Adj_20 + sel()$Adj_01 == 0, 1, 0)
        }  else {
          incon <- ifelse(sel()$Adj_10 + sel()$Adj_01 == 0, 1, 0)
        }
      })


      temp <- reactive({

        if(is.null(inputs_dots()$palette)){
          palette <- "Set2"
        } else {
          palette <- inputs_dots()$palette
        }
      })

      plots <- function(){

        plts <- plot(sel(), palette = temp(), ...)

        if(is.null(inputs_dots()$node_groups)){
          plt_h1 <- plts$plt +
            theme(legend.position = "none") +
            ggtitle("Conditional Dependence")
          plt_h0 <- plts$plt_null + theme(legend.position = "none") +
            ggtitle("Conditional Independence")

          plt_adj <- plot_adjacency(incon(), ... )

          cowplot::plot_grid(plt_h1, plt_h0, plt_adj, nrow = 1)




        } else {

          plt_h1 <- plts$plt +
            theme(legend.position = "none") +
            ggtitle("Conditional Dependence")
          plt_h0 <- plts$plt_null + theme(legend.position = "none") +
            ggtitle("Conditional Independence")

          plt_adj <- plot_adjacency(incon(), ... ) +
            theme(legend.position = "top",
                  legend.direction = "vertical") +
            scale_color_brewer(name = "",
                               palette = temp())

          leg <- cowplot::get_legend(plt_adj)

          plt_adj <- plot_adjacency(incon(),
                                    palette = temp(), ...) +
            theme(legend.position = "none") +
            ggtitle("Inconclusive")

          plts <- cowplot::plot_grid(plt_h1, plt_h0, plt_adj,
                                     rel_widths = c(1),
                                     nrow = 1)
          cowplot::plot_grid(plts,
                             leg,nrow = 1,
                             rel_widths = c(5, 1))

        }

      }

      output$plots <- renderPlot({
        req(plots())
        plots()

       })

      output$Download <- downloadHandler(
        filename = function(){
          paste('bggm_plot', '.pdf', sep = '')
        },
        content = function(file){
          req(plots())
          ggsave(file, plot = plots(), device = 'pdf', width = 10, height = 3.5)
        }
      )




    })
  shiny::shinyApp(ui = ui, server = server)
}



#' \code{shiny_bggm} Generic
#'
#' @param object object of class \code{explore}
#' @param ... currently ignored
#'
#' @export
shiny_bggm <- function(object,...){
  UseMethod("shiny_bggm", object)
}

