ui <- dashboardPage(
  dashboardHeader(title = "BGGM"),
  dashboardSidebar(),
  dashboardBody(fluidPage(
      # App title ----
      titlePanel("Bayesian Gaussian Graphical Models"),

      sidebarLayout(
        sidebarPanel(
          fileInput("file1", "Choose CSV File",
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv")
          ),

          checkboxInput("header", "Header", TRUE),
          tags$hr(),
          selectInput("region", "Method:",
                      choices= c("Estimate", "Hypothesis Testing", "GGM Compare") ),

        conditionalPanel(
          condition = "input.region == 'Estimate'",
          checkboxInput("header", "Analytic Solution", FALSE)
        ))))))


# ui <- fluidPage(
#   # App title ----
#   titlePanel("Bayesian Gaussian Graphical Models"),
#
#   sidebarLayout(
#     sidebarPanel(
#       fileInput("file1", "Choose CSV File",
#                 accept = c(
#                   "text/csv",
#                   "text/comma-separated-values,text/plain",
#                   ".csv")
#       ),
#
#       checkboxInput("header", "Header", TRUE),
#       tags$hr(),
#       selectInput("region", "Method:",
#                   choices= c("Estimate", "Hypothesis Testing", "GGM Compare") ),
#
#     conditionalPanel(
#       condition = "input.region == 'Estimate'",
#       checkboxInput("header", "Analytic Solution", FALSE)
#     )),
#
#
#     mainPanel(
#       tableOutput("contents")
#     )
#   ),
#
#   # # Sidebar panel for inputs ----
#   # sidebarPanel(
#   #   selectInput("region", "Method:",
#   #               choices= c("Estimate", "Hypothesis Testing", "GGM Compare"))
#   # ),
#
#   # Main panel for displaying outputs ----
#   mainPanel(
#     textOutput("test")
#   )
# )

library(shinydashboard)
