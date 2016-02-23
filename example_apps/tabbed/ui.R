
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)

shinyUI(navbarPage("Sample Shiny App",

  # Application title
#  titlePanel("Diamond selector"),
  tabPanel(title="Diamond Price Distribution",
          sidebarLayout(
            sidebarPanel(
              sliderInput("bins",
                          "Number of bins:",
                          min = 1,
                          max = 100,
                          value = 20),
              hr("Designed by Nick Waters 2016"),
              hr("(data from ggplot2's Diamonds dataset")
            ),
            # Show a plot of the generated distribution
            mainPanel(
              plotOutput("distPlot")
            )
          )
  ),
  tabPanel(title="Compare Features",
           sidebarLayout(
             sidebarPanel(
               fileInput('file1', 'Choose file to upload',
                         accept = c(
                           'text/csv',
                           'text/comma-separated-values',
                           'text/tab-separated-values',
                           'text/plain',
                           '.csv',
                           '.tsv'
                         )
               ),
               hr("Dataset"),
               selectizeInput("dataset", "Attribute, for X axis",
                              colnames(diamonds), 
                              selected = "price"),
               selectizeInput("x", "Attribute, for X axis",
                              colnames(diamonds), 
                              selected = "price"),
               selectizeInput("y", "Attribute, for X axis",
                              colnames(diamonds), 
                              selected = "caret"),
               selectizeInput("color", "Attribute, for X axis",
                              colnames(diamonds), 
                              selected = "clarity"),
               hr("Optional settings:").
               selecti
#               actionButton("goButton", "Go!")
               
               
             ),
             
             # Show a plot of the generated distribution
             mainPanel(
               plotOutput("compare")

             )
           )
      )
  )
)


#### Hosting
"
install.packages('devtools')

devtools::install_github('rstudio/rsconnect')


rsconnect::setAccountInfo(name='omicsguy-apps',
			  token='7A2FC72584D0E2259D09EBB343C0A570',
			  secret='<SECRET>')
library(rsconnect)
rsconnect::deployApp('~//to/your/app')
"
