
#test for installed packages
if(!"plotly" %in% rownames(installed.packages())){
  install.packages("plotly")
}
if(!"shiny" %in% rownames(installed.packages())){
  install.packages("shiny")
}
if(!"ggplot2" %in% rownames(installed.packages())){
  install.packages("ggplot2")
}
if(!"reshape2" %in% rownames(installed.packages())){
  install.packages("reshape2")
}
### load libraries

library(shiny)
library(ggplot2)
library(reshape2)
library(plotly)

shinyUI(navbarPage("Sample Shiny App",

  # Application title
#  titlePanel("Diamond selector"),
tabPanel(title="Plotly Demo with  Shiny",
         sidebarLayout(
           sidebarPanel(
             hr("integration of plotly and shiny using the iris dataset")
           ),
           # Show a plot of the generated distribution
           mainPanel(
             plotlyOutput("plotly")
           )
         )
),
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
               # fileInput('file1', 'Choose file to upload',
               #           accept = c(
               #             'text/csv',
               #             'text/comma-separated-values',
               #             'text/tab-separated-values',
               #             'text/plain',
               #             '.csv',
               #             '.tsv'
               #           )
               # ),
               # hr("Dataset"),
               # selectizeInput("dataset", "Attribute, for X axis",
               #                colnames(diamonds), 
               #                selected = "price"),
               selectizeInput("x", "Attribute, for X axis",
                              colnames(diamonds), 
                              selected = "price"),
               selectizeInput("y", "Attribute, for X axis",
                              colnames(diamonds), 
                              selected = "caret"),
               selectizeInput("color", "Attribute, for X axis",
                              colnames(diamonds), 
                              selected = "clarity")
               
             ),
             
             # Show a plot of the generated distribution
             mainPanel(
               plotOutput("compare")

             )
           )
      )
  )
)
