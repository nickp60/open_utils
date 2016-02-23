#Rscript

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

shinyServer(function(input, output) {
  output$plotly<-renderPlotly({
    data<- iris
    
    plot<-plot_ly(data, x = Sepal.Length, y = Sepal.Width, 
                  z = Petal.Length,
                  text = paste("Petal.Width: ",Petal.Width),
                  color= Species,
                  type = "scatter3d", mode = "markers")
    # set other attributes
    layout(plot, title = "Test 3d Scatter via Plotly")
  })
    
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- diamonds[, "price"]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white', xlab = "Price")
    
  })
  output$distPlot2 <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x2    <- diamonds[, "price"]
    bins2 <- seq(min(x2), max(x2), length.out = input$bins2 + 1)
    
    # draw the histogram with the specified number of bins
    hist(x2, breaks = bins2, col = 'darkgray', border = 'white', xlab = "Price")
    
  })
  ##################################################
  #  observeEvent(input$button, {
  #  input<-data.frame(fields=c("price", "carat", "x"), stringsAsFactors = F)
  output$compare <-renderPlot({
    ggplot(diamonds, aes_string(x=input$x, input$y, color=input$color))+geom_point()
  })

})

