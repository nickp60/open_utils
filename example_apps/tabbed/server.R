
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)
library(reshape2)

shinyServer(function(input, output) {

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
    })+
  # output$pc1.3 <-renderPlot({
  #   numeric_columns<-lapply(diamonds, class)
  #   numeric_columns<-colnames(diamonds)[numeric_columns=="numeric"]
  #   subset<-diamonds[,input$fields]
  #   pca <- prcomp(subset, scale=T)
  #   melted <- cbind(subset, melt(pca$rotation[,1:length(pca$sdev)]))
  #   scores <- cbind(melted, pca$x[,1:length(pca$sdev)])
  #   qplot(x=PC1, y=PC3, data=scores, colour=factor(price)) +
  #     theme(legend.position="none")})
  # output$pc2.3 <-renderPlot({
  #   numeric_columns<-lapply(diamonds, class)
  #   numeric_columns<-colnames(diamonds)[numeric_columns=="numeric"]
  #   subset<-diamonds[,input$fields]
  #   pca <- prcomp(subset, scale=T)
  #   melted <- cbind(subset, melt(pca$rotation[,1:length(pca$sdev)]))
  #   scores <- cbind(melted, pca$x[,1:length(pca$sdev)])
  #  qplot(x=PC2, y=PC3, data=scores, colour=factor(price)) +
  #     theme(legend.position="none")})
  #    
  

})
