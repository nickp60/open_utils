#Rscript

#test for installed packages
if(!"plotly" %in% rownames(installed.packages())){
  install.packages("plotly")
}
if(!"shiny" %in% rownames(installed.packages())){
  install.packages("shiny")
}
### load libraries
library(plotly)
library(shiny)
#input data

data<- iris

plot<-plot_ly(data, x = Sepal.Length, y = Sepal.Width, 
        z = Petal.Length,
        text = paste("Petal.Width: ",Petal.Width),
        color= Species,
        type = "scatter3d", mode = "markers")
# set other attributes
plot <- layout(plot,              
            title = "Test 3d Scatter via Plotly"
)
plot
