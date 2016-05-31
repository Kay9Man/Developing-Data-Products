#ui.R
library(shiny)
ui = shinyUI(pageWithSidebar(
        headerPanel('Comprehensive Learning Particle Swarm Optimization'),
        sidebarPanel(
                
                
                helpText("For this simple test, the parameters have been limited as follows:
                         Number of Iterations = Min: 100, Max: 2000
                         Number of Population = Min: 10, Max: 50
                         Number of Dimensions = Min: 2, Max: 50. As guessed,
                         the bigger the numbers, the bigger the problem!"
                ),

                    helpText("The Schwefel is a complex function is for a bigger dimensional space
                         (-500,500) as against the simple test dimension we are using (-5.12,5.150).
                             Hence it looks almost like poor performnance but it is not."   
                         
                         
                ),
                
                numericInput('noIter', 'Number Of Iterations',
                             value = 100,min = 100, max = 2000),
                sliderInput('noPop', 'Number of Population',
                            min = 10, max = 50,value = 10),
                sliderInput('noDim', 'Number of Dimensions',
                            min = 2, max = 50,value = 5),
                selectInput("testProb", "Test Problem",
                            c("deJong",
                              "Rastrigin",
                              "Schwefel",
                              "Michalewicz",
                              "griewank",
                              "parabola",
                              "rosenbrock",
                              "ackley")),
                submitButton("Run")
        ),
        mainPanel(
                h5("Current Parameters"),
                tableOutput("values"),
                h5("Initial Best Fitness"),
                verbatimTextOutput("inBestFit"),
                h5("Initial Best Solution"),
                verbatimTextOutput("inBestSoln"),
                h5("Global Best Fitness"),
                verbatimTextOutput("finBestFit"),
                h5("Global Best Solution"),
                verbatimTextOutput("finBestSoln"),
                plotOutput('plot1')
              
        )
        ))