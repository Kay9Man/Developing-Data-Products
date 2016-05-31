library(shiny)
library(ggplot2)
rm(list = ls())
# Define UI for random distribution application 
ui = shinyUI(pageWithSidebar(
        headerPanel('Comprehensive Learning Particle Swarm Optimization'),
        sidebarPanel(
                
                
                helpText("For this simple test, the parameters have been limited as follows:
                         Number of Iterations = Min: 100, Max: 2000
                         Number of Population = Min: 10, Max: 50
                         Number of Dimensions = Min: 2, Max: 50. As guessed,
                         the bigger the numbers, the bigger the problem!"
                         
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
#                 
#                 
#                
#                 plotOutput('plot1')
        )
))

deJong <- function(x){
        ans = sum(x^2)
        return(ans)
}

Rastrigin <- function(x) {
        n=length(x)
        ans2 = 10*n + sum(x^2 - 10*cos(2*pi*x))
        return(ans2)
}

Schwefel <- function(xx){
  d <- length(xx)
  
  sum <- sum(xx*sin(sqrt(abs(xx))))
  
  y <- 418.9829*d - sum
  return(y)
}


Michalewicz <- function(xx, m=10){ 
        
        ii <- c(1:length(xx))
        sum <- sum(sin(xx) * (sin(ii*xx^2/pi))^(2*m))
        
        y <- -sum
        return(y)
}

griewank <- function(xx) {
  ii <- c(1:length(xx))
  sum <- sum(xx^2/4000)
  prod <- prod(cos(xx/sqrt(ii)))
  
  y <- sum - prod + 1
  return(y)
}
parabola <- function(x) {
        return(sum(x*x))
}

rosenbrock <- function(xx) {
  d <- length(xx)
  xi <- xx[1:(d-1)]
  xnext <- xx[2:d]
  
  sum <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
  
  y <- sum
  return(y)
}
ackley <- function(xx, a=20, b=0.2, c=2*pi){
  d <- length(xx)
  
  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))
  
  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)
  
  y <- term1 + term2 + a + exp(1)
  return(y)
}


Toroidal <- function(Initial,a,b){
        #Initial <- matrix(runif(Initial),nrow=1,byrow = TRUE)
        Adjusted <- matrix(0,nrow=1,length(Initial))
        for (i in 1:length(Initial)){
                if (Initial[i]>=a && Initial[i]<=b){
                        Adjusted[i] <- Initial[i] 
                } else if (Initial[i]>b) {    
                        delta <-Initial[i]-b
                        Adjusted[i]<-a+delta
                }
                if (Adjusted[i]>b){
                        Adjusted[i] <- b
                } else if (Initial[i]<a) {
                        delta <- a-Initial[i]
                        Adjusted[i] <- a+delta
                }
                if(Adjusted[i]<a){
                        Adjusted[i] <- a
                } else if( Adjusted[i]> b) {
                        Adjusted[i] <- b
                }
        }
        return (Adjusted)
}



server = shinyServer(function(input,output,session) {
       

  testProbs <- reactive({
    data.frame(
      Problem = input$testProb,Population = input$noPop,
      Dimension = input$noDim, Iterations = input$noIter,
      stringsAsFactors = FALSE
    )
  })

  output$values <- renderTable({
    testProbs()
  })
  problem <- reactive({
    switch(input$testProb,
                              "deJong" = deJong,
                              "Rastrigin" = Rastrigin,
                              "Schwefel" = Schwefel,
                              "Michalewicz" = Michalewicz,
                              "griewank" = griewank,
                              "parabola" = parabola,
                              "rosenbrock" = rosenbrock,
                              "ackley" = ackley)
  })
        output$inBestFit <- renderPrint({
                func=problem()
                lower=-5.12; 
                upper=5.12;
                numPopulation=input$noPop
                dimPopulation=input$noDim
                VMax=(upper-lower)*0.25; 
                VMin=-VMax;
                Budget=input$noIter
                WMax=1;
                WMin=0.1;
                PCMax = 0.5
                PCMin = 0.05
                C1=0.4;
                C2=0.8;
                
                final <- numeric(Budget)
                xlb <- x2<- v <- x <- array(0,c(numPopulation,dimPopulation))
                
                #matrix(0,numPopulation,dimPopulation,byrow = T)
                Fitxlb <- Fitx  <- Fitx2 <- array(0,c(0,numPopulation))
                #matrix(0,0,numPopulation,byrow = T)
                #array(runif(dimPopulation),c(1,dimPopulation))
                i <- 1:numPopulation
                PC <- PCMin + (PCMax - PCMin) * ((exp(10*(i-1) / (numPopulation - 1))-1) / (exp(10)-1))
                
                for (i in 1:numPopulation){
                        #%*% (matrix(runif(M),nrow=1,byrow = T))
                        # (matrix(runif(dimPopulation),nrow=1,byrow=T))
                        x[i, ] <- lower + (matrix(runif(dimPopulation),nrow=1,byrow=T)) * (upper - lower) #position
                        v[i, ] <- VMin + (matrix(runif(dimPopulation),nrow=1,byrow=T)) * (VMax - VMin) # initialize velocity
                        
                        x2[i, ] <- v[i, ] + x[i, ]
                        x2[i, ] <- Toroidal(x2[i, ],upper,lower)
                        Fitx[i] <- func(x[i, ])
                        Fitx2[i] <- func(x2[i, ])
                        
                        if(Fitx2[i] < Fitx[i]){
                                xlb[i, ] <- x2[i, ]
                                Fitxlb[i] <- Fitx2[i]
                        }else{
                                xlb[i, ] <- x[i, ]
                                Fitxlb[i] <- Fitx[i]
                        }
                        
                }
                
                Fitxgb <- min(Fitx); gbIndex <- which.min(Fitx)
                xgb <- xlb[gbIndex, ]
                initialBestFit = Fitxgb
                initialBestSoln = xgb
                print(initialBestFit)
                for (k in 1:Budget){
                        final[k] <- Fitxgb
                      
                        w <- WMax - (WMax - WMin) * k/Budget
                        for (i in 1:numPopulation){ 
                                Index <- 1:numPopulation
                                y <- sample(Index[-i],2) # index of any two apart from the current one
                                y1<- y[1]; y2<- y[2]
                                ypc <- sample(dimPopulation,1)
                                Allsame = 0
                                for (ii in 1:dimPopulation){
                                        r1 <- (matrix(runif(1),nrow=1,byrow = T))
                                        if (r1 > PC[i] || ii==ypc) {
                                                Allsame = Allsame +1
                                                if (isTRUE(all.equal(Allsame,dimPopulation))) {  #ensure learning from itself
                                                        v[i,ii] <- w%*%v[i,ii] + C1%*% r1 %*%(xlb[i,ii] - x[i,ii]) + C2%*% r1 %*%(xgb[ii]-x[i,ii])
                                                        #else
                                                        if (Fitxlb[y1] < Fitxlb[y2]){
                                                                v[i,ii] <- w%*%v[i,ii] + r1%*%(xlb[y1,ii] - x[i,ii]) + C2%*%r1%*%(xgb[ii] - x[i,ii])
                                                        } else {  
                                                                v[i,ii] <- w%*%v[i,ii] + r1%*%(xlb[y2,ii] - x[i,ii]) + C2%*%r1%*%(xgb[ii] - x[i,ii])
                                                        }
                                                }
                                                else
                                                        if (Fitxlb[y1] < Fitxlb[y2]){
                                                                v[i,ii] <- w%*%v[i,ii] + r1%*%(xlb[y1,ii] - x[i,ii]) + C2%*%r1%*%(xgb[ii] - x[i,ii])
                                                        } else { 
                                                                v[i,ii] <- w%*%v[i,ii] + r1%*%(xlb[y2,ii] - x[i,ii]) + C2%*%r1%*%(xgb[ii] - x[i,ii])
                                                        }
                                        }
                                }
                                
                                for (isii in 1:dimPopulation){
                                        if (v[i,isii] < VMin){ 
                                                v[i,isii] <- VMin
                                        } else if(v[i,isii] > VMax) { 
                                                v[i,isii] <- VMax
                                        }
                                }
                                
                                x[i, ] <- x[i, ] + v[i, ] #some particles may fall outside the decision space
                                # if (sum(x[i, ] < lower) + sum(x[i, ] > upper)) == 0{
                                
                                if(all(x[i, ] >= lower) && all(x[i, ] <= upper)) { #computes within the decision space
                                        
                                        Fitx_1 <- func(x[i, ])
                                        
                                        if (Fitx_1 <= Fitxlb[i]){
                                                xlb[i, ] <- x[i, ]
                                                Fitxlb[i] <- Fitx_1
                                        }
                                        if (Fitxlb[i] <= Fitxgb){
                                                xgb <- x[i, ]
                                                Fitxgb <- Fitxlb[i]
                                        }
                                        
                                }
                        }
                        
                        
                }
                
               
                output$inBestSoln <- renderPrint({
                         print(initialBestSoln)
                                         })
                output$finBestFit <- renderPrint({
                        print(Fitxgb) 
                                 })
                output$finBestSoln <- renderPrint({
                        print(xgb)
                                 })
                output$plot1 <- renderPlot({
                                  df <- data.frame(Iterations = 1:Budget,BestValue = final)
                                  #prob <- paste("Minimization Results for ", problem())
                                  ggplot(df,aes(Iterations,BestValue))+
                                    geom_line(linetype="dashed",color = "red")+
                                    geom_point()+
                                    ggtitle("Minimization Results")
                                       #qplot(1:Budget,final,geom = "line")
                                 })
              
        })
                        
        
        
        
})

shinyApp(ui,server)
