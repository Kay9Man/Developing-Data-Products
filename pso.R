#Particle Swarm Optimization (PSO)
#rm(list=ls())
setwd("/Users/Oracle/Dropbox/Experiment5_ICAST")
source("Toroidal.R")
source("deJong.R")
source("Rastrigin.R")
source("Schwefel.R")
source("Michalewicz.R")
a=-5.12; b=5.12;
N=30
M=10
Vmax=0.25*(b-a); Vmin=-Vmax;
Kmax=1000;
wmax=1;
wmin=0.1;
c1=0.4;
c2=0.8;

final <- numeric(Kmax)
xlb <- x2<- v <- x <- matrix(0,N,M,byrow = T)
Fitxlb <- matrix(1,1,N,byrow = T)
for ( i in 1:N){
  
  x[i, ] <- a+(b-a) %*% (matrix(runif(M),nrow=1,byrow = T))
  v[i,  ] <- Vmin+(Vmax-Vmin) %*% (matrix(runif(M),nrow=1,byrow = T))
  x2[i,  ] <- (x[i, ]) + (v[i, ])
  x2[i,  ] <- Toroidal((x2[i, ]), a,b)
  Fitx <- Rastrigin(x[i,  ])
  Fitx2 <- Rastrigin(x2[i,  ])
  
  if (Fitx<Fitx2){
    xlb[i,  ] <- x[i,  ]
  } else {
    #Fitxlb[i] <- Fitx  ##
    xlb[i,  ] <- x2[i,  ]
  }
   Fitxlb[i] <- Rastrigin(xlb[i,  ]) ##
}
Bfit <- min(Fitxlb);BIndex <- which.min(Fitxlb)

xgb <- xlb[BIndex, ]
InitialBestFit <- min(Fitxlb)
#print(InitialBestFit)
cat("Intial Best Fit is: ",InitialBestFit,"\n")
cat("Intial Best solution is: ",xgb,"\n")

for (k in 1:Kmax) {
  #Update Local particle best position
  for (i in 1:N){
    Fxx <- Rastrigin(x[i, ])
    Flb <- Rastrigin(xlb[i, ])
    if (Fxx<Flb){
      xlb[i, ] <- x[i, ]
    }
    #update global best
    Flb <- Rastrigin(xlb[i, ])
    Fgb <- Rastrigin(xgb)
    
    #cat("Best Value: ",Fgb,"\n")
    if (Flb < Fgb){
      xgb <- xlb[i, ]
      
      #FinalBestFit <- Rastrigin(xgb)
      #cat("Best Value: ",FinalBestFit,"\n")
    }
  }
  #update particle velocity and position
  w <- wmax - (wmax-wmin) %*% (k/Kmax)
  for (i in 1:N){
    m2 <- (matrix(runif(1),nrow=1,byrow = T))
    v[i, ] <- w %*% v[i, ] + c1 %*% m2 %*% (xlb[i, ]-x[i, ]) + c2 %*% m2 %*% (xgb - x[i, ])
    v[i, ] <- Toroidal(v[i, ],Vmin,Vmax)
    x[i, ] <- x[i, ] + v[i, ]
    x[i, ] <- Toroidal(x[i, ],a,b)
  }
  final[k] <- Rastrigin(xgb)
}
FinalBestFit <- Rastrigin(xgb)
cat("Final Best Fit is: ",FinalBestFit,"\n")
#cat("Local Best is: ",xlb,"\n")
cat("Global Best is: ",xgb,"\n")
#print(FinalBestFit)

#x1 <-  (matrix(runif(20),=20))
#for (i in 1:20){
 
# d <- x1[NULL,i] <- 5 + x1
 #print(d)
#}

#x<-matrix(rep(0,1000),nrow=1000,byrow=T)
#for(i in 1:1000){
#  x[i,]<-rexp(1,rate=2/3)
#}