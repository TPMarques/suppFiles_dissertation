#convRJMCMC test script

#Uncomment the line bellow in case that devtools is not installed
#install.packages("devtools")

#Install package convRJMCMC
devtools::install_github("TPMarques/convRJMCMC")

#Load package convRJMCMC
library(convRJMCMC)

#Test plot for single parameter

set.seed(634963)
theta1<-rnorm(100000,mean=20,sd=1.5)
theta2<-rnorm(100000,mean=20,sd=1.5)
theta3<-rnorm(100000,mean=20,sd=1.5)
theta<-c(theta1,theta2,theta3)
chains<-c(rep(1,100000),rep(2,100000),rep(3,100000))
models<-factor(sample(c(1,2,3),300000,TRUE,c(0.25,0.5,0.25)))
mcmciterations<-rep(seq(1:100000),3)

#CZ_ANOVA object is saved on test_1
teste_1<-CZ_ANOVA(theta,chains,models,mcmciterations)

x11()
plot(teste_1)
#Note#
#You must press enter to see the other plots
###

#Test plot for multiple parameters

set.seed(500)

chains<-c(rep(1,100000),rep(2,100000),rep(3,100000))
models<-factor(sample(c(1,2,3),300000,TRUE,c(0.25,0.5,0.25)))
mcmciterations<-rep(seq(1:100000),3)

{
  theta1a<-rnorm(100000,mean=0,sd=1)
  theta2a<-rnorm(100000,mean=1,sd=2)
  theta3a<-rnorm(100000,mean=1,sd=3)
  thetaa<-cbind(theta1a,theta2a,theta3a)

  theta1b<-rnorm(100000,mean=0,sd=1)
  theta2b<-rnorm(100000,mean=1,sd=2)
  theta3b<-rnorm(100000,mean=1,sd=3)
  thetab<-cbind(theta1b,theta2b,theta3b)

  theta1c<-rnorm(100000,mean=0,sd=1)
  theta2c<-rnorm(100000,mean=1,sd=2)
  theta3c<-rnorm(100000,mean=1,sd=3)
  thetac<-cbind(theta1c,theta2c,theta3c)

  theta<-rbind(thetaa,thetab,thetac)
}

#CZ_MANOVA object is saved on test_2
teste_2<-CZ_MANOVA(theta,chains,models,mcmciterations)
x11()
plot(teste_2)

#Note#
#You must press enter to see the other plots
###

#See if the help text is ok
help("CZ_ANOVA")
help("CZ_MANOVA")
help("plot.CZ_ANOVA")
help("plot.CZ_MANOVA")

