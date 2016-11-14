library(deSolve)
library(maemod)

#maxtime <- 1000
maxtime <- 10000

out <- maemodrun("D:\\Dropbox\\IBM project_Sai\\ODE\\SIRSI.txt", timegrid = seq(0,maxtime,1)) #scenario2
#out <- maemodrun("SIRSI - Copy.txt", timegrid = seq(0,10000,1)) #scenario1
#out <- maemodrun("Scenario1.txt", timegrid = seq(0,10000,1))
head(out)
tail(out)

#humans
plot(out[,2]+out[,6], type='l', ylim=c(0,max(c(out[,2],out[,4]))), col="blue", main="Humans")
lines(out[,4], col="red")
lines(out[,7]+out[,8], col="purple")

#Mosquitos
plot(out[,9], type='l', ylim=c(0,max(c(out[,9],out[,10]))), col="blue", main="Mosquitos")
lines(out[,10], col="red")

#lambda
plot(lambda_H, type='l', ylim=c(0,max(c(lambda_H,lambda_M))), col="blue", main="lambda, blue=humans")
lines(lambda_M, col="red")

####################
####inside 1########
####################

#incidence #calculated from inside
plot(inc, type='l', col="blue", main="Incidence")

#monthly incidence #calculated from inside
period <- length(inc)/(maxtime/30)

inc_mnth <- unname(tapply(inc, (seq_along(inc)-1) %/% period, sum))

plot(inc_mnth, type='l', main = "Incidence per month (not per 1000)")


#prevalence (R_T+I_DA)/N #calculated from inside
plot(prev, type='l', col='blue', main="Prevalence (R_T+I_DA)/N")

#incidence vs prevalence #calculated from inside
plot(prev,inc, main="incidence vs prevalence")

####################
####outside1########
####################

#incidence calculated the normal way
#daily incidence
N_inc <- (out[,4]*.7*(1/20)) #*input$rep*input$gamma) 

#daily incidence
plot(N_inc, type = 'l', main = "daily incidence")

#changing daily incidence into monthly
N_inc_mnth <- unname(tapply(N_inc, (seq_along(N_inc)-1) %/% 30, sum))

#monthly incidence
plot(N_inc_mnth, type = 'l', main = "monthly incidence")

####################
####inside 2########
####################

#inside 2 is calcualted from ci, which is a state variable
ci <- out[,11] #cumulative incidence (per day)
inc2 <- ci[-1] - ci[-length(ci)] #incidence calc from inside 2 (daily)
#plot(inc2)
inc2_mnth <- unname(tapply(inc2, (seq_along(inc2)-1) %/% 30, sum)) #incidence calc from inside 2 (monthly)

#prevalence inside 2
prev2 <- (out[,5] + out[,8]) / (out[,2] + out[,3] +out[,4] +out[,5] +out[,6] +out[,7] +out[,8]) #incidence per timestep
#prev2 <-  (out[,3]) /(out[,2] + out[,3] +out[,4] +out[,5] +out[,6] +out[,7] +out[,8]) #incidence per timestep
prev2_mnth <- prev2[seq(1,length(prev2),30)] #incidence per month

#incidence 2 vs prevalence 2
PrevVsInc2 <- cbind(prev2_mnth,inc2_mnth)
#PrevVsInc2 <- PrevVsInc2[-(1:9),] #removing the noise from 1st epidemics

PVI2o <- PrevVsInc2[order(PrevVsInc2[,1]),] #ordered version
plot(PVI2o, type='l')

#prevelance 2 plot
plot(prev2_mnth, type='l')

#lambda_H
# (700/320)*3*.3*(200/700)
# 
# out1 <- maemodrun("ex1.txt", timegrid = seq(0,1000,1))
# out2 <- maemodrun("ex2.txt", timegrid = seq(0,1000,1))
# 
#out3 <- maemodrun("ex3.txt", timegrid = seq(0,800,1))
for(i in 2:10){
  print(min(out[,i])  )
}

