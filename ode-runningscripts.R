library(deSolve)
library(maemod)

maxtime <- 1000

out <- maemodrun("D:\\Dropbox\\IBM project_Sai\\ODE\\SIRSI.txt", timegrid = seq(0,maxtime,1)) #scenario2
#out <- maemodrun("SIRSI - Copy.txt", timegrid = seq(0,10000,1)) #scenario1
#out <- maemodrun("Scenario1.txt", timegrid = seq(0,10000,1))
head(out)
tail(out)

#humans
plot(out[,2], type='l', ylim=c(0,max(c(out[,2],out[,4]))), col="blue", main="Humans")
lines(out[,4], col="red")
lines(out[,7]+out[,8], col="purple")

#Mosquitos
plot(out[,9], type='l', ylim=c(0,max(c(out[,9],out[,10]))), col="blue", main="Mosquitos")
lines(out[,10], col="red")

#lambda
plot(lambda_H, type='l', ylim=c(0,max(c(lambda_H,lambda_M))), col="blue", main="lambda")
lines(lambda_M, col="red")

#incidence #calculated from inside
plot(inc, type='l', col="blue", main="Incidence")

#monthly incidence #calculated from inside
period <- length(inc)/(maxtime/30)

inc_mnth <- unname(tapply(inc, (seq_along(inc)-1) %/% period, sum))

plot(inc_mnth, type='l', main = "Incidence per month (not per 1000)")


#prevalence (R_T+I_DA)/N #calculated from inside
plot(prev, type='l', col='blue', main="Prevalence (R_T+I_DA)/N")

#incidence vs prevalence #calculated from inside
plot(inc,prev, main="incidence vs prevalence")


#incidence calculated the normal way
#daily incidence
N_inc <- (out[,4]*.7*(1/20)) #*input$rep*input$gamma) 

#daily incidence
plot(N_inc, type = 'l', main = "daily incidence")

#changing daily incidence into monthly
N_inc_mnth <- unname(tapply(N_inc, (seq_along(N_inc)-1) %/% 30, sum))

#monthly incidence
plot(N_inc_mnth, type = 'l', main = "monthly incidence")



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

