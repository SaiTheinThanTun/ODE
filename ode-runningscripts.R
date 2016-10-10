library(deSolve)
library(maemod)

out <- maemodrun("D:\\Dropbox\\IBM project_Sai\\ODE\\SIRSI.txt", timegrid = seq(0,10000,1)) #scenario2
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

