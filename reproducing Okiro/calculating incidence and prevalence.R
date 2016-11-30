#capture and plot incidence vs prevalence based on various parameters
#and see if the model conforms to the findings in Okiro paper

#R denotes range (has to use runif)
#V denotes vector (use sample)
library(deSolve)
rm(list = ls())

####parameters####
#parameters in the global environment
init.pop <- c(S=1000, E=0, I_S=1, R_T=0, S_I=0, I_UA=0, I_DA=0, D = 1000) 
maxtime <- 10000
mu_i = 1/(55*364) #birth rate of human
mu_S = 1/(55*364) #normal death rate of human

pmu <- rate2prob(mu_i) #birth/death probability of human, RUNS FASTER if defined here

omega_0 = 1/(5*364) #total loss of immunity (5 years), #then the model will need to run over 5 years
mu_E = 1 #same death rate as susceptibles
f = 1/10 #duration of E is approximately 10 days, Collins
mu_IS = 3 # increased death rate due to untreated malaria, 3 times higher chance
gamma = 1/(20) #recovery rate for both natural recovery and treatment
rep = .70 # treatment+report coverage
mu_RT = 1.0001 # increased death rate from treatment, 1.0001 times
sigma = 1/(6*30) # loss of prophylactic effect of drugs/treatment
mu_SI = 1 # same death rate as mu_S
xi1 = .3 #30% of reinfected individuals become sympatomatic
xi2 = .6 #60% of reinfected & Asymptomatic individuals become detectable with some Dx tests
omega_IUA = 1/365 # loss of all parasites at one year after becoming subpatent
omega_IDA = 1/100 # loss of detectability of parasites
mu_IUA = 1 #same death rate as mu_S
mu_IDA = 1 #same death rate as mu_S
timecounter = 0
#mosquitos
#init
S_M=500
I_M=200
#parameters
beta_ = 1/10 #birth rate of mosquitos
kappa_SM = 1/10 #death rate of mosquitos
kappa_IM = 1/10 #in search of equilibrium :) 1/8, #death rate of infected mosquitos

a = 3 #human feeding rate
b = .006 #.01, #probability of disease transmission per bite for humans
c_ = .02 #.03 #probability a mosquito becomes infected after biting an infected human
#seasonality
seas_switch <- 1 #logical switch for seasonality
amp <- .2 #.7 or 1 Lisa's advice
phi <- 0
magnitude <- 1

####ODE everything under here####
#parameterizing for ODE
times <- seq(0, maxtime, by = 1)
mosquitostates <- c(S_M=S_M,I_M=I_M, ci=0) #has to be EXACTLY like this!!!!
param_ode <- c(
  mu_i = mu_i,
  mu_S = mu_S,
  omega_0 = omega_0,
  mu_E = mu_E,
  f = f,
  mu_IS = mu_IS,
  gamma = gamma,
  rep = rep,
  mu_RT = mu_RT,
  sigma = sigma,
  mu_SI = mu_SI,
  xi1 = xi1,
  xi2 = xi2,
  omega_IUA = omega_IUA,
  omega_IDA = omega_IDA,
  mu_IUA = mu_IUA,
  mu_IDA = mu_IDA,
  beta_ = beta_,
  kappa_SM = kappa_SM,
  kappa_IM = kappa_IM,
  a = a,
  b = b,
  c_ = c_,
  seas_switch = seas_switch,
  amp = amp,
  phi = phi,
  magnitude = magnitude
)


#transient values in ODE
transient_ode <- c("N <- S + E + I_S + R_T + S_I + I_UA + I_DA",
                   "M <- S_M + I_M", 
                   "m <- M/N",
                   "x <- (I_S + I_UA + I_DA)/N",
                   "z <- I_M/M",
                   "seas <- ((amp*cos(2*pi*(t-phi)/365)+magnitude)*seas_switch)+(1-seas_switch)",
                   "lam_M <- a*c_*x*seas",
                   "lam_H <- m*a*b*z")

#calculation for prevalence and incidence can be found in the files with Sompob codes

#model structures
#ODE
rateofchange <- c("dS <- N*mu_i + S_I*omega_0 - S*mu_S - S*lam_H",
                  "dE <- lam_H*S - mu_S*mu_E*E - f*E",
                  "dI_S <- f*E + S_I*lam_H*xi1 - I_S*mu_S*mu_IS - I_S*gamma*rep - I_S*gamma*(1-rep) + lam_H*xi1*I_DA + lam_H*xi1*I_UA",
                  "dR_T <- I_S*gamma*rep - R_T*mu_S*mu_RT - R_T*sigma",
                  "dS_I <- R_T*sigma + omega_IUA*I_UA - S_I*mu_S*mu_SI - S_I*omega_0 - S_I*lam_H*xi1 - S_I*lam_H*xi2*(1-xi1) - S_I*lam_H*(1-xi2)*(1-xi1)",
                  "dI_UA <- omega_IDA*I_DA + S_I*lam_H*(1-xi2)*(1-xi1) - lam_H*(1-xi1)*I_UA - omega_IUA*I_UA - lam_H*xi1*I_UA - I_UA*mu_S*mu_IUA ",
                  "dI_DA <- I_S*gamma*(1-rep) + S_I*lam_H*xi2*(1-xi1) + lam_H*(1-xi1)*I_UA - lam_H*xi1*I_DA - omega_IDA*I_DA - I_DA*mu_S*mu_IDA",
                  "dS_M <- M*beta_ - S_M*kappa_SM - S_M*lam_M",
                  "dI_M <- S_M*lam_M - I_M*kappa_IM",
                  "ci <- rep*gamma*I_S"
)
#ODE function
testfun <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)),
       {
         # define variables
         eval(parse(text = transient_ode))
         
         # rate of change
         eval(parse(text = rateofchange))
         
         # return the rate of change
         list(c(dS, dE, dI_S, dR_T, dS_I, dI_UA, dI_DA, dS_M, dI_M, ci))
       })
}

out <-  ode(y = c(init.pop[-length(init.pop)], mosquitostates), times = times, func = testfun,parms = param_ode)


#inside 2 is calcualted from ci, which is a state variable
ci <- out[,11] #cumulative incidence (per day)
inc2 <- ci[-1] - ci[-length(ci)] #incidence calc from inside 2 (daily)
inc2_mnth <- unname(tapply(inc2, (seq_along(inc2)-1) %/% 30, sum)) #incidence calc from inside 2 (monthly)

#prevalence inside 2
prev2 <- (out[,5] + out[,8]) / (out[,2] + out[,3] +out[,4] +out[,5] +out[,6] +out[,7] +out[,8])  #incidence per timestep
prev2_mnth <- prev2[seq(1,length(prev2),30)] #incidence per month

#####plot ODE ALONE####
plot(out[,2]+out[,6], ylim=c(0, sum(init.pop[-length(init.pop)])), type='l', col="blue", main="ODE alone")
lines(out[,4], col="red")
lines(out[,7]+out[,8], col="purple")
legend("top",legend=c("Susceptibles","Infected","Asymptomatic Infected"),
       text.col=c("blue","red","purple"),pch= "__", col=c("blue", "red","purple"))

#####plot prev vs incidence####
plot(prev2_mnth,inc2_mnth, type='l', main="Monthly incidence vs prevalence \n Can't take last point since it maynot be a full month")
tail(prev2_mnth)
tail(inc2_mnth) ######Can't take last point since it maynot be a full month!!!!!!!!!!!!!!!!###########