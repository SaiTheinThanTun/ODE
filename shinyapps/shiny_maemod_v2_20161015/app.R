library(deSolve)
library(shiny)
library(maemod)
ui <- fluidPage(
  fluidRow(column(4, plotOutput(outputId = "humans")),
           column(4, plotOutput(outputId = "mosquitos")),
           column(4, plotOutput(outputId = "lambda"))
           ),
  fluidRow(column(4,
  tags$textarea(id="model", rows=30, cols=100, 
"
#SIRSI
!Inits
S=1000, E=0, I_S=1, R_T=0, S_I=0, I_UA=0, I_DA=0, S_M=500, I_M=200, ci=0

!Parameters
mu_i = 1/(55*364), #birth rate of human
mu_S = 1/(55*364), #normal death rate of human
omega_0 = 1/(5*364), #total loss of immunity (5 years), #then the model will need to run over 5 years
mu_E = 1, #same death rate as susceptibles
f = 1/10, #duration of E is approximately 10 days, Collins
mu_IS = 3, # increased death rate due to untreated malaria, 3 times higher chance
#gamma = 1/(10), # recovery with treatment
#gamma = 1/(20), # natural recovery, approximately from Collins
gamma = 1/(20), #recovery rate for both natural recovery and treatment
rep = .70, # treatment+report coverage
mu_RT = 1.0001, # increased death rate from treatment, 1.0001 times
sigma = 1/(6*30), # loss of prophylactic effect of drugs/treatment
mu_SI = 1, # same death rate as mu_S
xi1 = .3, #30% of reinfected individuals become sympatomatic
xi2 = .6, #60% of reinfected & Asymptomatic individuals become detectable with some Dx tests
#lam_SI = .60, #FOI to become truely sympatomatic infection is reduced by 40%
#lam_IUA = .85, #FOI of Undectable Asymptomatic infection becoming symptomatic again is decreased by 15% 
#lam_IDA = .75, #FOI of Dectable Asymptomatic infection becoming symptomatic again is decreased by 25%
omega_IUA = 1/365, # loss of all parasites at one year after becoming subpatent
omega_IDA = 1/100, # loss of detectability of parasites
mu_IUA = 1, #same death rate as mu_S
mu_IDA = 1, #same death rate as mu_S

beta = 1/10, #birth rate of mosquitos
kappa_SM = 1/10, #death rate of mosquitos
kappa_IM = 1/10, #in search of equilibrium :) 1/8, #death rate of infected mosquitos

a = 3, #human feeding rate
b = .006, #.01, #probability of disease transmission per bite for humans
c = .02 #.03 #probability a mosquito becomes infected after biting an infected human


!Equations
N <- S + E + I_S + R_T + S_I + I_UA + I_DA
M <- S_M + I_M
m <- M/N #ratio of mosquito to human
x <- (I_S + I_UA + I_DA)/N #ratio of infected humans
z <- I_M/M #ratio of infected mosquitos
lam_M <- a*c*x
lam_H <- m*a*b*z #FOI for humans

lambda_M <<- c(lambda_M,lam_M)
lambda_H <<- c(lambda_H,lam_H)

inc_i <- rep*gamma*I_S #incidence per timestep
inc <<- c(inc, inc_i) #incidence value written to outside vector 

ci <- rep*gamma*I_S #incidence calculated to the state variable

prev_i <- (R_T+I_DA)/N
prev <<- c(prev,prev_i)

dS <- N*mu_i + S_I*omega_0 - S*mu_S - S*lam_H
dE <- lam_H*S - mu_S*mu_E*E - f*E
dI_S <- f*E + S_I*lam_H*xi1 - I_S*mu_S*mu_IS - I_S*gamma*rep - I_S*gamma*(1-rep) + lam_H*xi1*I_DA + lam_H*xi1*I_UA
dR_T <- I_S*gamma*rep - R_T*mu_S*mu_RT - R_T*sigma
dS_I <- R_T*sigma + omega_IUA*I_UA - S_I*mu_S*mu_SI - S_I*omega_0 - S_I*lam_H*xi1 - S_I*lam_H*xi2*(1-xi1) - S_I*lam_H*(1-xi2)*(1-xi1)
dI_UA <- omega_IDA*I_DA + S_I*lam_H*(1-xi2)*(1-xi1) - lam_H*(1-xi1)*I_UA - omega_IUA*I_UA - lam_H*xi1*I_UA - I_UA*mu_S*mu_IUA 
dI_DA <- I_S*gamma*(1-rep) + S_I*lam_H*xi2*(1-xi1) + lam_H*(1-xi1)*I_UA - lam_H*xi1*I_DA - omega_IDA*I_DA - I_DA*mu_S*mu_IDA

dS_M <- M*beta - S_M*kappa_SM - S_M*lam_M
dI_M <- S_M*lam_M - I_M*kappa_IM


!Outputs
dS, dE, dI_S, dR_T, dS_I, dI_UA, dI_DA, dS_M, dI_M, ci

!ExtraFunctions
lambda_H <- c()
lambda_M <- c()
inc <- c()
prev <- c()


!MAEMOD_End
")),
###################
####comparison#####
###################
#column(4, plotOutput(outputId = "monthly_inc_compare")), #monthly incidence comparison
column(4, plotOutput(outputId = "incVsprev")),
#column(4, plotOutput(outputId = "incidence")), #incidence by timesteps 
column(4, plotOutput(outputId = "PrevVsInc2_plot"))),
fluidRow(column(4,""),
         column(4, plotOutput(outputId = "monthly_inc")),
         column(4, plotOutput(outputId = "prev2_plot"))),
fluidRow(column(4,
                numericInput(inputId = "maxtime", label= "Max time (days)", value=10000), #1080),
                numericInput(inputId = "rep", label = "Report/Rx rate", value = .70),
                numericInput(inputId = "gamma", label = "Recovery rate", value = 1/(20))
                )),
fluidRow(h3("backup graphs"),
         column(4, plotOutput(outputId = "prevalence")))

)

server <- function(input, output) {
  output$model <- renderText(input$model)
  out <- reactive(maemodrun(input.text=input$model,timegrid=seq(0,input$maxtime,1)))
  output$humans <- renderPlot({
    #humans
    plot(out()[,2]+out()[,6], type='l', ylim=c(0,max(c(out()[,2],out()[,4]))), col="blue", main="Humans")
    lines(out()[,4], col="red")
    lines(out()[,7]+out()[,8], col="purple")
  })
  output$mosquitos <- renderPlot({
    #Mosquitos
    plot(out()[,9], type='l', ylim=c(0,max(c(out()[,9],out()[,10]))), col="blue", main="Mosquitos")
    lines(out()[,10], col="red")
  })
  output$lambda <- renderPlot({
    #lambda
    plot(lambda_H, type='l', ylim=c(0,max(c(lambda_H,lambda_M))), col="blue", main="lambda, blue=humans")
    lines(lambda_M, col="red")
  })
  # output$incidence <- renderPlot({
  #   #incidence #calculated from inside
  #   plot(inc, type='l', col="blue", main="Incidence, inside1")
  # })
  ####################
  ####inside 2########
  ####################
  
  #inside 2 is calcualted from ci, which is a state variable
  ci <- reactive(out()[,11]) #cumulative incidence (per day)
  inc2 <- reactive({ci()[-1] - ci()[-length(ci())] })#incidence calc from inside 2 (daily)
  inc2_mnth <- reactive({unname(tapply(inc2(), (seq_along(inc2())-1) %/% 30, sum))}) #incidence calc from inside 2 (monthly)
  
  #prevalence inside 2
  prev2 <- reactive({(out()[,5] + out()[,8]) / (out()[,2] + out()[,3] +out()[,4] +out()[,5] +out()[,6] +out()[,7] +out()[,8])})  #incidence per timestep
  prev2_mnth <- reactive(prev2()[seq(1,length(prev2()),30)]) #incidence per month
  
  ###################
  ####comparison#####
  ###################
  # output$monthly_inc_compare <- renderPlot({
  #   #monthly incidence #calculated from inside
  #   period <- length(inc)/(input$maxtime/30)
  #   
  #   inc_mnth <- unname(tapply(inc, (seq_along(inc)-1) %/% period, sum))
  #   
  #   inc_mnth_n <- unname(tapply(N_inc_R(), (seq_along(N_inc_R())-1) %/% 30, sum))
  #   
  #   inc2_mnth <- unname(tapply(inc2(), (seq_along(inc2())-1) %/% 30, sum)) #incidence calc from inside 2 (monthly)
  #   
  #   
  #   plot(inc_mnth, type='p', main="Incidence per month (not per 1000), \nblue=inside1, red=outside, \n green=inside2", col="blue")
  #   lines(inc_mnth_n,type='p', col="red")
  #   lines(inc2_mnth, type='l', col="green")
  # })
  
  output$monthly_inc <- renderPlot({
    plot(inc2_mnth(), type='l', main="Incidence per month per 1000", col="blue")
  })
  output$prevalence <- renderPlot({
    #prevalence (R_T+I_DA)/N #calculated from inside1
    plot(prev, type='l', col='blue', main="Prevalence (R_T+I_DA)/N")
  })
  output$incVsprev <- renderPlot({
    #incidence vs prevalence #calculated from inside1
    plot(prev2()[-1],inc2(), main="incidence vs prevalence per day", type='l')
    #plot(prev2_mnth(),inc2_mnth(), main="incidence vs prevalence per mnth", type ='l')
  })
  
  #incidence calculated the normal way
  N_inc_R <- reactive(out()[,4]*input$rep*input$gamma) 
  

  output$PrevVsInc2_plot <- renderPlot({
    #incidence 2 vs prevalence 2
    PrevVsInc2 <- cbind(prev2_mnth(),inc2_mnth())
    #PrevVsInc2 <- PrevVsInc2[-(1:9),] #removing the noise from 1st epidemics
    
    PVI2o <- PrevVsInc2[order(PrevVsInc2[,1]),] #ordered version
    plot(PVI2o, type='l')
  })
  
  #prevelance 2 plot
  output$prev2_plot <- renderPlot({
    plot(prev2_mnth(), type='l')
  })
  
  

}

shinyApp(ui = ui, server = server)

