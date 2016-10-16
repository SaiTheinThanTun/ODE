
library(shiny)
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
S=300, E=0, I_S=1, R_T=0, S_I=0, I_UA=0, I_DA=0, S_M=500, I_M=200

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
kappa_IM = 1/8, #death rate of infected mosquitos

a = 3, #human feeding rate
b = .01, #probability of disease transmission per bite for humans
c = .03 #probability a mosquito becomes infected after biting an infected human


!Equations
N <- S + E + I_S + R_T + S_I + I_UA + I_DA
M <- S_M + I_M
m <- M/N #ratio of mosquito to human
x <- (I_S + I_UA + I_DA)/N #ratio of infected humans
z <- I_M/M #ratio of infected mosquitos
lam_M <- a*c*x
lam_H <- m*a*b*z

lambda_M <<- c(lambda_M,lam_M)
lambda_H <<- c(lambda_H,lam_H)

inc_i <- rep*gamma*I_S
inc <<- c(inc, inc_i)

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
dS, dE, dI_S, dR_T, dS_I, dI_UA, dI_DA, dS_M, dI_M

!ExtraFunctions
lambda_H <- c()
lambda_M <- c()
inc <- c()
prev <- c()


!MAEMOD_End
")),
column(4, plotOutput(outputId = "monthly_inc")), #monthly incidence
#column(4, plotOutput(outputId = "incidence")), #incidence by timesteps
column(4, plotOutput(outputId = "prevalence"))),
fluidRow(column(4,""),
         column(4, plotOutput(outputId = "incVsprev")),
         column(4, plotOutput(outputId = "N_inc_daily"))),
fluidRow(column(4,
                numericInput(inputId = "maxtime", label= "Max time (days)", value=1000),
                numericInput(inputId = "rep", label = "Report/Rx rate", value = .70),
                numericInput(inputId = "gamma", label = "Recovery rate", value = 1/(20))
                ))

)

server <- function(input, output) {
  output$model <- renderText(input$model)
  out <- reactive(maemodrun(input.text=input$model,timegrid=seq(0,input$maxtime,1)))
  output$humans <- renderPlot({
    #humans
    plot(out()[,2], type='l', ylim=c(0,max(c(out()[,2],out()[,4]))), col="blue", main="Humans")
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
    plot(lambda_H, type='l', ylim=c(0,max(c(lambda_H,lambda_M))), col="blue", main="lambda")
    lines(lambda_M, col="red")
  })
  output$incidence <- renderPlot({
    #incidence #calculated from inside
    plot(inc, type='l', col="blue", main="Incidence")
  })
  output$monthly_inc <- renderPlot({
    #monthly incidence #calculated from inside
    period <- length(inc)/(input$maxtime/30)
    
    inc_mnth <- unname(tapply(inc, (seq_along(inc)-1) %/% period, sum))
    
    plot(inc_mnth, type='l', main="Incidence per month (not per 1000)")
  })
  output$prevalence <- renderPlot({
    #prevalence (R_T+I_DA)/N #calculated from inside
    plot(prev, type='l', col='blue', main="Prevalence (R_T+I_DA)/N")
  })
  output$incVsprev <- renderPlot({
    #incidence vs prevalence #calculated from inside
    plot(inc,prev, main="incidence vs prevalence")
  })
  
  #incidence calculated the normal way
  N_inc_R <- reactive(out()[,4]*input$rep*input$gamma) 
  
  output$N_inc_daily <- renderPlot({
    #daily incidence
    plot(N_inc_R(), type = 'l', main = "daily incidence")
  })

}

shinyApp(ui = ui, server = server)

