library(deSolve)
library(shiny)

function(input, output) {
  parameters_r <- reactive({
    c(mui=input$mui,    # birth #lifespan of mosquito 10 days
      muo=input$muo,    # death
      #beta= #per capita effective contact with infected human per unit time
      #ce = (.3*.01), #probability of disease transmission per bite * biting rate
      a = input$a, #human blood feeding rate
      c = input$c, ##probability of disease transmission per bite for mosquitos
      #beta = input$beta, #probability of disease transmission per bite for mosquitos
      b = input$beta_h, #probability of disease transmission per bite for human
      durinf = input$durinf * 2, #duration of infection in days * 2, because of 1/2 day time step
      immunity = input$immunity * 2 #duration of immunity * 2, because of 1/2 day time step
    )})
  
  output$lam <- renderPrint({
    #parameter <- parameters_r()
    #lam <- parameter[3]*parameter[4]*(input$initIh/input$initPh)
    #lam
    print("xxx\n yyy")
  })
  
  ode_out <- reactive({
    # define the number of weeks to run the model
    times <- seq(0, input$weeks, by = (1/14)) #previously week was 52, 7 days/14 = 1/2 day timestep
    
    #MODEL PARAMETERS
    parameters <- parameters_r()
    
    # MODEL INITIAL CONDITIONS
    H <- input$initPh #total human population
    X <- input$initIh #no of infected people, ideally this should be changing either from the IBM or the ODE model itself
    #page 75 of A biologist's guide to mathematical modelling
    initRh <- 0 #no of recovered individuals
    initSh <- H-(X+initRh) #no of suscepitables
    
    initP<- input$initP 
    initI<-input$initI
    initS<-initP-initI
    initD <- 0
    
    state <- c(S = initS, Z = initI, D = initD, P = initP, Sh= initSh, X = X, Rh= initRh, Y=0)
    
    
    # set up a function to solve the model
    mosQ<-function(t, state, parameters) 
    {
      with(as.list(c(state, parameters)),
           {
             
             # define variables
             M <- (S+Z)
             H <- (Sh+X+Rh)
             m <- M/H #ratio of mosquitos to humans
             z <- Z/M #ratio of infectious mosquitos
             x <- X/H #ratio of infectious humans
             
             #seas<-1+amp*cos(2*pi*(Y-phi)/52)
             #beta<-R0*(muo+nui)*gamma/(muo+gamma)
             lam <- a*c*x
             lam_h <- m*a*b*z
             
             # rate of change for mosquitos
             dS <- mui*P-muo*S-lam*S
             dZ <- -muo*Z+lam*S
             dD <- muo*S+muo*Z
             dP <- 0
             
             # rate of change for humans
             dSh <- -lam_h*Sh+(1/immunity)*Rh
             dX <- lam_h*Sh-(1/durinf)*X
             dRh <- (1/durinf)*X-(1/immunity)*Rh
             dY <- 1
             
             # return the rate of change
             list(c(dS, dZ, dD,dP, dSh, dX, dRh, dY))
           }
      ) 
      
    }
    ode(y = state, times = times, func = mosQ, parms = parameters)
    
    #out <- ode(y = state, times = times, func = mosQ, parms = parameters)
    #out[,]
  })
  
  
  
  output$everything_mosq <- renderPlot({
    out <- ode_out()
    par(mar=c(5,4,4,4)) #default is par(mar=c(5,4,4,2))
    
    plot(out[,1],out[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main="mosq_pop")
    axis(2, ylim=c(0,17),col="blue") 
    mtext("Susceptible mosquitoes",side=2,line=2.5) 
    box()
    par(new=TRUE)
    plot(out[,1],out[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
    axis(4, ylim=c(0,17),col="red") 
    mtext("Infected mosquitoes",side=4,line=2.5)
    
    axis(1,pretty(range(out[,1]),10))
    mtext("Time (Weeks)",side=1,col="black",line=2.5)
    
    legend("top",legend=c("Susceptibles","Infected"),
           text.col=c("blue","red"),pch= "__", col=c("blue","red"))
  })
  
  output$human_pop <- renderPlot({
    out <- ode_out()
    par(mar=c(5,4,4,4)) #default is par(mar=c(5,4,4,2))
    
    plot(out[,1],out[,6], type="l", col="blue", axes=FALSE, xlab="", ylab="", main="human_pop")
    axis(2, ylim=c(0,17),col="blue") 
    mtext("Susceptible humans",side=2,line=2.5) 
    box()
    par(new=TRUE)
    plot(out[,1],out[,7], type="l", col="red", axes=FALSE, xlab="", ylab="")
    axis(4, ylim=c(0,17),col="red") 
    mtext("Infected humans",side=4,line=2.5)
    
    axis(1,pretty(range(out[,1]),10))
    mtext("Time (Weeks)",side=1,col="black",line=2.5)
    
    legend("top",legend=c("Susceptibles","Infected"),
           text.col=c("blue","red"),pch= "__", col=c("blue","red"))
  })
}