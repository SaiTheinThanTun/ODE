library(deSolve)
library(shiny)

function(input, output) {
  parameters_r <- reactive({
    c(mui=input$mui,    # birth #lifespan of mosquito 10 days
      muo=input$muo,    # death
      #beta= #per capita effective contact with infected human per unit time
      #ce = (.3*.01), #probability of disease transmission per bite * biting rate
      a = input$a, #human blood feeding rate
      b = input$beta_h, #probability of disease transmission per bite for human
      c = input$c, ##probability of disease transmission per bite for mosquitos
      #beta = input$beta, #probability of disease transmission per bite for mosquitos
      durinf = input$durinf * 2, #duration of infection in days * 2, because of 1/2 day time step
      immunity = input$immunity * 2 #duration of immunity * 2, because of 1/2 day time step
    )})
  
  lam_h <- reactive({
    parameter <- parameters_r()
    #lam_h <- m*a*b*z
    lam_h <- (input$M/input$H)*parameter[3]*parameter[4]*(input$Z/input$M)
  })
  
  ode_out <- reactive({
    # define the number of weeks to run the model
    times <- seq(0, input$days, by = (1/2)) #previously week was 52, 7 days/14 = 1/2 day timestep
    
    #MODEL PARAMETERS
    parameters <- parameters_r()
    
    # MODEL INITIAL CONDITIONS
    H <- input$H
    X <- input$X
    #H <- input$initPh #total human population
    #X <- input$initIh #no of infected people, ideally this should be changing either from the IBM or the ODE model itself
    #page 75 of A biologist's guide to mathematical modelling
    
    #
    #
    #initRh <- 0 #no of recovered individuals
    initSh <- H-X  #+initRh) #no of suscepitables
    
    M<- input$M 
    Z<-input$Z
    initS<-M-Z
    initD <- 0
    
    ##
    lam_h <- 0
    lam <- 0
    
    state <- c(Sh= initSh, X = X, lam_h = lam_h, S = initS, Z = Z, lam=lam,Y=0)
    
    
    # set up a function to solve the model
    mosQ<-function(t, state, parameters) 
    {
      with(as.list(c(state, parameters)),
           {
             
             # define variables
             M <- (S+Z)
             H <- (Sh+X)
             m <- M/H #ratio of mosquitos to humans
             z <- Z/M #ratio of infectious mosquitos
             x <- X/H #ratio of infectious humans
             
             #seas<-1+amp*cos(2*pi*(Y-phi)/52)
             #beta<-R0*(muo+nui)*gamma/(muo+gamma)
             lam <- a*c*x
             lam_h <- m*a*b*z
             
             # rate of change for mosquitos
             dS <- mui*M-muo*S-lam*S
             dZ <- -muo*Z+lam*S
             #dD <- muo*S+muo*Z #remove this!!!!
             #dM <- 0
             
             # rate of change for humans
             dSh <- -lam_h*Sh+(1/durinf)*X
             dX <- lam_h*Sh-(1/durinf)*X
             #dRh <- (1/durinf)*X-(1/immunity)*Rh
             dY <- 1
             
             # return the rate of change
             list(c(dSh, dX,lam_h, dS, dZ, lam, dY))
           }
      ) 
      
    }
    ode(y = state, times = times, func = mosQ, parms = parameters)
    
    #out <- ode(y = state, times = times, func = mosQ, parms = parameters)
    #out[,]
  })
  
  output$downloadODE <- downloadHandler(
    #out <- ode_out()
    filename= function(){paste('summary_ode_',Sys.Date(),'.csv',sep='')},
    content= function(file){
      write.csv(ode_out(),file)
    }
  )
  
  output$everything_mosq <- renderPlot({
    out <- ode_out()
    par(mar=c(5,4,4,4)) #default is par(mar=c(5,4,4,2))
    
    plot(out[,1],out[,5], type="l", col="blue", axes=FALSE, xlab="", ylab="", main="mosq_pop")
    axis(2, ylim=c(0,17),col="blue") 
    mtext("Susceptible mosquitoes",side=2,line=2.5) 
    box()
    par(new=TRUE)
    plot(out[,1],out[,6], type="l", col="red", axes=FALSE, xlab="", ylab="")
    axis(4, ylim=c(0,17),col="red") 
    mtext("Infected mosquitoes",side=4,line=2.5)
    
    axis(1,pretty(range(out[,1]),10))
    mtext("Time (days)",side=1,col="black",line=2.5)
    
    legend("top",legend=c("Susceptibles","Infected"),
           text.col=c("blue","red"),pch= "__", col=c("blue","red"))
  })
  
  output$human_pop <- renderPlot({
    out <- ode_out()
    lam_h <- lam_h()
    par(mar=c(5,4,4,4))
    plot(out[,1],out[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("human_pop with lambda ",lam_h()))
    axis(2, ylim=c(0,17),col="blue") 
    mtext("Susceptible humans",side=2,line=2.5) 
    
    box()
    par(new=TRUE)
    plot(out[,1],out[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
    axis(4, ylim=c(0,17),col="red") 
    mtext("Infected humans",side=4, line=2.5)
    
    axis(1,pretty(range(out[,1]),10))
    mtext("Time (days)",side=1,col="black",line=2.5)
    
    legend("top",legend=c("Susceptibles","Infected"),
           text.col=c("blue","red"),pch= "__", col=c("blue","red"))
  })
}