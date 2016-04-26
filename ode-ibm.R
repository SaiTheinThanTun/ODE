library(deSolve)
library(shiny)
library(ggplot2)

#required data
age_prob_0to97 <- read.csv("C:/wd/0to97_age_prob.csv", header=FALSE)
male_prob_0to97 <- read.csv("C:/wd/0to97_male_prob.csv", header=FALSE)
age <- 0:97

ui <- fluidPage(
  fluidRow(
    column(4,
           h3("Transmission parameters"),
           sliderInput(inputId = "mui",
                       label = "birth rate of mosquitos",
                       value = .05, min = 0, max = .1 #10 days survival= 20 half-days survival, therefore 1/20=.05
           ),
           sliderInput(inputId = "muo",
                       label = "death rate of mosquitos",
                       value = .05, min = 0, max = .1
           ),
           sliderInput(inputId = "durinf",
                       label = "duration of infection (days): durinf",
                       value = 30, min = 7, max = 365
           ) #, immunity is removed below
           #sliderInput(inputId = "immunity",
            #           label = "duration of immunity (days): immunity",
            #           value = 30, min = 7, max = 365
           #)
           ,numericInput(inputId = "no_sims",
                         label = "number of simulations",
                         value = 10)
    ),
    column(4,
           sliderInput(inputId = "a",
                       label = "human blood feeding rate: a",
                       value = .1, min = .01, max = .3
           ),
           #sliderInput(inputId= "beta",
            #           label = "probability of disease transmission per bite for mosquitos",
            #           value = .3, min = .01, max = 1
           #),
           sliderInput(inputId= "beta_h",
                       label = "probability of disease transmission per bite for humans: b",
                       value = .3, min = 0, max = 1 #change here
           ),
           sliderInput(inputId= "c",
                       label = "probability a mosquito becomes infected after biting an infected human: c",
                       value = .7, min = 0, max = 1 #change here
           )
           ),
    column(4,
           h3("Initial Conditions"),
           numericInput(inputId = "H",
                        label = "Human Population: H",
                        value = 80
           ),
           numericInput(inputId = "X",
                        label = "Infected Human Population: X",
                        value = 30
           ),
           numericInput(inputId = "M",
                        label = "Initial Mosquito Population: M",
                        value = 800
           ),
           sliderInput(inputId = "Z",
                       label = "Initial infected Mosquito Population: Z",
                       value = 200, min = 1, max=800
           ),
           numericInput(inputId= "days",
                        label= "No. of days for the model",
                        value = 28
           )
    )
  ),
  fluidRow(h2("ODE"),
    #column(5,
    #       h3("plot of inc2"),
    #       plotOutput(outputId = "inc2")
    #       ),
    column(5,
           h3("Human Population"),
           plotOutput(outputId="human_pop")
    ),
    column(5,
           h3("Mosquito Population"),
           plotOutput(outputId = "everything_mosq")
    ),
    column(2,h3("Download csv"),
#            p("m <- M/H #ratio of mosquitos to humans"),
#            p("z <- Z/M #ratio of infectious mosquitos"),
#            p("x <- X/H #ratio of infectious humans"),
#            p("lam <- a*c*x #FOI for mosquitos"),
#            p("lam_h <- m*a*b*z #FOI for humans"),
#            br(),
#            p("# rate of change for mosquitos"),
#            p("dS <- mui*M-muo*S-lam*S"),
#            p("dZ <- -muo*Z+lam*S"),
#            #p("dD <- muo*S+muo*Z"),
#            br(),
#            p("# rate of change for humans"),
#            p("dSh <- -lam_h*Sh+(1/durinf)*X"),
#            p("dX <- lam_h*Sh-(1/durinf)*X")
           #p("dRh <- (1/durinf)*X-(1/immunity)*Rh")
           downloadButton('downloadODE','Download (works only in browser)')
           #,textOutput(outputId = "lam")
    )
  ),
fluidRow(h2("IBM"),
         column(5,h3("Human population")
                #,plotOutput(outputId="human_pop")
         ),
         column(5,
                h3("Mosquito Population")
                #,plotOutput(outputId = "everything_mosq")
         ),
         column(2,h3("Download csv")
                ,downloadButton('downloadIBM','Download (works only in browser)')
                ))
  
  
)

server <- function(input, output) {
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
  #ibm outputs
  
  ibm_sims <- reactive({
    ####parameters####
    durinf <- 7 #duration of infection
    a <- .1 #human blood feeding rate
    b <- .3 #probability of disease transmission per bite for human
    c <- .7 #probability a mosquito becomes infected after biting an infected human
    muo <- .05 ##10 days survival= 20 half-days survival, therefore 1/20=.05
    moi <- .05
    
    H <- 80 #human population
    X <- 30 #infected humans
    M <- 800 #initial mosquito population
    Z <- 200 #initial infected mosquitos
    timesteps_days <- 28
    timesteps <- timesteps_days*2 #365*2 
    no_sims <- 10 #no. of simulations
    
    lci <- .05 #lowest confidence interval
    hci <- .95  #highest confidence interval
    
    #this also needs to be changed during the subsequent timesteps
    m <- M/H
    z <- Z/M
    x <- X/H #ratio of infectious humans
    
    lam_h <- m*a*b*z #lam_h <- 0.075 #lambda for humans
    lam <- a*c*x #lambda for mosquitos
    
    
    ####synthesizing age and gender####
    sim_age <- sample(age, H, replace=TRUE,prob=age_prob_0to97[,1])
    #hist(sim_age) 
    gender <- rep(NA,length(sim_age))
    
    for(i in 1:length(sim_age)){
      p <- male_prob_0to97[sim_age[i]+1,1] #male_prob is already arranged in ascending age
      gender[i] <- sample(2,1,prob=c(p,1-p))
    }
    
    ####testing the proportions####
    #tmp <- cbind(sim_age,gender)
    #tmp <- as.data.frame(tmp)
    #colnames(tmp) <- c('s_age','s_gender')
    #tmp$s_gender <- as.factor(tmp$s_gender)
    #qplot(s_age, data=tmp,fill=s_gender)
    
    ####synthesizing infected population####
    infected_h <- rep(NA,H)
    for(i in 1:H){
      infected_h[i] <- sample(c(0,1),1, prob=c(.8,.2))
    }
    
    tts <- rep(0, H) #time to become susceptable, 1/dur_inf in normal distribution
    
    random_no <- runif(H) # random uniform no. to decide the prob. of being infected if susceptable
    current <- rep(1, H) #infected in current timestep
    
    df <- cbind(sim_age,gender,infected_h, tts, random_no, current) #variable addition for populated dataframe
    
    ####codebook for df####
    #1. sim_age
    #2. gender
    #3. infected_h
    #4. tts #time to become susceptible again
    #5. random_n #random no. drawn from uniform distribution
    #6. current #a switch to detect if an individual is infected in current timestep or not
    
    
    ###initializing####
    for(i in 1:nrow(df)){
      if(df[i,3] && df[i,6]){ #if infected #at current timestep 
        
        
        df[i,4] <- rnorm(1,mean=1,sd=.2) * durinf #input into tts, time to susceptable
        
        
      }
      
      df[i,4] <- df[i,4]-.5 #tts-.5 per timestep
      df[i,6] <- 0 # resetting 'infected at current timestep'
    }
    
    #first row of the summary table
    #time 0
    X <- sum(df[,3]) #no. of infected humans
    x <- X/H #ratio of infectious humans
    #rate of change of Z from ODE
    lam <- a*c*x
    ###need to check this####
    #Z <- Z+lam*(M-Z) 
    
    #m <- M/H ###no. of mosquitos doesn't change FOR NOW
    z <- Z/M
    lam_h <- m*a*b*z
    
    time0 <- c(0, H-X, X, lam_h, M-Z, Z, lam) #variable addition for simulation table
    
    #######write an initialized file#####
    #write.csv(df, file='0.csv')
    
    
    #######Simulate Summary table function#####
    simulate_summ <- function(){#function for subsequent timesteps
      
      summ_tab <- matrix(NA, nrow=timesteps+1, ncol=7) # summary table for plotting, +1 because it starts from 0 #variable addition for simulation table
      colnames(summ_tab) <- c('timesteps','susceptables','infected', 'lam_h','S','Z','lam') #column names for the summary table
      #variable addition for simulation table
      summ_tab[1,] <- time0 #the first line of the table. the states at time0
      
      #there's an error which one to take as time 0 (or 0.5)
      summ_tab[,1] <- seq(0,timesteps_days,by=(1/2))
      
      for(j in 1:timesteps+1){ #this means 2:(timesteps+1)
        
        for(i in 1:nrow(df)){
          if(df[i,5]<=lam_h){ #if uniform random no. drawn for individual is <= prob of infected
            df[i,3] <- df[i,6] <- 1 #denoting this person is infected on this timestep
          }
          
          if(df[i,3]==1 && df[i,6]==1){ #if infected #at current timestep 
            
            df[i,4] <- rnorm(1,mean=1,sd=.2) * durinf #input into tts, time to become susceptable again
            
          }
          
          df[i,4] <- df[i,4]-.5 #tts-.5 per timestep
          
          if(df[i,4]<=0 && df[i,3]==1){ #currently infected, but durinf is over
            df[i,3] <- 0 #then he becomes suscepitable again
          }
          
          #resetting for the next round
          df[i,5] <- runif(1) #drawing random no. for each individual
          df[i,6] <- 0 # resetting 'infected at current timestep'
        }
        #at the end of big for loop
        #calculate summary variables and lam_h for the next timestep
        X <- sum(df[,3]) #no. of infected humans
        x <- X/H #ratio of infectious humans
        #rate of change of Z from ODE
        lam <- a*c*x #1-(1-(a*c))^x #a*c*x
        S_prev <- (M-Z)
        S <- S_prev+M*mui-muo*S_prev-lam*S_prev
        Z <- Z+lam*S_prev-muo*Z
        
        M <- S+Z #recalculating mosquito population
        #m <- M/H ###no. of mosquitos doesn't change FOR NOW
        z <- Z/M
        lam_h <- m*a*b*z #1-(1-(a*b*m))^z #m*a*b*z
        
        #writing a summary table
        #summ_tab[j,1] <- j
        summ_tab[j,2] <- H-X
        summ_tab[j,3] <- X
        summ_tab[j,4] <-lam_h
        summ_tab[j,5] <- S #############################
        summ_tab[j,6] <- Z #need to have some limitation on Z, infected mosquitos
        summ_tab[j,7] <- lam
        
        ######outputing csv of the simulation on each timestep#######
        #if(j<10 | j>(max(timesteps)-10)){
        #  write.csv(df, file=paste(j,".csv",sep=""))
        #}
      }
      summ_tab
    }
    ####plotting multiple simulation####
    
    
    #creating a list of simulations
    sims <- list()
    for(i in 1:no_sims){
      sims[[i]] <- simulate_summ()
    }
    sims
  })
  ibm_avg_out <- reactive({
    #averaging across the list
    
    #fetching variables
    sims <- ibm_sims()
    no_sims <- input$no_sims
    
    tmp_avg <- rep(NA,no_sims)
    avg_sims <- matrix(NA,nrow(sims[[1]]),ncol(sims[[1]])) #initializing a blank dataset of summary table
    for(i in 1:ncol(sims[[1]])){ #outer loop for the columns
      for(j in 1:nrow(sims[[1]])){ #inner loop for the rows
        for(k in 1:no_sims){#innermost loop for no. of simulations(3rd dimension)
          tmp_avg[k] <- sims[[k]][j,i]
        }
        avg_sims[j,i] <- mean(tmp_avg)
      }
    }
    colnames(avg_sims) <- c('timesteps','susceptables','infected', 'lam_h','S','Z','lam') #column names for the summary table #variable addition
    avg_sims
  })
  
  output$downloadIBM <- downloadHandler(
    filename= function(){paste('ibm_avg_',Sys.Date(),'.csv',sep='')},
    content= function(file){
      write.csv(ibm_avg_out(),file)
    })
}

shinyApp(ui = ui, server = server)

#output$inc2 <- renderPlot({
#  parameters <- parameters_r()
#  out <- ode_out()
#  #pop <- out[,2]+out[,3]

#  #inc <- parameters[3]*parameters[4]*(out[,7]/out[,5])*out[,3]*out[,2]/pop
#  inc2 <- parameters[3]*parameters[4]*(out[,7]/(out[,6]+out[,7]))*out[,3]
#  #plot(inc2, type="l")
#  plot(out[,3])
#})