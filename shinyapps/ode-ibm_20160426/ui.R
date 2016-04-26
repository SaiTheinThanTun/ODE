fluidPage(
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
           ,h3("Parameters for IBM"),
           numericInput(inputId = "no_sims",
                         label = "number of simulations",
                         value = 10),
           sliderInput(inputId = "lci",
                       label = "lower bound of CI : lci",
                       value = .05, min = 0, max = 1
           ),
           sliderInput(inputId = "hci",
                       label = "higher bound of CI : hci",
                       value = .95, min = 0, max = 1
           )
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
                  ,plotOutput(outputId="ibm_sims_plot")
           ),
           column(5,
                  h3("Mosquito Population")
                  ,plotOutput(outputId = "ibm_sims_plot_mosq")
           ),
           column(2,h3("Download csv")
                  ,downloadButton('downloadIBM','Download (works only in browser)')
           ))
  
  
)
