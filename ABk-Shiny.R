library(shiny)
library(psych)
library(shinyalert)

ui <- fluidPage(
  # useShinyalert(), README
  # Title
  titlePanel("(AB)k Power Calculator"),
  numericInput("k", label = "Enter number of phase repetitions (k): ", value = NULL),
  numericInput("n", label = "Enter number of observations per phase (n): ", value = NULL),
  numericInput("m", label = "Enter number of subjects (m): " , value = NULL),
  numericInput(
    "phi",
    label = "Enter autocorrelation (phi): ",
    value = NULL,
    min = -1,
    max = 1
  ),
  numericInput(
    "rho",
    label = "Enter intraclass correlation (rho): ",
    value = NULL,
    min = 0,
    max = 1
  ),
  numericInput("d", label =  "Enter effect size (d): ", value = NULL),
  numericInput(
    "alpha",
    label = "Enter type 1 error rate (alpha): ",
    value = NULL,
    min = 0,
    max = 1
  ),
  numericInput(
    "power",
    label = "Enter desired power: ",
    value = NULL,
    min = 0,
    max = 1
  ),
  actionButton("solve", "Solve for unknown variable"),
  verbatimTextOutput("outputTextBox"),
  tags$style(type = 'text/css', '#outputTextBox         {font-size: 25px; ; background-color: rgba(255,255,255,0.40); color: black; border-style: none;}'), 
)

server <- function(input, output, sessionl) {
  #### POWER CALCULATION FUNCTIONS ####
  h <- function (k, n, m, phi, tau2, s2) {
    h <-
      2 * B(k, n, m, phi, tau2, s2) ^ 2 / Cee(k, n, m, phi, tau2, s2)
    return(h)
  }
  
  
  COV <-
    function (n, phi, s2) {
      SIGMA <- matrix(nrow = n, ncol = n)
      for (i in 1:n) {
        for (j in 1:n) {
          SIGMA[i, j] <- phi ^ abs(i - j)
        }
      }
      SIGMA <- SIGMA * s2 / (1 - phi ^ 2)
      return(SIGMA)
    }
  
  
  VarD <- function (k, n, m, phi, s2) {
    N <- 2 * k * n
    Vee <- COV(N, phi, s2)
    Contrast <- matrix(nrow = 2 * k * n,
                       ncol = 1,
                       data = 1 / (k * n))
    for (i in 1:n) {
      for (j in 1:k) {
        Contrast[i + (j - 1) * 2 * n, 1] <-
          -1 / (k * n)
      }
    }
    Var <- t(Contrast) %*% Vee %*% Contrast
    Var <- Var / m
    return(Var)
  }
  
  
  EA <- function(k, n, m, phi, rho) {
    s2 = 1 - rho
    a <- VarD(k, n, m, phi, s2) / s2
    return(a)
  }
  
  
  ARow <- function (k, n, m, phi, tau2, s2, row) {
    one <- matrix(data = 1,
                  nrow = m,
                  ncol = 1)
    Ablock <- diag(m) - one %*% t(one) / m
    Oblock <- matrix(data = 0,
                     nrow = m,
                     ncol = m)
    BR <- Ablock
    limit <- 2 * k * n
    for (j in 1:limit) {
      if (j == row) {
        BR <- cbind(BR, Ablock)
      } else {
        BR <- cbind(BR, Oblock)
      }
    }
    LL <- m + 1
    UL <- 2 * k * n * m + m
    BR <- BR[, LL:UL]
    return(BR)
  }
  
  
  A <- function (k, n, m, phi, tau2, s2) {
    limit <- 2 * k * n
    A <- ARow(k, n, m, phi, tau2, s2, 1)
    for (i in 2:limit) {
      A <- rbind(A, ARow(k, n, m, phi, tau2, s2, i))
    }
    return(A)
  }
  
  
  CovYRow <- function (k, n, m, phi, tau2, s2, row) {
    Block <- diag(nrow = m)
    BR <- Block * s2 + Block * tau2
    limit <- 2 * k * n
    for (j in 1:limit) {
      BR <-
        cbind(BR, (s2 / (1 - phi ^ 2)) * phi ^ abs(j - row) * Block + Block *
                tau2) 
    }
    LL <- m + 1
    UL <- 2 * k * n * m + m
    BR <- BR[, LL:UL]
    return(BR)
  }
  
  
  CovY <- function (k, n, m, phi, tau2, s2) {
    limit <- 2 * k * n
    R <- CovYRow(k, n, m, phi, tau2, s2, 1)
    for (i in 2:limit) {
      R <- rbind(R, CovYRow(k, n, m, phi, tau2, s2, i))
    }
    return(R)
  }
  
  
  B <- function (k, n, m, phi, tau2, s2) {
    E <-
      tr(A(k, n, m, phi, tau2, s2) %*% CovY(k, n, m, phi, tau2, s2)) /
      (2 * k * n * (m - 1) * (s2 + tau2))
    return(E)
  }
  
  
  EB <- function(k, n, m, phi, rho) {
    s2 <- 1 - input$rho
    tau2 <- rhos <- 1 - rho
    b <- B(k, n, m, phi, tau2, s2) / (s2 + tau2)
    return(b)
  }
  
  
  Cee <- function (k, n, m, phi, tau2, s2) {
    V <-
      2 * tr(
        A(k, n, m, phi, tau2, s2) %*% CovY(k, n, m, phi, tau2, s2) %*% A(k, n, m, phi, tau2, s2) %*% CovY(k, n, m, phi, tau2, s2)
      ) /
      (2 * k * n * (m - 1) * (s2 + tau2)) ^ 2
    return(V)
  }
  
  
  EC <- function (k, n, m, phi, rho) {
    s2 <- 1 - input$rho
    tau2 <- rhos2 <- 1 - rho
    c <- Cee(k, n, m, phi, tau2, s2)
    return(c)
  }
  
  
  PowerNew <- function (k, n, m, phi, rho, d, alpha) {
    tau2 <- rho
    s2 <- 1 - rho
    df <- h(k, n, m, phi, tau2, s2)
    CV <- qf(1 - alpha, df1 = 1, df2 = df)
    a <- VarD(k, n, m, phi, s2) / s2
    b <- B(k, n, m, phi, tau2, s2)
    L <- (b / a) * d ^ 2
    p <- 1 - pf(CV,
                df1 = 1,
                df2 = df,
                ncp = L)
    return(p)
  }
  
  #### END ####
  
  
  clicked <<- 0
  
  output$outputTextBox <- renderPrint({
    ids <- c("k", "n", "m", "phi", "rho", "d", "alpha", "power")
    solvevar <- ""
    numna <- 0
    msg <- ""
    rangecheckpassed <- FALSE
    
    for (id in ids) {
      val <- input[[id]]
      if (is.na(val)) {
        numna <- numna + 1
        
      }
    }
    
    if (numna > 1) {
      print("Please enter exactly 7 input values")
    }
    if (numna == 0) {
      print("Please leave one input blank")
    }
    if (numna == 1) {
      for (id in ids) {
        if (is.na(input[[id]]))
          solvevar <- id
      }
    }
    
    
    numerrors <- 0
    if (numna == 1) {
      for (id in ids) {
        if (id != solvevar) {
          #### CHECK RANGES ####
          if (id == "k") {
            if (input$k <= 1 || input$k %% 1 != 0) {
              print("Please enter a valid number for k (integer value greater than 1)")
              numerrors <- numerrors + 1
            }
          }
          if (id == "n") {
            if (input$n <= 1 || input$n %% 1 != 0) {
              print("Please enter a valid number for n (integer value greater than 1)")
              numerrors <- numerrors + 1
            }
          }
          if (id == "m") {
            if (input$m <= 1 || input$m %% 1 != 0) {
              print("Please enter a valid number for m (integer value greater than 1)")
              numerrors <- numerrors + 1
            }
          }
          if (id == "phi") {
            if (input$phi <= -1 || input$phi >= 1) {
              print(
                "Please enter a valid number for phi (value between -0.99 and 0.99, inclusive)"
              )
              numerrors <- numerrors + 1
            }
          }
          if (id == "rho") {
            if (input$rho <= 0 || input$rho >= 1) {
              print(
                "Please enter a valid number for rho (value between 0 and 0.99, inclusive)"
              )
              numerrors <- numerrors + 1
            }
          }
          if (id == "d") {
            if (input$d < 0) {
              print("Please enter a valid number for d (value greater than or equal to 0)")
              numerrors <- numerrors + 1
            }
          }
          if (id == "alpha") {
            if (input$alpha <= 0 || input$alpha >= 1) {
              print(
                "Please enter a valid number for alpha (value greater than 0 and less than= 1)"
              )
              numerrors <- numerrors + 1
            }
          }
          if (id == "power") {
            if (input$power < 0 || input$power > 1) {
              print("Please enter a valid number for power (value between 0 and 1, inclusive)")
              numerrors <- numerrors + 1
            }
          }
          #### END ####
        }
      }
      if (numerrors == 0) {
        rangecheckpassed <- TRUE
      }
    }
    
    
    showcount <- 0
    if (showcount == 0) {
      show <- TRUE
    } else {
      showcount <- showcount + 1
      print(paste0(showcount))
    }
    if (input$solve > 0) {
      if (rangecheckpassed) {
        msg <- "PASSED CHECKS - BUTTON PRESSED"
        if (solvevar == "k") {
          if (input$n > 15 && input$m > 15) {
            shinyalert(
              "n and m have values of over 15, which may take over 30 seconds to compute k.",
              
            )
          }
          for (k in 2:9999) {
            tempPower = PowerNew(k,
                                 input$n,
                                 input$m,
                                 input$phi,
                                 input$rho,
                                 input$d,
                                 input$alpha)
            if (tempPower >= input$power) {
              print(paste0(solvevar, ": ", k))
              show <- TRUE
              break
            }
          }
        }
        
        if (solvevar == "n") {
          if (input$k > 15 && input$m > 15) {
            shinyalert(
              "k and m have values of over 15, which may take over 30 seconds to compute n.",
              
            )
          }
          for (n in 2:9999) {
            tempPower = PowerNew(input$k,
                                 n,
                                 input$m,
                                 input$phi,
                                 input$rho,
                                 input$d,
                                 input$alpha)
            if (tempPower >= input$power) {
              print(paste0(solvevar, ": ", n))
              show <- TRUE
              break
            }
          }
        }
        if (solvevar == "m") {
          if (input$k > 15 && input$n > 15) {
            shinyalert(
              "k and n have values of over 15, which may take over 30 seconds to compute m.",
              
            )
          }
          for (m in 2:9999) {
            tempPower = PowerNew(input$k,
                                 input$n,
                                 m,
                                 input$phi,
                                 input$rho,
                                 input$d,
                                 input$alpha)
            if (tempPower >= input$power) {
              print(paste0(solvevar, ": ", m))
              show <- TRUE
              break
            }
          }
        }
        if (solvevar == "phi") {
          closestPower = 1
          closestPhi = 0
          set <- seq(-0.99, 0.99, by = 0.01)
          for (phi in set) {
            tempPower = PowerNew(input$k,
                                 input$n,
                                 input$m,
                                 phi,
                                 input$rho,
                                 input$d,
                                 input$alpha)
            if ((tempPower > input$power) &&
                ((tempPower - input$power) < (closestPower - input$power))) {
              closestPhi = phi
              closestPower = tempPower
            }
          }
          print(paste0(solvevar, ": ", closestPhi))
          
        }
        if (solvevar == "rho") {
          closestPower = 1
          closestRho = 0
          set <- seq(0, 0.99, by = 0.01)
          for (rho in set) {
            tempPower = PowerNew(input$k,
                                 input$n,
                                 input$m,
                                 input$phi,
                                 rho,
                                 input$d,
                                 input$alpha)
            if ((tempPower > input$power) &&
                ((tempPower - input$power) < (closestPower - input$power))) {
              closestRho = rho
              closestPower = tempPower
            }
          }
          print(paste0(solvevar, ": ", closestRho))
        }
        if (solvevar == "d") {
          set <- seq(0, 10, by = 0.01)
          for (d in set) {
            tempPower = PowerNew(input$k,
                                 input$n,
                                 input$m,
                                 input$phi,
                                 input$rho,
                                 d,
                                 input$alpha)
            if (tempPower >= input$power) {
              print(paste0(solvevar, ": ", d))
              show <- TRUE
              break
            }
          }
        }
        if (solvevar == "alpha") {
          set <- seq(0.01, 1, by = 0.01)
          for (alpha in set) {
            tempPower = PowerNew(input$k,
                                 input$n,
                                 input$m,
                                 input$phi,
                                 input$rho,
                                 input$d,
                                 alpha)
            if (tempPower >= input$power) {
              print(paste0(solvevar, ": ", alpha))
              show <- TRUE
              break
            }
          }
        }
        if (solvevar == "power") {
          print(
            PowerNew(
              input$k,
              input$n,
              input$m,
              input$phi,
              input$rho,
              input$d,
              input$alpha
            )
          )
          show <- TRUE
        }
        
      }
      
      clicked <<- clicked + 1
    } else {
      msg <- "PASSED CHECKS - BUTTON NOT PRESSED"
      
    }
    
  })
}

shinyApp(ui = ui, server = server)