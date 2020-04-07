install.packages("deSolve")
library(deSolve)


## Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 0.25, gamma = 0.1)
## Time frame
times      <- seq(0, 500, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
head(out, 10)
tail(out)


## Plot
matplot(x = times, y = out, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:4)

## Add legend
legend(40, 0.7, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")


## SEIR model 
seir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -(1-delta * I) * beta * S * I
    dI <-  (1-delta * I) * beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}



## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter; delta: endogenous parameter
parameters <- c(beta = 0.25, gamma = 0.1, delta=8)
## Time frame
times      <- seq(0, 500, by = 1)


## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = seir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
head(out, 10)
tail(out, 10)
max(out$I)
which.max(out$I)

## Plot
matplot(x = times, y = out, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
        lwd = 2, lty = 1, bty = "l", col = 2:4)

## Add legend
legend(300, 0.7, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")



## SEIR model (two groups)
seir2 <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS1 <- -(1-alpha1 * n1 * I) * beta * S1 * I
    dS2 <- -(1-alpha2 * n2 * I) * beta * S2 * I
    dI <-  ((1-alpha1 * n1 * I) * beta * S1 + (1-alpha2 * n2 * I) * beta * S2) * I - gamma * I
    
    dD1 <- n1*(1-alpha1 * n1 * I) * beta * S1 * I
    dD2 <- n2*(1-alpha2 * n2 * I) * beta * S2 * I
    dR <-                 gamma * I
    
    return(list(c(dS1, dS2, dI, dR, dD1, dD2)))
  })
}




## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S1 = (1-1e-6)*0.5, S2 = (1-1e-6)*0.5, I = 1e-6, R = 0.0, D1=0.0, D2=0.0)
## beta: infection parameter; gamma: recovery parameter; alpha: endogenous parameter; eta: death rate
parameters <- c(beta = 0.25, gamma = 0.1, alpha1 = 1, alpha2=50, n1=0.001, n2=0.05)
## Time frame
times      <- seq(0, 500, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = seir2, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
out$D1 <- NULL
out$D2 <- NULL
## Show data

tail(out)
max(out$I)
which.max(out$I)



## Plot
matplot(x = times, y = out, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
        lwd = 2, lty = 1, bty = "l", col = 2:5)

## Add legend
legend(300, 0.7, c("Susceptible1", "Susceptible2", "Infected", "Recovered"), pch = 1, col = 2:5, bty = "n")



