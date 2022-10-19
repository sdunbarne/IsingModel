## WARNING: USE ONLY FOR N = 2, 3, 4

allBinarySeq <- function(N) { # USE ONLY FOR N <= 4!!
# list of N copies of -1, 1, simplify=FALSE returns lists, not vectors
# then all possible combinations of \pm 1 over N places
# resulting in data.frame of 2^N observations of N variables
# transpose to make N by 2^N matrix
# recast back to data.frame    
    as.data.frame(
        t(
            expand.grid(
                replicate(N, c(-1,1), simplify=FALSE)
            )
        )
    )
}

allConfigs <- function(N) { # USE ONLY FOR N <= 4
# returns LIST of all square configs
    X <- allBinarySeq(N*N)
    # data.frame of all 2^(N^2) observations over N^2 variables

    data = lapply(X, matrix, nrow=N, ncol=N)
}

downConfig <- function(N) { matrix(-1,N,N) } # 1 of 2 lowest energy configs
upConfig <- function(N) { matrix(1,N,N) }    # other lowest energy configs

energy  <-  function(config, J=1) {  #pass in a MATRIX
    ## This uses -J, same as both Richey and Schlusser
    ## J has units of energy
## Returns a potential energy, negative, with lowest energy -4N^2    
    N  <- NROW(config)

    below  <- config[ c(2:N, 1), ]
    left  <- config[ , c(N, 1:(N-1)) ]
    above  <- config[c(N, 1:(N-1)), ]
    right  <- config[ , c(2:N, 1) ]

    nearNeighProd <- config * (below + left + above + right)
    ## element-wise multiplication, note factorization of sum
    
    nrg  <- -J *sum(nearNeighProd)
}

allEnergies  <- function(configurations) { #list of configurations
    nrgs  <-  lapply(configurations, energy, J)
}

boltzWeight <- function(nrg, kTemp) { exp(- nrg/kTemp) }
## This has the negative sign, same as both Richey and Schlusser
## But this does not use J, already in energy, consistent with Richey (but not Schlusser)
## Uses kTemp, Boltzmann constant times temperature, with units of energy
## so the energy usits cancel.

netMag  <-  function(config) { abs( sum( config ) ) }

boltzDist  <- function(configurations, kTemperature) {
    bzs  <- unlist(
        lapply(allEnergies(configurations), boltzWeight, kTemperature)
    )
}

partitionFunc  <- function(configurations, kTemperature){ #list of energies
        Z  <- sum( boltzDist(configurations, kTemperature) )
}

netSpin <- function(configurations, N) {   #list of configurations
    lm <- lapply(configurations, matrix, N, N)
    sm  <- abs( unlist( lapply(lm, sum) ) ) #absolute value sum of entries in each matrix
    # unlist converts a list to a vector
}

magnetization  <- function(configurations, kTemperature) {
    Z <- partitionFunc(configurations, kTemperature)
    mag <- (1/Z) * sum(  netSpin(configurations, N) * boltzDist(configurations, kTemperature) )
}

J <- 1
N  <- 3

configs  <- allConfigs(N)

nKTemperatures  <-  11
Temperatures  <- seq(2,2.5, length=nKTemperatures)
# surrounding the critical value kT_c/J = 2.269 when J=1

m  <- rep(0, nKTemperatures)

for (i in 1:nKTemperatures) {
  m[i] = magnetization(configs, Temperatures[i])  
  cat("T = ", Temperatures[i], "mag = ",  m[i], "\n" )
}

plot(Temperatures, m)

## WARNING: USE ONLY FOR N = 2, 3, 4
## NAME: isingDist.R
## USAGE: within R, at interactive prompt
##        source("isingDist.R")
## REQUIRED ARGUMENTS: none
## OPTIONS: none
## DESCRIPTION: calculate magnetization directly from the Boltzmann
##              distribution of ALL configurations.
##              WARNING: USE ONLY FOR N = 2, 3, 4
## DIAGNOSTICS: none
## CONFIGURATION AND ENVIRONMENT: base R
## DEPENDENCIES:  base R
## INCOMPATIBILITIES: none known
## PROVENANCE: Steve Dunbar,
## BUGS AND LIMITATIONS: WARNING: USE ONLY FOR N = 2, 3, 4
## FEATURES AND POTENTIAL IMPROVEMENTS: 
## AUTHOR:  Steve Dunbar
## VERSION: Version 1.0 as of Tue 05 Jan 2021 07:35:03 AM CST
## KEYWORDS: Ising Model, statistical mechanic, magnetization,
##           Boltzmann distribution 

