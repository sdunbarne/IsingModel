randomConfig <- function(N) {
  config <- matrix(2*sample(0:1, N*N, replace=TRUE)-1, N,N)
}

boltzWeight <- function(nrg, kBoltzT) { exp(- nrg/kBoltzT) }

surroundNeighbors  <- function(site) {
    belowSite  <- if (site[1] == N) c(1, site[2]) else c(site[1]+1, site[2])
    leftSite   <- if (site[2] == 1) c(site[1], N) else c(site[1], site[2]-1)
    aboveSite  <- if (site[1] == 1) c(N, site[2]) else c(site[1]-1, site[2])
    rightSite  <- if (site[2] == N) c(site[1], 1) else c(site[1], site[2]+1)
    neighbors  <- c(belowSite, leftSite, aboveSite, rightSite)
}

sweepMetropolis <- function(config) {
  for (x in 1:N) {
      for (y in 1:N) {
          site <- c(x,y)
          configAtSite  <- config[ site[1], site[2] ]
          configPrimeAtSite  <- configAtSite * (-1)

          neigh  <- surroundNeighbors( site ) 
          deltaEnergy <- -2*J * (configPrimeAtSite - configAtSite)*
              (config[ neigh[1], neigh[2] ] +
               config[ neigh[3], neigh[4] ] +
               config[ neigh[5], neigh[6] ]+
               config[ neigh[7], neigh[8] ]
              )
          if ( runif(1) <= min(1, boltzWeight(deltaEnergy, kBoltzT)) )
        {
        config[ site[1], site[2] ]  <- configPrimeAtSite
        }   

    ## IF deltaEnergy < 0, so energy(configPrime) < energy(config)
    ## so boltzWeight(deltaEnergy, kBoltzT) > 1,
    ## THEN  runif(1) is certain to be less than p so certainly
    ## change config to configPrime
    ## ELSE deltaEnergy > 0, so 
    ## energy(configPrime) > energy(config), so 
    ## boltzWeight(deltaEnergy, kBoltzT) < 1,
    ## so change on chance that runif(1) < p
    }
  }
  return(config)
}

J  <- 1
N  <- 100
kBoltzT  <- 2.2
mcmcSteps  <- 50

config  <- randomConfig(N)

for (i in 1:mcmcSteps ) {
    config  <- sweepMetropolis(config)
}

image(1:N, 1:N, config[N:1, ], 
      col=gray.colors(2),
      xlab="", ylab="", xaxt="n", yaxt="n")
## config[N:1, ] so that the plot corresponds to matrix order
## i.e. first row of matrix is at top of plot, last row at bottom


## NAME: isingMetropolis.R
## USAGE: within R, at interactive prompt
##        source("isingMetropolis.R")
## REQUIRED ARGUMENTS: none
## OPTIONS: none
## DESCRIPTION: Find an steady state configuration in Ising model
##              for magnetization, using the Metropolis algorithm
## DIAGNOSTICS: none
## CONFIGURATION AND ENVIRONMENT: base R
## DEPENDENCIES:  base R
## INCOMPATIBILITIES: none known
## PROVENANCE: Steve Dunbar
## BUGS AND LIMITATIONS: 
## FEATURES AND POTENTIAL IMPROVEMENTS:
## AUTHOR:  Steve Dunbar
## VERSION: Version 1.0 as of Fri Oct 14 09:41:06 AM CDT 2022
## KEYWORDS: Ising Model, Metropolis algorithm

