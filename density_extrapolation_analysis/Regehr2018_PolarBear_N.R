## 
## 
## 
## Integrated Population Modeling Provides the First Empirical Estimates of Vital Rates
## and Abundance for Polar Bears in the Chukchi Sea
##
## E.V. Regehr, N.J. Hostetter, R.R. Wilson, K.D. Rode, M. St. Martin, and S.J. Converse
##
## Scientific Reports
## October 2018
## 
##
## SUPPLEMENT: DATA AND CODE FOR DENSITY EXTRAPOLATION
## 
## REQUIRED FILES:
##   1. Regehr2018_dextrp_data.txt                   # Auxiliary data (provided)
##   2. NStudyMean                                   # Sample of posterior iterations of mean study area abundance 
##                                                   # excluding AFC0 and C0, obtained from fitted IPM  
##
## DISCLAIMER: The authors provide no guarantee regarding the completeness or functionality of these data and programs, 
## and are not responsible for any consequences of their use. 

#clear all objects in workspace
rm(list=ls())

#load data
dat <- dget("Regehr2018_dextrp_data.txt")
mean.pred <- dat$mean.pred      #mean RSF habitat-quality metric
prop.in <- dat$prop.in          #proportion of study population inside study area at any given time during sampling period
grid.xy <- dat$grid.xy          #state space 25 x 25 km grid
p.AFC0 <- dat$p.AFC0            #approx. number of AFC0 divided by abundance excluding AFC0 and C0; depends on details of fitted IPM
C0.lit <- dat$C0.lit            #C0 litter size; depends on details of fitted IPM
n.bootstrap <- 30000

#function to estimate mode 
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)] }

#initialize results object
a1 <- c( "sampling", "pbsg"  )
results <- matrix( NA, ncol = length(a1), nrow = n.bootstrap )
colnames(results) <- a1; results <- as.data.frame(results)

for(i in 1:n.bootstrap) {
  
  #select values for parameters of interest
  n.study <- NStudyMean[[i]]
  n.sampling <- n.study * sample( x = prop.in, size = 1, replace = TRUE )
  which.rsf <- sample( x = seq( from = 3, to = 1002, by = 1 ), size = 1, replace = TRUE ) 
  tmp.rsf <- as.data.frame( mean.pred[ , c( 1, 2, which.rsf ) ] )
  names(tmp.rsf)[[3]] <- "rsf"
  tmp.p.AFC0 <- p.AFC0[[i]]
  tmp.C0.lit <- C0.lit[[i]] 
  
  #extrapolate abundance using habitat metric
  rsf.sampling <- tmp.rsf[ ( !is.na(grid.xy$study.area) & ( grid.xy$study.area == 1 ) ), ]
  rsf.pbsg <- tmp.rsf[ ( !is.na(grid.xy$PBSG) & ( grid.xy$PBSG == 1 ) ), ]
  N.pbsg <- sum( rsf.pbsg$rsf, na.rm = TRUE ) * ( n.sampling / sum( rsf.sampling$rsf, na.rm = TRUE ) )

  #adjust for AFC0 adn C0
  N.pbsg <- N.pbsg * ( 1 + tmp.p.AFC0 + tmp.p.AFC0 * tmp.C0.lit )

  #populate results
  results[ i, "sampling" ] <- n.sampling
  results[ i, "pbsg" ] <- N.pbsg  }  

#summarize the results
results <- as.matrix(results)
a1 <- c( "mode", "mean", "median", "lci", "uci", "cv" )
N.summary <- matrix( NA, nrow = length(a1), ncol = ncol(results) )
dimnames(N.summary) <- list( a1, colnames(results) )
N.summary[ "mode", ] <- apply( results, 2, estimate_mode )
N.summary[ "mean", ] <- colMeans(results)
N.summary[ "median", ] <- apply( results, 2, median )
a1 <- apply( results, 2, function(x) { xout <- quantile( x, probs = c( 0.025, 0.975 ) ); xout } )
N.summary[ "lci", ] <- a1[ "2.5%", ]
N.summary[ "uci", ] <- a1[ "97.5%", ]
N.summary[ "cv", ] <- round( apply( results, 2, sd ) / colMeans(results), digits = 2 )
N.summary


##
##
##
## END 
##
##
##
