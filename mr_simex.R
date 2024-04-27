library(simex)
set.seed(100)
mr_simex <- function(harmo_dat = harmo){
  dat_t <- harmo_dat
  #Rename required columns
  BetaXG <- dat_t$beta.exposure
  seBetaXG <- dat_t$se.exposure
  BetaYG <- dat_t$beta.outcome
  seBetaYG <- dat_t$se.outcome
  
  BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
  #BYG <- abs(BetaYG)
  BXG <- abs(BetaXG)         
  
  # MR-Egger regression (unweighted)
  Fit2 <- lm(BetaYG~BetaXG,x=TRUE,y=TRUE) 
  # Simulation extrapolation 
  mod.sim2 <- simex(Fit2, B=1000, measurement.error = seBetaYG, 
                    SIMEXvariable="BetaYG",fitting.method ="quad",asymptotic="FALSE") 
  mod2 <- summary(mod.sim2)
  
  #extract results in beta format
  beta2 <- mod2$coefficients$jackknife[2,1]
  se2 <- mod2$coefficients$jackknife[2,2]
  p2 <- mod2$coefficients$jackknife[2,4]
  
  l= cbind.data.frame("chl","CAD","outcome","exposure","MR Egger SIMEX", length(dat_t$SNP), 
                      beta2, se2, p2)
  colnames(l) <- c("id.exposure", "id.outcome" , 
                   "outcome","exposure","method","nsnp","b","se","pval" )
  
  return(l)
}
