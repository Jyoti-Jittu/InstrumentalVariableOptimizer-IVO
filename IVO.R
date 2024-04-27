IVO <- function(dat) # InstrumentalVariableOptimizer
{
 # Requirements:
  # beta AND standard error value for the given data should present
  
  if (length(dat)==0) 
     stop(paste("Exposure/Outcome data is empty"))
  
  if(!is.null(dat$beta.exposure))
  {
    dat$t_value <- dat$beta.exposure/dat$se.exposure # t value calculation
    t_avg_exp <- mean(abs(dat$t_value)) 
    exp_t <- dat[dat$t_value>t_avg_exp,] 
    return(exp_t)
  }
  
  if(!is.null(dat$beta.outcome))
  {
    dat$t_value <- dat$beta.outcome/dat$se.outcome # t value calculation
    t_avg_out <- mean(abs(dat$t_value)) 
    out_t <- dat[dat$t_value>t_avg_out,]
    return(out_t)
  }
  else
    stop(paste("beta value not found"))
}