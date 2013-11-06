#plot(ret)
# Pseudo-BLS:
# -- Fold time series at each candidate period
# -- Fit step function to folded series
# -- For each (i1,i2)
# -- Compute mean(folded_y[i1:i2])
# -- Compute var(folded_y[i1:i2])
# -- Minimum variance for which (min_transit_depth < mean < max_transit_depth) is top candidate
# -- Repeat for all periods, i1, i2
# -- Rank candidates by variances
"pbls" <- function(y, P_lo, P_hi, P_inc, a_alpha, b_alpha, a_t_d, b_t_d)
{
  y <- y-mean(y)
  P_candidates <- seq(P_lo,P_hi,by=P_inc)
  ret <- matrix(NA,nrow=length(P_candidates),ncol=5)
  colnames(ret) <- c("P","t_0","t_d","alpha","var")
  ret[,1] <- P_candidates
  for (i in 1:length(P_candidates)){
    P <- P_candidates[i]
    max_i1 <- as.integer(P-a_t_d)
    min_tvar <- Inf
    for (i1 in 1:max_i1){
      min_i2 <- i1+a_t_d
      max_i2 <- min(i1+b_t_d,n)
      for (i2 in min_i2:max_i2){
        
        tbar <- mean(y[i1:i2])
        tvar <- var(y[i1:i2])
        if (tbar<a_alpha || tbar>b_alpha){
        } else if (tvar < min_tvar){
          ret[i,2:5] <- c(i1,(i2-i1),tbar,tvar)
        }
      }
    }
    cat(paste0("Finished P = ",P,"...\n"))
  }
  return(ret)
}
#crude_pbls <- pbls(y=y$data$y,
#	P_lo=exp(a_logP),P_hi=exp(b_logP),P_inc=6000,
#	P_lo=as.integer(y$parameters$P),P_hi=as.integer(y$parameters$P)+1,P_inc=1,
#	a_alpha=a_alpha,b_alpha=b_alpha,
#	a_t_d=a_t_d,b_t_d=b_t_d)