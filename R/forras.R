library(goftest)

Kolmogorov_Smirnov_test<-function(xp, n_sim, alpha, beta, gamma, delta){

  simulated_samples<-rstable(length(xp)*n_sim, alpha=alpha, beta=beta,
                             gamma=gamma, delta=delta)
  simulated_samples<-matrix(simulated_samples, nrow=length(xp),ncol=n_sim)

  ks_test_statistics<-apply(simulated_samples, FUN=function(x){
    fitted_dist<-stableFit(x,type='q', doplot=FALSE)
    fitted_params<-fitted_dist@fit$estimate
    KS_test_stat<-ks.test(x,y='pstable',
                          alpha=fitted_params[1],
                          beta=fitted_params[2],
                          gamma=fitted_params[3],
                          delta=fitted_params[4])$statistic
    return(KS_test_stat)
  }, MARGIN=2)

  ks_stat_original_sample<-ks.test(xp,y='pstable',
                                   alpha=alpha,
                                   beta=beta,
                                   gamma=gamma,
                                   delta=delta)$statistic

  p_value<-sum(ks_test_statistics>ks_stat_original_sample)/n_sim

  return(p_value)
}

Chi_square_test<-function(xp, c, alpha, beta, gamma, delta){
  xx=sort(xp)
  p=rep(0,times=4)

  q=c(0:c)/c
  q[1]=0.00001
  q[c+1]=1-q[1]

  h=qstable(q,alpha,beta,gamma,delta)
  fr_emp=rep(1/c,times=c)
  fr_th=fr_emp
  n=length(xx)
  fr_th=fr_th*n

  for (i in 1:c)
  {fr_emp[i]=sum(xx<h[i+1] &xx>=h[i])}
  chi=sum((fr_emp-fr_th)^2/fr_th)
  p=1-pchisq(chi,df=c-1-4)  #p.value
  return(p)
}
