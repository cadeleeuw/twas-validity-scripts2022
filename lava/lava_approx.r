# run bivariate LAVA-rG analysis on locus object, using normal approximation
# assumes locus object contains only eQTL and single outcome phenotype (order doesn't matter)
# returns p-value
bivar.approx = function(locus, param.lim=1.25) {
	if (!any(locus$failed) && is.finite(locus$omega.cor[1,2]) && abs(locus$omega.cor[1,2]) <= param.lim) {
		return(approx.core(locus$omega, locus$sigma, locus$K))
	} else {
		return(NA)
	}
}


# core function based on sufficient statistics
approx.core = function(omega, sigma, K) {
  obs = omega[1,2] #observed covariance

  sigma = diag(sigma); omega = diag(omega) #reduce to variances
	if (all(c(sigma, omega) > 0)) {
		samp.var = (sigma[1]*sigma[2] + omega[1]*sigma[2] + omega[2]*sigma[1]) / K
	  pval = pnorm(abs(obs/sqrt(samp.var)), lower.tail=F)*2
  	return(pval)
	} else {
		return(NA)
	}
}

