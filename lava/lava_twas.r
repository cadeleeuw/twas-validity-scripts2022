# run LAVA-TWAS analysis on locus object
# assumes locus object contains only eQTL and single outcome phenotype, in that order
# returns p-value and variance explained in Y by E_hat
lava.twas = function(locus) {
	if (length(locus$phenos) != 2 || !(any(locus$phenos == "eqtl"))) return(rep(NA, 2))

	index = which(locus$phenos == "eqtl")
	index = c(index, 3 - index) # first eqtl index, then pheno index

	return(twas.core(locus$omega[index,index], locus$sigma[index,index], locus$K))
}

# core function based on sufficient statistics
# assumes order of: eqtl, then phenotype
twas.core = function(omega, sigma, K) {
	obs = omega[1,2]

  sigma = diag(sigma); omega = diag(omega) #reduce to variances
	samp.var = (sigma[1]*sigma[2] + omega[1]*sigma[2]) / K
  if (samp.var > 0) {
    stat = -abs(obs/sqrt(samp.var))
    p = pnorm(stat)*2
    r2 = K * obs^2 / (omega[1] + sigma[1])
    return(c(p, r2))
  } else {
    return(rep(NA, 2))
  }
}
