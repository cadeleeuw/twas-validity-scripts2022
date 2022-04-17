flags = commandArgs(T)

library(snpStats, quietly=T)
library(CoMM)
library(data.table)

cond = flags[1:2] # expression and phenotype conditions to analyse (in that order)
data.dir = flags[3]
input.dir = flags[4]
output.pref = flags[5]

# load in genotype data, and simulated E and Y values (block of 100 iterations)
N = gsub(".*N([0-9]*).*", "\\1", cond)	# extract sample size from condition specifiers
X = list(); Y = list(); covar = list()
for (i in 1:2) {
	curr.ids = read.table(paste0(input.dir, "/N", N[i], ".var.tpl"), header=T, stringsAsFactors=F) # template file that contains the individual IDs for the raw E and Y values (stored separately)
	curr.X = read.plink(paste0(data.dir, "/N", N[i]))
	if (!all(curr.ids$IID == curr.X$fam$member)) stop("ID mismatch")

	X[[i]] = as(curr.X$genotypes, "numeric")
	Y[[i]] = data.matrix(data.table::fread(paste0(input.dir, "/", cond[i], ".var"), data.table=F))
	covar[[i]] = matrix(1, nrow(X[[i]]), 1) # dummy covariate matrix
}


run.comm = function(X.expr, X.pheno, Y.expr, Y.pheno, cov.expr, cov.pheno) {
  h.alt = CoMM_covar_pxem(Y.expr, Y.pheno, X.expr, X.pheno, cov.expr, cov.pheno, constr = 0)
  h.null = CoMM_covar_pxem(Y.expr, Y.pheno, X.expr, X.pheno, cov.expr, cov.pheno, constr = 1)

  stat = 2*(max(h.alt$loglik, na.rm=T) - max(h.null$loglik, na.rm=T))
  pval = pchisq(stat, 1, lower.tail=F)
  return(pval)
}


pval = c()
for (i in 1:100) {
	pval[i] = run.comm(X[[1]], X[[2]], Y[[1]][index,], Y[[2]][index,], covar[[1]], covar[[2]])
}

write.table(pval, file=paste0(output.pref, ".pval"), row.names=F, col.names=F, quote=F, sep="\t")

