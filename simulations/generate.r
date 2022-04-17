flags = commandArgs(T)

library(snpStats, quietly=T)
library(mvtnorm)


data.id = as.numeric(flags[1]) # from 1 to 10, corresponding to SNP blocks taken from the middle of chromosomes 1 through 10
iter.block = as.numeric(flags[2]) # unique ID to generate additional iterations for the same settings
N.pheno = as.numeric(flags[3]) # sample size for phenotype sample, must be at least as big as largest N.expr
K.snp = as.numeric(flags[4]) # number of SNPs in the block

store.raw = (iter.block == 1) #store expression and phenotype vectors for FUSION/CoMM
full.range = (iter.block <= 10) #reduced parameter range for iterations beyond first 10,000

no.iter = 100  # number of iterations per condition
prune.thresh = 0.99  # pruning threshold during PC projection

numstr = function(x) {format(x, scientific=F, trim=T)}
make.label = function(type, N, K, h2) {paste0(type, "_N", numstr(N), "_K", K, "_h", numstr(h2))}
expr.label = function(Np, Ne, K, h2) {make.label(paste0("expr_p", numstr(Np)), Ne, K, h2)}
pheno.label = function(N, K, h2) {make.label("pheno", N, K, h2)}

if (full.range) {
	N.expr = c(250, 1000)  # sample sizes for gene expression sample
} else {
	N.expr = 1000
}
N.used = sort(unique(c(N.expr, N.pheno)))
N.sim = max(N.used)  # use largest N as simulation base
n.index = function(N) {which(N.used == N)}


# h2 values are defined in percentages, then scaled up to whole numbers for easier file naming/processing later
if (full.range) {
	h2.expr = c(1, 5, 10)
	h2.pheno = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)
} else {
	h2.expr = 5
	h2.pheno = c(0.001, 0.01, 0.1, 1)
}

perc.scale = 1000
h2.scale = perc.scale * 100 # constant to divide by to get explained variance fractions
h2.expr = perc.scale * h2.expr
h2.pheno = perc.scale * h2.pheno


prefix.base = paste0("block", data.id)

root = "..."  # root directory
sim.root = paste0(root, "/simulations/generate")  # simulation directory
data.root = paste0(root, "/data/ukb") # data directory
block.root = paste0(data.root, "/blocks/", prefix.base) # genotype data directory
out.root = paste0(sim.root, "/stats")


# set up output storage and subindices
stats = list(); vars = list()
for (h2 in h2.pheno) {
	stats[[pheno.label(N.pheno, K.snp, h2)]] = matrix(NA, nrow=no.iter, ncol=K.snp)
	if (store.raw) vars[[pheno.label(N.pheno, K.snp, h2)]] = matrix(NA, nrow=no.iter, ncol=N.pheno)
}
for (Ne in N.expr) {
	for (h2 in h2.expr) {
		stats[[expr.label(N.pheno, Ne, K.snp, h2)]] = matrix(NA, nrow=no.iter, ncol=K.snp)
		if (store.raw) vars[[expr.label(N.pheno, Ne, K.snp, h2)]] = matrix(NA, nrow=no.iter, ncol=Ne)
	}
}

#set and check output files
out.files = list(stats=list(), var=list())
for (cond in names(stats)) {
  out.dir = paste0(out.root, "/", cond)
  if (!file.exists(out.dir)) dir.create(out.dir, recursive=T, showWarnings=F)
	out.prefix = paste0(out.dir, "/", prefix.base, "-", iter.block, "_", cond)
	out.files$stats[[cond]] = paste0(out.prefix, ".stats")
	if (store.raw) out.files$var[[cond]] = paste0(out.prefix, ".var")
}
if (all(file.exists(unlist(out.files)))) quit("no") #if all output files already exist, abort (for easy resumption of partially completed simulation runs)


# read genotype data, make index of individuals for each block
# - samples are nested, so largest sample used will contain all smaller samples as subsets
# - family and indivual ID are identical in input PLINK files for all subjects
# - data has been processed in advance to remove MAC=0 SNPs and mean-impute missing genotype values
geno = read.plink(paste0(block.root, "/", prefix.base, "_N", numstr(N.sim), "_K", numstr(K.snp)))

indiv.index = list()
id.files = paste0(data.root, "/subj_N", numstr(N.used), ".ids")
if (!all(file.exists(id.files))) stop("missing subject ID files")
for (i in seq_along(N.used)) {
  matched = which(geno$fam$member %in% read.table(id.files[i], stringsAsFactors=F)[,1])
  if (length(matched) == N.used[i]) {
    indiv.index[[i]] = matched
  } else {
    stop(paste0("sample for N = ", N.used[i], " has incorrect size"))
  }
}


# create template file if not yet exists
tpl.file = paste0(out.root, "/tpl/", prefix.base, ".K", K.snp, ".stats.tpl")
if (!file.exists(tpl.file)) {
  sum.stats = data.frame(CHR=geno$map$chromosome, BP=geno$map$position, SNP=geno$map$snp.name, A1=geno$map$allele.1, A2=geno$map$allele.2)
  write.table(sum.stats, file=tpl.file, row.names=F, quote=F, sep="\t")
}

# process genotype data for simulation
X = scale(as(geno$genotypes, "numeric")); rm(geno)
no.snps = ncol(X)
if (no.snps != K.snp) stop("incorrect number of SNPs in data")

# create subsamples of genotypes for analysis later
X.sub = list()
for (i in seq_along(N.used)) {
  if (i < length(N.used)) {
    X.sub[[i]] = scale(X[indiv.index[[i]],])
  } else {
    X.sub[[i]] = X
  }
  tpl.file = paste0(out.root, "/tpl/", prefix.base, ".N", numstr(N.used[i]), ".var.tpl")
	if (!file.exists(tpl.file)) {
		var.tpl = data.frame(FID=rownames(X.sub[[i]]), IID=rownames(X.sub[[i]]), stringsAsFactors=F)
		write.table(var.tpl, file=tpl.file, row.names=F, quote=F, sep="\t")
	}
}

# compute decomposition, and project genotypes onto principal component matrix (standardized) W
decomp = eigen(cor(X))
lambda = decomp$values
comp.use =  1:min(which(cumsum(lambda / sum(lambda)) >= prune.thresh))
no.comp = length(comp.use)

R = decomp$vectors[,comp.use] %*% diag(1/sqrt(lambda[comp.use]))
W = X %*% R


# generates outcome, assuming standardized G
make.Y = function(G, h2) {
  sigma = (1 - h2) / h2
  err = scale(rnorm(length(G))) * as.numeric(sqrt(sigma))
  return(scale(G+err))
}

# computes sumstats (simple regression) for each SNP, assuming standardized X and Y; outputs t-statistic
compute.stats = function(X, Y) {
  N = length(Y)
  beta = t(X) %*% Y / (N-1)
  svar = (1-beta^2) / (N-2) # since residual variance = (1-beta^2) * (N-1) / (N-2), and sampling variance = residual variance / N-1
  stat = beta / sqrt(svar)
  return(stat)
}


# main simulation loop
for (i in 1:no.iter) {
  # create true effect sizes for principal components under null of complete independence
  delta.e = scale(rnorm(no.comp))
  delta.y = scale(lm(rnorm(no.comp)~delta.e)$residuals)

  # create corresponding genetic effect vectors
  G.expr = scale(W %*% delta.e)
  G.pheno = scale(W %*% delta.y)

  # create full variable Y and corresponding association statistics for all h2 for phenotype
  for (h2.curr in h2.pheno) {
    Y = make.Y(G.pheno, h2.curr/h2.scale)
    stats[[pheno.label(N.pheno, K.snp, h2.curr)]][i,] = compute.stats(X, Y)
		if (store.raw) vars[[pheno.label(N.pheno, K.snp, h2.curr)]][i,] = Y
  }

  # create full variable Y and corresponding association statistics for all N and h2 for expression
  for (j.expr in seq_along(N.expr)) {
    Ne = N.expr[j.expr]
    index.curr = indiv.index[[n.index(Ne)]]
    X.curr = X.sub[[n.index(Ne)]]

    G.curr = scale(G.expr[index.curr])
    for (h2.curr in h2.expr) {
      Y = make.Y(G.curr, h2.curr/h2.scale)
      stats[[expr.label(N.pheno, Ne, K.snp, h2.curr)]][i,] = compute.stats(X.curr, Y)
      if (store.raw) vars[[expr.label(N.pheno, Ne, K.snp, h2.curr)]][i,] = Y
    }
  }
}


#save output (skip existing)
for (cond in names(stats)) {
	if (!file.exists(out.files$stats[[cond]])) write.table(stats[[cond]], file=out.files$stats[[cond]], row.names=F, col.names=F, quote=F, sep="\t")
	if (store.raw && !file.exists(out.files$var[[cond]])) write.table(vars[[cond]], file=out.files$var[[cond]], row.names=F, col.names=F, quote=F, sep="\t")
}





