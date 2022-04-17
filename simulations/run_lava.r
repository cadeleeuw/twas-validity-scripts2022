library(LAVA)

flags = commandArgs(T)

data.id = flags[1]
iter.block = flags[2]
conditions = flags[-(1:2)] # any number of condition identifiers, runs all expr-pheno combinations (must all have the same number of SNPs)

numstr = function(x) {format(x, scientific=F, trim=T)}

N.ref = 1000

no.snps = as.numeric(gsub(".*K([0-9]*).*", "\\1", conditions))
if (!all(no.snps == no.snps[1])) stop("inconsistent K values")
no.snps = no.snps[1]

prefix.base = paste0("block", data.id)
prefix.iter = paste0(prefix.base, "-", iter.block)

root = "..."  # root directory
lava.dir = paste0(root, "/prog/lava") # directory with additional LAVA script files
data.dir = paste0(root, "/data/ukb/blocks/", prefix.base) # directory with reference genotype data
stat.root = paste0(root, "/simulations/generate/stats")
out.root = paste0(root, "/simulations/analyse/lava/output")

ref.data = paste0(data.dir, "/", prefix.base, "_N", numstr(N.ref), "_K", no.snps)
tpl.file = paste0(stat.root, "/tpl/", prefix.base, ".K", no.snps, ".stats.tpl") # template file containing the SNP IDs and alleles corresponding to the simulated summary statistics (stored separately)
if (!file.exists(tpl.file)) stop("missing .tpl file")


# parse conditions to analyse (expression condition specifiers are prefixed by 'expr')
# - running rG only for subset of conditions to reduce computational burden, specified via conditions$do.rg
# - iter.block == 1 corresponds to first 1,000 iterations, all subsequent iterations (up to 100,000) are only analyses for LAVA-TWAS
conditions = data.frame(condition=conditions, stringsAsFactors=F)
conditions$is.expr = grepl("expr", conditions$condition)
conditions$N.stats = as.numeric(gsub(".*N([0-9]*)_.*", "\\1", conditions$condition))
conditions$N.pheno = suppressWarnings(as.numeric(gsub(".*p([0-9]*)_.*", "\\1", conditions$condition)))
conditions$h2 = as.numeric(gsub(".*h([0-9]*)", "\\1", conditions$condition))
conditions$do.rg = (iter.block == 1) & (no.snps == 1000) & (conditions$N.stats != 250)
conditions$file = paste0(stat.root, "/", conditions$condition, "/", prefix.iter, "_", conditions$condition, ".stats")
conditions = conditions[file.exists(conditions$file),]
no.cond = nrow(conditions)


condition.pairs = NULL
for (i.e in which(conditions$is.expr)) {
	condition.pairs = rbind(
		condition.pairs,
		cbind(i.e, which(!conditions$is.expr & conditions$N.stats == conditions$N.pheno[i.e]))
	)
}

condition.pairs = data.frame(
	index.expr=condition.pairs[,1],
	index.pheno=condition.pairs[,2],
	stringsAsFactors=F
)
condition.pairs$do.rg = conditions$do.rg[condition.pairs$index.expr] & conditions$do.rg[condition.pairs$index.pheno]
condition.pairs$label = paste0("K", numstr(no.snps),
	"_N",	numstr(conditions$N.stats[condition.pairs$index.expr]),
			"-", numstr(conditions$N.stats[condition.pairs$index.pheno]),
	"_h", numstr(conditions$h2[condition.pairs$index.expr]),
			"-", numstr(conditions$h2[condition.pairs$index.pheno])
)
no.pairs = nrow(condition.pairs)


# load additional LAVA scripts
source(paste0(lava.dir, "/lava_twas.r"))
source(paste0(lava.dir, "/lava_approx.r"))

# read reference data for alignment and checking
bim = read.bim.custom(ref.data, as.env=T)
if (no.snps != nrow(bim$bim)) stop("wrong number of SNPs")

# read summary statistics
stats = list()
for (i in 1:nrow(conditions)) stats[[i]] = data.matrix(read.table(conditions$file[i]))
if (!all(sapply(stats, ncol) == no.snps)) stop("sumstat files contain incorrect number of summary statistics")
if (length(unique(sapply(stats, nrow))) > 1) stop("inconsistent number of iterations in sumstat files")
no.iter = nrow(stats[[1]])

stat.info = read.table(tpl.file, header=T, stringsAsFactors=F); stat.info$SNP = tolower(stat.info$SNP)
if (nrow(stat.info) != no.snps || !all(stat.info$SNP == bim$bim$snp.name)) stop("mismatch between reference data and sumstat files")


# check and correct alignment
mismatch = !apply(bim$bim[,c("allele.1", "allele.2")] == stat.info[,c("A1", "A2")], 1, all)
if (any(mismatch)) {
  if (!all(bim$bim[mismatch,c("allele.1", "allele.2")] == stat.info[mismatch,c("A2", "A1")])) stop("encountered inconsistent alleles")
  for (i in 1:2) stats[[i]][,mismatch] = -1*stats[[i]][,mismatch]
}


# manually construct input environment and locus definition object
input = new.env(parent=globalenv())
input$info = data.frame(phenotypes=conditions$condition, binary=F, stringsAsFactors=F)
input$ref = bim
input$ref.prefix = ref.data
input$analysis.snps = stat.info$SNP
input$sum.stats = list()
for (i in 1:no.cond) input$sum.stats[[input$info$phenotypes[i]]] = data.frame(SNP=stat.info$SNP, N=conditions$N.stats[i], STAT=NA, stringsAsFactors=F)

locus = data.frame(LOC="locus", SNPS=paste(stat.info$SNP, collapse=";"), stringsAsFactors=F)


out.tpl.twas = data.frame(univ.p.expr=rep(NA, no.iter), univ.p.pheno=NA, twas.p=NA, locus.failed=F)
out.tpl.full = data.frame(univ.p.expr=rep(NA, no.iter), univ.p.pheno=NA, bivar.approx.p=NA, bivar.full.p=NA, twas.p=NA, locus.failed=F)

out = list()
for (i in 1:no.pairs) {
	if (condition.pairs$do.rg[i]) {
  	out[[i]] = out.tpl.full
	} else {
  	out[[i]] = out.tpl.twas
	}
}

for (i in 1:no.iter) {
  for (j in 1:no.cond) input$sum.stats[[j]]$STAT = stats[[j]][i,]
  locus.data = process.locus(locus, input, drop.failed=F) # loads each condition as a separate phenotype into the locus object

  univ = run.univ(locus.data) # run univariate test
	for (j in 1:no.pairs) {
		ind = c(expr=condition.pairs$index.expr[j], pheno=condition.pairs$index.pheno[j])

    out[[j]][i,c("univ.p.expr", "univ.p.pheno")] = univ$p[ind]

		# run LAVA-TWAS
		out[[j]]$twas.p[i] = twas.core(locus.data$omega[ind,ind], locus.data$sigma[ind,ind], locus.data$K)[1]

		# run LAVA-rG
		if(condition.pairs$do.rg[j] && !any(locus.data$failed[ind])) {
			bivar = try(run.bivar(locus.data, phenos=locus.data$phenos[ind], adap.thresh=NULL, CIs=F), silent=T)
			if (class(bivar) != "try-error") {
				out[[j]]$bivar.full.p[i] = bivar$p
				out[[j]]$bivar.approx.p[i] = approx.core(locus.data$omega[ind,ind], locus.data$sigma[ind,ind], locus.data$K)[1]
			}
		} else {
	    out[[j]]$locus.failed[i] = T
		}
	}
}

for (i in 1:no.pairs) {
	out.dir = paste0(out.root, "/", condition.pairs$label[i])
	if (!file.exists(out.dir)) dir.create(out.dir)
	out.file = paste0(out.dir, "/", prefix.iter, "_", condition.pairs$label[i], ".res")
	write.table(out[[i]], file=out.file, row.names=F, quote=F, sep="\t")
}



