library(LAVA)

flags = commandArgs(T)

pheno = flags[1]
chromosome = flags[2]

tissue = "whole_blood"


root = "..."  # root directory
rev.root = paste0(root, "/revision")

eqtl.dir = paste0(rev.root, "/data/gtex_v8")  # directory with eQTL input files formatted for LAVA
prog.dir = paste0(rev.root, "/prog/lava") # directory with additional LAVA scripts
pheno.dir = paste0(root, "/data/pheno") # directory with GWAS summary statistics

out.dir = paste0(rev.root, "/output/bivar/", pheno)
if (!file.exists(out.dir)) dir.create(out.dir, recursive=T, showWarnings=F)

g1000.prefix = paste0(rev.root, "/data/ref/g1000_maf005") # 1,000 Genomes reference data
pheno.file = paste0(pheno.dir, "/phenotypes.info")  # phenotype info file
out.file = paste0(out.dir, "/", pheno, "-", tissue, "-chr", chromosome, ".res")
if (file.exists(out.file)) stop("output file already exists")

source(paste0(prog.dir, "/lava_twas.r"))
source(paste0(prog.dir, "/lava_approx.r"))

# load in phenotype GWAS summary stats
gwas.input = process.input(input.info.file=pheno.file, ref.prefix=g1000.prefix, sample.overlap.file=NULL, input.dir=pheno.dir, phenos=pheno)

# append eQTL sumstats and annotation
process.eqtl.input(gwas.input, eqtl.dir, tissue=tissue, chromosomes=chromosome)


tmp.file = Sys.glob(paste0(out.file, ".tmp*"))

out = data.frame(
	gene=gwas.input$current.genes, chromosome=chromosome,
	no.snps=NA, K=NA, N.pheno=NA, N.tissue=NA,
	failed.pheno=T, failed.tissue=T,
	omega.pheno=NA, omega.tissue=NA, omega.cov=NA,
	sigma.pheno=NA, sigma.tissue=NA,
	univ.pheno.h2=NA, univ.tissue.h2=NA,
	univ.pheno.p=NA, univ.tissue.p=NA,
	rho=NA, rho.lower=NA, rho.upper=NA,
	r2=NA, r2.lower=NA, r2.upper=NA,
	bivar.p=NA, bivar.approx.p=NA,
	twas.p=NA, twas.r2=NA,
	stringsAsFactors=F
)


bivar.names = c("rho", "rho.lower", "rho.upper", "r2", "r2.lower", "r2.upper")
i.range = which(is.na(out$no.snps))
for (i in i.range) {
	locus = process.eqtl.locus(i, gwas.input, drop.failed=F)  # create locus object with processed eQTL and phenotype data
	if (length(locus) > 1) {
		#NB: for extracting/storing values, order in locus object is always phenotype, eqtl
		out$no.snps[i] = locus$n.snps
		out$K[i] = locus$K
		out[i,c("N.pheno", "N.tissue")] = round(locus$N)
		out[i,c("omega.pheno", "omega.tissue")] = diag(locus$omega)
		out$omega.cov[i] = locus$omega[1,2]
		out[i,c("sigma.pheno", "sigma.tissue")] = diag(locus$sigma)
		out[i,c("failed.pheno", "failed.tissue")] = locus$failed

		univ = run.univ(locus)  # run univariate tests; filtering on univariate p-value is done post-hoc, when processing results
		out[i,c("univ.pheno.h2", "univ.tissue.h2")] = univ$h2.obs
		out[i,c("univ.pheno.p", "univ.tissue.p")] = univ$p

		twas = lava.twas(locus) # run LAVA-TWAS
		out$twas.p[i] = twas[1]
		out$twas.r2[i] = twas[2]

		if (!any(locus$failed)) { # run LAVA-rG, original version and normal approximation
			bivar = try(run.bivar(locus), silent=T)
			if (class(bivar) != "try-error") {
				out[i,c(bivar.names,"bivar.p")] = bivar[,c(bivar.names,"p")]
				out$bivar.approx.p[i] = bivar.approx(locus)
			}
		}
	}
}

write.table(out, file=out.file, row.names=F, quote=F, sep="\t")



