flags = commandArgs(T)

prefix = flags[1]
tpl.file = flags[2]

stats = data.matrix(read.table(paste0(prefix, ".stats"), stringsAsFactors=F))
tpl = read.table(tpl.file, stringsAsFactors=F, header=T)[,c("SNP", "A1", "A2")]
tpl$Z = NA

for (i in 1:100) {
  tpl$Z = stats[i,]
  write.table(tpl, file=paste0(prefix, ".iter", i, ".stats"), row.names=F, quote=F, sep="\t")
}


