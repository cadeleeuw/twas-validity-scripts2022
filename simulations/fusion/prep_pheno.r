flags = commandArgs(T)

prefix = flags[1]
in.dir = flags[2]
tpl.file = flags[3]

in.file = paste0(in.dir, "/", prefix, ".var")

expr = read.table(tpl.file, header=T, stringsAsFactors=F)
vars = data.matrix(read.table(in.file))

for (i in 1:nrow(vars)) {
	expr$GEX = vars[i,]
	write.table(expr, file=paste0(in.dir, "/", prefix, ".iter", i, ".var"), row.names=F, col.names=F, quote=F, sep="\t")
}
