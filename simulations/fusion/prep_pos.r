flags = commandArgs(T)

prefix = flags[1]
bim.file = flags[2]
N = flags[3]

bim = read.table(bim.file, stringsAsFactors=F)
coords = range(bim[,4])
pos = data.frame(PANEL="sim", WGT="", ID="SIM", CHR=bim[1,1], P0=coords[1]-1, P1=coords[2]+1, N=N, stringsAsFactors=F)

wgt.files = paste0(prefix, ".iter", 1:100, ".wgt.RDat")
out.files = paste0(prefix, ".iter", 1:100, ".pos")

for (i in 1:100) {
  pos$WGT = wgt.files[i]
  write.table(pos, file=out.files[i], row.names=F, quote=F, sep="\t")
}





