## useDevel() if need be
library(methylumi)
## biocLite('IlluminaHumanMethylation27k.db') if need be
mset <- normalizeMethyLumiSet(methylumi.bgcorr(methylumIDAT(getBarcodes())))
show(mset)
sessionInfo()
