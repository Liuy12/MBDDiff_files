library(doParallel)
library(dplyr)
library(MBDDiff)
library(data.table)
library(Biostrings)
library(pryr)
library(doMC)

registerDoMC(40)

### Generate background annotation file for mm9
Promoter_anno <- GetPromoterAnno('mm9', save = T, Dir = './')
t1 <- Sys.time()
bed_100bp_bg <- IdentifyBackground(organism = 'mm9', bed_path = './', binsize = 100, promo_bed = Promoter_anno, fa ='/illumina/iGenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa')
t2 <- Sys.time()
t2 - t1
write.table(bed_100bp_bg, 'mm9Bg.bed', sep = '\t', col.names = F, quote =F)

### Generate background annotation file for mm10
Promoter_anno <- GetPromoterAnno('mm10', save = T, Dir = './')
t1 <- Sys.time()
bed_100bp_bg <- IdentifyBackground(organism = 'mm10', bed_path = './', binsize = 100, promo_bed = Promoter_anno, fa ='/illumina/iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa')
t2 <- Sys.time()
t2 - t1
write.table(bed_100bp_bg, 'mm10Bg.bed', sep = '\t', col.names = F, quote =F)

### Generate background annotation file for hg18
Promoter_anno <- GetPromoterAnno('hg18', save = T, Dir = './')
t1 <- Sys.time()
bed_100bp_bg <- IdentifyBackground(organism = 'hg18', bed_path = './', binsize = 100, promo_bed = Promoter_anno, fa ='/illumina/iGenomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa')
t2 <- Sys.time()
t2 - t1
write.table(bed_100bp_bg, 'hg1Bg.bed', sep = '\t', col.names = F, quote =F)

### Generate background annotation file for hg19
Promoter_anno <- GetPromoterAnno('hg19', save = T, Dir = './')
t1 <- Sys.time()
bed_100bp_bg <- IdentifyBackground(organism = 'hg19', bed_path = './', binsize = 100, promo_bed = Promoter_anno, fa ='/illumina/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa')
t2 <- Sys.time()
t2 - t1
write.table(bed_100bp_bg, 'hg19Bg.bed', sep = '\t', col.names = F, quote =F)

