library(doParallel)
library(dplyr)
library(MBDDiff)
library(data.table)
library(Biostrings)
library(pryr)
library(doMC)

registerDoMC(40)
Promoter_anno <- GetPromoterAnno('mm9', save = T, Dir = './')
t1 <- Sys.time()
bed_100bp_bg <- IdentifyBackground(organism = 'mm9', bed_path = './', binsize = 100, promo_bed = Promoter_anno, fa ='/illumina/iGenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa')
t2 <- Sys.time()
t2 - t1





#Promoter_anno <- GetPromoterAnno('hg18', save = T, Dir = './')
Promoter_anno <- fread('Promoter_Anno.bed', data.table= F)
#bed_100bp_bg <- IdentifyBackground(organism = 'hg19', bed_path = './', binsize = 100, promo_bed = Promoter_anno, fa ='/illumina/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa' , cores = 20)


# cl <- makeCluster(40, outfile = "Process.txt")
# registerDoParallel(cl)
registerDoMC(40)

chromlen <- as.data.frame(GetChromLength('hg18'))
chromlen <- chromlen[which(sapply(strsplit(chromlen$chrom, "_"), length) ==1),]
bed_100bp <- as.data.frame(CreateWindows(chromlen, 100))
rm(chromlen)
bed_path <- paste('./', 'bed_100bp.bed', sep = '')
storage.mode(bed_100bp[,2]) <-
  storage.mode(bed_100bp[,3]) <- 'integer'
cat("Write bed file to disk", "\t")
write.table(
  bed_100bp, bed_path, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE
)
Exclude_regions <-
  as.data.frame(ConstructExRegions('hg18', Promoter_anno))
Getfasta('/illumina/iGenomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa', bed_path)
Exclude_index <- ExcludeIntersection(bed_100bp, Exclude_regions)
rm(Exclude_regions)
letter_freq <- CountFreqency('./bed_100bp.fa')
Exclude_index1 <- which(letter_freq$ATCG != 1)
Exclude_index_all <-
  Reduce(dplyr::union, list(Exclude_index, Exclude_index1))
rm(list = c("Exclude_index", "Exclude_index1"))
bed_100bp_filtered <- 
  bind_cols(bed_100bp[-Exclude_index_all,], letter_freq[-Exclude_index_all,][1])
rm(list = c("Exclude_index_all", "bed_100bp", "letter_freq"))
colnames(bed_100bp_filtered) <-
  c('chrom', 'chromStart', 'chromEnd', 'GC')

t1 <- Sys.time()
Background_region <-
  foreach (
    i = 1:nrow(Promoter_anno), .options.multicore = list(preschedule = FALSE),
    .combine = bind_rows, .inorder = TRUE, .verbose =
      TRUE,
    .errorhandling = 'stop', .multicombine =
      TRUE
  ) %dopar% {
    cat('promoter', i, '\n')
    preset <-
      filter(bed_100bp_filtered, chrom == Promoter_anno$chrom[i])
    preset <-
      mutate(
        preset, proximity = abs(chromStart - Promoter_anno$chromStart[i]),
        geneName = Promoter_anno$name[i]
      )
    preset <-
      arrange(preset, proximity)
    k <- 0
    bgregion <- c()
    for (j in 1:nrow(preset)) {
      if (preset[j,4] < 0.4) {
        k <- k + 1
        bgregion <-
          bind_rows(bgregion, preset[j,])
      }
      if (k == 80)
        break
    }
    bgregion
  }
t2 <- Sys.time()
t2 - t1


rm(list = c("bed_100bp_filtered", "Promoter_anno"))




