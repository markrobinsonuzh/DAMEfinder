#plots for README.md

load("data/tupleASM_fullCancer.RData")
load("data/derASM_fullcancer2.RData")
metadata <- read.table("data/fullCancerSamples.txt", stringsAsFactors = FALSE)
colData(ASM)$group <- metadata$V2
colData(derASM)$group <- metadata$V2
colData(ASM)$samples <- colnames(ASM)

DAME <- GRanges(9, IRanges(99984206,99984364))
myColor <- RColorBrewer::brewer.pal(9, "Set1")


##SNP ASM
dame_track(DAME, derASM = derASM[,seq(2,12,2)], colvec = myColor, window = 3)
ggsave("../DAMEfinder_git/vignettes/DAME_snp_allsamps.png")
dame_track_median(DAME, derASM = derASM[,seq(2,12,2)], colvec = myColor, window = 3)
ggsave("../DAMEfinder_git/vignettes/DAME_snp_allsamps_median.png")

reference_file <- "/home/Shared_taupo/data/annotation/Human/Ensembl_GRCh37.91/GRCh37.91.fa"
metadata <- read.table("data/fullCancerSamples.txt", stringsAsFactors = FALSE)
sample_names <- metadata$V1
path <- "/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/"
snp <- GRanges(9, IRanges(99984349, width = 1))

allps <- mapply(methyl_circle_plot,
                vcfFile = gsub("/home/",path, metadata$V4)[c(2,8)],
                bamFile = gsub("/home/",path,metadata$V3)[c(2,8)],
                sampleName = sample_names[c(2,8)],
                MoreArgs=list(
                  snp = snp,
                  refFile = gsub("/home/",path,reference_file),
                  dame = DAME,
                  sampleReads = TRUE,
                  numReads = 15,
                  pointSize = 2,
                  letterSize = 3
                ),
                SIMPLIFY = FALSE)

cowplot::plot_grid(plotlist = allps, nrow = 2, ncol = 1)
ggsave("../DAMEfinder_git/vignettes/DAME_snp_sampledreads.png")

##Tuple ASM
dame_track(DAME, ASM = ASM[,seq(2,12,2)], colvec = myColor, window = 3)
ggsave("../DAMEfinder_git/vignettes/DAME_tuple_allsamps.png")

dame_track_median(DAME, ASM = ASM[,seq(2,12,2)], colvec = myColor, window = 3)
ggsave("../DAMEfinder_git/vignettes/DAME_tuple_allsamps_median.png")

allps <- mapply(methyl_circle_plotCpG,
                #vcfFile = gsub("/home/",path, metadata$V4)[c(2,8)],
                bamFile = gsub("/home/",path,metadata$V3)[c(2,8)],
                sampleName = sample_names[c(2,8)],
                MoreArgs=list(
                  cpgsite = DAME,
                  #snp = snp,
                  refFile = gsub("/home/",path,reference_file),
                  dame = DAME,
                  sampleReads = TRUE,
                  numReads = 30,
                  pointSize = 2
                ),
                SIMPLIFY = FALSE)

cowplot::plot_grid(plotlist = allps, nrow = 2, ncol = 1)
ggsave("../DAMEfinder_git/vignettes/DAME_only30reads.png")

