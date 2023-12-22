library(getopt)
suppressMessages(library(devtools))

args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args), collapse = " ")
listoptions <- unlist(strsplit(hh, "--"))[-1]
options.args <- sapply(listoptions, function(x) {
    unlist(strsplit(x, " "))[-1]
})
options.names <- sapply(listoptions, function(x) {
    option <- unlist(strsplit(x, " "))[1]
})
names(options.args) <- unlist(options.names)
fragfile <- options.args[1]
outdir <- options.args[2]
binfile <- options.args[3]
targetfile <- options.args[4]

if (targetfile == "nova") {
    library(PlasmaToolsNovaseq.hg19)
    targetfile <- target20
} else if (targetfile == "hi") {
    library(PlasmaToolsHiseq.hg19)
    targetfile <- target55
}


sample <- gsub(".rds", "", basename(fragfile))
print(sample)

library(data.table)
suppressMessages(library(devtools))
suppressMessages(library(GenomicRanges))
suppressMessages(library(PlasmaTools.lucas))

setDTthreads(1)

target <- as.data.table(targetfile)
bins <- fread(binfile)
bins <- bins[map >= 0.90 & gc >= 0.30]

fragments <- readRDS(fragfile)

fragments <- keepSeqlevels(fragments, paste0("chr", c(1:22, "X")),
                           pruning.mode="coarse")
fragments <- as.data.table(fragments)
setnames(fragments,  "seqnames", "chr")
#
### Filter blacklisted regions
filters <- as.data.table(filters.hg19)
setnames(filters, "seqnames",  "chr")
setkey(filters, chr, start, end)
fragdt <- foverlaps(fragments, filters, type="any")
fragdt <- fragdt[is.na(start)][,`:=`(start=NULL, end=NULL, width=NULL,
                                     strand=NULL, name=NULL, score=NULL)]

setnames(fragdt, c("i.start", "i.end", "i.width",  "i.strand"),
         c("start", "end", "width", "strand"))
##########
setnames(target, c("seqnames", "gcmed"), c("chr", "target"))
fragdt <- gcCorrectTarget(fragments, target)

bins <- bins[,chr:=factor(chr, paste0("chr", c(1:22, "X")))]
fragdt <- fragdt[,chr:=factor(chr, paste0("chr", c(1:22, "X")))]

setnames(fragdt, "gc", "fraggc")
setkey(fragdt, chr, start, end)
setkey(bins, chr, start, end)

bins[,bin:=1:.N]
bins2 <- binFrags(fragdt, bins)

bins2[,id:=gsub("_downsamp", "", sample)]
setcolorder(bins2, c("id", "chr", "start", "end", "bin"))


filename <- file.path(outdir, paste0(sample, "_5mb.csv"))

fwrite(bins2, filename)
q('no')
