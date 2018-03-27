library(data.table)
library(assertthat)
library(pheatmap)
source("utils/sv_classifier_probs.R")
source("utils/sv_classifier_counts.R")



# counts = fread(paste("zcat","counts/simulation23-100000/100000_fixed.txt.gz"))
# info   = fread("counts/simulation23-100000/100000_fixed.info")
# strand = fread("strand_states/simulation23-100000/final.txt")
# segs   = fread("segmentation2/simulation23-100000/100000_fixed.few.txt") 

# counts = fread(paste("zcat","counts/simulation7-50000/50000_fixed.txt.gz"))
# info   = fread("counts/simulation7-50000/50000_fixed.info")
# strand = fread("strand_states/simulation7-50000/final.txt")
# segs   = fread("segmentation2/simulation7-50000/50000_fixed.few.txt")

counts = fread(paste("zcat",snakemake@input[["counts"]]))
info   = fread(snakemake@input[["info"]])
strand = fread(snakemake@input[["states"]])
segs   = fread(snakemake@input[["bp"]])






################################################################################
# Check input data 
#
# counts
assert_that("chrom" %in% colnames(counts),
            "start" %in% colnames(counts),
            "end"   %in% colnames(counts),
            "sample"%in% colnames(counts),
            "cell"  %in% colnames(counts))
counts <- counts[order(sample,cell,chrom,start,end),]
setkey(counts,sample,cell)

# info
assert_that("sample"%in% colnames(info),
            "cell"  %in% colnames(info),
            "nb_p"  %in% colnames(info),
            "nb_r"  %in% colnames(info),
            "nb_a"  %in% colnames(info),
            "pass1" %in% colnames(info))
info <- info[order(sample,cell),]
setkey(info,sample,cell)

# strand
assert_that("sample"%in% colnames(strand),
            "cell"  %in% colnames(strand),
            "chrom" %in% colnames(strand),
            "start" %in% colnames(strand),
            "end"   %in% colnames(strand),
            "class" %in% colnames(strand))
strand <- strand[order(sample,cell,chrom,start,end),]

# segs
assert_that("chrom" %in% colnames(segs),
            "bps"   %in% colnames(segs))
segs <- segs[order(chrom, bps),]

message("[SV classifier] Problem size: ", nrow(info), " cells x ", nrow(segs), " segments.")
################################################################################






# Kick out non-PASS cells     # To test, use something like `info[seq(1,nrow(info),3)]$pass1 = 0`
if (nrow(info[pass1 != 1])> 0) message("[SV classifier] Kicking out ", nrow(info[pass1 != 1]), " low quality cells. ", nrow(info[pass1 == 1]), " remain.")
info <- info[pass1 == 1,]
counts <- counts[ paste(sample,cell) %in% info[,paste(sample,cell)] ]
assert_that(all(unique(counts[,.(sample,cell)]) == unique(info[,.(sample,cell)])))


# Get mean and median from count data
# When calculated, ignore "None" bins and also use the trimmed mean!
info <- merge(info, counts[, .(mean = mean((w+c)[class != "None"], trim = 0.05),
                               median = median((w+c)[class != "None"])), by = .(sample, cell)],
              by = c("sample","cell"))



# Prepare function to assess the strand_state of a cell/chrom/region
bins <- unique(counts[, .(chrom, start, end)])
get_strand_state <- function(sample_, cell_, chrom_, from_, to_) {
    x = strand[sample == sample_ & cell == cell_ & chrom == chrom_]
    assert_that(nrow(x)>0)
    if (nrow(x) == 1) return (x$class)
    min_pos = bins[chrom == chrom_]$start[from_]
    max_pos = bins[chrom == chrom_]$end[to_]
    x = x[start <= min_pos & end >= max_pos]
    if (nrow(x) == 1) return (x$class)
    return ("sce")
}


# Expand segments into {chrom, [from, to]}
segs[, from := shift(bps,fill = 0) + 1, by = chrom]
segs[, `:=`(to = bps, bps = NULL, k = NULL)]




# Do all the work:
message("[SV classifier] Preparing large table [segments + cells x features]")
probs <- segs[,cbind(.SD,info[,.(sample,cell,nb_p,nb_r,medbin,mean)]), by = .(chrom,from)]

message("[SV classifier] Annotating strand-state (slow)")
probs[, state := get_strand_state(sample, cell, chrom, from, to), by = .(sample, cell, chrom, from, to)]

message("[SV classifier] Annotating expected coverage")
probs[, expected := (to - from +1)*mean, by = .(sample, cell, chrom, from, to)]

message("[SV classifier] Annotating observed W/C counts")
probs <- add_seg_counts(probs, counts)
probs[, scalar := 1]

message("[SV classifier] Annotating NB probabilities")
probs <- add_NB_probs(probs)

message("[SV classifier] Post-processing NB probabilities, e.g. adding a prior")
probs[state == "sce", `:=`(p_ref     = -1000, 
                           p_hetInv  = -1000, 
                           p_homInv  = -1000, 
                           p_hetDup  = -1000,
                           p_homDup  = -1000,
                           p_hetDel  = -1000, 
                           p_homDel  = -1000,
                           p_hetIdup = -1000)]
PRIORS = c(ref     = 0.77, 
           homInv  = 0.05, 
           hetInv  = 0.05, 
           hetDup  = 0.05,
           homDup  = 0.01,
           hetDel  = 0.05, 
           homDel  = 0.01,
           hetIdup = 0.01)
probs[, `:=`(p_ref    = p_ref    + log(PRIORS["ref"]),
             p_homInv = p_homInv + log(PRIORS["homInv"]),
             p_hetInv = p_hetInv + log(PRIORS["hetInv"]),
             p_hetDup = p_hetDup + log(PRIORS["hetDup"]),
             p_homDel = p_homDel + log(PRIORS["homDel"]),
             p_hetDel = p_hetDel + log(PRIORS["hetDel"]))]
# Add log likelihood ratio log( p(SV) / p(REF) )
probs[,max_loklikratio := pmax(p_hetInv, p_hetDel, p_hetDup, p_homDel, p_homDup, p_homInv) - p_ref]

probs[, `:=`(start = bins[from]$start, end = bins[to]$end)]
saveRDS(probs, snakemake@output[[1]])


