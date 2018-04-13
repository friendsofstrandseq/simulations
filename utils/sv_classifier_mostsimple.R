library(assertthat)
library(data.table)
library(dplyr)

LLR = as.numeric(snakemake@params[["llr_cutoff"]])
prob = readRDS(snakemake@input[[1]]) # prob = readRDS("sv_probabilities/simulation5-50000/50000_fixed.few/raw_probabilities.Rdata")


#message("[SV classifier] Post-processing NB probabilities, e.g. adding a prior")
#prob[state == "sce", `:=`(p_ref     = -1000,
#                           p_hetInv  = -1000,
#                           p_homInv  = -1000,
#                           p_hetDup  = -1000,
#                           p_homDup  = -1000,
#                           p_hetDel  = -1000,
#                           p_homDel  = -1000,
#                           p_hetIdup = -1000)]
#PRIORS = c(ref     = 0.77,
#           homInv  = 0.05,
#           hetInv  = 0.05,
#           hetDup  = 0.05,
#           homDup  = 0.01,
#           hetDel  = 0.05,
#           homDel  = 0.01,
#           hetIdup = 0.01)
#prob[, `:=` (p_ref     = p_ref     + log(PRIORS["ref"]),
#             p_homInv  = p_homInv  + log(PRIORS["homInv"]),
#             p_hetInv  = p_hetInv  + log(PRIORS["hetInv"]),
#             p_hetDup  = p_hetDup  + log(PRIORS["hetDup"]),
#             p_homDup  = p_homDup  + log(PRIORS["homDup"]),
#             p_homDel  = p_homDel  + log(PRIORS["homDel"]),
#             p_hetDel  = p_hetDel  + log(PRIORS["hetDel"]),
#             p_hetIdup = p_hetIdup + log(PRIORS["hetIdup"]))]


x = melt(prob, id.vars = c("chrom","start","end","sample","cell","p_ref"), 
     measure.vars = c("p_homInv", "p_hetInv", "p_hetDel","p_homDel", "p_homDup", "p_hetDup", "p_hetIdup"),
     variable.name = "SV_class",
     value.name    = "loglik")
x = x[, .SD[loglik == max(loglik)][1], by = .(chrom, start, end, sample, cell)]

# Filter by log likelihood ratio
x = x[loglik - p_ref > LLR]

# Rename SV classes
rename_svs = data.table(SV_class = c("p_homInv", "p_hetInv", "p_hetDel", "p_homDel", "p_homDup", "p_hetDup", "p_hetIdup"),
                        SV_class_renamed =      c("inv_hom",  "inv_h1",  "del_h1",  "del_hom",  "dup_hom",  "dup_h1",  "idup_h1"))
assert_that(all(x$SV_class %in% rename_svs$SV_class)) %>% invisible
x = merge(x, rename_svs, by = "SV_class")

write.table(x[, .(chrom, start, end, sample, cell, SV_class = SV_class_renamed, loglikratio = loglik - p_ref)],
            file = snakemake@output[[1]], quote=F, sep="\t", row.names = F)

