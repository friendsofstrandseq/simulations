library(assertthat)
library(data.table)
library(dplyr)

LLR = 1

prob = readRDS(snakemake@input[[1]]) # prob = readRDS("sv_probabilities/simulation5-50000/50000_fixed.few/raw_probabilities.Rdata")

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

