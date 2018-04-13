library(assertthat)
library(data.table)
library(dplyr)

MIN_CELLS = as.integer(snakemake@params[["mincells"]])
prob = readRDS(snakemake@input[[1]]) # prob = readRDS("sv_probabilities/simulation5-50000/50000_fixed.fraction12/raw_probabilities.Rdata")


#message("[SV classifier] Adding a prior")
#prob[state == "sce", `:=`(p_ref     = -1000,
#                          p_hetInv  = -1000,
#                          p_homInv  = -1000,
#                          p_hetDup  = -1000,
#                          p_homDup  = -1000,
#                          p_hetDel  = -1000,
#                          p_homDel  = -1000,
#                          p_hetIdup = -1000)]
#PRIORS = c(ref     = 0.3,
#           homInv  = 0.1,
#           hetInv  = 0.1,
#           hetDup  = 0.1,
#           homDup  = 0.1,
#           hetDel  = 0.1,
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



message("[SV classifier] Model across cells. Allow only one SV class, seen in at least ", MIN_CELLS, " cells")
mod = prob[, data.table( model = c("ref","hetDel","homDel","hetInv","homInv","hetDup","homDup","hetIdup"),
                         loglik = c(sum(p_ref),
                                    sum(log(exp(p_ref) + exp(p_hetDel))),
                                    sum(log(exp(p_ref) + exp(p_homDel))),
                                    sum(log(exp(p_ref) + exp(p_hetInv))),
                                    sum(log(exp(p_ref) + exp(p_homInv))),
                                    sum(log(exp(p_ref) + exp(p_hetDup))),
                                    sum(log(exp(p_ref) + exp(p_homDup))),
                                    sum(log(exp(p_ref) + exp(p_hetIdup)))), 
                         num    = c(sum(p_ref == pmax(p_ref, p_hetDel, p_homDel, p_homInv, p_hetInv, p_hetDup, p_homDup, p_hetIdup)),
                                    sum(p_hetDel  > p_ref),
                                    sum(p_homDel  > p_ref),
                                    sum(p_hetInv  > p_ref),
                                    sum(p_homInv  > p_ref),
                                    sum(p_hetDup  > p_ref),
                                    sum(p_homDup  > p_ref),
                                    sum(p_hetIdup > p_ref)) )[order(loglik, decreasing = T)],
          by = .(chrom, from, to)]

# Pick model with highest likelihood for each segment
mod = mod[, .SD[num>=MIN_CELLS][1,], by = .(chrom, from, to)]



# Apply the best model to the prob by overwriting probabilities
newprobs = merge(prob, mod[, .(chrom, from, to, model)], by = c("chrom","from","to"))
newprobs[model == "ref", p_ref := 0]
newprobs[model == "hetDup" & p_hetDup > p_ref,         p_hetDup := 0]
newprobs[model == "hetDup" & p_hetDup <= p_ref,        p_ref    := 0]
newprobs[model == "homDup" & p_homDup > p_ref,         p_homDup := 0]
newprobs[model == "homDup" & p_homDup <= p_ref,        p_ref    := 0]
newprobs[model == "hetDel" & p_hetDel > p_ref,         p_hetDel := 0]
newprobs[model == "hetDel" & p_hetDel <= p_ref,        p_ref    := 0]
newprobs[model == "hetInv" & p_hetInv > p_ref,         p_hetInv := 0]
newprobs[model == "hetInv" & p_hetInv <= p_ref,        p_ref    := 0]
newprobs[model == "homInv" & p_homInv > p_ref + 0.01,  p_homInv := 0]
newprobs[model == "homInv" & p_homInv <= p_ref + 0.01, p_ref := 0]
newprobs[model == "homDel" & p_homDel > p_ref,         p_homDel := 0]
newprobs[model == "homDel" & p_homDel <= p_ref,        p_ref    := 0]
newprobs[model == "hetIdup" & p_hetIdup > p_ref,       p_hetIdup := 0]
newprobs[model == "hetIdup" & p_hetIdup <= p_ref,      p_ref    := 0]

newprobs[model != "ref"]

x = melt(newprobs[model!="ref"], 
         id.vars = c("chrom","start","end","sample","cell", "state"), 
         measure.vars = c("p_ref","p_homInv", "p_hetInv", "p_hetDel","p_homDel", "p_homDup", "p_hetDup", "p_hetIdup"),
         variable.name = "SV_class",
         value.name    = "loglik")

setkey(x, chrom, start, end, sample, cell)
x = x[, .SD[order(loglik, decreasing=T)][1], by = .(chrom, start, end, sample, cell)]
x = x[SV_class != "p_ref"]


# Rename SV classes
rename_svs = data.table(SV_class =              c("p_homInv", "p_hetInv", "p_hetDel", "p_homDel", "p_homDup", "p_hetDup", "p_hetIdup"),
                        SV_class_renamed =      c("inv_hom",  "inv_h1",   "del_h1",   "del_hom",  "dup_hom",  "dup_h1",   "idup_h1"))
assert_that(all(x$SV_class %in% rename_svs$SV_class)) %>% invisible
x = merge(x, rename_svs, by = "SV_class")

write.table(x[, .(chrom, start, end, sample, cell, SV_class = SV_class_renamed, loglikratio = NA)],
            file = snakemake@output[[1]], quote=F, sep="\t", row.names = F)




