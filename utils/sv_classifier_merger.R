library(dplyr)
library(data.table)
library(assertthat)
source("utils/sv_classifier_merger_utils.R")
source("utils/sv_classifier_probs.R")


pseudo_count = 1e-10
SV_states = c("ref","hetDel","homDel","hetDup","homDup","hetInv","homInv","hetIdup")




# 1) Read input data (prob table + LLR value)
print(snakemake@input)
print(snakemake@params)

LLR  = as.numeric(snakemake@params[["llr_cutoff"]])
prob = readRDS(snakemake@input[["prob"]]) 
# prob = readRDS("../run_20180401/sv_probabilities/simulation5-50000/50000_fixed.fraction15/raw_probabilities.Rdata")


# 2) Annotate prob table with consecutive seg_id
segment_ids = unique(prob[,.(chrom, from, to)])
segment_ids[, seg_id     := 1:.N, by = chrom]
segment_ids[, seg_id_max := .N,   by = chrom]
prob = merge(prob, segment_ids, by = c("chrom", "from", "to"))


# 3) Group SV states into gain, inv, loss
prob[, `:=`(P_hetDel  = p_hetDel,
            P_homDel  = p_homDel,
            P_hetDup  = p_hetDup,
            P_homDup  = p_homDup,
            P_hetInv  = p_hetInv,
            P_homInv  = p_homInv,
            P_hetIdup = p_hetIdup,
            p_gain    = log(exp(p_hetDup) + exp(p_homDup)),
            p_loss    = log(exp(p_homDel) + exp(p_hetDel)),
            p_inv     = log(exp(p_hetInv) + exp(p_homInv) + exp(p_hetIdup)),
            p_hetDel  = NULL,
            p_homDel  = NULL,
            p_hetDup  = NULL,
            p_homDup  = NULL,
            p_hetInv  = NULL,
            p_homInv  = NULL,
            p_hetIdup = NULL)]


# 4) Normalize SV probabilities per segment and cell
prob = unlog_and_normalize_probs(prob, pseudo_count)


# 5) calculate jump probabilities for all cells, using the SV classes gain, loss, inv
calc_jump_probabilities_byref(prob, c("ref", "inv","loss","gain"))



# 6) Apply three criteria in order to merge segments:
#    1) The segment or its left neighbor must be small enough
#    2) The jump criteria (mean_pr_merge > mean_pr_separate) must be given
#    3) A minimum number of cells with alternative state must be given
MAX_MERGE_SIZE = 5e5
MIN_ALT_CELLS  = 3
XXX = unique(prob[,.(chrom, start, seg_id, end, mean_pr_merge, mean_pr_separate, mean_pr_n)])
XXX[, `:=`(jump_criteria = !is.na(mean_pr_merge) & mean_pr_merge > mean_pr_separate,
           size_criteria = end - start <= MAX_MERGE_SIZE | c(F, end[1:(.N-1)] - start[1:(.N-1)] <= MAX_MERGE_SIZE),
           num_criteria  = !is.na(mean_pr_n) & (mean_pr_n >= MIN_ALT_CELLS)), 
    by = .(chrom)]

# 7) Group into consecutive segments to be merged
YYY = XXX[, .(seg_id, 
              start, 
              end, 
              criteria = jump_criteria & size_criteria & num_criteria,
              grouping = cumsum(abs(sign(c(0,diff(jump_criteria & size_criteria & num_criteria)))))), 
          by = chrom]
YYY[, new_seg_id := seg_id]
YYY[criteria == TRUE, new_seg_id := as.integer(min(new_seg_id - 1)), by = .(chrom, grouping)]



# 8) Report about the merging
YYY[, { if (nrow(.SD) > 1) {
        message("[SV merger] Merging ", nrow(.SD), " segments on chromosome ", chrom, ": ", 
                format_Mb(min(start)), " - ", format_Mb(max(start)) ) } }, 
    by = .(chrom, new_seg_id)]


# 9) do the final merging
prob = merge(prob, YYY[, .(chrom, seg_id, new_seg_id)], by = c("chrom", "seg_id"))
new_prob = prob[,
                .(start = min(start), 
                  end = max(end), 
                  nb_p = nb_p[1], 
                  state = paste(unique(state), collapse = ":"),
                  expected = sum(expected), 
                  W = sum(W), 
                  C = sum(C)), 
                by = .(chrom, new_seg_id, sample, cell)]

# 10) Check if segments with different state have been merged --> set those to state "SCE"
new_prob[grepl(":", state, fixed=T), state := "sce"]
new_prob[, scalar := 1]



# 11) Now, recompute probabilities for these segments
message("[SV merger] Annotating NB probabilities")
new_prob <- add_NB_probs(new_prob)



# 12) Extract SV call per locus
x = melt(new_prob, 
         id.vars = c("chrom","start","end","sample","cell","p_ref"), 
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







### Stop here

### TEST AREA - plots for a particular locus
if (FALSE) {
  
  carrier_cells = fread("simulation/variants/genome5-50000.txt")[chrom == "chr1" & start >= 78e6 & end <= 80e6]$cell
  
  plot_change_probs_for_single_segment(prob, "chr1", 7, c("inv","loss","gain","ref")) + 
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 78.8 Mb-78.9 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 7, .(V1 = W+C, segment = "chr1: 78.8 Mb-78.9 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0)
  
  
  plot_change_probs_for_single_segment(prob, "chr1", 8,  c("inv","loss","gain","ref")) + 
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 78.8 Mb-78.9 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 78.9 Mb-79 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 7, .(V1 = W+C, segment = "chr1: 78.8 Mb-78.9 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 8, .(V1 = W+C, segment = "chr1: 78.9 Mb-79 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0)
  
  
  plot_change_probs_for_single_segment(prob, "chr1", 9,  c("inv","loss","gain","ref")) + 
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 79 Mb-79.1 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 78.9 Mb-79 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 9, .(V1 = W+C, segment = "chr1: 79 Mb-79.1 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 8, .(V1 = W+C, segment = "chr1: 78.9 Mb-79 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0)
  
  
  plot_change_probs_for_single_segment(prob, "chr1", 10,  c("inv","loss","gain","ref")) + 
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 79 Mb-79.1 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 79.1 Mb-79.2 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 9, .(V1 = W+C, segment = "chr1: 79 Mb-79.1 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 10, .(V1 = W+C, segment = "chr1: 79.1 Mb-79.2 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0)
  
  
  plot_change_probs_for_single_segment(prob, "chr1", 11,  c("inv","loss","gain","ref")) + 
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 79.2 Mb-79.3 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 79.1 Mb-79.2 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 11, .(V1 = W+C, segment = "chr1: 79.2 Mb-79.3 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 10, .(V1 = W+C, segment = "chr1: 79.1 Mb-79.2 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0)
  
  
  plot_change_probs_for_single_segment(prob, "chr1", 12,  c("inv","loss","gain","ref")) + 
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 79.2 Mb-79.3 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_point(data = data.table(x = carrier_cells, segment = "chr1: 79.3 Mb-81.5 Mb"), 
               aes(x = x), inherit.aes = F, color = "black", y = 1, size = 2) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 11, .(V1 = W+C, segment = "chr1: 79.2 Mb-79.3 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0) +
    geom_text(data = prob[chrom == "chr1" & seg_id == 12, .(V1 = W+C, segment = "chr1: 79.3 Mb-81.5 Mb"), by = cell], 
              aes(x = cell, y = 0, label = V1), angle = 90, inherit.aes = F, hjust = 0)
  
  
  
  
  XXX = unique(prob[grepl('^chr1[0-7]$', chrom),.(chrom, start, seg_id, end, mean_pr_merge, mean_pr_separate, mean_pr_n)])
  ggplot(XXX) + 
    #geom_segment(aes(x = start, xend = end, y = mean_pr_separate, yend = mean_pr_separate), col = "dodgerblue3") + 
    #geom_segment(aes(x = start, xend = end, y = mean_pr_merge,    yend = mean_pr_merge), col = "red") + 
    geom_segment(aes(x = seg_id-1, xend = seg_id, y = mean_pr_separate, yend = mean_pr_separate), col = "dodgerblue3") +
    geom_segment(aes(x = seg_id-1, xend = seg_id, y = mean_pr_merge,    yend = mean_pr_merge), col = "red") +
    geom_text(aes(x = seg_id-1, y = 0, label = format_Mb(end-start)), hjust=0,vjust=0,size=3) +
    geom_text(aes(x = seg_id-1, y = 1, label = paste0("n=", mean_pr_n)), hjust=0,vjust=1,size=3) +
    facet_grid(chrom ~ .) +
    theme_minimal()
}


