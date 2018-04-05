library(data.table)
library(dplyr)
library(stringr)
library(assertthat)
library(ggplot2)
source("utils/evaluation.utils.R")

input.truth = c("simulation/variants/genome5-50000.txt", 
                "simulation/variants/genome6-50000.txt", 
                "simulation/variants/genome7-50000.txt", 
                "simulation/variants/genome8-50000.txt")
input.calls = c("sv_calls/simulation5-50000/50000_fixed.fraction02/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction02/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction02/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction02/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction05/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction05/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction05/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction05/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction08/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction08/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction08/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction08/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction12/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction12/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction12/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction12/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction15/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction15/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction15/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction15/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction18/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction18/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction18/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction18/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction20/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction20/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction20/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction20/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction25/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction25/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction25/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction25/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction30/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction30/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction30/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction30/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction35/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction35/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction35/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction35/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction45/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction45/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction45/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction45/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction55/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction55/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction55/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction55/simple.txt", 
                "sv_calls/simulation5-50000/50000_fixed.fraction65/simple.txt", 
                "sv_calls/simulation6-50000/50000_fixed.fraction65/simple.txt", 
                "sv_calls/simulation7-50000/50000_fixed.fraction65/simple.txt", 
                "sv_calls/simulation8-50000/50000_fixed.fraction65/simple.txt")

input.truth = snakemake@input[["truth"]]
input.calls = snakemake@input[["calls"]]



### Read simulated SVs
regex = "simulation/variants/genome(\\d+)-(\\d+)\\.txt"
truth = NULL
for (f in input.truth) {
  vars = str_match(f, regex)
  t = fread(f)
  assert_that("chrom"   %in% colnames(t),
              "start"   %in% colnames(t),
              "end"     %in% colnames(t),
              "SV_type" %in% colnames(t),
              "sample"  %in% colnames(t),
              "cell"    %in% colnames(t),
              all(t$SV_type %in% c("hom_dup","het_dup","hom_inv","het_inv",
                                   "hom_del","het_del","inv_dup","false_del"))) %>%
              invisible
  truth = rbind(truth, cbind(t, SIMUL_simulation_id = as.integer(vars[2]), SIMUL_window_size = as.integer(vars[3])))
}




### Read SV calls one by one
Recall = NULL
Precision = NULL
regex = "sv_calls/simulation(\\d+)-(\\d+)/(\\d+)_fixed\\.fraction(\\d+)/([a-zA-Z0-9_-]+)\\.txt"
for (f in input.calls) {
  message("[Evaluation] Reading ", f)
  vars = str_match(f, regex)
  d = fread(f)
  assert_that("chrom"    %in% colnames(d),
              "start"    %in% colnames(d),
              "end"      %in% colnames(d),
              "sample"   %in% colnames(d),
              "cell"     %in% colnames(d),
              "SV_class" %in% colnames(d)) %>% invisible
  
  # Get recall and precision
  truth_ = truth[SIMUL_simulation_id == as.integer(vars[2])]
  rp = recall_precision(truth_, d)
  recall = rp[["recall"]]
  precision = rp[["precision"]]
  
  # Major conceptual caveat
  # =======================
  # The categorization (in SV size and VAF) for recall is based on the true SVs, whereas
  # it is based on the predicted SVs for precision.
  # This means for example that precision at high VAFs will always be 100% (because it
  # would require many cells to be consistently called wrong, which is very unlikely
  # in simulated data)
  recall = categorization(recall, n_cells = max(max(recall$SV_vaf,100)))
  precision = categorization(precision, n_cells = max(max(recall$SV_vaf,100)))

  recall[, `:=`(SIMUL_simulation_id = as.integer(vars[2]),
                SIMUL_window_size   = as.integer(vars[3]),
                SIMUL_segmentation  = as.integer(vars[5]),
                method = vars[6])]
  precision[, `:=`(SIMUL_simulation_id = as.integer(vars[2]),
                SIMUL_window_size   = as.integer(vars[3]),
                SIMUL_segmentation  = as.integer(vars[5]),
                method = vars[6])]

  Recall = rbind(Recall, recall)
  Precision = rbind(Precision, precision)
}

assert_that(all(!is.na(Precision$SV_size_factor)),
            all(!is.na(Precision$SV_vaf_factor))) %>% invisible


# Summarize results:
xxx = Recall[, .(N            = .N,
                 N_           = nrow(unique(.SD[,.(chrom, start, end)])),
                 matches_call = sum(matches_call)/.N,
                 correct_gt   = mean(correct_gt/SV_vaf),
                 correct_sv   = mean(correct_sv/SV_vaf)),
              by = .(SV_size_factor, SV_vaf_factor, SIMUL_segmentation, method)]
assert_that(xxx[, all(N == N_)]) %>% invisible


yyy = Precision[, .(N            = .N,
                    matches_SV   = sum(matches_SV)/.N,
                    correct_gt   = mean(correct_gt/SV_vaf),
                    correct_sv   = mean(correct_sv/SV_vaf)),
                by = .(SV_size_factor, SV_vaf_factor, SIMUL_segmentation, method)]


zzz = merge(xxx[,.(SIMUL_segmentation, method, SV_size_factor, SV_vaf_factor,
                   N, recall1 = matches_call, recall2 = correct_sv, recall3 = correct_gt)],
            yyy[,.(SIMUL_segmentation, method, SV_size_factor, SV_vaf_factor,
                   N, precision1 = matches_SV, precision2 = correct_sv, precision3 = correct_gt)],
            by = c("SIMUL_segmentation", "method", "SV_size_factor", "SV_vaf_factor"))


cairo_pdf(file = snakemake@output[[1]], width=16, height = 12, onefile=T)
  p <- ggplot(zzz) +
    geom_line(aes(recall2, precision2), linetype = "dashed", color = "darkgrey") +
    geom_point(aes(recall2, precision2, col = SIMUL_segmentation), shape = 17) +
    geom_line(aes(recall1, precision1)) +
    geom_point(aes(recall1, precision1, col = SIMUL_segmentation)) +
    geom_text(data = zzz[, .(min = min(N.y), max = max(N.y)),
                         by = .(SV_size_factor, SV_vaf_factor)],
              x=0,y=0,hjust=0,vjust=0, aes(label = paste("Calls =",min,"-",max))) +
    geom_text(data = zzz[, mean(N.x),
                         by = .(SV_size_factor, SV_vaf_factor)],
              x=0,y=0.1,hjust=0,vjust=0, aes(label = paste("True SVs =",V1))) +
    facet_grid(SV_vaf_factor ~ SV_size_factor) +
    theme_bw() + theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0,1), xlim = c(0,1))
  print(p)
dev.off()
