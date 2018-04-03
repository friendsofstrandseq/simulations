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
Results = NULL
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
  
  recall = categorization(rp[["recall"]])
  recall[, `:=`(SIMUL_simulation_id = as.integer(vars[2]),
                SIMUL_window_size   = as.integer(vars[3]),
                SIMUL_segmentation  = as.integer(vars[5]),
                precision = rp[["precision"]],
                method = vars[6])]
  
  Results = rbind(Results, recall)
}

message("[Evaluation] Number of simulated cases:")
table(Results[SIMUL_segmentation==2]$SV_vaf_factor, Results[SIMUL_segmentation==2]$SV_size_factor) %>% print


# Summarize recall:
xxx = Results[, .(N  = .N,
                  N_ = nrow(unique(.SD[,.(chrom, start, end)])),
                  breakpoint_found = sum(breakpoint_found)/.N,
                  correct_gt       = mean(correct_gt/SV_vaf),
                  correct_gt_q1    = quantile(correct_gt/SV_vaf, 0.25),
                  correct_gt_q2    = quantile(correct_gt/SV_vaf, 0.75),
                  correct_sv       = mean(correct_sv/SV_vaf),
                  precision = mean(unique(precision))),
              by = .(SV_real, SV_vaf_factor, SV_size_factor, SIMUL_segmentation)]
assert_that(xxx[, all(N == N_)]) %>% invisible


cairo_pdf(file = snakemake@output[[1]], width=16, height = 14, onefile=T)
for (svclass in c("het_del","hom_del","het_dup","hom_dup","het_inv","hom_inv","inv_dup")) {
  message("[Evaluation] Plotting ", svclass)
  p <- ggplot(xxx[SV_real == svclass]) + 
    aes(breakpoint_found, precision) +
    geom_line() +
    facet_grid(SV_size_factor ~ SV_vaf_factor) + 
    geom_point(aes(col = SIMUL_segmentation)) +
    geom_label(x=0, y=0, aes(label = paste0("N=",N)), hjust = 0, vjust = 0) +
    theme_bw() + theme(legend.position = "bottom") +
    ggtitle(svclass)
  print(p)
}
dev.off()
