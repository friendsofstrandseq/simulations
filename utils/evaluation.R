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
  
  # Major issue is that categorization does NOT work for precision. Only for recall
  #recall = categorization(rp[["recall"]])

  # Hence I go without categorization for now
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

Results

# Summarize recall:
xxx = Results[, .(N  = .N,
                  N_ = nrow(unique(.SD[,.(chrom, start, end)])),
                  breakpoint_found = sum(breakpoint_found)/.N,
                  correct_gt       = mean(correct_gt/SV_vaf),
                  correct_sv       = mean(correct_sv/SV_vaf)),
              by = .(SV_real, SIMUL_segmentation, method)]
assert_that(xxx[, all(N == N_)]) %>% invisible
xxx


yyy = Precision[, .(loci_80  = .SD[type_of_precision == "Loci at 80% overlap", sum(value_of_precision)],
              loci_any = .SD[type_of_precision == "Loci at any overlap", sum(value_of_precision)],
              loci_all = .SD[type_of_precision == "all loci", sum(value_of_precision)]), 
          by = SIMUL_segmentation]
yyy[,`:=`(precision_80  = loci_80 / loci_all,
          precision_any = loci_any / loci_all)]

zzz = merge(xxx, yyy[,.(SIMUL_segmentation, precision_80, precision_any)], by = "SIMUL_segmentation")

cairo_pdf(file = snakemake@output[[1]], width=16, height = 14, onefile=T)
p <- ggplot(zzz) +
    geom_line(aes(breakpoint_found, precision_80, col = SV_real)) +
    geom_line(aes(correct_sv, precision_80, col = SV_real), linetype = "dashed") +
    #geom_point(aes(breakpoint_found, precision, col = SIMUL_segmentation)) +
    #geom_point(aes(breakpoint_found, correct_sv, col = SIMUL_segmentation)) +
    #geom_label(x=0, y=0, aes(label = paste0("N=",N)), hjust = 0, vjust = 0) +
    theme_bw() + theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0,1), xlim = c(0,1))
  print(p)
}
dev.off()
