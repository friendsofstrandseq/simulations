library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(assertthat)
source("utils/evaluation.utils.R")


#### START

input.truth = snakemake@input[["truth"]]
input.calls = snakemake@input[["calls"]]




### Read simulated SVs
message("[Evaluation] Reading ", length(input.truth), " simulated SV sets ...")
regex = "simulation_new/seed(\\d+)_size(\\d+)-(\\d+)_vaf(\\d+)-(\\d+)/variants-(\\d+).txt"
Truth = NULL
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
  t[, `:=`(SIMUL_id = as.integer(vars[2]), 
           SIMUL_minsize = as.integer(vars[3]),
           SIMUL_maxsize = as.integer(vars[4]),
           SIMUL_minvaf  = as.integer(vars[5]),
           SIMUL_maxvaf  = as.integer(vars[6]),
           SIMUL_binsize = as.integer(vars[7]))]
  Truth = rbind(Truth, t)
}



### Read SV calls
message("[Evaluation] Reading ", length(input.calls), " data sets ...")
regex = "sv_calls/seed(\\d+)_size(\\d+)-(\\d+)_vaf(\\d+)-(\\d+)-(\\d+)/(\\d+)_fixed\\.fraction(\\d+)/([a-zA-Z0-9_-]+)\\.txt"
Recall = NULL
Precision = NULL
for (f in input.calls) {
  
  vars = str_match(f, regex)
  d = fread(f)
  assert_that("chrom"    %in% colnames(d),
              "start"    %in% colnames(d),
              "end"      %in% colnames(d),
              "sample"   %in% colnames(d),
              "cell"     %in% colnames(d),
              "SV_class" %in% colnames(d)) %>% invisible

  # Get corresponding true SV calls
  truth = Truth[SIMUL_id      == as.integer(vars[2]) &
                SIMUL_minsize == as.integer(vars[3]) &
                SIMUL_maxsize == as.integer(vars[4]) & 
                SIMUL_minvaf  == as.integer(vars[5]) &
                SIMUL_maxvaf  == as.integer(vars[6]), ]
  assert_that(nrow(truth) > 0) %>% invisible
  
  
  if (nrow(truth)==0) {
    message("Skipping file because there are no simulated SVs: ", f)
    next
  }
  
  # Calculate precision and recall
  rp = recall_precision(truth, d)

  recall = rp[["recall"]]
  recall[, `:=`(SIMUL_id        = as.integer(vars[2]),
                SIMUL_minsize   = as.integer(vars[3]),
                SIMUL_maxsize   = as.integer(vars[4]),
                SIMUL_minvaf    = as.integer(vars[5]),
                SIMUL_maxvaf    = as.integer(vars[6]),
                SIMUL_binsize   = as.integer(vars[7]),
                SIMUL_fraction  = as.integer(vars[9]),
                SIMUL_method    = vars[10]) ]
  Recall = rbind(Recall, recall)
  
  precision = rp[["precision"]]
  if (nrow(precision) == 0) {
    precision[, `:=`(SIMUL_id        = integer(),
                     SIMUL_minsize   = integer(),
                     SIMUL_maxsize   = integer(),
                     SIMUL_minvaf    = integer(),
                     SIMUL_maxvaf    = integer(),
                     SIMUL_binsize   = integer(),
                     SIMUL_fraction  = integer(),
                     SIMUL_method    = character()) ]
  } else {
    precision[, `:=`(SIMUL_id        = as.integer(vars[2]),
                     SIMUL_minsize   = as.integer(vars[3]),
                     SIMUL_maxsize   = as.integer(vars[4]),
                     SIMUL_minvaf    = as.integer(vars[5]),
                     SIMUL_maxvaf    = as.integer(vars[6]),
                     SIMUL_binsize   = as.integer(vars[7]),
                     SIMUL_fraction  = as.integer(vars[9]),
                     SIMUL_method    = vars[10]) ]
  }
  Precision = rbind(Precision, precision)
}


 ### Summarize across all the different simulations
xxx = Recall[, .(N        = .N,
                 N_       = nrow(unique(.SD[,.(chrom, start, end)])),
                 bp       = sum(matches_call)/.N,
                 bp_sv    = mean(correct_gt/SV_vaf),
                 bp_sv_gt = mean(correct_sv/SV_vaf)),
             by = .(SIMUL_minsize,
                    SIMUL_maxsize,
                    SIMUL_minvaf,
                    SIMUL_maxvaf,
                    SIMUL_fraction)]
assert_that(xxx[, all(N == N_)]) %>% invisible
yyy = Precision[, .(N        = .N,
                    bp       = sum(matches_SV)/.N,
                    bp_sv    = mean(correct_sv/SV_vaf),
                    bp_sv_gt = mean(correct_gt/SV_vaf)),
                by = .(SIMUL_minsize,
                       SIMUL_maxsize,
                       SIMUL_minvaf,
                       SIMUL_maxvaf,
                       SIMUL_fraction)]
zzz = merge(xxx,yyy, suffixes = c(".recall", ".precision"),
            by = c("SIMUL_minsize", "SIMUL_maxsize", "SIMUL_minvaf", "SIMUL_maxvaf", "SIMUL_fraction"))


# Beatify data set prior to plotting
format_Mb = function(x) {
  ifelse(x>=1e6, 
         paste(round(x/1e6,1),"Mb"),
         ifelse(x>=1e3,
                paste(round(x/1e3,1),"kb"),
                paste(x, "bp")))
  return (x)
}
zzz[, SV_size := paste0(format_Mb(SIMUL_minsize),"-",format_Mb(SIMUL_maxsize))]
zzz[, SV_vaf  := factor(paste0(SIMUL_minvaf,"-",SIMUL_maxvaf,"%"),
                        levels = unique(paste0(SIMUL_minvaf,"-",SIMUL_maxvaf,"%"))[order(unique(SIMUL_minvaf))],
                        ordered = T)]
zzz = zzz[order(SIMUL_fraction, SIMUL_minvaf, SIMUL_minsize)]

p <- ggplot(zzz) +
  geom_path(aes(bp.recall, bp.precision), 
            linetype = "solid", color = "darkgrey") +
  geom_point(aes(bp.recall, bp.precision, col = SIMUL_fraction)) +
  geom_path(aes(bp_sv.recall, bp_sv.precision), 
            linetype = "dashed") +
  geom_point(aes(bp_sv.recall, bp_sv.precision, col = SIMUL_fraction), 
             shape = 17) +
  geom_text(x = 0, y = 0, hjust = 0, vjust = 0, aes(label = paste("SVs =",V1)),
            data = zzz[, N.recall[1], by = .(SV_size, SV_vaf)]) +
  geom_text(x = 0, y = 0.1, hjust = 0, vjust = 0, aes(label = paste("Calls =",V1, "-", V2)),
            data = zzz[, .(min(N.precision),max(N.precision)), by = .(SV_size, SV_vaf)]) +
  facet_grid(SV_size ~ SV_vaf) +
  theme_bw() + theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0,1), xlim = c(0,1))

cairo_pdf(file = snakemake@output[[1]], width=21, height = 14, onefile=T)
print(p)
dev.off()

