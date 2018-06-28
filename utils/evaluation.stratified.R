library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(assertthat)
source("utils/evaluation.utils.R")

# SETTINGS:
# minimum overlap require to say that a predicted SV matches a true call.
min_bp_overlap   = 1
# minimum fraction of a true SV to be coverd by SV prediction in order to be considered "found"
min_frac_covered = 0.8


#### START

input.truth = snakemake@input[["truth"]]
input.calls = snakemake@input[["calls"]]




### Read simulated SVs
message("[Evaluation] Reading ", length(input.truth), " simulated SV sets ...")
regex = "simulation_new/seed(\\d+)_size(\\d+)-(\\d+)_vaf(\\d+)-(\\d+)_([a-z_]+)/variants-(\\d+).txt"
Truth = NULL
num_truth = 0
for (f in input.truth) {
  vars = str_match(f, regex)
  t = fread(f)

  if (nrow(t) == 0) {
    message("No SV calls in ", t)
    next
  }

  num_truth = num_truth + 1
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
           SIMUL_svclass = vars[7],
           SIMUL_binsize = as.integer(vars[8]))]
  Truth = rbind(Truth, t)
}
message(num_truth, "/", length(input.truth), " Simulated SV call sets were not empty")



### Read SV calls
message("[Evaluation] Reading ", length(input.calls), " data sets ...")
regex = "sv_calls/seed(\\d+)_size(\\d+)-(\\d+)_vaf(\\d+)-(\\d+)_([a-z_]+)-(\\d+)/(\\d+)_fixed\\.fraction(\\d+)/([a-zA-Z0-9_-]+)\\.txt"
LOCI_SUMMARY = NULL

for (f in input.calls) {

  vars = str_match(f, regex)
  message(" * Reading ", f)
  d = fread(f)
  assert_that("chrom"    %in% colnames(d),
              "start"    %in% colnames(d),
              "end"      %in% colnames(d),
              "sample"   %in% colnames(d),
              "cell"     %in% colnames(d),
              "sv_call_name" %in% colnames(d)) %>% invisible

  # Get corresponding true SV calls
  truth = Truth[SIMUL_id      == as.integer(vars[2]) &
                SIMUL_minsize == as.integer(vars[3]) &
                SIMUL_maxsize == as.integer(vars[4]) & 
                SIMUL_minvaf  == as.integer(vars[5]) &
                SIMUL_maxvaf  == as.integer(vars[6]) &
                SIMUL_svclass == vars[7], ]

  if (nrow(truth)==0) {
    message("No simulated SV calls for call set. Ignore ", f)
    next
  }
  
  # Calculate precision and recall just based on the coordinates
  # (ignoring single cells)
  if (nrow(d) > 0) {
    rp = find_overlapping_calls(truth, d, min_bp_overlap)

    loci_summary = data.table(recall_base    = nrow(rp$true_loci),
                              recall         = nrow(rp$true_loci[covered > min_frac_covered]),
                              precision_base = nrow(rp$called_loci),
                              precision      = nrow(rp$called_loci[match_true_SVs==T]),
                              SIMUL_id        = as.integer(vars[2]),
                              SIMUL_minsize   = as.integer(vars[3]),
                              SIMUL_maxsize   = as.integer(vars[4]),
                              SIMUL_minvaf    = as.integer(vars[5]),
                              SIMUL_maxvaf    = as.integer(vars[6]),
                              SIMUL_svclass   = vars[7],
                              SIMUL_binsize   = as.integer(vars[8]),
                              SIMUL_fraction  = as.integer(vars[10]))
  } else {
    loci_summary = data.table(recall_base    = nrow(unique(truth[, .(chrom, start, end)])),
                              recall         = 0,
                              precision_base = 0,
                              precision      = 0,
                              SIMUL_id        = as.integer(vars[2]),
                              SIMUL_minsize   = as.integer(vars[3]),
                              SIMUL_maxsize   = as.integer(vars[4]),
                              SIMUL_minvaf    = as.integer(vars[5]),
                              SIMUL_maxvaf    = as.integer(vars[6]),
                              SIMUL_svclass   = vars[7],
                              SIMUL_binsize   = as.integer(vars[8]),
                              SIMUL_fraction  = as.integer(vars[10]))
  }
  LOCI_SUMMARY = rbind(LOCI_SUMMARY, loci_summary)
}




### Make labels a bit nicer.
format_Mb = function(x) {
  y = ifelse(x>=1e6, 
             paste(round(x/1e6,1),"Mb"),
             ifelse(x>=1e3,
                    paste(round(x/1e3,1),"kb"),
                    paste(x, "bp")))
  return (y)
}
LOCI_SUMMARY[, SV_size := factor(paste0(format_Mb(SIMUL_minsize),"-",format_Mb(SIMUL_maxsize)),
                        levels = unique(paste0(format_Mb(SIMUL_minsize),"-",format_Mb(SIMUL_maxsize)))[order(unique(SIMUL_minsize))],
                        ordered = T)]
LOCI_SUMMARY[, SV_vaf  := factor(paste0(SIMUL_minvaf,"-",SIMUL_maxvaf,"%"),
                        levels = unique(paste0(SIMUL_minvaf,"-",SIMUL_maxvaf,"%"))[order(unique(SIMUL_minvaf))],
                        ordered = T)]

### Summarize across all the different simulations
zzz = LOCI_SUMMARY[, .(.N,
                       num_SVs   = sum(recall_base),
                       num_calls = sum(precision_base),
                       recall    = sum(recall)/sum(recall_base),
                       precision = sum(precision)/sum(precision_base)),
                   by = .(SV_size, SV_vaf, SV_class = SIMUL_svclass, segmentation = SIMUL_fraction) ]
zzz = zzz[order(segmentation)]



write.table(zzz, file = snakemake@output[["table"]],
            quote=F, sep = "\t", row.names=F, col.names = T)

cairo_pdf(file = snakemake@output[["supplement"]], width=16, height = 12, onefile=T)
for (sv_class in unique(zzz$SV_class)) {

    p <- ggplot(zzz[SV_class == sv_class]) +
      geom_path(aes(recall, precision),
                linetype = "solid", color = "darkgrey") +
      geom_point(aes(recall, precision, col = segmentation)) +
      geom_text(x = 0, y = 0, hjust = 0, vjust = 0, aes(label = paste("SVs =",V1)),
                data = zzz[, num_SVs[1], by = .(SV_size, SV_vaf)]) +
      geom_text(x = 0, y = 0.1, hjust = 0, vjust = 0, aes(label = paste("Calls =",V1, "-", V2)),
                data = zzz[, .(min(num_calls),max(num_calls)), by = .(SV_size, SV_vaf)]) +
      facet_grid(SV_size ~ SV_vaf) +
      scale_color_gradientn(colours = c("darkgrey","firebrick", "gold","olivedrab2","dodgerblue3")) +
      theme_bw() +
      theme(legend.position = "bottom") +
      coord_cartesian(ylim = c(0,1), xlim = c(0,1)) +
      ggtitle(sv_class)

    print(p)
}
dev.off()



# More summarized version
num_segments = zzz[,.N,by = .(SV_size, SV_class, SV_vaf)][1,N]
plt <- ggplot(zzz) +
  geom_path(aes(recall, precision, col = SV_vaf), size = 0.8) +
  #geom_point(aes(recall,precision, col = SV_vaf), shape = 5) +
  facet_grid(SV_size ~ SV_class) +
  theme_bw() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim=c(0,1), xlim = c(0,1)) +
  scale_color_manual(values = grey.colors(num_segments, start = 0.8, end = 0),
                    name = "Clonal fraction of SVs") +
  ggtitle("Precision/recall curves for Strand-seq simulations")

cairo_pdf(file = snakemake@output[["figure"]], width=8, height = 6)
print(plt)
dev.off()



# test
plt <- ggplot(zzz[segmentation == 30]) +
  geom_line(aes(recall, precision, col = SV_vaf, group = SV_vaf)) + 
  geom_point(aes(recall, precision, col = SV_vaf, alpha = SV_size), size = 2) +
  facet_wrap(~ SV_class) +
  theme_bw() + 
  theme(legend.position = "bottom")
cairo_pdf(file = snakemake@output[["jan"]], width=8, height = 4)
print(plt)
dev.off()

