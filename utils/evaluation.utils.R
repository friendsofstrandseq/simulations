library(data.table)
library(dplyr)

find_overlapping_calls <- function(truth, calls, rec_ovl = 0.8) {
  
  assert_that(is.data.table(truth),
              "chrom"   %in% colnames(truth),
              "start"   %in% colnames(truth),
              "end"     %in% colnames(truth),
              "SV_type" %in% colnames(truth)) %>% invisible
  assert_that(is.data.table(calls),
              "chrom"   %in% colnames(calls),
              "start"   %in% colnames(calls),
              "end"     %in% colnames(calls),
              "SV_class" %in% colnames(calls)) %>% invisible
  
  
  
  ### 1) Get unique list of loci
  
  # list of simulated loci (only chrom, start, end)
  true_loci = unique(truth[, .(chrom, start, end)])
  setkey(true_loci, chrom ,start, end)
  true_loci[, truth_id := 1:.N]
  
  # list of predicted loci (chrom, start, end) --> only use real SV calls here!
  called_loci = unique(calls[, .(chrom, start, end)])
  setkey(called_loci, chrom, start, end)
  called_loci[, called_id := 1:.N]
  
  
  
  ### 2) Find overlap by generating all combinations of loci
  
  # Merge all combinations of loci by chromosoems ...
  combined = merge(true_loci, called_loci,
                   by = "chrom",
                   all = T,
                   suffixes = c("",".call"),
                   allow.cartesian = T)
  
  # Mark those entries that overlap by more than xx %
  combined$match = FALSE
  combined[!is.na(start) & end.call > start & start.call < end & (pmin(end, end.call) - pmax(start, start.call)) / (pmax(end, end.call) - pmin(start, start.call)) >= rec_ovl,
           match := TRUE]
  
  # Moreover, mark those entries that at least touch SVs
  combined[, overlap := rep(0, .N)]
  combined[!is.na(start) & end.call > start & start.call < end,
           overlap := (pmin(end, end.call) - pmax(start, start.call))/(end.call - start.call)]
  assert_that(all(combined[match == TRUE, overlap >= rec_ovl])) %>% invisible
  
  
  
  ### 3) Annotate overlapping SVs
  
  # First, check that matches are OK
  assert_that(all(unique(combined[, .(chrom, start = start.call, end = end.call)])[!is.na(start), ] == called_loci[, .(chrom, start, end)]))  %>% invisible
  assert_that(all(unique(combined[, .(chrom, start, end)])[!is.na(start)] == true_loci[, .(chrom, start, end)])) %>% invisible
  
  # Next, annotate "true_loci" with "called_id"
  mult_ovl_warn = function(chrom, start, end, SD) {
      if (nrow(SD) > 1) {
          message("[Evaluation] Warning: ", nrow(SD), 
                  " overlapping SV calls were found for SV ", chrom, ":", 
                  format(start, big.mark = ","), "-", 
                  format(end, big.mark = ","), 
                  ". I chose only one of them")
      }
  }
  combined[match == TRUE,         # Warning of more than 1 call overlap same SV
           mult_ovl_warn(chrom, start, end, .SD), 
           by = .(chrom, start, end)] %>% invisible
  true_loci <- merge(true_loci,
                     combined[match == TRUE,
                        .(called_id = called_id[1]),
                        by = .(chrom, start, end)],
                     by = c("chrom", "start", "end"),
                     all.x = T)
  
  # Then, annotate "called_loci" with "truth_id"
  mult_ovl_warn = function(chrom, start, end, SD) {
    if (nrow(SD) > 1) {
      message("[Evaluation] Warning: ", nrow(SD), 
              " true SVs were found to overlap the predicted SV ", 
              chrom, ":", format(start, big.mark = ","), "-", 
              format(end, big.mark = ","), ". I chose only one of them")
    }
  }
  combined[match == TRUE, 
           mult_ovl_warn(chrom, start.call, end.call, .SD), 
           by = .(chrom, start.call, end.call)] %>% invisible
  called_loci <- merge(called_loci,
                       combined[match == TRUE, 
                                .(truth_id = truth_id[1]),
                                by = .(chrom, start.call, end.call)] %>%
                                .[,.(chrom, start = start.call, end = end.call, truth_id)],
                       by = c("chrom", "start", "end"),
                       all.x = T)
  
  # At last, also annotate "called_loci" with whether they partly overlap an SV
  called_loci <- merge(called_loci,
                       combined[,.(max_overlap = max(overlap)), by = .(chrom, start = start.call, end = end.call)],
                       by = c("chrom", "start", "end"),
                       all.x = T)
                       
  
  
  # Return these loci lists
  return(list(true_loci = true_loci,
              called_loci = called_loci))
}


rename_SV_classes <- function(calls) {
  rename_svs = data.table(SV_class =         c("dup_h1",  "dup_h2",  "dup_hom", "del_h1",  "del_h2",  "del_hom", "inv_h1",  "inv_h2",  "inv_hom", "idup_h1", "idup_h2"),
                          SV_class_renamed = c("het_dup", "het_dup", "hom_dup", "het_del", "het_del", "hom_del", "het_inv", "het_inv", "hom_inv", "inv_dup", "inv_dup"))
  assert_that(all(calls$SV_class %in% rename_svs$SV_class)) %>% invisible
  calls = merge(calls, rename_svs, by = "SV_class")
  calls[, `:=`(SV_class = SV_class_renamed, SV_class_renamed = NULL)][]
  return(calls)
}


recall_precision <- function(truth, calls, rec_ovl = 0.8) {
  
  # Find overlapping SV calls
  L = find_overlapping_calls(truth, calls, rec_ovl)
  
  # Prepare table of True SVs
  x = merge(truth, 
            L[["true_loci"]], 
            by = c("chrom", "start", "end"))
  
  # Prepare table of called SVs - this requires renaming of SV_classes
  y = merge(calls,  
            L[["called_loci"]], 
            by = c("chrom", "start", "end"))
  y = rename_SV_classes(y)
  
  
  # View centered on true SV calls (recall)
  recall = merge(x[, .(truth_id, sample, cell, chrom, start, end, called_id, SV_real = SV_type)],
                 y[, .(truth_id, sample, cell, SV_found = SV_class)],
                 by = c("truth_id", "sample", "cell"),
                 all.x = T) %>%
    .[, .(SV_size = (end - start + 1)[1],
          SV_vaf = .N,
          breakpoint_found = ifelse(any(is.na(called_id)), FALSE, TRUE),
          correct_gt = sum(!is.na(SV_found) & SV_real == SV_found),
          correct_sv = sum(!is.na(SV_found) & substr(SV_real,5,nchar(SV_real)) == substr(SV_found,5,nchar(SV_found)))), 
      by = .(chrom, start, end, truth_id, SV_real)]
  
  
  # View centered on called SVs (precision)
  precision = 3
  y[is.na(truth_id)] %>% .[,.(max_overlap[1]),.(chrom,start,end)] %>% unique
  
  message("@work: precision")
  
  return(list(recall = recall,
              precision = precision))
}



# Categorize SVs based on SV size ("start", "end") and VAF "SV_vaf"
categorization <- function(d, n_cells = 100) {
  assert_that(is.data.table(d),
              "chrom"   %in% colnames(d),
              "start"   %in% colnames(d),
              "end"     %in% colnames(d),
              "SV_vaf"  %in% colnames(d),
              "chrom" %in% colnames(d),
              "chrom" %in% colnames(d)) %>% invisible
  assert_that(n_cells >= max(d$SV_vaf)) %>% invisible
  d[, SV_size_factor := factor(cut(d$end - d$start, c(1e5, 2e5, 4e5, 8e5, 1.6e6, 3.2e6, 250e6)),
                              labels = c("100-200 kb", "200-400 kb", "400-800 kb", "0.8-1.6 Mb", "1.6-3.2 Mb",">3.2Mb"),
                              ordered = T)]
  vaf_boarders = c(0,0.05,0.1,0.2,0.4,0.75,1)
  d[, SV_vaf_factor := factor(cut(SV_vaf/n_cells, vaf_boarders),
                              labels = paste0(vaf_boarders[1:6]*100, "-", vaf_boarders[2:7]*100, "%"),
                              ordered = T)]
  d
}
