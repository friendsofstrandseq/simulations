library(data.table)
library(dplyr)
library(GenomicRanges)


rename_SV_classes <- function(calls) {
  if(nrow(calls)==0) {
      return(calls)
  }
  rename_svs = data.table(sv_call_name =         c("dup_h1",  "dup_h2",  "dup_hom", "del_h1",  "del_h2",  "del_hom", "inv_h1",  "inv_h2",  "inv_hom", "idup_h1", "idup_h2", "complex"),
                          sv_call_renamed = c("het_dup", "het_dup", "hom_dup", "het_del", "het_del", "hom_del", "het_inv", "het_inv", "hom_inv", "inv_dup", "inv_dup", "complex"))
  assert_that(all(calls$sv_call_name %in% rename_svs$sv_call_name)) %>% invisible
  calls = merge(calls, rename_svs, by = "sv_call_name")
  calls[, `:=`(sv_call_name = sv_call_renamed, sv_call_renamed = NULL)][]
  return(calls)
}




recall_precision <- function(truth, calls, rec_ovl = 0.8) {
  
  # Only allow if both tables are non-empty
  assert_that(is.data.table(truth),
              nrow(truth)>0,
              "chrom"   %in% colnames(truth),
              "start"   %in% colnames(truth),
              "end"     %in% colnames(truth),
              "SV_type" %in% colnames(truth)) %>% invisible
  assert_that(is.data.table(calls),
              # nrow(calls) > 0,
              "chrom"   %in% colnames(calls),
              "start"   %in% colnames(calls),
              "end"     %in% colnames(calls),
              "sv_call_name" %in% colnames(calls)) %>% invisible
  
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



  # Special case when there are no calls
  if (nrow(y)==0) {
    
      recall = x[, .(SV_size = (end - start + 1)[1],
                     SV_vaf = .N,
                     matches_call = FALSE,
                     correct_gt = 0,
                     correct_sv = 0),
             by = .(chrom, start, end, truth_id, SV_real = SV_type)]

      precision = data.table(chrom      = character(),
                             start      = integer(),
                             end        = integer(),
                             called_id  = integer(),
                             SV_size    = integer(),
                             SV_vaf     = integer(),
                             SV_found   = character(),
                             matches_SV = logical(),
                             correct_gt = numeric(),
                             correct_sv = numeric())
      return(list(recall = recall,
                  precision = precision,
                  sv_call_set_empty = TRUE))

  # Normal mode
  } else {

    # View centered on true SV calls (recall)
    recall = merge(x[, .(truth_id, sample, cell, chrom, start, end, called_id, SV_real = SV_type)],
                   y[, .(truth_id, sample, cell, SV_found = sv_call_name)],
                   by = c("truth_id", "sample", "cell"),
                   all.x = T) %>%
      .[, .(SV_size = (end - start + 1)[1],
            SV_vaf = .N,
            matches_call = ifelse(any(is.na(called_id)), FALSE, TRUE),
            correct_gt = sum(!is.na(SV_found) & SV_real == SV_found),
            correct_sv = sum(!is.na(SV_found) & substr(SV_real,5,nchar(SV_real)) == substr(SV_found,5,nchar(SV_found)))), 
        by = .(chrom, start, end, truth_id, SV_real)]


    # View centered on called SVs (precision)
    # Note that SV size and VAF are now defined based on the CALLED SVs !!! this can be slighlty counter-intuitive.
    precision = merge(y[, .(called_id, sample, cell, chrom, start, end, SV_found = sv_call_name)],
          x[, .(called_id, sample, cell, truth_id, SV_real = SV_type)],
          by = c("called_id", "sample", "cell"),
          all.x = T) %>%
      .[, .(SV_size    = (end - start + 1)[1],
            SV_vaf     = .N,
            SV_found   = names(table(SV_found))[which.max(table(SV_found))],
            matches_SV = ifelse(all(is.na(truth_id)), FALSE, TRUE),
            correct_gt = sum(!is.na(SV_real) & SV_real == SV_found),
            correct_sv = sum(!is.na(SV_real) & substr(SV_real,5,nchar(SV_real)) == substr(SV_found,5,nchar(SV_found)))),
        by = .(chrom, start, end, called_id)]
  }

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

  d[, SV_size_factor := factor(cut(end - start, c(0, 2.5e5, 5e5, 1e6, 300e6)),
                               levels = c("(0,2.5e+05]", "(2.5e+05,5e+05]", "(5e+05,1e+06]", "(1e+06,3e+08]"),
                              labels = c("<250 kb", "250-500 kb", "0.5-1 Mb", ">1 Mb"),
                              ordered = T)]

  vaf_boarders = c(0, 0.05, 0.1, 0.2, 0.5, 1)
  d[, SV_vaf_factor := factor(cut(SV_vaf/n_cells, vaf_boarders),
                              levels = c("(0,0.05]", "(0.05,0.1]", "(0.1,0.2]", "(0.2,0.5]", "(0.5,1]"),
                              labels = paste0(vaf_boarders[1:5]*100, "-", vaf_boarders[2:6]*100, "%"),
                              ordered = T)]
  d
}










# New precision / recall. Use just 1bp overlap
find_overlapping_calls <- function(truth, calls, min_bp_overlap = 1) {

  assert_that(is.data.table(truth),
              "chrom"   %in% colnames(truth),
              "start"   %in% colnames(truth),
              "end"     %in% colnames(truth)) %>% invisible
  assert_that(nrow(truth)>0) %>% invisible
  assert_that(is.data.table(calls),
              "chrom"   %in% colnames(calls),
              "start"   %in% colnames(calls),
              "end"     %in% colnames(calls)) %>% invisible


  # list of simulated loci (only chrom, start, end)
  true_loci = unique(truth[, .(chrom, start, end)])
  setkey(true_loci, chrom ,start, end)
  true_loci.gr = makeGRangesFromDataFrame(true_loci, seqinfo = c(paste0("chr",1:22),"chrX","chrY"))

  # list of predicted loci (chrom, start, end) --> only use real SV calls here!
  called_loci = unique(calls[, .(chrom, start, end)])
  setkey(called_loci, chrom, start, end)
  called_loci.gr = makeGRangesFromDataFrame(called_loci, seqinfo = c(paste0("chr",1:22),"chrX","chrY"))

  # Check that the GR is still sorted the same way
  assert_that(all(called_loci$start == start(called_loci.gr))) %>% invisible
  assert_that(all(true_loci$start == start(true_loci.gr))) %>% invisible

  # Annotate whether calls overlap true SVs
  called_loci[, match_true_SVs := overlapsAny(called_loci.gr, true_loci.gr, minoverlap = min_bp_overlap)]

  # Annotate whether true SV is covered by calls
  ovlp <- as.data.table(findOverlaps(called_loci.gr, true_loci.gr))
  ovlp <- ovlp[,
       .(covered = sum(width(intersect(called_loci.gr[queryHits], true_loci.gr[subjectHits]))) / width(true_loci.gr[subjectHits])),
       by = subjectHits]

  true_loci[, covered := 0.0]
  true_loci[ovlp$subjectHits, covered := ovlp$covered]

  return( list(true_loci   = true_loci,
               called_loci = called_loci) )
}
