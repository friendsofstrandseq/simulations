library(dplyr)
library(data.table)
library(assertthat)



################################################################################
### Function section for sv_classifier_merger.R
################################################################################


format_Mb = function(x) {
  y = ifelse(x>=1e6, 
             paste(round(x/1e6,1),"Mb"),
             ifelse(x>=1e3,
                    paste(round(x/1e3,1),"kb"),
                    paste(x, "bp")))
   return (y)
}


# normalize_probs
# ---------------
# transform probabilities from raw logged values into normalized (rowSum = 1)
# probabilities for the different SV states.
# First, normalize to sum up to 1, then add a pseudo_count (aka Dirichle prior)
# to all values, and normalize again.
unlog_and_normalize_probs <- function(d, pseudo_count = 1e-10) {
  
  d <- copy(d)
  assert_that(is.data.table(d),
              "p_ref" %in% colnames(d) || "p_ref" %in% SV_classes,
              all(d[, p_ref <= 0]))
  
  # These are the columns with (log) probabilities
  col_idx = grep('^p_', colnames(d))
  
  # Transform to non-log Probabilities
  for (j in col_idx) set(d, j = j, value = exp(d[[j]]))
  
  # Normalize to a rowSum of 1
  # Then add pseudo_count and normalize again
  norm = d[, rowSums(.SD), .SDcols = col_idx]
  for (j in col_idx) set(d, j = j, value = d[[j]] / norm + pseudo_count)
  
  norm = d[, rowSums(.SD), .SDcols = col_idx]
  for (j in col_idx) set(d, j = j, value = d[[j]] / norm)
  
  d
}

# calc_jump_probabilities_segment
# -------------------------------
# Calculate jump probabilities between two segments
# `SV_states`` defines the SV classes that are present as columns in `left` 
# and `right` --> They must contain "ref".
# Then, multiply the probabilities of all combinations of SV (or non-SV) states 
# between the two segments to get the joint probability of them being the same
# or not. For example:
# Pr(merge)     = Pr(ref left) * Pr(ref right)  + Pr(gain left) * Pr(gain right) + ...
# Pr(separate)  = Pr(ref left) * Pr(gain right) + Pr(ref left) * Pr(loss right)  + ...
calc_jump_probabilities_segment <- function(left, 
                                    right, 
                                    pseudo_count, 
                                    SV_classes) {
  
  # Check input:
  # SV_classes must start be "ref", "hetInv", a.s.o and there must be respective columns with a "p_" prefix
  # in both left and right.
  assert_that(is.character(SV_classes),
              "ref" %in% SV_classes || "p_ref" %in% SV_classes)
  
  if (!all(grepl("^p_", SV_classes))) {
    SV_classes = paste0("p_", SV_classes)
  }
  
  assert_that(is.data.table(left),
              "chrom"  %in% colnames(left),
              "seg_id" %in% colnames(left),
              "p_ref"  %in% SV_classes,
              all(SV_classes %in% colnames(left)))
  assert_that(is.data.table(right),
              "chrom"  %in% colnames(right),
              "seg_id" %in% colnames(right),
              "p_ref"  %in% SV_classes,
              all(SV_classes %in% colnames(right)))
  assert_that(nrow(left) == nrow(right))
  # assert that probabilities are between zero and one (not log!)
  assert_that(left[, all(.SD >= 0) && all(.SD <= 1), .SDcols = SV_classes])
  assert_that(right[, all(.SD >= 0) && all(.SD <= 1), .SDcols = SV_classes])
  
  
  # Calculate jump proabilities
  p_not_matching = rep(0,nrow(left))
  p_matching     = rep(0,nrow(left))
  for (j in SV_classes) {
    for (k in SV_classes) {
      if (j == k) {
        p_matching = p_matching + left[,j,with=F] * right[,k,with=F] 
      } else {
        p_not_matching = p_not_matching + left[,j,with=F] * right[,k,with=F]
      }
    }
  }
  
  # Look only at non-REF cells:
  MIN_LLR = 0.5
  select_alt_cells = left[, do.call(pmax, .SD) / p_ref > exp(MIN_LLR), .SDcols = SV_classes[!grepl("p_ref",SV_classes)]] | right[, do.call(pmax, .SD) / p_ref > exp(MIN_LLR), .SDcols = SV_classes[!grepl("p_ref",SV_classes)]]
  
  if (sum(select_alt_cells) > 1) {
    mean_pr_separate = mean(p_not_matching[[1]][select_alt_cells])
    mean_pr_merge    = mean(p_matching[[1]][select_alt_cells])
    mean_pr_n        = sum(select_alt_cells)
  } else {
    mean_pr_separate = 0
    mean_pr_merge    = 0
    mean_pr_n        = 0
  }
  
  return(list(pr_merge      = p_matching[[1]],
              pr_separate   = p_not_matching[[1]],
              mean_pr_separate = mean_pr_separate,
              mean_pr_merge = mean_pr_merge,
              mean_pr_n     = mean_pr_n))
}



# calc_jump_probabilities_byref
# -----------------------------
# For each segment in `prob` (except for the first segment in each chromosome)
# calculate the jump probabilities compared to the previous segment.
# `pr_merge` is the probability to merge the segment with its left neighbour per cell,
# `pr_separate` the opposite. These probabilities are based on the reduced SV states
# "ref", "gain", "loss", and "inv".
# Then, I accumulate them into one value which is the mean of these values among the
# cells that show evidence for an alternative allele --> leading to `mean_pr_separate`,
# `mean_pr_merge`, and `mean_pr_n`.
# For every neighboring pair of segments, `calc_jump_probabilities_segment` gets called
# internally.
calc_jump_probabilities_byref <- function(prob, SV_classes, quiet = F) {
  
  assert_that(is.character(SV_classes),
              "ref" %in% SV_classes)
  
  if (!all(grepl("^p_", SV_classes))) {
    SV_classes = paste0("p_", SV_classes)
  }
  
  assert_that(is.data.table(prob),
              "chrom"   %in% colnames(prob),
              "seg_id"  %in% colnames(prob),
              all(SV_classes %in% colnames(prob)))
  
  # check that values are between 0 and 1 (not log scale)
  assert_that(prob[, all(.SD >= 0) && all(.SD <= 1), .SDcols = SV_classes])
  
  prob[, `:=`(pr_merge         = as.numeric(NA),
              pr_separate      = as.numeric(NA),
              mean_pr_separate = as.numeric(NA),
              mean_pr_merge    = as.numeric(NA),
              mean_pr_n        = as.integer(NA))]
  for (chrom_ in unique(prob$chrom)) {
    
    if (!quiet) message("[SV merger] Calculating jump proabilities for ", chrom_)
    
    # first segment gets a default prob (NA).
    # for every other segment, I calculate the jump/stay probabilities to the previous segment
    for (seg_id_ in 2:prob[chrom == chrom_]$seg_id_max[1]) {
      
      jump_probs = calc_jump_probabilities_segment(
            prob[chrom == chrom_ & seg_id == (seg_id_ - 1)], 
            prob[chrom == chrom_ & seg_id == seg_id_], 
            pseudo_count,
            SV_classes )
  
      set(prob, 
          i = which(prob[, chrom == chrom_ & seg_id == seg_id_]), 
          j = "pr_merge",
          jump_probs$pr_merge)
      set(prob, 
          i = which(prob[, chrom == chrom_ & seg_id == seg_id_]), 
          j = "pr_separate",
          jump_probs$pr_separate)
      set(prob,
          i = which(prob[, chrom == chrom_ & seg_id == seg_id_]), 
          j = "mean_pr_separate",
          jump_probs$mean_pr_separate)
      set(prob,
          i = which(prob[, chrom == chrom_ & seg_id == seg_id_]), 
          j = "mean_pr_merge",
          jump_probs$mean_pr_merge)
      set(prob,
          i = which(prob[, chrom == chrom_ & seg_id == seg_id_]), 
          j = "mean_pr_n",
          jump_probs$mean_pr_n)
    }
  }
}





# Inspect jump probabilities in two neighboring segments
plot_change_probs_for_single_segment <- function(prob, chrom_, seg_id_, SV_classes) {
  
  assert_that(is.character(SV_classes),
              "ref" %in% SV_classes)
  
  if (!all(grepl("^p_", SV_classes))) {
    SV_classes = paste0("p_", SV_classes)
  }
  
  assert_that(is.data.table(prob),
              "chrom" %in% colnames(prob),
              "from" %in% colnames(prob),
              "to" %in% colnames(prob),
              "seg_id" %in% colnames(prob),
              "sample" %in% colnames(prob),
              "cell" %in% colnames(prob),
              "p_ref" %in% colnames(prob),
              "pr_merge" %in% colnames(prob),
              "pr_separate" %in% colnames(prob),
              all(SV_classes %in% colnames(prob)))
  
  left  = prob[chrom == chrom_ & seg_id == seg_id_ - 1]
  right = prob[chrom == chrom_ & seg_id == seg_id_]
  
  assert_that(all(left$cell == right$cell))
  
  xxx = melt(rbind(left,right), 
             id.vars = c("sample","cell","chrom","start","end"), 
             measure.vars = SV_classes)
  xxx[, `:=`(segment = paste0(chrom, ": ", format_Mb(start), "-", format_Mb(end)),
             chrom = NULL,
             start = NULL,
             end  = NULL)]
  xxx = rbind(xxx, data.table(sample = right$sample,
                              cell   = right$cell,
                              segment = "Change probabilities",
                              variable = "pr_merge",
                              value  = right$pr_merge))
  xxx = rbind(xxx, data.table(sample = right$sample,
                              cell   = right$cell,
                              segment = "Change probabilities",
                              variable = "pr_separate",
                              value  = right$pr_separate))
  
  ggplot(xxx) + 
    aes(cell,value,fill=variable) + 
    geom_bar(stat = "identity") + 
    facet_grid(segment ~ .) + 
    theme(legend.position = "bottom", axis.text.x = element_text(angle=60,hjust=1)) + 
    scale_fill_brewer(palette = 8, type = "qual") +
    ggtitle(paste(xxx$sample[1],"-",nrow(xxx), "cells"))
}

