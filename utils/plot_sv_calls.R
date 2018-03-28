#
# Copyright (C) 2017 Sascha Meiers
# Distributed under the MIT software license, see the accompanying
# file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.

library(dplyr)
library(data.table)
library(ggplot2)
library(scales)
library(assertthat)


zcat_command = "zcat"
format_Mb   <- function(x) {paste(comma(x/1e6), "Mb")}
cells_per_page <- 8


# f_counts   = "counts/simulation5-50000/50000_fixed.txt.gz"
# f_segments = "segmentation2/simulation5-50000/50000_fixed.few.txt"
# f_sv_calls = "sv_calls/simulation5-50000/50000_fixed.few/maryam.txt"
# f_sv_simul = "simulation/variants/genome5-50000.txt"
# f_out      = "sv_calls/simulation5-50000/50000_fixed.few/maryam"

f_counts   = snakemake@input[["counts"]]
f_segments = snakemake@input[["segs"]]
f_sv_calls = snakemake@input[["svs"]]
f_sv_simul = snakemake@input[["simul"]]
f_out      = snakemake@output[[1]]
f_out      = sub('\\.[chr0-9_]+\\.pdf$','', f_out)




### Colors for background
manual_colors = c(  # duplications
                  simul_hom_dup   = "firebrick4",
                  dup_hom         = muted("firebrick4", 70, 50),
                  simul_het_dup   = "firebrick2",
                  dup_h1          = muted("firebrick2", 90, 30),
                  dup_h2          = muted("firebrick2", 90, 30),
                    # deletions
                  simul_hom_del   = "dodgerblue4",
                  del_hom         = muted("dodgerblue4", 50, 60),
                  simul_het_del   = "dodgerblue2",
                  del_h1          = muted("dodgerblue2", 80, 50),
                  del_h2          = muted("dodgerblue2", 80, 50),
                    # inversions
                  simul_hom_inv   = "chartreuse4",
                  inv_hom         = muted("chartreuse4", 80, 50),
                  simul_het_inv   = "chartreuse2",
                  inv_h1          = muted("chartreuse2", 100, 60),
                  inv_h2          = muted("chartreuse2", 100, 60),
                    # other SVs
                  simul_false_del = "darkgrey",
                  simul_inv_dup   = "darkgoldenrod2",
                  idup_h1         = muted("darkgoldenrod2", 80, 70),
                  idup_h2         = muted("darkgoldenrod2", 80, 70),
                    # background
                  bg1 = "#ffffff",
                  bg2 = "khaki2")





#############
# Read counts
message(" * Reading count data ", f_counts, "...")
if (grepl('\\.gz$',f_counts)) {
    counts = fread(paste(zcat_command, f_counts))
} else {
    counts = fread(f_counts)
}

assert_that("chrom"  %in% colnames(counts),
            "start"  %in% colnames(counts),
            "end"    %in% colnames(counts),
            "class"  %in% colnames(counts),
            "sample" %in% colnames(counts),
            "cell"   %in% colnames(counts),
            "w"      %in% colnames(counts),
            "c"      %in% colnames(counts)) %>% invisible
counts[, sample_cell := paste(sample, "-", cell)]
setkey(counts, sample_cell, chrom)



### Read SV call table
message(" * Reading SV calls from ", f_sv_calls, "...")
svs = fread(f_sv_calls)
assert_that("chrom"    %in% colnames(svs),
            "start"    %in% colnames(svs),
            "end"      %in% colnames(svs),
            "sample"   %in% colnames(svs),
            "cell"     %in% colnames(svs),
            "SV_class" %in% colnames(svs)) %>% invisible
assert_that(all(svs$SV_class %in% names(manual_colors))) %>% invisible
svs[, sample_cell := paste(sample, "-", cell)]
assert_that(all(unique(svs$sample_cell) %in% unique(counts$sample_cell))) %>% invisible



### Read simulated variants
message(" * Reading simulated variants from ", f_sv_simul, "...")
simul = fread(f_sv_simul)
assert_that("chrom"   %in% colnames(simul),
            "start"   %in% colnames(simul),
            "end"     %in% colnames(simul),
            "sample"  %in% colnames(simul),
            "cell"    %in% colnames(simul),
            "SV_type" %in% colnames(simul)) %>% invisible
simul[, `:=`(SV_class = paste0("simul_",SV_type), SV_type = NULL, sample_cell = paste(sample, "-", cell))]



### Read segment file
message(" * Reading segmentation file from ", f_segments, "...")
seg = fread(f_segments)
assert_that("chrom" %in% colnames(seg),
            "bps"   %in% colnames(seg)) %>% invisible
bins = unique(counts[, .(chrom, start, end)])
seg = merge(seg, bins[,.N, by =chrom][, .(chrom, N = c(0,cumsum(N)[1:(.N-1)]))], by = "chrom")
seg[, `:=`(from = c(1,bps[1:(.N-1)]+1), to = bps), by = chrom]
seg[, `:=`(start = bins[from + N]$start,
           end   = bins[to   + N]$end  )]
seg[, SV_class := rep(c("bg1","bg2"),.N)[1:.N], by = chrom]






##########################
# Plot, one file per chrom
y_lim = 3 * counts[,median(w+c)]
n_cells = length(unique(counts[,sample_cell]))

for (CHROM in unique(counts[, chrom])) {

    out = paste0(f_out,".",CHROM,".pdf")
    message(" * Plotting ", CHROM, " (", out, ")")

    cairo_pdf(out, width=14, height=10, onefile = T)

    # Plot only a couple of cells per page
    i = 1
    while (i <= n_cells) {

        # Subset to this set of cells:
        CELLS = unique(counts[,.(sample_cell)])[i:(min(i+cells_per_page-1,n_cells))]
        setkey(CELLS, sample_cell)
        local_counts = counts[CELLS, on = .(sample_cell), nomatch = 0][chrom == CHROM]

        local_svs = svs[CELLS, on = .(sample_cell), nomatch = 0][chrom == CHROM]
        local_sim = simul[CELLS, on = .(sample_cell), nomatch = 0][chrom == CHROM]

        # Segments need to be multiplied by "CELLS"
        local_seg = seg[chrom == CHROM]
        local_seg = CELLS[, cbind(local_seg, sample_cell), by = sample_cell]


        # Start major plot
        plt <- ggplot(local_counts)

        # Add background colors:
        if (nrow(local_seg)>0) {
            plt <- plt +
                geom_rect(data = local_seg, alpha = 0.4,
                          aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = SV_class)) 
        }
        if (nrow(local_svs)>0) {
          plt <- plt +
            geom_rect(data = local_svs, alpha = 1,
                      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = SV_class))
        }
        if (nrow(local_sim)>0) {
          plt <- plt +
            geom_rect(data = local_sim,
                      aes(xmin = start, xmax = end, ymin = y_lim, ymax = Inf, fill = SV_class))
        }

        plt <- plt +
            geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax = -w), fill='sandybrown') +
            geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax =  c), fill='paleturquoise4') +
            facet_wrap(~ sample_cell, ncol = 1) +
            ylab("Watson | Crick") + xlab(NULL) +
            scale_x_continuous(breaks = pretty_breaks(12), labels = format_Mb) +
            scale_y_continuous(breaks = pretty_breaks(3)) +
            coord_cartesian(ylim = c(-y_lim, y_lim)) +
            scale_fill_manual(values = manual_colors) +
            theme_minimal() +
            theme(panel.spacing = unit(0, "lines"),
                  axis.ticks.x = element_blank(),
                  strip.background = element_rect(color = "#eeeeee", fill = "#eeeeee"),
                  strip.text = element_text(size = 5),
                  legend.position = "bottom") +
            ggtitle(paste("data:", basename(f_counts), "chromosome:", CHROM))

        print(plt)
        i = i + cells_per_page
    } # while
    dev.off()
} # for CHROM

