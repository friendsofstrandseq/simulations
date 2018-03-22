#
# Copyright (C) 2017 Sascha Meiers
# Distributed under the MIT software license, see the accompanying
# file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.

library(data.table)
library(ggplot2)
library(scales)
library(assertthat)


zcat_command = "zcat"
format_Mb   <- function(x) {paste(comma(x/1e6), "Mb")}
cells_per_page <- 8


f_counts   = "counts/simulation5-50000/50000_fixed.txt.gz"
f_segments = "segmentation2/..."
f_sv_calls = "sv_calls/..."
f_out      = "sv_calls/..."

f_counts   = snakemake@input[["counts"]]
f_segments = snakemake@input[["segs"]]
f_sv_calls = snakemake@input[["svs"]]
f_out      = snakemake@output[[1]]



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
            "c"      %in% colnames(counts))
counts[, sample_cell := paste(sample, "-", cell)]
setkey(counts, sample_cell, chrom)



### Read segment file
seg = fread(f_segments)
assert_that("chrom"      %in% colnames(seg),
            "breakpoint" %in% colnames(seg),
            "start"      %in% colnames(seg),
            "end"        %in% colnames(seg))

### Read SV call table
svs = fread(f_sv_calls)
assert_that("chrom"  %in% colnames(svs),
            "start"  %in% colnames(svs),
            "end"    %in% colnames(svs),
            "sample" %in% colnames(svs),
            "cell"   %in% colnames(svs),
            "sv"     %in% colnames(svs))




# @work!

background_colors = data.table()
manual_colors = c(ref = "#bbbbbb",
                  del_hom = "dodgerblue4", del_h1 = "dodgerblue1", del_h2 = "dodgerblue3",
                  inv_hom = "seagreen",    inv_h1 = "seagreen1",   inv_h2 = "seagreen3",
                  dup_hom = "red4",        dup_h1 = "red1",        dup_h2 = "red3",
                  idup_h1 = "yellow1",    idup_h2 = "yellow3",
                  seg1 = "#edf8b1", seg2 = "#7fcdbb", seg3 = "#2c7fb8")




    # Read 'other' table - only then we can find out which type of data it is
    message(" * Reading addtional data ", f_other, "...")
    if (grepl('\\.gz$',f_other)) {
        other = fread(paste(zcat_command, f_other))
    } else {
        other = fread(f_other)
    }
    assert_that(is.data.table(other),
                "chrom" %in% colnames(other),
                "start" %in% colnames(other),
                "end"   %in% colnames(other))



    ##################
    # Switch data type
    #
    ### SV prob file
    if (all(c("sample","cell","p_ref") %in% colnames(other))) {

            if (length(splitted)>1) {
                message("You specified an SV prob file, but additional '='. This confuses me.")
                stop()
            }

            message("   -> detected SV probability file")
            prob = other

            prob[, sample_cell := paste(sample, "-", cell)]

            # Remove CN columns
            if (length(grep("^p_cn", colnames(prob)))>0) {
                prob[, grep("^p_cn", colnames(prob)) := NULL]
            }

            # Add prior probabilities: this is not the correct way to do it! [TODO]
            prob$p_ref = prob$p_ref * 1.1

            # Check that cells are the same as in "counts"
            c1 = unique(counts[,.(sample_cell)])
            c2 = unique(prob[,.(sample_cell)])
            c3 = c1[c2,, on = .(sample_cell), nomatch=0]

            if(nrow(c1) != nrow(c3)) {
                message("WARNING: Not all the cells from the count table are also covered in the probabilities")
                message("Cells in the count table:")
                print(paste(c1$sample_cell))
                message("Cells in the probability table:")
                print(paste(c2$sample_cell))
                message("Cells of the count table that are also covered by probs:")
                print(paste(c3$sample_cell))
            }

            #######################
            # Reformat prob as long
            background_colors <- melt(prob,
                              id.vars = c("chrom","start","end","sample_cell"),
                              measure.vars = colnames(prob)[grepl('^p_',colnames(prob))],
                              variable.name = "SV_class",
                              value.name    = "probability")

            background_colors[, SV_class := substr(SV_class,3,nchar(as.character(SV_class)))]


            # Select
            background_colors <- background_colors[, .SD[probability == max(probability),][order(probability),][1,], by = .(chrom,start,end,sample_cell)]
            setkey(background_colors, chrom, sample_cell)



    ### Seg file
    } else if (all(c("k","breakpoint") %in% colnames(other))) {


            message("   -> detected segmentation file")
            segments = other[, .SD[k == quantile(1:max(k), f_seg_quantile, type = 3)], by = chrom]
            assert_that(all(segments[, .N == k[1], by = chrom]$V1))
            segments[, SV_class := rep(paste0("seg",1:3), ceiling(.N/3))[1:.N], by = chrom]


            background_colors <- NULL
            for ( samcell in unique(counts[,sample_cell]))
                background_colors <- rbind(background_colors, cbind(segments, sample_cell = samcell))


    ### other, unknown input
    } else {
            message("I don't understand the file ", f_other)
            message("Segmentation files need to contain the columns k, breakpoint, chorm, start, end.")
            message("SV prob files need to have columns chrom, start, end, sample, cell and then columns")
            message("for the different SV probabilities, called p_ref, p_inv_het, p_del_hom, ...")
            print_usage()
            stop()
    }
} # !is.null(f_other)







##########################
# Plot, one file per chrom
y_lim = 3 * counts[,median(w+c)]
n_cells = length(unique(counts[,sample_cell]))

for (CHROM in unique(counts[, chrom])) {

    out = paste0(f_out,".",CHROM,".pdf")
    message(" * Plotting ", CHROM, " (", out, ")")
    if (nrow(background_colors)>0 && "seg1" %in% background_colors$SV_class)
        message("  -> Chose ", background_colors[chrom == CHROM, k][1], " segments (", f_seg_quantile, " quantile)")


    cairo_pdf(out, width=14, height=10, onefile = T)

    i = 1
    while (i <= n_cells) {


        # Subset to this set of cells:
        CELLS = unique(counts[,.(sample_cell)])[i:(min(i+cells_per_page-1,n_cells))]
        setkey(CELLS, sample_cell)
        local_counts = counts[CELLS, on = .(sample_cell), nomatch = 0][chrom == CHROM]

        # Also subset the background coloring, if available at all
        local_background_colors = data.table()
        if (nrow(background_colors)>0) {
            local_background_colors = background_colors[CELLS, on = .(sample_cell), nomatch=0][chrom == CHROM,]
        }

        plt <- ggplot(local_counts)

        # Add GT colors:
        if(nrow(local_background_colors)>0) {
            plt <- plt +
                geom_rect(data = local_background_colors, alpha = 0.5, size = 0.1,
                          aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = SV_class)) +
                scale_fill_manual(values = manual_colors)
        }

        plt <- plt +
            geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax = -w), fill='sandybrown') +
            geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax =  c), fill='paleturquoise4') +
            facet_wrap(~ sample_cell, ncol = 1) +
            ylab("Watson | Crick") + xlab(NULL) +
            scale_x_continuous(breaks = pretty_breaks(12), labels = format_Mb) +
            scale_y_continuous(breaks = pretty_breaks(3)) +
            coord_cartesian(ylim = c(-y_lim, y_lim)) +
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

