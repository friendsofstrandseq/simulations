library(data.table)
d = fread(snakemake@input[[1]])
minsize = as.integer(snakemake@params[["minsvsize"]])
mincells = as.integer(snakemake@params[["minvaf"]]) / 100 * as.integer(snakemake@params[["numcells"]])

d = merge(d, d[, .N, by = .(chrom, start, end)], by = c("chrom", "start", "end"))
write.table(d[end-start >= minsize & N >= mincells],
            file = snakemake@output[[1]],
            quote=F, row.names = F, sep="\t")
