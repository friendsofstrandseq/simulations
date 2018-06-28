library(data.table)
library(ggplot2)


# Beatify data set prior to plotting
format_Mb = function(x) {
  y = ifelse(x>=1e6, 
             paste(round(x/1e6,1),"Mb"),
             ifelse(x>=1e3,
                    paste(round(x/1e3,1),"kb"),
                    paste(x, "bp")))
  return (y)
}



D = NULL
for (f in snakemake@input) {
  D = rbind(D, cbind(fread(f), method = sub(".pdf.txt","",basename(f))))
}
D

cairo_pdf(snakemake@output[[1]], width = 12, height = 8, onefile = T)
for (SV_class_ in unique(D$SV_class)) {
  for (SV_size_ in unique(D$SV_size)) {
    for (SV_vaf_ in unique(D$SV_vaf) ) {
      E = D[SV_size == SV_size_ & SV_class == SV_class_ & SV_vaf == SV_vaf_]
      p <- ggplot(E) +
        geom_path(aes(recall, precision, col = method)) +
        geom_point(aes(recall, precision, col = method)) +
        theme_bw() + 
        theme(legend.position = "bottom") +
        coord_cartesian(ylim = c(0,1), xlim = c(0,1)) + 
        ggtitle(paste0(SV_class_, "s with a size of ", SV_size_, " and a clonal fraction of ", SV_vaf_))
      print(p)
    }
  }
}
dev.off()
