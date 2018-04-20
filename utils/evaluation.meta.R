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

cairo_pdf(snakemke@output[[1]], width = 12, height = 8, onefile = T)
for (sv_size_ in unique(D$SIMUL_minsize)) {
  for (vaf_ in unique(D$SIMUL_minvaf) ) {
    E = D[SIMUL_minsize == sv_size_ & SIMUL_minvaf == vaf_]
    size_string = paste(format_Mb(unique(E$SIMUL_minsize)), "-", format_Mb(unique(E$SIMUL_maxsize)))
    vaf_string  = paste0(unique(E$SIMUL_minvaf),"-",unique(E$SIMUL_maxvaf), " %")
    p <- ggplot(E) +
      geom_path(aes(bp.recall, bp.precision, col = method)) +
      geom_point(aes(bp.recall, bp.precision, col = method)) +
      theme_bw() + 
      theme(legend.position = "bottom") +
      coord_cartesian(ylim = c(0,1), xlim = c(0,1)) + 
      ggtitle(paste0("SV size ", size_string, " with a VAF of ", vaf_string))
    print(p)
  }
}
dev.off()