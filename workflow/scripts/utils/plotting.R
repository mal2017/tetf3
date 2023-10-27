# func to summarize by fixed length windows
import_and_summarize <- function(x, tiles, which) {
  
  gr <- rtracklayer::import.bw(x,which=which)
  
  gr$score <- replace_na(gr$score, replace = 0)
  
  ba <- binnedAverage(bins=tiles,
                      numvar=mcolAsRleList(gr, varname = "score")[seqlevels(tiles)],
                      varname="score", na.rm=T)
}


# func to make plot of binned granges genome wide, uses ggbio. It returns an s3 object, one slot of which
# holds a plain ggplot2 / grob objec
plot_genome_signal <-  function(x) {
  g <- autoplot(x, aes(y=score), 
                geom="area", 
                chr.weight=c("2L"=1/6,"2R"=1/6,"3R"=1/6,"3L"=1/6,"4"=1/6,"X"=1/6)) +
    #theme_clear() +
    theme_genome() +
    facet_grid(grl_name~seqnames, scales="free",space = "free_x") +
    ylab("log2(IP/input)") +
    xlab("chromosomal position")
  
  return(g)
}
