#!/usr/bin/env Rscript
root <- normalizePath(if (length(commandArgs(trailingOnly=TRUE))) commandArgs(trailingOnly=TRUE)[1] else ".")
out <- file.path(root, "data", "burgold_analysis", "anchored_comparison")
d <- read.delim(file.path(out, "plasmid_comparison.tsv"), stringsAsFactors=FALSE)
draw <- function() {
  par(mfrow=c(2,2), mar=c(3.5,4.5,3,1), oma=c(4,0,4,0), family="sans")
  labels <- c(guide1="gRNA1", guide2="gRNA2", both="gRNA1 + gRNA2", swap="SWAP")
  colors <- c("#999999", "#1677B8", "#D55E00")
  for (metric in names(labels)) {
    r <- d[d$metric==metric,]
    values <- c(r$paper_percent[1], r$dnacomb_percent[r$policy=="baseline_control"], r$dnacomb_percent[r$policy=="agreement"])
    lim <- if (metric=="swap") 12 else 108
    bp <- barplot(values, col=colors, border=NA, ylim=c(0,lim), axes=FALSE,
                  ylab="Read pairs (%)", names.arg=c("Paper", "Original", "Anchored\n+ agreement"),
                  cex.names=0.9, main=labels[metric])
    axis(2, las=1, at=if(metric=="swap") seq(0,12,2) else seq(0,100,20))
    text(bp,values+lim*0.035,labels=sprintf("%.2f%%",values),cex=0.95)
    box(bty="l")
  }
  mtext("Burgold Figure 1F: corrected paired-read extraction", outer=TRUE, side=3, line=1.5, cex=1.25, font=2)
  mtext("Pilot plasmid library | 1,963,508 read pairs | Exact library matching", outer=TRUE, side=3, line=0, cex=0.9)
  mtext("Anchored extraction recognizes complete guides at read boundaries; agreement requires both mates to support the same guide.", outer=TRUE, side=1, line=1, cex=0.76)
  mtext("Both-guide definitions and scaffold handling retain small differences. Structural single-guide bars require separate validation.", outer=TRUE, side=1, line=2.4, cex=0.76)
}
png(file.path(out,"figure1f_corrected.png"),width=1800,height=1350,res=180)
draw(); dev.off()
pdf(file.path(out,"figure1f_corrected.pdf"),width=10,height=7.5)
draw(); dev.off()
