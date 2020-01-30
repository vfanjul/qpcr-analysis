time0 = proc.time()

project = "Heart progerin" # mir29   Heart progerin   Heart lamin  ISO Challenge  PCTX Treatment  Zmpste-Rankl  HGPS Amanda  DBU Alberto

## Import data
if (.Platform$OS.type == "unix") setwd("/Volumes/Victor/") else setwd("S:/LAB_VA/LAB/Victor/")
baseroute = paste0(project, " project/")

source("Methodology/utility.r")
# source("S:/LAB_VA/LAB/Victor/Methodology/utility.r")
# source("/Volumes/Victor/Methodology/utility.r")

rconfig = data.frame(fread(paste0(baseroute, "Design/rconfig.txt"), encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "")))
rconfig = rconfig[rconfig$Run == 1 & rconfig$technique == "qPCR",-1]
invisible(sapply(names(rconfig), function (x) if (class(rconfig[,x]) == "integer") rconfig[,x] <<- as.logical(rconfig[,x])))

exp = 3
for (exp in 1:nrow(rconfig)) {
  for(opt in 1:ncol(rconfig)) assign(names(rconfig)[opt], rconfig[exp,opt])
  tryCatch({
    route = paste0(baseroute, "Results/", study, "/", technique)    
    expdata = data.frame(fread(paste0(route, "/idmeanexpdata.xls"), encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "", "nan")))
    if (!onlymeanplots) {
      allexpdata = data.frame(fread(paste0(route, "/expdata.xls"), encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "", "nan")))
      expdata = allexpdata[,names(expdata)]
    }
    expstats = data.frame(fread(paste0(route, "/expstats.xls"), encoding = "Latin-1", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "", "nan")))
    expstats$Condition = gsub("vs[[:print:]]*", "", expstats$Comparison)
    expstats$Gene = gsub(" relative expression", "", expstats$Parameter)

    genecols = grep("relative.expression", names(expdata), value = T)
    plotdata = c()
    for (i in genecols) {
      meanctrl = mean(expdata[expdata$Condition.Id == 1,i], na.rm = T)
      if (is.na(meanctrl) | meanctrl == 0) {
        genecols = genecols[-which(genecols == i)]
      } else plotdata = rbind.data.frame(plotdata, cbind.data.frame(expdata[!is.na(expdata[,i]),1:11],
                                                                    Gene = gsub("\\.relative.expression", "", i), 
                                                                    Relative.expression = expdata[!is.na(expdata[,i]),i]/meanctrl, 
                                                                    stringsAsFactors = F), stringsAsFactors = F)
    }
    if (length(plotdata$Color) == 0) plotdata$Color = plotdata$Condition.Id
    if (study %in%  c("HMosloiWATdiff", "HMoslopWATdiff")) plotdata = plotdata[plotdata$Gene %in% c("Adipoq", "Ccl2", "Cebpb", "Fabp4", "Il6", "Pgc1a", "Pparg"),]
    if (study == "HMpWAT16") plotdata = plotdata[plotdata$Gene != "Ucp1",]
    expstats = expstats[expstats$Condition %in% unique(plotdata$Condition) & expstats$Gene %in% unique(plotdata$Gene),c("Gene", "Condition", "p.value", "Sig")]
    if (onlyplotsigcomps) expstats$Sig[expstats$p.value > 0.05] = ""

    genes = sort(unique(plotdata$Gene))
    positions = cumsum(rep(c(1.5, 1), length(genes))) - 1.5
    # xlabels = parse(text = c(paste0("atop(italic(", genes, "))")))
    xlabels1 = parse(text = c(paste0("atop(italic(", genes[seq(1, by = 2, length.out = length(genes)/2)], "))")))
    xlabels2 = parse(text = c(paste0("atop(,italic(", genes[seq(2, by = 2, length.out = length(genes)/2)], "))")))
    statpos = c()
    for (i in genes) statpos = c(statpos, max(log(plotdata$Relative.expression[plotdata$Gene == i & plotdata$Condition.Id == 2])) + (max(log(plotdata$Relative.expression)) - min(log(plotdata$Relative.expression)))*0.1)
    statlabels = merge(data.frame(Gene = genes), expstats)$Sig
    
    if (study %in%  c("HMosloiWATdiff", "HMoslopWATdiff")) {
      pdf(paste0(route, "/Relative.expression.pdf"), 1.6 + 1.7*round(length(genes)/2)*1.5, 6, pointsize = 24, useDingbats = F)
    } else pdf(paste0(route, "/Relative.expression.pdf"), 1.6 + 1.7*round(length(genes)/2), 6, pointsize = 24, useDingbats = F)
    par(bty = "l", cex.axis = 0.75, mar = c(2,3,2,1), mgp = c(1.5,0.25,0), tck = - 0.02, xaxs = "i")
    boxplot(log(Relative.expression) ~ Condition.Id + Gene, plotdata, pch = 20, xaxt = "n", xlab ="", ylab = "Log relative mRNA expression", border = 0, at = positions)
    # for (i in seq(0, by = 5, length.out = trunc(length(genes)/2))) rect(2.75 + i,log(min(plotdata$Relative.expression))*1.1, 5.25 + i,log(max(plotdata$Relative.expression))*1.1, col = adjustcolor(1, alpha.f = 0.05),  border = NA)
    for (i in seq(0, by = 5, length.out = trunc(length(genes)/2))) rect(1.75 + i,log(min(plotdata$Relative.expression))*1.1, 4.25 + i,log(max(plotdata$Relative.expression))*1.1, col = adjustcolor(1, alpha.f = 0.05),  border = NA)
    abline(h = 0, col = 8, lty = 1)
    boxplot(log(Relative.expression) ~ Condition.Id + Gene, plotdata, pch = 20, xaxt = "n", xlab ="", ylab = "", add = T, boxlwd = 2.5, border = unique(plotdata[,c("Gene", "Color")])[[2]], at = positions)
    beeswarm(log(Relative.expression) ~ Condition.Id + Gene, plotdata, pch = 20, xaxt = "n", xlab ="", ylab = "", add = T, col = adjustcolor(unique(plotdata[,c("Gene", "Color")])[[2]], alpha.f = 0.2), at = positions, corral = "wrap")
    # mtext(xlabels, side = 1, line = 1.5, at = positions[seq(1, by = 2, length.out = length(positions)/2)] + 0.5, cex = 0.75)
    mtext(xlabels1, side = 1, line = 1, at = positions[seq(1, by = 4, length.out = length(positions)/4)] + 0.5, cex = 0.9)
    mtext(xlabels2, side = 1, line = 1, at = positions[seq(3, by = 4, length.out = length(positions)/4)] + 0.5, cex = 0.9)
    # text(positions[seq(2, by = 2, length.out = length(statpos))], statpos, labels = statlabels, xpd = NA, cex = 0.75)
    # text(positions[seq(2, by = 2, length.out = length(statpos))], statpos, labels = statlabels, xpd = NA, cex = 1.2)
    mtext(statlabels, side = 3, line = - 0.5, at = positions[seq(2, by = 2, length.out = length(statpos))] - 0.5, cex = 1.2)
    dev.off()
  }, error = function (e) cat("ERROR :", study, " ", technique, " ", conditionMessage(e), "\n"))
}
graphics.off()
time1 = (proc.time() - time0)[[3]]
paste0(round(round(round(time1/60)/60)/24), "d ", round(round(time1/60)/60)%%24, "h ", round(time1/60)%%60, "min ", round(time1%%60), "s")
