time0 = proc.time()

###----- 1. Set variables
project = "Heart progerin" # mir29   Heart progerin   Heart lamin
# techniques = 
technique = "qPCR"
analyzeall = F # Default is F
housekeeping = c("Arbp", "Tbp", "Hprt", "Ywhaz", "Lmnb", "Ube2b")
maxctdif = 0.5 # Oslo 1

packages = c("data.table", "matrixStats", "zoo", "tidyr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) install.packages(setdiff(packages, rownames(installed.packages())))  
for (i in packages) library(i, character.only = T)

if (.Platform$OS.type == "unix") setwd("/Volumes/Victor/") else setwd("S:/LAB_VA/LAB/Victor/")

datefun = function (x) {
  if (any(grepl("/",x))) {
    date1 = sub("/.*", "", x)
    date2 = sub("/.*", "", sub("[[:digit:]]*/", "", x))
    date3 = sub(".*/", "", x)
  } else if (any(grepl("-",x))) {
    date1 = sub("-.*", "", x)
    date2 = sub("-.*", "", sub("[[:digit:]]*/", "", x))
    date3 = sub(".*-", "", x)
  }
  if (max(nchar(date3)) == 4) {
    if (max(as.numeric(date2), na.rm = T) > 12 & max(as.numeric(date1), na.rm = T) <= 12) {
      as.Date(paste0(date3, "/", date1, "/", date2))
    } else if (max(as.numeric(date1), na.rm = T) > 12 & max(as.numeric(date2)) <= 12){
      as.Date(paste0(date3, "/", date2, "/", date1))
    } else as.Date(paste0(date3, "/", date2, "/", date1)) # Might mistake days for months
  } else if (max(nchar(date1)) == 4) {
    as.Date(x)
  } else if (max(as.numeric(date1), na.rm = T) > 31 | max(as.numeric(date3), na.rm = T) > as.integer(format(Sys.Date(),"%Y"))-2000) {
    as.Date(paste0(20,x))
  } else if (max(as.numeric(date3), na.rm = T) > 31 | max(as.numeric(date1), na.rm = T) > as.integer(format(Sys.Date(),"%Y"))-2000) {
    if (max(as.numeric(date2), na.rm = T) > 12 & max(as.numeric(date1), na.rm = T) <= 12) {
      as.Date(paste0(20, date3, "/", date1, "/", date2))
    } else if (max(as.numeric(date1), na.rm = T) > 12 & max(as.numeric(date2), na.rm = T) <= 12){
      as.Date(paste0(20, date3, "/", date2, "/", date1))
    } else as.Date(paste0(20, date3, "/", date2, "/", date1)) # Might mistake days for months
  } else if (max(as.numeric(date2), na.rm = T) > 12) {
    as.Date(paste0(20, date3, "/", date1, "/", date2))
  } else as.Date(paste0(20, date3, "/", date2, "/", date1)) # Might mistake between days, months and years
}
# outliers = function (x) {
#   if (length(x[!is.na(x)]) == 1) {
#     rep(NA, length(x))
#   } else if (sd(x, na.rm = T) >= maxctdif/sqrt(3)) {
#     y = which(abs(x-mean(x)) == max(abs(x-mean(x))))
#     x[y] = NA
#     if (sd(x, na.rm = T) >= maxctdif/sqrt(2)) rep(NA, length(x)) else x
#   } else x
# }
outliers = function (x) {
  if (length(x[!is.na(x)]) <= 1) {
    rep(NA, length(x))
  } else if (sd(x, na.rm = T) >= maxctdif/sqrt(3)) {
    y = which(abs(x-mean(x)) == max(abs(x-mean(x))))
    x[y] = NA
    if (length(x[!is.na(x)]) <= 1) {
      rep(NA, length(x))
    } else if (sd(x, na.rm = T) >= maxctdif/sqrt(2)) rep(NA, length(x)) else x
  } else x
}

baseroute = paste0(project, " project/Raw data/", technique, "/")
route = paste0(baseroute, "Raw txt/")


###----- 2. Import data
previous = c()
if (file.exists(paste0(baseroute, "normalized.txt")) & analyzeall == F) previous = data.frame(fread(paste0(baseroute, "normalized.txt"), header = F, encoding = "Latin-1", sep = "\t", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "")))[,1]

files = setdiff(grep("Run information", dir(route), value = T, ignore.case = T), previous)
normalized = c()
info = c()
i = files[1]
for (i in files) {
  importdata = data.frame(t(fread(paste0(route, i), header = F, encoding = "Latin-1", sep = "\t", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "NaN", ""))), stringsAsFactors = F)
  names(importdata) = importdata[1,]
  importdata = importdata[-1,]
  importdata$`Run Started` = as.POSIXct(importdata$`Run Started`, format = "%m/%d/%Y %H:%M:%S")
  importdata$`Run Ended` = as.POSIXct(importdata$`Run Ended`, format = "%m/%d/%Y %H:%M:%S")
  info = rbind.data.frame(info, importdata)
  normalized = c(normalized,i)
}
files = paste0(setdiff(gsub(" -  Quantification Summary_0.txt", "", grep("Quantification Summary_0", dir(route), value = T, ignore.case = T)), gsub("_Run information.txt", "",ignore.case = T, previous)), " -  Quantification Summary_0.txt")
rawdata_exp = c()
i = files[1]
for (i in files) {
  importdata = data.frame(fread(paste0(route, i), header = T, encoding = "Latin-1", sep = "\t", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "NaN", "NeuN", "")))
  importdata$Cq = as.numeric(gsub(",", ".", importdata$Cq))
  importdata = importdata[, c("Well", "Target", "Sample", "Cq")]
  importdata$Target = tolower(importdata$Target)
  importdata$Target = paste0(toupper(substring(importdata$Target, 1,1)), substring(importdata$Target, 2))
  importdata$Sample = toupper(importdata$Sample)
  for (target in unique(importdata$Target)) for (sample in unique(importdata$Sample)) importdata$Well[importdata$Target == target & importdata$Sample == sample] = 1:nrow(importdata[importdata$Target == target & importdata$Sample == sample,])

  # importdata$Cq2 = importdata$Cq
  contamination = data.frame(table(importdata$Target[importdata$Sample == "AGUA" & !is.na(importdata$Cq[importdata$Sample == "AGUA"])]))
  contamination = contamination[contamination$Freq > 1,-2]
  importdata$Cq[importdata$Target %in% contamination] = NA
  importdata = importdata[importdata$Sample != "AGUA",]
  for (target in unique(importdata$Target[!is.na(importdata$Cq)])) {
    for (sample in unique(importdata$Sample[!is.na(importdata$Cq) & importdata$Target == target])) {
      importdata$Cq[importdata$Target == target & importdata$Sample == sample] = outliers(importdata$Cq[importdata$Target == target & importdata$Sample == sample])
    }
  }
  importdata_hk = importdata[!duplicated(importdata[,c("Target", "Sample")]) & importdata$Target %in% housekeeping,]
  for (target in unique(importdata_hk$Target)) for (sample in unique(importdata_hk$Sample)) importdata_hk$Cq[importdata_hk$Target == target & importdata_hk$Sample == sample] = median(importdata$Cq[importdata$Target == target & importdata$Sample == sample], na.rm = T)
  for (target in unique(importdata_hk$Target)) {
    importdata = merge(importdata, importdata_hk[importdata_hk$Target == target, c("Sample", "Cq")], by = "Sample", suffixes = c(".Target", target), all.x = T)
    names(importdata)[ncol(importdata)] = paste0("Cq.", target)
    importdata = importdata[importdata$Target != target,]
  }
  importdata[,"Relexp"] = 2^(rowMeans(importdata[,grep("Target", grep("Cq.", names(importdata), value = T), value = T, invert = T)], na.rm = T) - importdata$Cq.Target)
  
  Experiment = gsub(" -  Quantification Summary_0.txt", "", i)
  importdata = cbind.data.frame(Experiment = Experiment, importdata[,c(1:3,ncol(importdata))])
  rawdata_exp = rbind.data.frame(rawdata_exp, importdata)
}

###----- 3. Consolidate data
rawdata_exp = merge(cbind.data.frame(gsub(".pcrd", "", info[,1]), as.Date(info[,5])), rawdata_exp, by = 1, all = T)
rawdata_exp = rawdata_exp[,c(1,3,2,4:ncol(rawdata_exp))]
names(rawdata_exp)[1:3] = c("Experiment", "Id", "Date")

targets = rawdata_exp[!duplicated(rawdata_exp[,c(1,3,5)]),c(1,3,5)]

rawdata_exp = rawdata_exp[,-1]
rawdata_exp$Time = 0
rawdata_exp$Time[grep("_", rawdata_exp$Id)] = as.numeric(gsub(".*_", "", grep("_", rawdata_exp$Id, value = T)))
rawdata_exp$Id[grep("_", rawdata_exp$Id)] = gsub("_.*", "", grep("_", rawdata_exp$Id, value = T))



###----- 4. Export data
if (length(files) > 0) {
  if (!file.exists(paste0(baseroute, "rawdatalong.txt")) | analyzeall == T) write.table(t(data.frame(names(rawdata_exp))), file = paste0(baseroute, "rawdatalong.txt"), row.names = F, col.names = F, sep = "\t", append = F, quote = F)
  fwrite(rawdata_exp, file = paste0(baseroute, "rawdatalong.txt"), row.names = F, col.names = F, sep = "\t", append = T)
  write.table(normalized, file = paste0(baseroute, "normalized.txt"), row.names = F, col.names = F, sep = "\t", append = !analyzeall, quote = F, eol = "\r\n")
}

rawdata = data.frame(fread(paste0(baseroute, "rawdatalong.txt"), header = T, encoding = "Latin-1", sep = "\t", stringsAsFactors = F, na.strings = c("N/A", "NA", "#N/A", "")))
rawdata = spread(rawdata, Target, Relexp)[,-3]
names(rawdata)[-1:-3] = paste0(names(rawdata)[-1:-3], " relative expression")
fwrite(rawdata, file = paste0(baseroute, "rawdata.txt"), row.names = F, col.names = T, sep = "\t")


# option to exclude experiments, samples and genes
# normalize to level 1 in stats
# collapse plots in one

time1 = (proc.time() - time0)[[3]]
paste0(round(round(round(time1/60)/60)/24), "d ", round(round(time1/60)/60)%%24, "h ", round(time1/60)%%60, "min ", round(time1%%60), "s")
