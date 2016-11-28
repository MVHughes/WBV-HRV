if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, stats, scales, AICcmodavg, arm, faraway, lme4, astsa, RHRV, pbkrtest, MASS, ggplot2, ggthemes, dplyr, timeSeries, timeDate, TSA, memisc, lubridate, forecast, zoo, Hmisc, reshape2, RColorBrewer, caTools, xts, scatterplot3d, lattice, latticeExtra, gridExtra, car)

options(max.print = 222)
#FileLocation = "G:/HrvValuesOnly.csv"
#FileLocation = "~/Documents/WBV-HRV"
#RecordLocation = "G:/"
#RecordLocation = "~/Documents/WBV-HRV"
FileLocation = "D:/Backup/WIP/Maggie/Thesis/9. Lotta-Bose - Repeated Measures Study of Whole Body Vibration Exposure and Heart Rate Variability in Truck Drivers/Data/Intermediate Data/WBV-HRV"
RecordLocation = "D:/Backup/WIP/Maggie/Thesis/9. Lotta-Bose - Repeated Measures Study of Whole Body Vibration Exposure and Heart Rate Variability in Truck Drivers/Data/Intermediate Data/WBV-HRV"
#setwd("~/Documents/WBV-HRV")
setwd("D:/Backup/WIP/Maggie/Thesis/9. Lotta-Bose - Repeated Measures Study of Whole Body Vibration Exposure and Heart Rate Variability in Truck Drivers/Data/Intermediate Data/WBV-HRV")


rm(RRSubsets)
rm(HrvValuesOnly)
rm(SubjSubsets)

par(mfrow=c(1,1))
HrvValuesOnly <- read.csv("HrvValuesOnly.csv", stringsAsFactors=FALSE)
attach(HrvValuesOnly)


HrvValuesOnly$Subject <- factor(HrvValuesOnly$Subject)
HrvValuesOnly$Period <- factor(HrvValuesOnly$Period)
HrvValuesOnly$Condition <- factor(HrvValuesOnly$Condition)
HrvValuesOnly$HRV.Segment <- factor(HrvValuesOnly$HRV.Segment)
HrvValuesOnly$id <- factor(HrvValuesOnly$id)
levels(HrvValuesOnly$id)


HrvValuesOnly$DateTime <- as.POSIXct(HrvValuesOnly$DateTime, format = "%m/%d/%Y %H:%M")


#Reformat to how RHRV expects date times to look
#HrvValuesOnly$DateTime <- format(HrvValuesOnly$DateTime, "%d/%m/%Y %H:%M:%S")

HrvValuesOnly <- subset(HrvValuesOnly, HrvValuesOnly$id != "S10NULL5M2" & 
                          HrvValuesOnly$id != "S10R3M1" & 
                          HrvValuesOnly$id != "S10WBV1M7" & 
                          HrvValuesOnly$id != "S10WBV3M7" & 
                          HrvValuesOnly$id != "S12WBV2M7"  & 
                          HrvValuesOnly$id != "S12WBV3M7" &
                          HrvValuesOnly$id != "S15NULL2M2" &
                          HrvValuesOnly$id != "S16NULL2M2" &
                          HrvValuesOnly$id != "S17WBV1M7" & 
                          HrvValuesOnly$id != "S18WBV2M7" &
                          HrvValuesOnly$id != "S1NULL1M9" &
                          HrvValuesOnly$id != "S4NULL4M2" &
                          HrvValuesOnly$id != "S4NULL6M3" &
                          HrvValuesOnly$id != "S4NULL6M4" &
                          HrvValuesOnly$id != "S4NULL6M5" &
                          HrvValuesOnly$id != "S4NULL6M6" &
                          HrvValuesOnly$id != "S4NULL6M7" &
                          HrvValuesOnly$id != "S6WBV3M7" &
                          HrvValuesOnly$id != "S7WBV1M7" &
                          HrvValuesOnly$id != "S8NULL0M2"
                        )

#Create additinal
#HrvValuesOnly$SubjCond <- mutate(HrvValuesOnly, concated_column = paste(Subject, Condition, sep = '_')) 

RRSubsets <- split(HrvValuesOnly, HrvValuesOnly$id, drop=TRUE)

ExportFileTitle <- "RRStatistics"
ExportFileName = paste(ExportFileTitle, "csv", sep = ".")
statistics <- as.matrix(t(data_frame(c("id", "Time", "SDNN", "pNN50", "rMSSD", "LFnu", "HFnu", "LFHF"))))
write.table(statistics, ExportFileName, sep=",", col.names=FALSE, row.names = FALSE, append = FALSE)


SubjSubsets <- split(HrvValuesOnly, HrvValuesOnly$Subject, drop=TRUE)

#SubjCondSubsets <- split(HrvValuesOnly, HrvValuesOnly$SubjCond, drop=TRUE)

#how to get list of first DateTime per each id value?
#tapply(HrvValuesOnly$DateTime,HrvValuesOnly$id,)

options(warn=-1)
for (i in 1:795) {

BegTime <- format(strptime(as.character(RRSubsets[[i]]$DateTime[1]), "%Y-%m-%d %H:%M:%S"), "%d/%m/%Y %H:%M:%S")

export <- data.frame(RRSubsets[[i]]$RR)
FileName = paste(RRSubsets[[i]]$id[1], "csv", sep = ".")
write.table(export, FileName, sep=",", col.names=FALSE, row.names = FALSE)

md <- CreateHRVData(Verbose = TRUE)
md <- LoadBeatRR(md, FileName, RecordPath = ".", scale = 0.001, datetime = BegTime)

md <- BuildNIHR(md)
md <- FilterNIHR(md)

md <- InterpolateNIHR(md, freqhr = 4, method = "linear")

MainTitle = paste("Non Interpolated Heart Rate During ", RRSubsets[[i]]$Condition[1], " Condition for Subject ", RRSubsets[[i]]$Subject[1], sep = " ")
XTitle = "Time"
YTitle = "Milliseconds"
filename <-  paste("NIHR", RRSubsets[[i]]$id[1], "jpg", sep=".")
jpeg(file=filename)
PlotNIHR(md)
dev.off()

MainTitle = paste("Interpolated Heart Rate During ", RRSubsets[[i]]$Condition[1], " Condition for Subject ", RRSubsets[[i]]$Subject[1], sep = " ")
XTitle = "Time"
YTitle = "Milliseconds"
filename <-  paste("HR", RRSubsets[[i]]$id[1], "jpg", sep=".")
jpeg(file=filename)
PlotHR(md)
dev.off()


md <- CreateTimeAnalysis(md, size = 300, interval = 7.8125)


#Current Dumpster Fire -- export summary files so can write to a file
#md <- CreateFreqAnalysis(md)
#md <- CalculatePowerBand(md, indexFreqAnalysis= 1,
#                         type = "wavelet", wavelet = "d4", bandtolerance = 0.01)

#md <- CreateFreqAnalysis(md)
#md <- CalculatePowerBand(md, indexFreqAnalysis= 2,
#                         type = "fourier", shift = 300, size = 300)

#PlotPowerBand(md, indexFreqAnalysis = 1, ymax = 700, ymaxratio = 50)

#HRVData$TimeAnalysis[[num+1]]$SDNN=sd(HRVData$Beat$RR)
SDNN <- md$TimeAnalysis[[1]]$SDNN


#RRDiffs = diff(HRVData$Beat$RR)
#RRDiffs50=RRDiffs[abs(RRDiffs)>50]
#HRVData$TimeAnalysis[[num+1]]$pNN50=100.0*length(RRDiffs50)/length(RRDiffs)
pNN50 <- md$TimeAnalysis[[1]]$pNN50

#HRVData$TimeAnalysis[[num+1]]$rMSSD=sqrt(mean(RRDiffs^2))
rMSSD <- md$TimeAnalysis[[1]]$rMSSD


#nu <- md$FreqAnalysis[[1]]$LF[1] + md$FreqAnalysis[[1]]$HF[1]
#LFnu <- md$FreqAnalysis[[1]]$LF[1] / nu
#HFnu <- md$FreqAnalysis[[1]]$HF[1] / nu
#LFHF <- md$FreqAnalysis[[1]]$LFHF[1]

id <- as.character.factor(RRSubsets[[i]]$id[1])
statistics <- data_frame(id, BegTime, SDNN, pNN50, rMSSD)

write.table(statistics, file = ExportFileName, append = TRUE, quote = FALSE, sep = ",", eol = "\n", row.names = FALSE, col.names = FALSE)

Subject <- RRSubsets[[i]]$Subject
DateTime <- RRSubsets[[i]]$DateTime
RR <- RRSubsets[[i]]$RR
Condition <- RRSubsets[[i]]$Condition
tmp <- data_frame(Subject, DateTime, RR, Condition)

tmp$DateTime <- as.POSIXct(tmp$DateTime, format = "%m/%d/%Y %H:%M:%S")


MainTitle = paste("R-R Durations Over ", RRSubsets[[i]]$Condition[1], " Condition for Subject ", RRSubsets[[i]]$Subject[1], sep = " ")
XTitle = "Time"
YTitle = "Milliseconds"
filename <-  paste("Time Distribution of ", MainTitle, "jpg", sep=".")
pdf(filename, width=6, height=4, paper='special')
ggplot(data=tmp, aes(y=RR, x=DateTime, colour = Condition)) + 
  geom_point(na.rm=FALSE) +
  xlab(XTitle) + ylab(YTitle) + ggtitle(MainTitle) +
  scale_x_datetime(labels = scales::date_format("%H:%M"))
dev.off()
}


SubjExportFileTitle <- "SubStatistics"
SubjExportFileName = paste(SubjExportFileTitle, "csv", sep = ".")
SubjStatistics <- as.matrix(t(data_frame(c("Subject", "SDNN", "pNN50", "rMSSD", "LFnu", "HFnu", "LFHF"))))
write.table(SubjStatistics, ExportFileName, sep=",", col.names=FALSE, row.names = FALSE, append = FALSE)


SubjSubsets <- split(HrvValuesOnly, HrvValuesOnly$Subject, drop=TRUE)
for (i in 1:18) {
 
  
  BegTime <- format(strptime(as.character(SubjSubsets[[i]]$DateTime[1]), "%Y-%m-%d %H:%M:%S"), "%d/%m/%Y %H:%M:%S")
  

  
  export <- data.frame(SubjSubsets[[i]]$RR)
  FileName = paste(SubjSubsets[[i]]$Subject[1], "csv", sep = ".")
  write.table(export, FileName, sep=",", col.names=FALSE, row.names = FALSE)
  
  md <- CreateHRVData(Verbose = TRUE)
  md <- LoadBeatRR(md, FileName, RecordPath = ".", scale = 0.001, datetime = BegTime)
  
  md <- BuildNIHR(md)
  md <- FilterNIHR(md)
  
  md <- InterpolateNIHR(md, freqhr = 4, method = "linear")
  
  
  MainTitle = paste("Non Interpolated Heart Rate for Subject ", SubjSubsets[[i]]$Subject[1], sep = " ")
  XTitle = "Time"
  YTitle = "Milliseconds"
  filename <-  paste("NIHR", SubjSubsets[[i]]$Subject[1], "jpg", sep=".")
  jpeg(file=filename)
  PlotNIHR(md)
  dev.off()
  
  MainTitle = paste("Interpolated Heart Rate for Subject ", SubjSubsets[[i]]$Subject[1], sep = " ")
  XTitle = "Time"
  YTitle = "Milliseconds"
  filename <-  paste("HR", SubjSubsets[[i]]$Subject[1], "jpg", sep=".")
  jpeg(file=filename)
  PlotHR(md)
  dev.off()
  
  md <- CreateTimeAnalysis(md, size = 300, interval = 7.8125)
  
  
  #Current Dumpster Fire -- export summary files so can write to a file
  #export <- as.character.factor(RRSubsets[[i]]$id[1])
  
  
  #md <- CreateFreqAnalysis(md)
  #md <- CalculatePowerBand(md, indexFreqAnalysis= 1,
  #                         type = "wavelet", wavelet = "d4", bandtolerance = 0.01)
  
  #md <- CreateFreqAnalysis(md)
  #md <- CalculatePowerBand(md, indexFreqAnalysis= 2,
  #                         type = "fourier", shift = 300, size = 300)
  
  #PlotPowerBand(md, indexFreqAnalysis = 1, ymax = 700, ymaxratio = 50)
  
  #HRVData$TimeAnalysis[[num+1]]$SDNN=sd(HRVData$Beat$RR)
  SDNN <- md$TimeAnalysis[[1]]$SDNN
  
  
  #RRDiffs = diff(HRVData$Beat$RR)
  #RRDiffs50=RRDiffs[abs(RRDiffs)>50]
  #HRVData$TimeAnalysis[[num+1]]$pNN50=100.0*length(RRDiffs50)/length(RRDiffs)
  pNN50 <- md$TimeAnalysis[[1]]$pNN50
  
  #HRVData$TimeAnalysis[[num+1]]$rMSSD=sqrt(mean(RRDiffs^2))
  rMSSD <- md$TimeAnalysis[[1]]$rMSSD
  
  
  #nu <- md$FreqAnalysis[[1]]$LF[1] + md$FreqAnalysis[[1]]$HF[1]
  #LFnu <- md$FreqAnalysis[[1]]$LF[1] / nu
  #HFnu <- md$FreqAnalysis[[1]]$HF[1] / nu
  #LFHF <- md$FreqAnalysis[[1]]$LFHF[1]
  
  id <- as.character.factor(SubjSubsets[[i]]$Subject[1])
  
  SubjStatistics <- data_frame(id, BegTime, SDNN, pNN50, rMSSD)
  
  write.table(SubjStatistics, file = SubjExportFileName, append = TRUE, quote = FALSE, sep = ",", eol = "\n", row.names = FALSE, col.names = FALSE)
 Subject <- SubjSubsets[[i]]$Subject
 DateTime <- SubjSubsets[[i]]$DateTime
 RR <- SubjSubsets[[i]]$RR
 Condition <- SubjSubsets[[i]]$Condition
  tmp <- data_frame(Subject, DateTime, RR, Condition)

  tmp$DateTime <- as.POSIXct(tmp$DateTime, format = "%m/%d/%Y %H:%M:%S tz=-5")
  
  

  MainTitle = paste("R-R Durations Over Entire Study for Subject", SubjSubsets[[i]]$Subject[1], sep = " ")
  filename <-  paste(MainTitle, "pdf", sep=".")
   XTitle = "Time"
    YTitle = "Milliseconds"
    pdf(filename, width=6, height=4, paper='special')
    ggplot(data=tmp, aes(y=RR, x=DateTime, colour = Condition)) + 
      geom_point(na.rm=FALSE) +
      xlab(XTitle) + ylab(YTitle) + ggtitle(MainTitle) +
      scale_x_datetime(labels = scales::date_format("%m/%d %I:%M"))
  dev.off()
  
}


SubjSubsets[[18]]$DateTime[1]


#To Do With Vibration 
#p <- ggplot() + 
#  geom_line(data = SubjStatistics[[i]], aes(x = SubjSubsets[[i]]$DateTime, y = SubjSubsets[[i]]$RR, color = "red")) +
#  geom_line(data = jobsAFAM2, aes(x = data_date, y = Percent.Change, color = "blue"))  +
#  xlab('data_date') +
#  ylab('percent.change') 