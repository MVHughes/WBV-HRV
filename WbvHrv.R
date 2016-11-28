if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, stats, AICcmodavg, arm, faraway, lme4, astsa, RHRV, pbkrtest, MASS, ggplot2, ggthemes, dplyr, timeSeries, timeDate, TSA, memisc, lubridate, forecast, zoo, Hmisc, reshape2, RColorBrewer, caTools, xts, scatterplot3d, lattice, latticeExtra, gridExtra, car)

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



#RRSubsets[[i]]$DateTime <- format(RRSubsets[[i]]$DateTime, "%m/%d/%Y %H:%M:%S")

BegTime <- RRSubsets[[i]]$DateTime[1]

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
PlotNIHR(md, main=MainTitle, xlab=XTitle, ylab=YTitle)
dev.off()

MainTitle = paste("Interpolated Heart Rate During ", RRSubsets[[i]]$Condition[1], " Condition for Subject ", RRSubsets[[i]]$Subject[1], sep = " ")
XTitle = "Time"
YTitle = "Milliseconds"
filename <-  paste("HR", RRSubsets[[i]]$id[1], "jpg", sep=".")
jpeg(file=filename)
PlotHR(md, main=MainTitle, xlab=XTitle, ylab=YTitle)
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
jpeg(file=filename)
ggplot(data=tmp, aes(y=RR, x=DateTime, colour = Condition)) + 
  geom_point(na.rm=FALSE) +
  xlab(XTitle) + ylab(YTitle) + ggtitle(MainTitle) +
  scale_x_datetime(labels = date_format("%H:%M"))
dev.off()
}


SubjExportFileTitle <- "SubStatistics"
SubjExportFileName = paste(SubjExportFileTitle, "csv", sep = ".")
SubjStatistics <- as.matrix(t(data_frame(c("Subject", "SDNN", "pNN50", "rMSSD", "LFnu", "HFnu", "LFHF"))))
write.table(SubjStatistics, ExportFileName, sep=",", col.names=FALSE, row.names = FALSE, append = FALSE)


SubjSubsets <- split(HrvValuesOnly, HrvValuesOnly$Subject, drop=TRUE)
for (i in 1:18) {
  id <- as.character.factor(SubjSubsets[[i]]$Subject[1])
  
  SubjSubsets[[i]]$DateTime <- as.POSIXct(SubjSubsets[[i]]$DateTime, format = "%d/%m/%Y %H:%M:%S")
  SubjSubsets[[i]]$DateTime <- format(SubjSubsets[[i]]$DateTime, "%m/%d/%Y %H:%M:%S")
  
  BegTime <- SubjSubsets[[i]]$DateTime[1]
  
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
  PlotNIHR(md, main=MainTitle, xlab=XTitle, ylab=YTitle)
  dev.off()
  
  MainTitle = paste("Interpolated Heart Rate for Subject ", SubjSubsets[[i]]$Subject[1], sep = " ")
  XTitle = "Time"
  YTitle = "Milliseconds"
  filename <-  paste("HR", SubjSubsets[[i]]$Subject[1], "jpg", sep=".")
  jpeg(file=filename)
  PlotHR(md, main=MainTitle, xlab=XTitle, ylab=YTitle)
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
  
  SubjStatistics <- data_frame(id, BegTime, SDNN, pNN50, rMSSD)
  
  write.table(SubjStatistics, file = SubjExportFileName, append = TRUE, quote = FALSE, sep = ",", eol = "\n", row.names = FALSE, col.names = FALSE)
 Subject <- SubjSubsets[[i]]$Subject
 DateTime <- SubjSubsets[[i]]$DateTime
 RR <- SubjSubsets[[i]]$RR
 Condition <- SubjSubsets[[i]]$Condition
  tmp <- data_frame(Subject, DateTime, RR, Condition)

  tmp$DateTime <- as.POSIXct(tmp$DateTime, format = "%m/%d/%Y %H:%M:%S")
  
  
  MainTitle = paste("Interpolated Heart Rate During ", RRSubsets[[i]]$Condition[1], " Condition for Subject ", RRSubsets[[i]]$Subject[1], sep = " ")
  filename <-  paste(MainTitle, "jpg", sep=".")
  jpeg(file=filename)
  PlotHR(md, main=MainTitle, xlab=XTitle, ylab=YTitle)
  dev.off()
  MainTitle = paste("R-R Durations Over Entire Study for Subject ", SubjSubsets[[i]]$Subject[1], sep = " ")
  XTitle = "Time"
    YTitle = "Milliseconds"
 
    ggplot(data=tmp, aes(y=RR, x=DateTime, colour = Condition)) + 
      geom_point(na.rm=FALSE) +
      xlab(XTitle) + ylab(YTitle) + ggtitle(MainTitle) +
      scale_x_datetime(labels = date_format("%H:%M"))
  dev.off()
  
}
)



#Steaming Dumpster Fire
#To Do With Vibration 
#p <- ggplot() + 
#  geom_line(data = SubjStatistics[[i]], aes(x = SubjSubsets[[i]]$DateTime, y = SubjSubsets[[i]]$RR, color = "red")) +
#  geom_line(data = jobsAFAM2, aes(x = data_date, y = Percent.Change, color = "blue"))  +
#  xlab('data_date') +
#  ylab('percent.change') 
#



#Giant Dumpster Fire
plot.xts(tmp3, main="Time Series Occupational Exposures and Heart Rate Variability")  #using all the defaults
#png("chartTimeSeries.png",width=640,height=567,units="px")
plot.xts(tmp3, screens=1, #screens=1 probably most appropriate for this application
         lwd = 2,
         legend.loc = "bottomright", auto.legend=TRUE,
          )
#dev.off()
#png("chartTimeSeries with extra.png",width=640,height=567,units="px")

plot(as.zoo(hear), las=1, xlab="", ylab="R-R Standard Deviation (ms)", col=2)
par(new=TRUE)               
plot(as.zoo(part),
     bty='n',               
     xaxt="n",               
     yaxt="n",              
     xlab="Time", ylab="", main="Time Series of Log PM2.5 Concentrations and SDNN")

axis(4, las=1)


#######TO FIX###################
sdnn64<-0
hrv <- CreateHRVData()
hrv <- SetVerbose(hrv, TRUE)
hrv <- BuildNIHR(hrv.64)
hrv = CreateFreqAnalysis(hrv)
hrv = CalculatePowerBand(hrv, indexFreqAnalysis = 1, size = 120, 
                        shift = 10, sizesp = 1024)


results.df = data.frame(cbind(In=results$resultIn,Out=results$resultOut))
results =  CreateTimeAnalysisByEpisodes(hrv, Tag="11")
sdnn64 <- c(sdnn64, results.df$In$SDNN)

segment = SplitPowerBandByEpisodes(hrv, indexFreqAnalysis = 1, Tag="1")
mean(segment$InEpisodes$LF, na.rm=TRUE)
mean(segment$InEpisodes$HF, na.rm=TRUE)

#########################################################################################################################################
pro$Measurement <- as.factor(pro$Measurement)
pro$HourDay <- as.factor(pro$HourDay)
pro$DateTime <- as.POSIXct(pro$DateTime, format = "%m/%d/%Y %H:%M")
pro$Subject<-factor(pro$Subject,levels=c(1,2,3,4,5,6,7,8),
               labels=c("Subject 1", "Subject 2", "Subject 3", "Subject 4", "Subject 5", "Subject 6", "Subject 7", "Subject 8")) 

tapply(pro$ZsAeq, pro$Seat, mean)
tapply(pro$ZsAeq, pro$Seat, sd)



##########################################################################################################################################


xyplot(Air ~ Time.Index | Subject, data = pro, type = c("b"), 
       main=expression('Time Series of  PM'[2.5]*' Exposure by Subject'), par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3}~" "))
main = "Time Series of PM2.5 Exposure by Subject"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()




mod <- lmList(SDANN ~ log(Air) | Category, data = pro, na.action = na.pass)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[,4,2]
rp1 = vector('expression',2)
rp1[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2[1],dig=3)))[2]
rp1[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p[1], digits = 2)))[2]
rp2 = vector('expression',2)
rp2[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2[2],dig=3)))[2]
rp2[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p[2], digits = 2)))[2]
rp3 = vector('expression',2)
rp3[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2[3],dig=3)))[2]
rp3[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p[3], digits = 2)))[2]
xyplot(SDANN ~ Air | Category, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distributions of Log Fine Particulate Matter and SDANN by Vehicle Speed", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "SDNN (ms)", par.settings = theEconomist.theme()
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################
main = "Distribution of Sound Exposures"
bwplot(pro$Subject ~ pro$Sound, 
       panel = function(x,y, ...) { panel.bwplot(x,y,...)}, main=main, xlab="dBA", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Time Series of Sound Exposures by Subject"
xyplot(Sound ~ Time.Index | Subject, data = pro, type = c("b"), 
       main=main, par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "dBA")
dev.copy(png, paste(main=main, "png", sep="."))
dev.off()
dev.off()


main = "Distribution of Z-Axis VDV(8) Exposure by Subject"
x = pro$ZsVDV8
y = pro$Subject
bwplot(y ~ x, data=pro, main=main, xlab=expression("m/" ~ s^{1.75}), ylab="Subject", par.settings = theEconomist.theme(), xlim=c(0, 3.5))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

summary(x)
sd(x)


###############################################################################################################################################
###############################################################################################################################################

###############################################################################################################################################
###############################################################################################################################################

main = "Distribution of Vector Sum VDV(8) Exposures and rMSSD"
y = pro$RMSSD
ylab = "rMSSD (ms)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
rp1 = vector('expression',1)
rp1[1] = substitute(expression(VDV(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]

mod2 <- lm(y ~ x2, data = pro, na.action = na.omit)
modsum2 <- summary(mod2)
r22 = modsum2$adj.r.squared
my.p2 = modsum2$coefficients[2,1]
rp2 = vector('expression',1)
rp2[1] = substitute(expression(Lag~VDV(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r22,dig=1)))[2]


mod3 <- lm(y ~ x3, data = pro, na.action = na.omit)
modsum3 <- summary(mod3)
r23 = modsum3$adj.r.squared
my.p3 = modsum3$coefficients[2,1]
rp3 = vector('expression',1)
rp3[1] = substitute(expression(Cum.~VDV(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r23,dig=1)))[2]


xyplot(y ~ x, data = pro,  type = c("p", "r"), main=main, xlab = xlab, ylab = ylab, par.settings = theEconomist.theme(), col="navy") + 
  #  xyplot(y ~ x2, data = na.omit(pro),  type = c("p", "r"), par.settings = theEconomist.theme(), col="grey") +
  xyplot(y ~ x3, data = na.omit(pro),  type = c("p", "r"), par.settings = theEconomist.theme(), col="black")

trellis.focus("toplevel")
panel.text(0.3, 0.87, rp1, cex = 1.2, font = 2, col="navy")
#panel.text(0.3, 0.84, rp2, cex = 1.2, font = 2, col="grey")
panel.text(0.3, 0.81, rp3, cex = 1.2, font = 2, col="black")
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

###############################################################################################################################################


###############################################################################################################################################

main = "Distribution of Z-Axis Mean Amplitude Exposures and rMSSD"
y = pro$RMSSD
ylab = "rMSSD (ms)"
mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
rp1 = vector('expression',1)
rp1[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]

xyplot(y ~ x, data = pro,  type = c("p", "r"), main=main, xlab = xlab, ylab = ylab, par.settings = theEconomist.theme(), col="navy")


trellis.focus("toplevel")
panel.text(0.3, 0.87, rp1, cex = 1.2, font = 2, col="navy")
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (FreqPeakZs) | Category, data = pro, na.action = na.pass)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[,4,2]
rp1 = vector('expression',2)
rp1[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2[1],dig=3)))[2]
rp1[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p[1], digits = 2)))[2]
rp2 = vector('expression',2)
rp2[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2[2],dig=3)))[2]
rp2[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p[2], digits = 2)))[2]
rp3 = vector('expression',2)
rp3[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2[3],dig=3)))[2]
rp3[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p[3], digits = 2)))[2]
xyplot(LF.HF ~ FreqPeakZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Frequency Peak (Hz)and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################


xyplot(
  SDANN ~ Sound,
  data = pro, par.settings = theEconomist.theme(),
  main="Sound Exposure and SDANN by Subject", xlab="Sound (dBA)", ylab="SDNN (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)



mod <- lm(RMSSD ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and RMSSD", xlab = "Mean Amplitude", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


mod.dat <- groupedData(Hfnu ~ FreqAmpZs | Subject, data=pro) 
mod0 <- lme(Hfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Hfnu ~ FreqAmpZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Hfnu ~ FreqAmpZs, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1())
mod3 <- lme(Hfnu ~ FreqAmpZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1())


names <- c("mod1", "mod2", "mod3") 
aictab(cand.set=list(mod1, mod2, mod3), modnames=names, second.ord=TRUE, sort=TRUE)

summary(mod3)


##########################################################################################################################################################

##########################################################################################################################################################
###################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ MeanPowerXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(MeanPowerXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ MeanPowerXs, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

opar <- par(mfrow = c(2,2), oma = c(0,0,1.1,0))
plot(m, las = 1)

par(mfrow = c(1,1))

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
pro[c(1, 26), 1:3]
d1 <- cooks.distance(m)
r <- stdres(m)
a <- cbind(pro, d1, r)
a[d1 > 4/78, ]
rabs <- abs(r)
a <- cbind(pro, d1, r, rabs)
asorted <- a[order(-rabs), ]
asorted[1:10, ]

mod.dat <- groupedData(LF.HF ~ MeanPowerXs | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ MeanPowerXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ MeanPowerXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ MeanPowerXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))
