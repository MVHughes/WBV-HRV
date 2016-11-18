if (!require("pacman")) install.packages("pacman")
pacman::p_load(stats, AICcmodavg, arm, faraway, lme4, astsa, RHRV, pbkrtest, MASS, ggplot2, ggthemes, dplyr, timeSeries, timeDate, TSA, memisc, lubridate, forecast, zoo, stats, Hmisc, reshape2, RColorBrewer, caTools, xts, scatterplot3d, lattice, latticeExtra, gridExtra, car)

options(max.print = 22)

par(mfrow=c(1,1))
attach(HrvValuesOnly)

proposal <- read.delim("~/Dropbox/thesis/HRVPilotTidyData_Cuts_Static.txt", header=FALSE, stringsAsFactors=FALSE)
proposal <- rename(proposal, Measurement=V1, DateTime=V2, Speed=V3, Category=V4, Air=V5, Big_Air=V6, Sound=V7, A8=V8, VSA8=V9, VDV8=V10, VSVDV8=V11, SDANN=V12, RMSSD=V13, R_RMSSD=V14, pNN50=V15, R50=V16, pNN20=V17, R20=V18)
proposal <- proposal[-1,]
proposal$DateTime <- as.POSIXct(proposal$DateTime, format="%m/%d/%y %H:%M")

#Create table of subject measurement durations, 
sub59 <- subset(proposal, Measurement == 59)


hrvall66 <- read.csv("~/Dropbox/hrvall66.csv", stringsAsFactors=FALSE)
hrvall66$Date<-seq.POSIXt(c(ISOdatetime(2015, 8, 6, 13, 25, 23)), by="min", length.out=118)
SDNN <- as.xts(hrvall66$SDNN, hrvall66$Date)
xts::plot.xts(SDNN)
xts::plot.xts(as.xts(hrvall66$Pm2.5, hrvall66$Date))
xts::plot.xts(as.xts(hrvall66$Noise, hrvall66$Date))
xts::plot.xts(as.xts(hrvall66$RMS.VDV, hrvall66$Date))
datepm <- na.omit(hrvall66$Pm2.5)




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

tapply(pro$VSAeq, pro$Seat, mean)
tapply(pro$VSAeq, pro$Seat, sd)

tapply(pro$ZsVDV8, pro$Seat, mean)
tapply(pro$ZsVDV8, pro$Seat, sd)

tapply(pro$VSVDV8, pro$Seat, mean)
tapply(pro$VSVDV8, pro$Seat, sd)

##########################################################################################################################################

main="Distribution of PM2.5 Exposures by Subject"
bwplot(pro$Subject ~ pro$Air, 
       panel = function(x,y, ...) { panel.bwplot(x,y,...)}, main=expression('Distribution of PM'[2.5]*' Exposures by Subjects'), 
       xlab=expression(mu~"g/" ~ m^{3}), par.settings = theEconomist.theme(), na.action = na.omit)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distribution of PM2.5 Exposures by Subject with Logs"
bwplot(pro$Subject ~ pro$Air, panel = function(x,y, ...) { panel.bwplot(x,y,...)}, main=expression('Distribution of PM'[2.5]*' Exposures by Subjects'), 
       xlab=expression("log "~mu~"g/"~ m^{3}), par.settings = theEconomist.theme(), na.action = na.omit,scales=list(x=list(log=2)))
  dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

bwplot(pro$Subject ~ log(pro$Air), 
       panel = function(x,y, ...) { panel.bwplot(x,y,...)}, main=expression('Distribution of PM'[2.5]*' Exposures by Subject'), xlab=expression("log " ~ mu~"g/" ~ m^{3}), par.settings = theEconomist.theme(), na.action = na.omit)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(Air ~ Time.Index | Subject, data = pro, type = c("b"), 
       main=expression('Time Series of  PM'[2.5]*' Exposure by Subject'), par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3}~" ")))
main = "Time Series of PM2.5 Exposure by Subject"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(Air, RMSSD, lag.max=10, main=expression('Lag of  PM'[2.5]*' Exposure and rMSSD'),  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
main="Lag Correlation Fine Particulate Matter and RMSSD"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(Air, SDANN, lag.max=10, main="Lag Correlation Fine Particulate Matter and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(Air, pNN50, lag.max=10, main="Lag Correlation Fine Particulate Matter and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(Air, Lfnu, lag.max=10, main=expression('Lag of  PM'[2.5]*' Exposure and LFnu'), xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
main= "Lag Correlation Fine Particulate Matter and LFnu"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(Air, Hfnu, lag.max=10, main="Lag Correlation Fine Particulate Matter and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(Air, LF.HF, lag.max=10, main="Lag Correlation Fine Particulate Matter and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(AirLag, RMSSD, lag.max=10, na.action = na.pass, main=expression('Lag Correlation Lagged  PM'[2.5]*' Exposure and rMSSD'), xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
main="Lag Correlation Lagged PM2.5 Exposure and rMSSD"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(AirLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Fine Particulate Matter and SDNN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
main="Lag Correlation Lagged Fine Particulate Matter and SDANN"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(AirLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Fine Particulate Matter and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(AirLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Fine Particulate Matter and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(AirLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Fine Particulate Matter and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(AirLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Fine Particulate Matter and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ log(Air), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ Air, data = pro, scales = list (log = TRUE), type = c("p", "r"), main=expression('Distributionn  PM'[2.5]*' Exposure and rMSSD'), 
       xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "rMSSD (ms)", par.settings = theEconomist.theme()) + 
  xyplot(RMSSD ~ AirLag, data = pro, scales = list (log = TRUE), type = c("p", "r"), main=expression('Distributionn  PM'[2.5]*' Exposure and rMSSD'), 
         xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "rMSSD (ms)", par.settings = theEconomist.theme(), col="grey")
mod1 <- lm(RMSSD ~ log(AirLag), data = pro)
modsum1 <- summary(mod1)
r2 = modsum1$adj.r.squared
my.p = modsum1$coefficients[2,4]
rp2 = vector('expression',2)
rp2[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp2[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]


main = "Distribution of Log Fine Particulate Matter Exposures and RMSSD"
trellis.focus("toplevel")
panel.text(0.3, 0.87, format(round(rp, 3), nsmall = 3), cex = 1.2, font = 2, col="navy")
panel.text(0.3, 0.84, format(round(rp2, 3), nsmall = 3) cex = 1.2, font = 2, col="grey")
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ log(Air), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ Air, data = pro, scales = list (log = TRUE), type = c("p", "r"), main=expression('Distributionn  PM'[2.5]*' Exposure and SDNN'), 
       xlab = expression(Log ~ of ~ PM[2.5]* ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "SDNN (ms)", par.settings = theEconomist.theme()) + 
  xyplot(SDANN ~ AirLag, data = pro, scales = list (log = TRUE), type = c("p", "r"),  
         par.settings = theEconomist.theme(), col="grey")
mod1 <- lm(SDANN ~ log(AirLag), data = pro)
modsum1 <- summary(mod1)
r2 = modsum1$adj.r.squared
my.p = modsum1$coefficients[2,4]
rp2 = vector('expression',2)
rp2[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2,dig=3)))[2]
rp2[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p, digits = 2)))[2]

main = "Distribution of Log Fine Particulate Matter Exposures and SDNN"
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2, col="navy")
panel.text(0.3, 0.84, rp2, cex = 1.2, font = 2, col="grey")
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


###proposal 
options(digits=2)
x <- VSVDV8
x2 <- VSVDV8Cum
y <- pNN50

mod <- lm(y ~ x, data = pro)
modsum <- summary(mod)
r2 = round(modsum$adj.r.squared, 2)
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(VDV(8) ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(y ~ x, data = na.omit(pro), type = c("p", "r"), main=expression('Vector Sum VDV(8) Exposure and pNN50'), 
       xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme()) + 
  xyplot(y ~ x2, data = na.omit(pro), type = c("p", "r"), col="grey")

mod1 <- lm(y ~ x2, data = pro)
modsum1 <- summary(mod1)
r2 = round(modsum1$adj.r.squared, 2)
my.p = modsum1$coefficients[2,4]
rp2 = vector('expression',2)
rp2[1] = substitute(expression(Cum.~VDV(8) ~ italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r2,dig=3)))[2]
rp2[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                    list(MYOTHERVALUE = format(my.p, digits = 2)))[2]


main = "Distribution of Vector Sum VDV(8) Exposures and pNN50"
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2, col="navy")
panel.text(0.3, 0.84, rp2, cex = 1.2, font = 2, col="grey")
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()



mod <- lm(Lfnu ~ log(Air), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ Air, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distribution of Log Fine Particulate Matter Exposures and LFnu", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.88, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ log(Air), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ Air, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distribution of Log Fine Particulate Matter Exposures and HFnu", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ log(Air), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ Air, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distribution of Log Fine Particulate Matter Exposures and LF:HF Ratio", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ log(Air) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ Air | Category, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distributions of Log Fine Particulate Matter and RMSSD by Vehicle Speed", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
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
xyplot(SDANN ~ Air | Category, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distributions of Log Fine Particulate Matter and SDANN by Vehicle Speed", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ log(Air) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ Air | Category, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distributions of Log Fine Particulate Matter and pNN50 by Vehicle Speed", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ log(Air) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ Air | Category, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distributions of Log Fine Particulate Matter and LFnu by Vehicle Speed", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ log(Air) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ Air | Category, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distributions of Log Fine Particulate Matter and HFnu by Vehicle Speed", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ log(Air) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ Air | Category, data = pro, scales = list (log = TRUE), type = c("p", "r"), main="Distributions of Log Fine Particulate Matter and LF:HF by Vehicle Speed", xlab = expression(Log ~ of ~ Fine ~ Particulate ~ Matter ~ Exposure ~ (mu ~ "g/" ~ m^{3})), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
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

main="Lag Correlation Sound and rMSSD"
ccf(Sound, RMSSD, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Lag Correlation Sound and SDNN"
ccf(Sound, SDANN, lag.max=10, main=main, xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Lag Correlation Sound and pNN50"
ccf(Sound, pNN50, lag.max=10, main=main, xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Lag Correlation Sound and LFnu"
ccf(Sound, Lfnu, lag.max=10, main=main, xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Lag Correlation Sound and HFnu"
ccf(Sound, Hfnu, lag.max=10, main=main, xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Lag Correlation Sound and LF-HF"
ccf(Sound, LF.HF, lag.max=10, main=main, xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Lag Correlation Lagged Sound and rMSSD"
ccf(SoundLag, RMSSD, lag.max=10, na.action = na.pass, main=main, xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(SoundLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Sound and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(SoundLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Sound and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(SoundLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Sound and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(SoundLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Sound and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Lag Correlation Lagged Sound and LF-HF"
ccf(SoundLag, LF.HF, lag.max=10, na.action = na.pass, main=main, xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distribution of Sound Exposures and rMSSD"
mod <- lm(RMSSD ~ (Sound), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ Sound, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "rMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distribution of Sound Exposures and SDNN"
mod <- lm(pro$SDANN ~ pro$Sound, data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pro$SDANN ~ pro$Sound, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distribution of Sound Exposures and pNN50"
mod <- lm(pNN50 ~ (Sound), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ Sound, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distribution of Sound Exposures and LFnu"
mod <- lm(Lfnu ~ (Sound), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ Sound, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.88, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distribution of Sound Exposures and HFnu"
mod <- lm(Hfnu ~ (Sound), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ Sound, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distribution of Sound Exposures and LF-HF Ratio"
mod <- lm(LF.HF ~ Sound, data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ Sound, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distributions of Sound and rMSSD by Vehicle Speed"
mod <- lmList(RMSSD ~ (Sound) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ Sound | Category, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "rMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distributions of Sound and SDNN by Vehicle Speed"
mod <- lmList(SDANN ~ (Sound) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ Sound | Category, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distributions of Sound and pNN50 by Vehicle Speed"
mod <- lmList(pNN50 ~ (Sound) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ Sound | Category, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distributions of Sound and LFnu by Vehicle Speed"
mod <- lmList(Lfnu ~ (Sound) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ Sound | Category, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distributions of Sound and HFnu by Vehicle Speed"
mod <- lmList(Hfnu ~ (Sound) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ Sound | Category, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distributions of Sound and LF-HF by Vehicle Speed"
mod <- lmList(LF.HF ~ (Sound) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ Sound | Category, data = pro,  type = c("p", "r"), main=main, xlab = "dBA", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

attach(pro)
###########################################################################################################################################################
#pro$XsAeq[pro$XsAeq>2] <- .165

##########################################################################################################################################################
#proposal
##########################################################################################################################################################
####pro$YsAeq[pro$YsAeq>0.4]<-0.2
####pro$ZsAeq[pro$ZsAeq>1.0]<-0.4
#####pro$MeanAmpXs[pro$MeanAmpXs>1000]
#####pro$MeanPowerXs[pro$MeanPowerXs>3000]<-1571
#####pro$VSAeq <- (((1.4*(pro$XsAeq*pro$XsAeq))+((1.4*(pro$YsAeq * pro$YsAeq))+(pro$ZsAeq*pro$ZsAeq)))^(0.5))
setwd("C:\\users\\mhughes\\desktop")

#main = list(
#  label = expression('Distribution of Log PM'[2.5]*' Exposure by Subject'),
#  label = "Distribution of Noise Exposure by Subject (dBA)"
#  cex = 1.5)
x = pro$VSVDV8
y = pro$Subject
#  label=expression(paste("PM"[2.5], " ", mu, "g/m"^3)),
                  
bwplot(y ~ x, data=pro, main=" ", xlab=list(
  label = expression("Vector Sum Vibration Dose Value [VDV(8)] m/s"^1.75),
  cex=1.5),
  ylab=list(label=" ", cex=1.5), par.settings = theEconomist.theme(), 
       scales=list(tck=c(1,0), x=list(cex=1.5), y=list(cex=1.5)))
dev.copy(png, paste("Distribution of Vector Sum VDV(8) Exposure by Subject", "png", sep="."))
dev.off()
dev.off()


#xlab=expression("m/" ~ s^{1.75})



summary(x)
sd(x)
####pro$XsVDV8[pro$YsVDV8>10]<-1.4
####pro$YsVDV8[pro$YsVDV8>2]<-1.2
####pro$ZsVDV8[pro$ZsVDV8>10.0]<-1.7
##pro$VSVDV8 <- (((1.4*(pro$XsVDV8*pro$XsVDV8))+((1.4*(pro$YsVDV8 * pro$YsVDV8))+(pro$ZsVDV8*pro$ZsVDV8)))^(0.5))



x = pro$YsVDV8
x2 = pro$YsVDV8Cum
y = pro$pNN50


#x2 = pro$VSAeqLag
#x3 = pro$VSAeqCum
#xlab = expression(A(8)~"m/" ~ s^{2})

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
rp2[1] = substitute(expression(Cum~VDV(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r22,dig=1)))[2]


#mod3 <- lm(y ~ x3, data = pro, na.action = na.omit)
#modsum3 <- summary(mod3)
#r23 = modsum3$adj.r.squared
#my.p3 = modsum3$coefficients[2,1]
#rp3 = vector('expression',1)
#rp3[1] = substitute(expression(Cum.~VDV(8)~italic(R)^2 == MYVALUE), 
#                    list(MYVALUE = format(r23,dig=1)))[2]


xyplot(y ~ x, data = na.omit(pro),  type = c("p", "r"), main = list(
  label = expression(paste("Correlation Y-Axis VDV(8) Exposure and pNN50")),
  cex = 1.5), 
  xlab= list(
    label = expression(paste("VDV(8) m/s"^1.75)),
    cex= 1.5), 
  ylab=list(
    label = "Successive Beats Differing by More Than 50 ms (%)",
    cex = 1.5), 
  par.settings = theEconomist.theme(), col="navy", scales=list(cex=1.5)) + 
  xyplot(y ~ x2, data = na.omit(pro),  type = c("p", "r"), main = list(
    label = expression(paste("Correlation Y-Axis VDV(8) Exposure and pNN50")),
    cex = 1.5), 
    par.settings = theEconomist.theme(), col="grey") #+
#  xyplot(y ~ x3, data = na.omit(pro),  type = c("p", "r"), par.settings = theEconomist.theme(), col="black")

trellis.focus("toplevel")
panel.text(0.28, 0.86, rp1, cex = 1.2, font = 2, col="navy")
panel.text(0.30, 0.82, rp2, cex = 1.2, font = 2, col="grey")
#panel.text(0.3, 0.81, rp3, cex = 1.2, font = 2, col="black")
trellis.unfocus()
dev.copy(png, paste("Correlation Y-Axis VDV(8) and pNN50", "png", sep="."))
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






main = "Distribution of Vector Sum A(8) Exposure by Subject"
x = pro$VSAeq
y = pro$Subject
bwplot(y ~ x, data=pro, main=main, xlab=expression("m/" ~ s^{2}), ylab="Subject", par.settings = theEconomist.theme(), xlim=c(0, 0.7))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

summary(x)
sd(x)
####pro$XsVDV8[pro$YsVDV8>10]<-1.4
####pro$YsVDV8[pro$YsVDV8>2]<-1.2
####pro$ZsVDV8[pro$ZsVDV8>10.0]<-1.7
##pro$VSVDV8 <- (((1.4*(pro$XsVDV8*pro$XsVDV8))+((1.4*(pro$YsVDV8 * pro$YsVDV8))+(pro$ZsVDV8*pro$ZsVDV8)))^(0.5))

main = "Distribution of Z-Axis VDV(8) Exposure by Subject"
x = pro$ZsVDV8
y = pro$Subject
bwplot(y ~ x, data=pro, main=main, xlab=expression("m/" ~ s^{1.75}), ylab="Subject", par.settings = theEconomist.theme(), xlim=c(0, 3.5))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

summary(x)
sd(x)

####pro$MeanPowerZs[pro$MeanPowerZs>3500]<-396



main = "Distribution of Z-Axis Mean Power Exposure by Subject"
x = pro$MeanPowerZs
y = pro$Subject
xlab = "Mean Power Exposure"
ylab = "Subject"
xlim=c(0, 380)


bwplot(y ~ x, data=pro, main=main, xlab=xlab, ylab=ylab, par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

summary(x)
sd(x)

##pro$FreqPeakXs[pro$FreqPeakXs>3100]<-396
##pro$FreqPeakYs[pro$FreqPeakYs>3000]<-750
main = "Distribution of Z-Axis Mean Frequency Exposure by Subject"
x = pro$FreqPeakZs
y = pro$Subject
xlab = "Frequency (Hz)"
ylab = "Subject"
xlim=c(0, 2800)

bwplot(y ~ x, data=pro, main=main, xlab=xlab, ylab=ylab, par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

summary(x)
sd(x)


##pro$FreqAmpYs[pro$FreqAmpYs>1000]<-0.460
##pro$FreqAmpYs[pro$FreqAmpYs>1.5]<-0.460
##pro$FreqAmpZs[pro$FreqAmpZs>6]<-0.976
main = "Distribution of Z-Axis Mean Frequency Amplitude Exposure by Subject"
x = pro$FreqAmpZs
y = pro$Subject
xlab = "Mean Frequency Amplitude"
ylab = "Subject"
xlim=c(0, 4.2)

bwplot(y ~ x, data=pro, main=main, xlab=xlab, ylab=ylab, par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

summary(x)
sd(x)















main = "Lag Correlation Vector Sum A(8) and SDNN"
y = pro$SDANN
x = pro$VSAeq

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum A(8) and rMSSD"
y = pro$RMSSD

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum A(8) and pNN50"
y = pro$pNN50

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum A(8) and LF-HF"
y = pro$LF.HF

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum A(8) and LFnu"
y = pro$Lfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum A(8) and HFnu"
y = pro$Hfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


###############################################################################################################################################
###############################################################################################################################################


main = "Lag Correlation Vector Sum VDV(8) and SDNN"
y = pro$SDANN
x = pro$VSVDV8

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum VDV(8) and rMSSD"
y = pro$RMSSD

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum VDV(8) and pNN50"
y = pro$pNN50

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum VDV(8) and LF-HF"
y = pro$LF.HF

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum VDV(8) and LFnu"
y = pro$Lfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Vector Sum VDV(8) and HFnu"
y = pro$Hfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


###############################################################################################################################################
###############################################################################################################################################


main = "Lag Correlation Z-Axis Average Power and SDNN"
y = pro$SDANN
x = pro$MeanPowerZs

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Power and rMSSD"
y = pro$RMSSD

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Power and pNN50"
y = pro$pNN50

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Power and LF-HF"
y = pro$LF.HF

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Power and LFnu"
y = pro$Lfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Power and HFnu"
y = pro$Hfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


###############################################################################################################################################
###############################################################################################################################################


main = "Lag Correlation Z-Axis Average Frequency and SDNN"
y = pro$SDANN
x = pro$FreqPeakZs

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Frequency and rMSSD"
y = pro$RMSSD

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Frequency and pNN50"
y = pro$pNN50

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Frequency and LF-HF"
y = pro$LF.HF

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Frequency and LFnu"
y = pro$Lfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Frequency and HFnu"
y = pro$Hfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


###############################################################################################################################################
###############################################################################################################################################



main = "Lag Correlation Z-Axis Average Amplitude and SDNN"
y = pro$SDANN
x = pro$FreqAmpZs

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Amplitude and rMSSD"
y = pro$RMSSD

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Amplitude and pNN50"
y = pro$pNN50

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Amplitude and LF-HF"
y = pro$LF.HF

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Amplitude and LFnu"
y = pro$Lfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main = "Lag Correlation Z-Axis Average Amplitude and HFnu"
y = pro$Hfnu

ccf(x, y, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

###############################################################################################################################################
###############################################################################################################################################

###############################################################################################################################################
###############################################################################################################################################

###############################################################################################################################################
###############################################################################################################################################

###############################################################################################################################################
###############################################################################################################################################

###############################################################################################################################################
###############################################################################################################################################

###############################################################################################################################################
###############################################################################################################################################
ggplot(m72, aes(DateTime, SDANN)) + geom_line() + 
  xlab("Time of Day (24-Hour : Minute)") + ylab("SDNN (ms)") + ggtitle("Time Series of SDNN") + theme_economist() + scale_colour_economist() 
dev.copy(png, paste("Time Series of SDNN", "png", sep="."))
dev.off()
dev.off()
###############################################################################################################################################
###############################################################################################################################################
#proposal
###############################################################################################################################################
###############################################################################################################################################

main = "Distribution of Vehicle Speed and SDNN"
y = pro$SDANN
x = pro$Speed
#x2 = pro$VSAeqLag
#x3 = pro$VSAeqCum
#xlab = expression(A(8)~"m/" ~ s^{2})
xlab = "Vehicle Speed (km/h)"
ylab = "SDNN (ms)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
rp1 = vector('expression',1)
rp1[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]

#mod2 <- lm(y ~ x2, data = pro, na.action = na.omit)
#modsum2 <- summary(mod2)
#r22 = modsum2$adj.r.squared
#my.p2 = modsum2$coefficients[2,1]
#rp2 = vector('expression',1)
#rp2[1] = substitute(expression(Lag~A(8)~italic(R)^2 == MYVALUE), 
#                    list(MYVALUE = format(r22,dig=1)))[2]


#mod3 <- lm(y ~ x3, data = pro, na.action = na.omit)
#modsum3 <- summary(mod3)
#r23 = modsum3$adj.r.squared
#my.p3 = modsum3$coefficients[2,1]
#rp3 = vector('expression',1)
#rp3[1] = substitute(expression(Cum.~A(8)~italic(R)^2 == MYVALUE), 
#                    list(MYVALUE = format(r23,dig=1)))[2]


xyplot(y ~ x, data = pro,  type = c("p", "r"), main=main, xlab = xlab, ylab = ylab, par.settings = theEconomist.theme(), col="navy") #+ 
#  xyplot(y ~ x2, data = na.omit(pro),  type = c("p", "r"), par.settings = theEconomist.theme(), col="grey") +
#  xyplot(y ~ x3, data = na.omit(pro),  type = c("p", "r"), par.settings = theEconomist.theme(), col="black")
  
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp1, cex = 1.2, font = 2, col="navy")
#panel.text(0.3, 0.84, rp2, cex = 1.2, font = 2, col="grey")
#panel.text(0.3, 0.81, rp3, cex = 1.2, font = 2, col="black")
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()



###############################################################################################################################################

main = "Distribution of Vector Sum A(8) Exposures and rMSSD"
y = pro$RMSSD
ylab = "rMSSD (ms)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
rp1 = vector('expression',1)
rp1[1] = substitute(expression(A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]

mod2 <- lm(y ~ x2, data = pro, na.action = na.omit)
modsum2 <- summary(mod2)
r22 = modsum2$adj.r.squared
my.p2 = modsum2$coefficients[2,1]
rp2 = vector('expression',1)
rp2[1] = substitute(expression(Lag~A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r22,dig=1)))[2]


mod3 <- lm(y ~ x3, data = pro, na.action = na.omit)
modsum3 <- summary(mod3)
r23 = modsum3$adj.r.squared
my.p3 = modsum3$coefficients[2,1]
rp3 = vector('expression',1)
rp3[1] = substitute(expression(Cum.~A(8)~italic(R)^2 == MYVALUE), 
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

main = "Distribution of Vector Sum A(8) Exposures and pNN50"
y = pro$pNN50
ylab = "pNN50 (%)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
rp1 = vector('expression',1)
rp1[1] = substitute(expression(A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]

mod2 <- lm(y ~ x2, data = pro, na.action = na.omit)
modsum2 <- summary(mod2)
r22 = modsum2$adj.r.squared
my.p2 = modsum2$coefficients[2,1]
rp2 = vector('expression',1)
rp2[1] = substitute(expression(Lag~A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r22,dig=1)))[2]


mod3 <- lm(y ~ x3, data = pro, na.action = na.omit)
modsum3 <- summary(mod3)
r23 = modsum3$adj.r.squared
my.p3 = modsum3$coefficients[2,1]
rp3 = vector('expression',1)
rp3[1] = substitute(expression(Cum.~A(8)~italic(R)^2 == MYVALUE), 
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

main = "Distribution of Vector Sum A(8) Exposures and LF-HF"
y = pro$LF.HF
ylab = "LF-HF"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
rp1 = vector('expression',1)
rp1[1] = substitute(expression(A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]

mod2 <- lm(y ~ x2, data = pro, na.action = na.omit)
modsum2 <- summary(mod2)
r22 = modsum2$adj.r.squared
my.p2 = modsum2$coefficients[2,1]
rp2 = vector('expression',1)
rp2[1] = substitute(expression(Lag~A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r22,dig=1)))[2]


mod3 <- lm(y ~ x3, data = pro, na.action = na.omit)
modsum3 <- summary(mod3)
r23 = modsum3$adj.r.squared
my.p3 = modsum3$coefficients[2,1]
rp3 = vector('expression',1)
rp3[1] = substitute(expression(Cum.~A(8)~italic(R)^2 == MYVALUE), 
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


main = "Distribution of Vector Sum A(8) Exposures and LFnu"
y = pro$Lfnu
ylab = "LFnu"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
rp1 = vector('expression',1)
rp1[1] = substitute(expression(A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]

mod2 <- lm(y ~ x2, data = pro, na.action = na.omit)
modsum2 <- summary(mod2)
r22 = modsum2$adj.r.squared
my.p2 = modsum2$coefficients[2,1]
rp2 = vector('expression',1)
rp2[1] = substitute(expression(Lag~A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r22,dig=1)))[2]


mod3 <- lm(y ~ x3, data = pro, na.action = na.omit)
modsum3 <- summary(mod3)
r23 = modsum3$adj.r.squared
my.p3 = modsum3$coefficients[2,1]
rp3 = vector('expression',1)
rp3[1] = substitute(expression(Cum.~A(8)~italic(R)^2 == MYVALUE), 
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


main = "Distribution of Vector Sum A(8) Exposures and HFnu"
y = pro$Hfnu
ylab = "HFnu"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
rp1 = vector('expression',1)
rp1[1] = substitute(expression(A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]

mod2 <- lm(y ~ x2, data = pro, na.action = na.omit)
modsum2 <- summary(mod2)
r22 = modsum2$adj.r.squared
my.p2 = modsum2$coefficients[2,1]
rp2 = vector('expression',1)
rp2[1] = substitute(expression(Lag~A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r22,dig=1)))[2]


mod3 <- lm(y ~ x3, data = pro, na.action = na.omit)
modsum3 <- summary(mod3)
r23 = modsum3$adj.r.squared
my.p3 = modsum3$coefficients[2,1]
rp3 = vector('expression',1)
rp3[1] = substitute(expression(Cum.~A(8)~italic(R)^2 == MYVALUE), 
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

###############################################################################################################################################
###############################################################################################################################################

main = "Distribution of Vector Sum VDV(8) Exposures and SDNN"
y = pro$SDANN
x = pro$VSVDV8
x2 = pro$VSVDV8Lag
x3 = pro$VSVDV8Cum
xlab = expression(VDV(8)~"m/" ~ s^{1.75})
ylab = "SDNN (ms)"

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




main = "Distribution of Vector Sum VDV(8) Exposures and pNN50"
y = pro$pNN50
ylab = "pNN50 (%)"

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


main = "Distribution of Vector Sum VDV(8) Exposures and LF-HF"
y = pro$LF.HF
ylab = "LF-HF"

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




main = "Distribution of Vector Sum VDV(8) Exposures and LFnu"
y = pro$Lfnu
ylab = "LFnu"

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



main = "Distribution of Vector Sum VDV(8) Exposures and HFnu"
y = pro$Hfnu
ylab = "HFnu"

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

###############################################################################################################################################
###############################################################################################################################################

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
x = pro$FreqAmpZs
xlab = "Mean Frequency (Hz)"

main = "Distribution of Z-Axis Mean Amplitude Exposures and SDNN"
y = pro$SDANN
ylab = "SDNN (ms)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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


###############################################################################################################################################


main = "Distribution of Z-Axis Mean Amplitude Exposures and pNN50"
y = pro$pNN50
ylab = "pNN50 (%)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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

###############################################################################################################################################

main = "Distribution of Z-Axis Mean Amplitude Exposures and LF-HF"
y = pro$LF.HF
ylab = "LF-HF"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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

###############################################################################################################################################

main = "Distribution of Z-Axis Mean Amplitude Exposures and LFnu"
y = pro$Lfnu
ylab = "LFnu"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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


###############################################################################################################################################

main = "Distribution of Z-Axis Mean Amplitude Exposures and HFnu"
y = pro$Hfnu
ylab = "HFnu"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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


###############################################################################################################################################
###############################################################################################################################################

x = pro$FreqPeakZs
xlab = "Mean Frequency (Hz)"

main = "Distribution of Z-Axis Mean Frequency Exposure and SDNN"
y = pro$SDANN
ylab = "SDNN (ms)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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



###############################################################################################################################################



main = "Distribution of Z-Axis Mean Frequency Exposures and rMSSD"
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


###############################################################################################################################################




main = "Distribution of Z-Axis Mean Frequency Exposures and pNN50"
y = pro$pNN50
ylab = "pNN50 (%)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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

###############################################################################################################################################




main = "Distribution of Z-Axis Mean Frequency Exposures and LF-HF"
y = pro$LF.HF
ylab = "LF-HF"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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

###############################################################################################################################################




main = "Distribution of Z-Axis Mean Frequency Exposures and LFnu"
y = pro$Lfnu
ylab = "LFnu"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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


###############################################################################################################################################

main = "Distribution of Z-Axis Mean Frequency Exposures and HFnu"
y = pro$Hfnu
ylab = "HFnu"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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
###############################################################################################################################################
###############################################################################################################################################
main = "Distribution of Z-Axis Mean Power Exposures and SDNN"
y = pro$SDANN
x = pro$MeanPowerZs
xlab = "Mean Power"
ylab = "SDNN (ms)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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



###############################################################################################################################################



main = "Distribution of Z-Axis Mean Power Exposures and rMSSD"
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


###############################################################################################################################################




main = "Distribution of Z-Axis Mean Power Exposures and pNN50"
y = pro$pNN50
ylab = "pNN50 (%)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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

###############################################################################################################################################




main = "Distribution of Z-Axis Mean Power Exposures and LF-HF"
y = pro$LF.HF
ylab = "LF-HF"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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

###############################################################################################################################################




main = "Distribution of Z-Axis Mean Power Exposures and LFnu"
y = pro$Lfnu
ylab = "LFnu"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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


###############################################################################################################################################

main = "Distribution of Z-Axis Mean Power Exposures and HFnu"
y = pro$Hfnu
ylab = "HFnu"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
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

###############################################################################################################################################
###############################################################################################################################################


###############################################################################################################################################
###############################################################################################################################################



main = "Distribution of Vector Sum A(8) Exposures and HFnu"
y = pro$Hfnu
ylab = "HFnu"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
rp1 = vector('expression',1)
rp1[1] = substitute(expression(A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]

mod2 <- lm(y ~ x2, data = pro, na.action = na.omit)
modsum2 <- summary(mod2)
r22 = modsum2$adj.r.squared
my.p2 = modsum2$coefficients[2,1]
rp2 = vector('expression',1)
rp2[1] = substitute(expression(Lag~A(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r22,dig=1)))[2]


mod3 <- lm(y ~ x3, data = pro, na.action = na.omit)
modsum3 <- summary(mod3)
r23 = modsum3$adj.r.squared
my.p3 = modsum3$coefficients[2,1]
rp3 = vector('expression',1)
rp3[1] = substitute(expression(Cum.~A(8)~italic(R)^2 == MYVALUE), 
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
#Proposal Exposure Correlations
###############################################################################################################################################
###############################################################################################################################################

main = "Distribution of Log PM2.5 and Noise Exposures"
y = pro$Sound
x = log(pro$Air)

#xlab = expression(VDV(8)~"m/" ~ s^{1.75})
xlab = expression(Log~PM2.5~(~mu~"g/" ~ m^{3}))
ylab = "Sound (dBA)"

mod1 <- lm(y ~ x, data = pro, na.action = na.pass)
modsum1 <- summary(mod1)
r21 = modsum1$adj.r.squared
my.p1 = modsum1$coefficients[2,1]
rp1 = vector('expression',1)
rp1[1] = substitute(expression(VDV(8)~italic(R)^2 == MYVALUE), 
                    list(MYVALUE = format(r21,dig=1)))[2]


xyplot(y ~ x, data = pro,  type = c("p", "r"), main=main, xlab = xlab, ylab = ylab, par.settings = theEconomist.theme(), col="navy") 

trellis.focus("toplevel")
panel.text(0.3, 0.87, rp1, cex = 1.2, font = 2, col="navy")
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

###############################################################################################################################################
###############################################################################################################################################


###############################################################################################################################################
###############################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ YsAeq, data=pro, main="Distribution of Y-Aeq Exposures", xlab=expression("m/" ~ s^{2}), ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(YsAeq ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Y-Aeq Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("Average Weighted Vibration (m/" ~ s^{2}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeq, RMSSD, lag.max=10, main="Lag Correlation Y-Aeq and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeq, SDANN, lag.max=10, main="Lag Correlation Y-Aeq and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeq, pNN50, lag.max=10, main="Lag Correlation Y-Aeq and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeq, Lfnu, lag.max=10, main="Lag Correlation Y-Aeq and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeq, Hfnu, lag.max=10, main="Lag Correlation Y-Aeq and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeq, LF.HF, lag.max=10, main="Lag Correlation Y-Aeq and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeqLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-Aeq and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeqLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-Aeq and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeqLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-Aeq and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeqLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-Aeq and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeqLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-Aeq and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsAeqLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-Aeq and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (YsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ YsAeq, data = pro,  type = c("p", "r"), main="Distribution of Y-Aeq Exposures and RMSSD", xlab = expression("m/" ~ s^{2}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (YsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ YsAeq, data = pro,  type = c("p", "r"), main="Distribution of Y-Aeq Exposures and SDANN", xlab = expression("m/" ~ s^{2}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (YsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ YsAeq, data = pro,  type = c("p", "r"), main="Distribution of Y-Aeq Exposures and pNN50", xlab = expression("m/" ~ s^{2}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (YsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ YsAeq, data = pro,  type = c("p", "r"), main="Distribution of Y-Aeq Exposures and LFnu", xlab = expression("m/" ~ s^{2}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (YsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ YsAeq, data = pro,  type = c("p", "r"), main="Distribution of Y-Aeq Exposures and HFnu", xlab = expression("m/" ~ s^{2}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (YsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ YsAeq, data = pro,  type = c("p", "r"), main="Distribution of Y-Aeq Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{2}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (YsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ YsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-Aeq and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (YsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ YsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-Aeq and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (YsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ YsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-Aeq and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (YsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ YsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-Aeq and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (YsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ YsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-Aeq and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (YsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ YsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-Aeq and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

#proposal vibration
##########################################################################################################################################################


bwplot(Subject ~ ZsAeq, data=pro, main="Distribution of Z-Aeq Exposures", xlab=expression("m/" ~ s^{2}), par.settings = theEconomist.theme(), xlim=c(0, 2))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Distribution of  X-Axis VDV(8) Exposures"
bwplot(Subject ~ XsVDV8, data=pro, main=main, xlab=expression("m/" ~ s^{1.75}), par.settings = theEconomist.theme(), xlim=c(0, 4))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Time Series of Vector Sum VDV(8) Exposures by Subject"
xyplot(VSVDV8 ~ Time.Index | Subject, data = pro, type = c("b"), 
       main=main, par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("VDV(8) (m/" ~ s^{1.75}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Lag Correlation Vector Sum VDV(8) and pNN50"
ccf(pro$VSVDV8, pro$pNN50, lag.max=10, main=main,  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeq, SDANN, lag.max=10, main="Lag Correlation Z-Aeq and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeq, pNN50, lag.max=10, main="Lag Correlation Z-Aeq and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeq, Lfnu, lag.max=10, main="Lag Correlation Z-Aeq and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeq, Hfnu, lag.max=10, main="Lag Correlation Z-Aeq and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeq, LF.HF, lag.max=10, main="Lag Correlation Z-Aeq and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeqLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-Aeq and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeqLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-Aeq and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeqLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-Aeq and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeqLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-Aeq and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeqLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-Aeq and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsAeqLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-Aeq and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (ZsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ ZsAeq, data = pro,  type = c("p", "r"), main="Distribution of Z-Aeq Exposures and RMSSD", xlab = expression("m/" ~ s^{2}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (ZsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ ZsAeq, data = pro,  type = c("p", "r"), main="Distribution of Z-Aeq Exposures and SDANN", xlab = expression("m/" ~ s^{2}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (ZsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ ZsAeq, data = pro,  type = c("p", "r"), main="Distribution of Z-Aeq Exposures and pNN50", xlab = expression("m/" ~ s^{2}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (ZsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ ZsAeq, data = pro,  type = c("p", "r"), main="Distribution of Z-Aeq Exposures and LFnu", xlab = expression("m/" ~ s^{2}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (ZsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ ZsAeq, data = pro,  type = c("p", "r"), main="Distribution of Z-Aeq Exposures and HFnu", xlab = expression("m/" ~ s^{2}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (ZsAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ ZsAeq, data = pro,  type = c("p", "r"), main="Distribution of Z-Aeq Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{2}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (ZsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ ZsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-Aeq and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (ZsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ ZsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-Aeq and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (ZsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ ZsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-Aeq and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (ZsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ ZsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-Aeq and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (ZsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ ZsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-Aeq and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (ZsAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ ZsAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-Aeq and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

##########################################################################################################################################################

##########################################################################################################################################################
##pro$VSAeq[pro$VSAeq>4]<-.5

bwplot(Subject ~ VSAeq, data=pro, main="Distribution of Vector Sum Aeq Exposures", xlab=expression("m/" ~ s^{2}), ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(VSAeq ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Vector Sum Aeq Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("Cumulative Impulsive Vibration (m/" ~ s^{2}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeq, RMSSD, lag.max=10, main="Lag Correlation Vector Sum Aeq and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeq, SDANN, lag.max=10, main="Lag Correlation Vector Sum Aeq and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeq, pNN50, lag.max=10, main="Lag Correlation Vector Sum Aeq and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeq, Lfnu, lag.max=10, main="Lag Correlation Vector Sum Aeq and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeq, Hfnu, lag.max=10, main="Lag Correlation Vector Sum Aeq and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeq, LF.HF, lag.max=10, main="Lag Correlation Vector Sum Aeq and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeqLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum Aeq and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeqLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum Aeq and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeqLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum Aeq and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeqLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum Aeq and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeqLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum Aeq and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSAeqLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum Aeq and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (VSAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ VSAeq, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum Aeq Exposures and RMSSD", xlab = expression("m/" ~ s^{2}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (VSAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ VSAeq, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum Aeq Exposures and SDANN", xlab = expression("m/" ~ s^{2}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (VSAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ VSAeq, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum Aeq Exposures and pNN50", xlab = expression("m/" ~ s^{2}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (VSAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ VSAeq, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum Aeq Exposures and LFnu", xlab = expression("m/" ~ s^{2}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (VSAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ VSAeq, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum Aeq Exposures and HFnu", xlab = expression("m/" ~ s^{2}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (VSAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ VSAeq, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum Aeq Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{2}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (VSAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ VSAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum Aeq and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (VSAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ VSAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum Aeq and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (VSAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ VSAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum Aeq and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (VSAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ VSAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum Aeq and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (VSAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ VSAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum Aeq and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (VSAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ VSAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum Aeq and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ ZfAeq, data=pro, main="Distribution of Zf-Aeq Exposures", xlab=expression("m/" ~ s^{2}), ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(ZfAeq ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Zf-Aeq Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("Average Weighted Vibration (m/" ~ s^{2}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeq, RMSSD, lag.max=10, main="Lag Correlation Zf-Aeq and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeq, SDANN, lag.max=10, main="Lag Correlation Zf-Aeq and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeq, pNN50, lag.max=10, main="Lag Correlation Zf-Aeq and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeq, Lfnu, lag.max=10, main="Lag Correlation Zf-Aeq and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeq, Hfnu, lag.max=10, main="Lag Correlation Zf-Aeq and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeq, LF.HF, lag.max=10, main="Lag Correlation Zf-Aeq and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeqLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-Aeq and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeqLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-Aeq and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeqLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-Aeq and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeqLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-Aeq and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeqLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-Aeq and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfAeqLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-Aeq and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (ZfAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ ZfAeq, data = pro,  type = c("p", "r"), main="Distribution of Zf-Aeq Exposures and RMSSD", xlab = expression("m/" ~ s^{2}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (ZfAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ ZfAeq, data = pro,  type = c("p", "r"), main="Distribution of Zf-Aeq Exposures and SDANN", xlab = expression("m/" ~ s^{2}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (ZfAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ ZfAeq, data = pro,  type = c("p", "r"), main="Distribution of Zf-Aeq Exposures and pNN50", xlab = expression("m/" ~ s^{2}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (ZfAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ ZfAeq, data = pro,  type = c("p", "r"), main="Distribution of Zf-Aeq Exposures and LFnu", xlab = expression("m/" ~ s^{2}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (ZfAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ ZfAeq, data = pro,  type = c("p", "r"), main="Distribution of Zf-Aeq Exposures and HFnu", xlab = expression("m/" ~ s^{2}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (ZfAeq), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ ZfAeq, data = pro,  type = c("p", "r"), main="Distribution of Zf-Aeq Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{2}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (ZfAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ ZfAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-Aeq and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (ZfAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ ZfAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-Aeq and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (ZfAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ ZfAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-Aeq and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (ZfAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ ZfAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-Aeq and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (ZfAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ ZfAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-Aeq and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (ZfAeq) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ ZfAeq | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-Aeq and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{2}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################



bwplot(Subject ~ MeanPowerXs, data=pro, main="Distribution of X-VDV8 Exposures", xlab=expression("m/" ~ s^{1.75}), ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(MeanPowerXs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of X-VDV8 Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("Cumulative Impulsive Vibration (m/" ~ s^{1.75}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, RMSSD, lag.max=10, main="Lag Correlation X-VDV8 and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8, SDANN, lag.max=10, main="Lag Correlation X-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8, pNN50, lag.max=10, main="Lag Correlation X-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8, Lfnu, lag.max=10, main="Lag Correlation X-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8, Hfnu, lag.max=10, main="Lag Correlation X-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8, LF.HF, lag.max=10, main="Lag Correlation X-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8Lag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged X-VDV8 and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8Lag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged X-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8Lag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged X-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8Lag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged X-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8Lag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged X-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(XsVDV8Lag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged X-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (XsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ XsVDV8, data = pro,  type = c("p", "r"), main="Distribution of X-VDV8 Exposures and RMSSD", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (XsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ XsVDV8, data = pro,  type = c("p", "r"), main="Distribution of X-VDV8 Exposures and SDANN", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (XsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ XsVDV8, data = pro,  type = c("p", "r"), main="Distribution of X-VDV8 Exposures and pNN50", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (XsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ XsVDV8, data = pro,  type = c("p", "r"), main="Distribution of X-VDV8 Exposures and LFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (XsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ XsVDV8, data = pro,  type = c("p", "r"), main="Distribution of X-VDV8 Exposures and HFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (XsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ XsVDV8, data = pro,  type = c("p", "r"), main="Distribution of X-VDV8 Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (XsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ XsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of X-VDV8 and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (XsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ XsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of X-VDV8 and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (XsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ XsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of X-VDV8 and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (XsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ XsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of X-VDV8 and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (XsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ XsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of X-VDV8 and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (XsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ XsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of X-VDV8 and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ YsVDV8, data=pro, main="Distribution of Y-VDV8 Exposures", xlab=expression("m/" ~ s^{1.75}), ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(YsVDV8 ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Y-VDV8 Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("Cumulative Impulsive Vibration (m/" ~ s^{1.75}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8, RMSSD, lag.max=10, main="Lag Correlation Y-VDV8 and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8, SDANN, lag.max=10, main="Lag Correlation Y-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8, pNN50, lag.max=10, main="Lag Correlation Y-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8, Lfnu, lag.max=10, main="Lag Correlation Y-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8, Hfnu, lag.max=10, main="Lag Correlation Y-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8, LF.HF, lag.max=10, main="Lag Correlation Y-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8Lag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-VDV8 and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8Lag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8Lag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8Lag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8Lag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(YsVDV8Lag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Y-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (YsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ YsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Y-VDV8 Exposures and RMSSD", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (YsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ YsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Y-VDV8 Exposures and SDANN", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (YsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ YsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Y-VDV8 Exposures and pNN50", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (YsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ YsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Y-VDV8 Exposures and LFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (YsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ YsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Y-VDV8 Exposures and HFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (YsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ YsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Y-VDV8 Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (YsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ YsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-VDV8 and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (YsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ YsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-VDV8 and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (YsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ YsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-VDV8 and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (YsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ YsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-VDV8 and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (YsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ YsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-VDV8 and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (YsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ YsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Y-VDV8 and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ ZsVDV8, data=pro, main="Distribution of Z-VDV8 Exposures", xlab=expression("m/" ~ s^{1.75}), par.settings = theEconomist.theme(), xlim=c(0, 3.5))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(ZsVDV8 ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Z-VDV8 Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("Cumulative Impulsive Vibration (m/" ~ s^{1.75}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8, RMSSD, lag.max=10, main="Lag Correlation Z-VDV8 and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8, SDANN, lag.max=10, main="Lag Correlation Z-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8, pNN50, lag.max=10, main="Lag Correlation Z-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8, Lfnu, lag.max=10, main="Lag Correlation Z-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8, Hfnu, lag.max=10, main="Lag Correlation Z-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8, LF.HF, lag.max=10, main="Lag Correlation Z-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8Lag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-VDV8 and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8Lag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8Lag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8Lag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8Lag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZsVDV8Lag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Z-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (ZsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ ZsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Z-VDV8 Exposures and RMSSD", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (ZsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ ZsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Z-VDV8 Exposures and SDANN", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (ZsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ ZsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Z-VDV8 Exposures and pNN50", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (ZsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ ZsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Z-VDV8 Exposures and LFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (ZsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ ZsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Z-VDV8 Exposures and HFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (ZsVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ ZsVDV8, data = pro,  type = c("p", "r"), main="Distribution of Z-VDV8 Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (ZsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ ZsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-VDV8 and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (ZsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ ZsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-VDV8 and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (ZsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ ZsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-VDV8 and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (ZsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ ZsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-VDV8 and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (ZsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ ZsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-VDV8 and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (ZsVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ ZsVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Z-VDV8 and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################
##pro$VSVDV8[pro$VSVDV8>20]<-2.1

bwplot(Subject ~ VSVDV8, data=pro, main="Distribution of Vector Sum VDV8 Exposures", xlab=expression("m/" ~ s^{1.75}), par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(VSVDV8 ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Vector Sum VDV8 Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("Cumulative Impulsive Vibration (m/" ~ s^{1.75}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8, RMSSD, lag.max=10, main="Lag Correlation Vector Sum VDV8 and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8, SDANN, lag.max=10, main="Lag Correlation Vector Sum VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8, pNN50, lag.max=10, main="Lag Correlation Vector Sum VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8, Lfnu, lag.max=10, main="Lag Correlation Vector Sum VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8, Hfnu, lag.max=10, main="Lag Correlation Vector Sum VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8, LF.HF, lag.max=10, main="Lag Correlation Vector Sum VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8Lag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum VDV8 and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8Lag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8Lag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8Lag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8Lag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(VSVDV8Lag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Vector Sum VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (VSVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ VSVDV8, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum VDV8 Exposures and RMSSD", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (VSVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ VSVDV8, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum VDV8 Exposures and SDANN", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (VSVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ VSVDV8, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum VDV8 Exposures and pNN50", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (VSVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ VSVDV8, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum VDV8 Exposures and LFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (VSVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ VSVDV8, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum VDV8 Exposures and HFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (VSVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ VSVDV8, data = pro,  type = c("p", "r"), main="Distribution of Vector Sum VDV8 Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (VSVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ VSVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum VDV8 and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (VSVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ VSVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum VDV8 and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (VSVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ VSVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum VDV8 and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (VSVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ VSVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum VDV8 and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (VSVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ VSVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum VDV8 and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (VSVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ VSVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Vector Sum VDV8 and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################



bwplot(Subject ~ ZfVDV8, data=pro, main="Distribution of Zf-VDV8 Exposures", xlab=expression("m/" ~ s^{1.75}), ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(ZfVDV8 ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Zf-VDV8 Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("Cumulative Impulsive Vibration (m/" ~ s^{1.75}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, RMSSD, lag.max=10, main="Lag Correlation Zf-VDV8 and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, SDANN, lag.max=10, main="Lag Correlation Zf-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, pNN50, lag.max=10, main="Lag Correlation Zf-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, Lfnu, lag.max=10, main="Lag Correlation Zf-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, Hfnu, lag.max=10, main="Lag Correlation Zf-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, LF.HF, lag.max=10, main="Lag Correlation Zf-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and RMSSD", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and SDANN", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and pNN50", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and LFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and HFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()
bwplot(Subject ~ ZfVDV8, data=pro, main="Distribution of Zf-VDV8 Exposures", xlab=expression("m/" ~ s^{1.75}), ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(ZfVDV8 ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Zf-VDV8 Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = expression("Cumulative Impulsive Vibration (m/" ~ s^{1.75}~")"))
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, RMSSD, lag.max=10, main="Lag Correlation Zf-VDV8 and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, SDANN, lag.max=10, main="Lag Correlation Zf-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, pNN50, lag.max=10, main="Lag Correlation Zf-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, Lfnu, lag.max=10, main="Lag Correlation Zf-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, Hfnu, lag.max=10, main="Lag Correlation Zf-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8, LF.HF, lag.max=10, main="Lag Correlation Zf-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(ZfVDV8Lag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Zf-VDV8 and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and RMSSD", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and SDANN", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and pNN50", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and LFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and HFnu", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (ZfVDV8), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ ZfVDV8, data = pro,  type = c("p", "r"), main="Distribution of Zf-VDV8 Exposures and LF:HF Ratio", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and RMSSD by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and SDANN by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and pNN50 by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and LFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and HFnu by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (ZfVDV8) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ ZfVDV8 | Category, data = pro,  type = c("p", "r"), main="Distributions of Zf-VDV8 and LF:HF by Vehicle Speed", xlab = expression("m/" ~ s^{1.75}), ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################



bwplot(Subject ~ MeanPowerXs, data=pro, main="Distribution of Mean X Power Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(MeanPowerXs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean X Power Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, RMSSD, lag.max=10, main="Lag Correlation Mean X Power and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, SDANN, lag.max=10, main="Lag Correlation Mean X Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, pNN50, lag.max=10, main="Lag Correlation Mean X Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, Lfnu, lag.max=10, main="Lag Correlation Mean X Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, Hfnu, lag.max=10, main="Lag Correlation Mean X Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, LF.HF, lag.max=10, main="Lag Correlation Mean X Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ MeanPowerYs, data=pro, main="Distribution of Mean Y Frequency Peak (Hz) Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(MeanPowerYs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Y Power Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, RMSSD, lag.max=10, main="Lag Correlation Mean Y Power and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, SDANN, lag.max=10, main="Lag Correlation Mean Y Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, pNN50, lag.max=10, main="Lag Correlation Mean Y Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, Lfnu, lag.max=10, main="Lag Correlation Mean Y Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, Hfnu, lag.max=10, main="Lag Correlation Mean Y Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, LF.HF, lag.max=10, main="Lag Correlation Mean Y Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ MeanPowerZs, data=pro, main="Distribution of Mean Z Power Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(MeanPowerZs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Z Power Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, RMSSD, lag.max=10, main="Lag Correlation Mean Z Power and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, SDANN, lag.max=10, main="Lag Correlation Mean Z Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, pNN50, lag.max=10, main="Lag Correlation Mean Z Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, Lfnu, lag.max=10, main="Lag Correlation Mean Z Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, Hfnu, lag.max=10, main="Lag Correlation Mean Z Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, LF.HF, lag.max=10, main="Lag Correlation Mean Z Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################
###pro$MeanPowerZf[pro$MeanPowerZf>20000]<-925


bwplot(Subject ~ MeanPowerZf, data=pro, main="Distribution of Mean Floor Z Power Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(MeanPowerZf ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Floor Z Power Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, RMSSD, lag.max=10, main="Lag Correlation Mean Floor Z Power and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, SDANN, lag.max=10, main="Lag Correlation Mean Floor Z Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, pNN50, lag.max=10, main="Lag Correlation Mean Floor Z Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, Lfnu, lag.max=10, main="Lag Correlation Mean Floor Z Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, Hfnu, lag.max=10, main="Lag Correlation Mean Floor Z Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, LF.HF, lag.max=10, main="Lag Correlation Mean Floor Z Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()
bwplot(Subject ~ MeanPowerZf, data=pro, main="Distribution of Mean Floor Z Power Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################





##########################################################################################################################################################

##########################################################################################################################################################



bwplot(Subject ~ FreqPeakXs, data=pro, main="Distribution of Mean X Frequency Peak (Hz) Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(FreqPeakXs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean X Frequency Peak (Hz) Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXs, RMSSD, lag.max=10, main="Lag Correlation Mean X Frequency Peak (Hz)and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXs, SDANN, lag.max=10, main="Lag Correlation Mean X Frequency Peak (Hz)and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXs, pNN50, lag.max=10, main="Lag Correlation Mean X Frequency Peak (Hz)and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXs, Lfnu, lag.max=10, main="Lag Correlation Mean X Frequency Peak (Hz)and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXs, Hfnu, lag.max=10, main="Lag Correlation Mean X Frequency Peak (Hz)and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXs, LF.HF, lag.max=10, main="Lag Correlation Mean X Frequency Peak (Hz)and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Frequency Peak (Hz)and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Frequency Peak (Hz)and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Frequency Peak (Hz)and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Frequency Peak (Hz)and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Frequency Peak (Hz)and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakXsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Frequency Peak (Hz)and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (FreqPeakXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ FreqPeakXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Frequency Peak (Hz) Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (FreqPeakXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ FreqPeakXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Frequency Peak (Hz) Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (FreqPeakXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ FreqPeakXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Frequency Peak (Hz) Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (FreqPeakXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ FreqPeakXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Frequency Peak (Hz) Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (FreqPeakXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ FreqPeakXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Frequency Peak (Hz) Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (FreqPeakXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ FreqPeakXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Frequency Peak (Hz) Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (FreqPeakXs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ FreqPeakXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Frequency Peak (Hz)and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (FreqPeakXs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ FreqPeakXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Frequency Peak (Hz)and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (FreqPeakXs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ FreqPeakXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Frequency Peak (Hz)and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (FreqPeakXs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ FreqPeakXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Frequency Peak (Hz)and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (FreqPeakXs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ FreqPeakXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Frequency Peak (Hz)and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (FreqPeakXs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ FreqPeakXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Frequency Peak (Hz)and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ FreqPeakYs, data=pro, main="Distribution of Mean Y Frequency Peak (Hz) Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(FreqPeakYs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Y Frequency Peak (Hz) Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYs, RMSSD, lag.max=10, main="Lag Correlation Mean Y Frequency Peak (Hz)and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYs, SDANN, lag.max=10, main="Lag Correlation Mean Y Frequency Peak (Hz)and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYs, pNN50, lag.max=10, main="Lag Correlation Mean Y Frequency Peak (Hz)and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYs, Lfnu, lag.max=10, main="Lag Correlation Mean Y Frequency Peak (Hz)and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYs, Hfnu, lag.max=10, main="Lag Correlation Mean Y Frequency Peak (Hz)and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYs, LF.HF, lag.max=10, main="Lag Correlation Mean Y Frequency Peak (Hz)and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Frequency Peak (Hz)and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Frequency Peak (Hz)and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Frequency Peak (Hz)and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Frequency Peak (Hz)and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Frequency Peak (Hz)and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakYsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Frequency Peak (Hz)and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (FreqPeakYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ FreqPeakYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Frequency Peak (Hz) Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (FreqPeakYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ FreqPeakYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Frequency Peak (Hz) Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (FreqPeakYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ FreqPeakYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Frequency Peak (Hz) Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (FreqPeakYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ FreqPeakYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Frequency Peak (Hz) Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (FreqPeakYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ FreqPeakYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Frequency Peak (Hz) Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (FreqPeakYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ FreqPeakYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Frequency Peak (Hz) Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (FreqPeakYs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ FreqPeakYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Frequency Peak (Hz)and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (FreqPeakYs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ FreqPeakYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Frequency Peak (Hz)and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (FreqPeakYs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ FreqPeakYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Frequency Peak (Hz)and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (FreqPeakYs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ FreqPeakYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Frequency Peak (Hz)and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (FreqPeakYs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ FreqPeakYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Frequency Peak (Hz)and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (FreqPeakYs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ FreqPeakYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Frequency Peak (Hz)and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ FreqPeakZs, data=pro, main="Distribution of Mean Z Frequency Peak (Hz) Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(FreqPeakZs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Z Frequency Peak (Hz) Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZs, RMSSD, lag.max=10, main="Lag Correlation Mean Z Frequency Peak (Hz)and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZs, SDANN, lag.max=10, main="Lag Correlation Mean Z Frequency Peak (Hz)and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZs, pNN50, lag.max=10, main="Lag Correlation Mean Z Frequency Peak (Hz)and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZs, Lfnu, lag.max=10, main="Lag Correlation Mean Z Frequency Peak (Hz)and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZs, Hfnu, lag.max=10, main="Lag Correlation Mean Z Frequency Peak (Hz)and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZs, LF.HF, lag.max=10, main="Lag Correlation Mean Z Frequency Peak (Hz)and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Frequency Peak (Hz)and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Frequency Peak (Hz)and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Frequency Peak (Hz)and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Frequency Peak (Hz)and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Frequency Peak (Hz)and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Frequency Peak (Hz)and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (FreqPeakZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ FreqPeakZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Frequency Peak (Hz) Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (FreqPeakZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ FreqPeakZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Frequency Peak (Hz) Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (FreqPeakZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ FreqPeakZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Frequency Peak (Hz) Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (FreqPeakZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ FreqPeakZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Frequency Peak (Hz) Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (FreqPeakZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ FreqPeakZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Frequency Peak (Hz) Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (FreqPeakZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ FreqPeakZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Frequency Peak (Hz) Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (FreqPeakZs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ FreqPeakZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Frequency Peak (Hz)and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (FreqPeakZs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ FreqPeakZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Frequency Peak (Hz)and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (FreqPeakZs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ FreqPeakZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Frequency Peak (Hz)and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (FreqPeakZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ FreqPeakZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Frequency Peak (Hz)and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (FreqPeakZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ FreqPeakZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Frequency Peak (Hz)and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
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

##########################################################################################################################################################
###pro$FreqPeakZf[pro$FreqPeakZf>20000]<-925


bwplot(Subject ~ FreqPeakZf, data=pro, main="Distribution of Mean Floor Z Frequency Peak (Hz) Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(FreqPeakZf ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Floor Z Frequency Peak (Hz) Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZf, RMSSD, lag.max=10, main="Lag Correlation Mean Floor Z Frequency Peak (Hz)and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZf, SDANN, lag.max=10, main="Lag Correlation Mean Floor Z Frequency Peak (Hz)and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZf, pNN50, lag.max=10, main="Lag Correlation Mean Floor Z Frequency Peak (Hz)and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZf, Lfnu, lag.max=10, main="Lag Correlation Mean Floor Z Frequency Peak (Hz)and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZf, Hfnu, lag.max=10, main="Lag Correlation Mean Floor Z Frequency Peak (Hz)and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZf, LF.HF, lag.max=10, main="Lag Correlation Mean Floor Z Frequency Peak (Hz)and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZfLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Frequency Peak (Hz)and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZfLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Frequency Peak (Hz)and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZfLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Frequency Peak (Hz)and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZfLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Frequency Peak (Hz)and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZfLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Frequency Peak (Hz)and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqPeakZfLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Frequency Peak (Hz)and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (FreqPeakZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ FreqPeakZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Frequency Peak (Hz) Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (FreqPeakZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ FreqPeakZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Frequency Peak (Hz) Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (FreqPeakZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ FreqPeakZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Frequency Peak (Hz) Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (FreqPeakZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ FreqPeakZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Frequency Peak (Hz) Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (FreqPeakZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ FreqPeakZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Frequency Peak (Hz) Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (FreqPeakZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ FreqPeakZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Frequency Peak (Hz) Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (FreqPeakZf) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ FreqPeakZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Frequency Peak (Hz)and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (FreqPeakZf) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ FreqPeakZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Frequency Peak (Hz)and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (FreqPeakZf) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ FreqPeakZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Frequency Peak (Hz)and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (FreqPeakZf) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ FreqPeakZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Frequency Peak (Hz)and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (FreqPeakZf) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ FreqPeakZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Frequency Peak (Hz)and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (FreqPeakZf) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ FreqPeakZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Frequency Peak (Hz)and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()
bwplot(Subject ~ FreqPeakZf, data=pro, main="Distribution of Mean Floor Z Frequency Peak (Hz) Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################



bwplot(Subject ~ MeanPowerXs, data=pro, main="Distribution of Mean X Power Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(MeanPowerXs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean X Power Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, RMSSD, lag.max=10, main="Lag Correlation Mean X Power and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, SDANN, lag.max=10, main="Lag Correlation Mean X Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, pNN50, lag.max=10, main="Lag Correlation Mean X Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, Lfnu, lag.max=10, main="Lag Correlation Mean X Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, Hfnu, lag.max=10, main="Lag Correlation Mean X Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXs, LF.HF, lag.max=10, main="Lag Correlation Mean X Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerXsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (MeanPowerXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ MeanPowerXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Power Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (MeanPowerXs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ MeanPowerXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Power and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ MeanPowerYs, data=pro, main="Distribution of Mean Y Power Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(MeanPowerYs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Y Power Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, RMSSD, lag.max=10, main="Lag Correlation Mean Y Power and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, SDANN, lag.max=10, main="Lag Correlation Mean Y Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, pNN50, lag.max=10, main="Lag Correlation Mean Y Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, Lfnu, lag.max=10, main="Lag Correlation Mean Y Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, Hfnu, lag.max=10, main="Lag Correlation Mean Y Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYs, LF.HF, lag.max=10, main="Lag Correlation Mean Y Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerYsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (MeanPowerYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ MeanPowerYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Power Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (MeanPowerYs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ MeanPowerYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Power and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ MeanPowerZs, data=pro, main="Distribution of Mean Z Power Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(MeanPowerZs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Z Power Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, RMSSD, lag.max=10, main="Lag Correlation Mean Z Power and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, SDANN, lag.max=10, main="Lag Correlation Mean Z Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, pNN50, lag.max=10, main="Lag Correlation Mean Z Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, Lfnu, lag.max=10, main="Lag Correlation Mean Z Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, Hfnu, lag.max=10, main="Lag Correlation Mean Z Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZs, LF.HF, lag.max=10, main="Lag Correlation Mean Z Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (MeanPowerZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ MeanPowerZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Power Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (MeanPowerZs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ MeanPowerZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Power and LF:HF by Vehicle Speed", 
       xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################
###pro$MeanPowerZf[pro$MeanPowerZf>20000]<-925


bwplot(Subject ~ MeanPowerZf, data=pro, main="Distribution of Mean Floor Z Power Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(MeanPowerZf ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Floor Z Power Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Power"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, RMSSD, lag.max=10, main="Lag Correlation Mean Floor Z Power and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, SDANN, lag.max=10, main="Lag Correlation Mean Floor Z Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, pNN50, lag.max=10, main="Lag Correlation Mean Floor Z Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, Lfnu, lag.max=10, main="Lag Correlation Mean Floor Z Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, Hfnu, lag.max=10, main="Lag Correlation Mean Floor Z Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZf, LF.HF, lag.max=10, main="Lag Correlation Mean Floor Z Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(MeanPowerZfLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Power and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and RMSSD", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and SDANN", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and pNN50", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and LFnu", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and HFnu", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (MeanPowerZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ MeanPowerZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Power Exposures and LF:HF Ratio", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and RMSSD by Vehicle Speed", xlab = "Mean Power", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and SDANN by Vehicle Speed", xlab = "Mean Power", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and pNN50 by Vehicle Speed", xlab = "Mean Power", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and LFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and HFnu by Vehicle Speed", xlab = "Mean Power", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (MeanPowerZf) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ MeanPowerZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Power and LF:HF by Vehicle Speed", xlab = "Mean Power", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()
bwplot(Subject ~ MeanPowerZf, data=pro, main="Distribution of Mean Floor Z Power Exposures", xlab="Mean Power", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################



bwplot(Subject ~ FreqAmpXs, data=pro, main="Distribution of Mean X Amplitude Exposures", xlab="Mean Amplitude", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(FreqAmpXs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean X Amplitude Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Amplitude"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXs, RMSSD, lag.max=10, main="Lag Correlation Mean X Amplitude and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXs, SDANN, lag.max=10, main="Lag Correlation Mean X Amplitude and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXs, pNN50, lag.max=10, main="Lag Correlation Mean X Amplitude and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXs, Lfnu, lag.max=10, main="Lag Correlation Mean X Amplitude and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXs, Hfnu, lag.max=10, main="Lag Correlation Mean X Amplitude and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXs, LF.HF, lag.max=10, main="Lag Correlation Mean X Amplitude and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Amplitude and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Amplitude and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Amplitude and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Amplitude and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Amplitude and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpXsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean X Amplitude and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (FreqAmpXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ FreqAmpXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Amplitude Exposures and RMSSD", xlab = "Mean Amplitude", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (FreqAmpXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ FreqAmpXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Amplitude Exposures and SDANN", xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (FreqAmpXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ FreqAmpXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Amplitude Exposures and pNN50", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (FreqAmpXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ FreqAmpXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Amplitude Exposures and LFnu", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (FreqAmpXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ FreqAmpXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Amplitude Exposures and HFnu", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (FreqAmpXs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ FreqAmpXs, data = pro,  type = c("p", "r"), main="Distribution of Mean X Amplitude Exposures and LF:HF Ratio", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (FreqAmpXs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ FreqAmpXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Amplitude and RMSSD by Vehicle Speed", xlab = "Mean Amplitude", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (FreqAmpXs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ FreqAmpXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Amplitude and SDANN by Vehicle Speed", 
       xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (FreqAmpXs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ FreqAmpXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Amplitude and pNN50 by Vehicle Speed", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (FreqAmpXs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ FreqAmpXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Amplitude and LFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (FreqAmpXs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ FreqAmpXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Amplitude and HFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (FreqAmpXs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ FreqAmpXs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean X Amplitude and LF:HF by Vehicle Speed", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ FreqAmpYs, data=pro, main="Distribution of Mean Y Amplitude Exposures", xlab="Mean Amplitude", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(FreqAmpYs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Y Amplitude Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Amplitude"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYs, RMSSD, lag.max=10, main="Lag Correlation Mean Y Amplitude and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYs, SDANN, lag.max=10, main="Lag Correlation Mean Y Amplitude and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYs, pNN50, lag.max=10, main="Lag Correlation Mean Y Amplitude and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYs, Lfnu, lag.max=10, main="Lag Correlation Mean Y Amplitude and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYs, Hfnu, lag.max=10, main="Lag Correlation Mean Y Amplitude and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYs, LF.HF, lag.max=10, main="Lag Correlation Mean Y Amplitude and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Amplitude and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Amplitude and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Amplitude and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Amplitude and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Amplitude and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpYsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Y Amplitude and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (FreqAmpYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ FreqAmpYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Amplitude Exposures and RMSSD", xlab = "Mean Amplitude", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (FreqAmpYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ FreqAmpYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Amplitude Exposures and SDANN", xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (FreqAmpYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ FreqAmpYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Amplitude Exposures and pNN50", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (FreqAmpYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ FreqAmpYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Amplitude Exposures and LFnu", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (FreqAmpYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ FreqAmpYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Amplitude Exposures and HFnu", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (FreqAmpYs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ FreqAmpYs, data = pro,  type = c("p", "r"), main="Distribution of Mean Y Amplitude Exposures and LF:HF Ratio", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (FreqAmpYs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ FreqAmpYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Amplitude and RMSSD by Vehicle Speed", xlab = "Mean Amplitude", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (FreqAmpYs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ FreqAmpYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Amplitude and SDANN by Vehicle Speed", xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (FreqAmpYs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ FreqAmpYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Amplitude and pNN50 by Vehicle Speed", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (FreqAmpYs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ FreqAmpYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Amplitude and LFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (FreqAmpYs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ FreqAmpYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Amplitude and HFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (FreqAmpYs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ FreqAmpYs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Y Amplitude and LF:HF by Vehicle Speed", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################


bwplot(Subject ~ FreqAmpZs, data=pro, main="Distribution of Mean Z Amplitude Exposures", xlab="Mean Amplitude", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(FreqAmpZs ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Z Amplitude Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Amplitude"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZs, RMSSD, lag.max=10, main="Lag Correlation Mean Z Amplitude and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Lag Correlation Mean Z Amplitude and SDANN"
ccf(FreqAmpZs, SDANN, lag.max=10, main=main, xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


ccf(FreqAmpZs, pNN50, lag.max=10, main="Lag Correlation Mean Z Amplitude and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZs, Lfnu, lag.max=10, main="Lag Correlation Mean Z Amplitude and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZs, Hfnu, lag.max=10, main="Lag Correlation Mean Z Amplitude and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZs, LF.HF, lag.max=10, main="Lag Correlation Mean Z Amplitude and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZsLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Amplitude and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZsLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Amplitude and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZsLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Amplitude and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZsLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Amplitude and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZsLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Amplitude and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZsLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Z Amplitude and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

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

mod <- lm(SDANN ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and SDANN", xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and pNN50", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and LFnu", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and HFnu", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and LF:HF Ratio", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and RMSSD by Vehicle Speed", xlab = "Mean Amplitude", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and SDANN by Vehicle Speed", xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and pNN50 by Vehicle Speed", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and LFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and HFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and LF:HF by Vehicle Speed", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


##########################################################################################################################################################

##########################################################################################################################################################
###pro$FreqAmpZf[pro$FreqAmpZf>20000]<-925


bwplot(Subject ~ FreqAmpZf, data=pro, main="Distribution of Mean Floor Z Amplitude Exposures", xlab="Mean Amplitude", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


xyplot(FreqAmpZf ~ Time.Index | Subject, data = pro, type = c("b"), 
       main="Time Series of Mean Floor Z Amplitude Exposures by Subject", par.settings = theEconomist.theme(),
       xlab = "Time (5-Minute Intervals)", ylab = "Mean Amplitude"
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZf, RMSSD, lag.max=10, main="Lag Correlation Mean Floor Z Amplitude and RMSSD",  
    xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZf, SDANN, lag.max=10, main="Lag Correlation Mean Floor Z Amplitude and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZf, pNN50, lag.max=10, main="Lag Correlation Mean Floor Z Amplitude and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZf, Lfnu, lag.max=10, main="Lag Correlation Mean Floor Z Amplitude and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZf, Hfnu, lag.max=10, main="Lag Correlation Mean Floor Z Amplitude and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZf, LF.HF, lag.max=10, main="Lag Correlation Mean Floor Z Amplitude and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZfLag, RMSSD, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Amplitude and RMSSD", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZfLag, SDANN, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Amplitude and SDANN", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZfLag, pNN50, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Amplitude and pNN50", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZfLag, Lfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Amplitude and LFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZfLag, Hfnu, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Amplitude and HFnu", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

ccf(FreqAmpZfLag, LF.HF, lag.max=10, na.action = na.pass, main="Lag Correlation Lagged Mean Floor Z Amplitude and LF:HF", xlab="Lag (5-Minute Multiples)", ylab="Correlation Significance", cex.main=0.8, cex = 0.8, cex.lab=0.9, cex.axis=0.8, family="serif")
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(RMSSD ~ (FreqAmpZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(RMSSD ~ FreqAmpZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Amplitude Exposures and RMSSD", xlab = "Mean Amplitude", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.28, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(SDANN ~ (FreqAmpZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ FreqAmpZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Amplitude Exposures and SDANN", xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (FreqAmpZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ FreqAmpZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Amplitude Exposures and pNN50", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (FreqAmpZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ FreqAmpZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Amplitude Exposures and LFnu", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (FreqAmpZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ FreqAmpZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Amplitude Exposures and HFnu", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (FreqAmpZf), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ FreqAmpZf, data = pro,  type = c("p", "r"), main="Distribution of Mean Floor Z Amplitude Exposures and LF:HF Ratio", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (FreqAmpZf) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ FreqAmpZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Amplitude and RMSSD by Vehicle Speed", xlab = "Mean Amplitude", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (FreqAmpZf) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ FreqAmpZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Amplitude and SDANN by Vehicle Speed", xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (FreqAmpZf) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ FreqAmpZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Amplitude and pNN50 by Vehicle Speed", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (FreqAmpZf) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ FreqAmpZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Amplitude and LFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (FreqAmpZf) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ FreqAmpZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Amplitude and HFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (FreqAmpZf) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ FreqAmpZf | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Floor Z Amplitude and LF:HF by Vehicle Speed", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()
bwplot(Subject ~ FreqAmpZf, data=pro, main="Distribution of Mean Floor Z Amplitude Exposures", xlab="Mean Amplitude", ylab="Subject", par.settings = theEconomist.theme())
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

#############################################################################################################################################################
####################################################################################################################################################
######################################################################################################################################################
####################################################################################################################################################################
###################################################################################################################################################
########################################################################################################################################################
##########################################################################################################################################################








#####################################################################################################################################################################
###################################################################################################################################################
#############################################################################################################################################################
############################################################################################################################################################
########################################################################################################################################################
###########################################################################################################################################################
#####pro$SDANN[pro$SDANN>200] <- 175

xyplot(
  SDANN ~ log(Air),
  data = pro, par.settings = theEconomist.theme(),
  main="Log Fine Particulate Matter and SDANN by Subject", xlab="Log of Fine Particulate Matter", ylab="SDNN (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  RMSSD ~ log(Air),
  data = pro, par.settings = theEconomist.theme(),
  main="Log Fine Particulate Matter and RMSSD by Subject", xlab="Log of Fine Particulate Matter", ylab="RMSSD (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  pNN50 ~ log(Air),
  data = pro, par.settings = theEconomist.theme(),
  main="Log Fine Particulate Matter and pNN50 by Subject", xlab="Log of Fine Particulate Matter", ylab="pNN50",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  Lfnu ~ log(Air),
  data = pro, par.settings = theEconomist.theme(),
  main="Log Fine Particulate Matter and LFnu by Subject", xlab="Log of Fine Particulate Matter", ylab="LFnu",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  LF.HF ~ log(Air),
  data = pro, par.settings = theEconomist.theme(),
  main="Log Fine Particulate Matter and LF:HF by Subject", xlab="Log of Fine Particulate Matter", ylab="LF:HF",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)





xyplot(
  SDANN ~ Sound,
  data = pro, par.settings = theEconomist.theme(),
  main="Sound Exposure and SDANN by Subject", xlab="Sound (dBA)", ylab="SDNN (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  RMSSD ~ Sound,
  data = pro, par.settings = theEconomist.theme(),
  main="Sound Exposure and RMSSD by Subject", xlab="Sound (dBA)", ylab="RMSSD (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  pNN50 ~ Sound,
  data = pro, par.settings = theEconomist.theme(),
  main="Sound Exposure and pNN50 by Subject", xlab="Sound (dBA)", ylab="pNN50",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  Lfnu ~ Sound,
  data = pro, par.settings = theEconomist.theme(),
  main="Sound Exposure and LFnu by Subject", xlab="Sound (dBA)", ylab="LFnu",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  LF.HF ~ Sound,
  data = pro, par.settings = theEconomist.theme(),
  main="Sound Exposure and LF:HF by Subject", xlab="Sound (dBA)", ylab="LF:HF",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)




xyplot(
  SDANN ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and SDANN by Subject", xlab="VSAeq", ylab="SDNN (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  RMSSD ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and RMSSD by Subject", xlab="VSAeq", ylab="RMSSD (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  pNN50 ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and pNN50 by Subject", xlab="VSAeq", ylab="pNN50",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  Lfnu ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and LFnu by Subject", xlab="VSAeq", ylab="LFnu",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  LF.HF ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and LF:HF by Subject", xlab="VSAeq", ylab="LF:HF",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)



xyplot(
  SDANN ~ MeanPowerZf,
  data = pro, par.settings = theEconomist.theme(),
  main="MeanPowerZf Exposure and SDANN by Subject", xlab="MeanPowerZf", ylab="SDNN (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  RMSSD ~ MeanPowerZf,
  data = pro, par.settings = theEconomist.theme(),
  main="MeanPowerZf Exposure and RMSSD by Subject", xlab="MeanPowerZf", ylab="RMSSD (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  pNN50 ~ MeanPowerZf,
  data = pro, par.settings = theEconomist.theme(),
  main="MeanPowerZf Exposure and pNN50 by Subject", xlab="MeanPowerZf", ylab="pNN50",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  Lfnu ~ MeanPowerZf,
  data = pro, par.settings = theEconomist.theme(),
  main="MeanPowerZf Exposure and LFnu by Subject", xlab="MeanPowerZf", ylab="LFnu",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  LF.HF ~ MeanPowerZf,
  data = pro, par.settings = theEconomist.theme(),
  main="MeanPowerZf Exposure and LF:HF by Subject", xlab="MeanPowerZf", ylab="LF:HF",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)


xyplot(
  SDANN ~ ZsAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsAeq Exposure and SDANN by Subject", xlab="ZsAeq", ylab="SDNN (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  RMSSD ~ ZsAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsAeq Exposure and RMSSD by Subject", xlab="ZsAeq", ylab="RMSSD (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  pNN50 ~ ZsAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsAeq Exposure and pNN50 by Subject", xlab="ZsAeq", ylab="pNN50",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  Lfnu ~ ZsAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsAeq Exposure and LFnu by Subject", xlab="ZsAeq", ylab="LFnu",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  LF.HF ~ ZsAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsAeq Exposure and LF:HF by Subject", xlab="ZsAeq", ylab="LF:HF",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)



xyplot(
  SDANN ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and SDANN by Subject", xlab="VSAeq", ylab="SDNN (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  RMSSD ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and RMSSD by Subject", xlab="VSAeq", ylab="RMSSD (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  pNN50 ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and pNN50 by Subject", xlab="VSAeq", ylab="pNN50",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  Lfnu ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and LFnu by Subject", xlab="VSAeq", ylab="LFnu",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  LF.HF ~ VSAeq,
  data = pro, par.settings = theEconomist.theme(),
  main="VSAeq Exposure and LF:HF by Subject", xlab="VSAeq", ylab="LF:HF",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)


xyplot(
  SDANN ~ VSVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="VSVDV8 Exposure and SDANN by Subject", xlab="VSVDV8", ylab="SDNN (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  RMSSD ~ VSVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="VSVDV8 Exposure and RMSSD by Subject", xlab="VSVDV8", ylab="RMSSD (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  pNN50 ~ VSVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="VSVDV8 Exposure and pNN50 by Subject", xlab="VSVDV8", ylab="pNN50",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  Lfnu ~ VSVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="VSVDV8 Exposure and LFnu by Subject", xlab="VSVDV8", ylab="LFnu",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  LF.HF ~ VSVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="VSVDV8 Exposure and LF:HF by Subject", xlab="VSVDV8", ylab="LF:HF",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)



xyplot(
  SDANN ~ ZsVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsVDV8 Exposure and SDANN by Subject", xlab="ZsVDV8", ylab="SDNN (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  RMSSD ~ ZsVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsVDV8 Exposure and RMSSD by Subject", xlab="ZsVDV8", ylab="RMSSD (ms)",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  pNN50 ~ ZsVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsVDV8 Exposure and pNN50 by Subject", xlab="ZsVDV8", ylab="pNN50",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  Lfnu ~ ZsVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsVDV8 Exposure and LFnu by Subject", xlab="ZsVDV8", ylab="LFnu",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)

xyplot(
  LF.HF ~ ZsVDV8,
  data = pro, par.settings = theEconomist.theme(),
  main="ZsVDV8 Exposure and LF:HF by Subject", xlab="ZsVDV8", ylab="LF:HF",
  groups = Subject, type=c('r', 'g', 'p'), auto.key = list(title='Subject', space='right')
)


##########################################################################################################################################################
 
##########################################################################################################################################################
main="Distributions of Sound and Zs VDV8 by Vehicle Speed"
mod <- lmList(Sound ~ ZsVDV8 | Category, data = pro, na.action = na.pass)
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
xyplot(Sound ~ ZsVDV8 | Category, data = pro,  type = c("p", "r"), main=main, xlab = expression("m/" ~ s^{1.75}), ylab = "Sound (dBA)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

##########FIND HERE LATER############

##########################################################################################################################################################

##########################################################################################################################################################


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

mod <- lm(SDANN ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(SDANN ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and SDANN", xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(pNN50 ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(pNN50 ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and pNN50", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"), na.action = na.pass)
trellis.focus("toplevel")
panel.text(0.3, 0.87, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Lfnu ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Lfnu ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and LFnu", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.89, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(Hfnu ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(Hfnu ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and HFnu", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lm(LF.HF ~ (FreqAmpZs), data = pro)
modsum <- summary(mod)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
xyplot(LF.HF ~ FreqAmpZs, data = pro,  type = c("p", "r"), main="Distribution of Mean Z Amplitude Exposures and LF:HF Ratio", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.86, rp, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(RMSSD ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(RMSSD ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and RMSSD by Vehicle Speed", xlab = "Mean Amplitude", ylab = "RMSSD (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(SDANN ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(SDANN ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and SDANN by Vehicle Speed", xlab = "Mean Amplitude", ylab = "SDNN (ms)", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(pNN50 ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(pNN50 ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and pNN50 by Vehicle Speed", xlab = "Mean Amplitude", ylab = "pNN50 (%)", par.settings = theEconomist.theme(), lattice.options = list(panel.error = "warning"))
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Lfnu ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Lfnu ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and LFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Lfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.4, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.4, 0.458, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.455, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(Hfnu ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(Hfnu ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and HFnu by Vehicle Speed", xlab = "Mean Amplitude", ylab = "Hfnu", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.3, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.3, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

mod <- lmList(LF.HF ~ (FreqAmpZs) | Category, data = pro, na.action = na.pass)
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
xyplot(LF.HF ~ FreqAmpZs | Category, data = pro,  type = c("p", "r"), main="Distributions of Mean Z Amplitude and LF:HF by Vehicle Speed", xlab = "Mean Amplitude", ylab = "LF:HF", par.settings = theEconomist.theme())
trellis.focus("toplevel")
panel.text(0.354, 0.84, rp3, cex = 1.2, font = 2)
panel.text(0.354, 0.45, rp1, cex = 1.2, font = 2)
panel.text(0.7, 0.45, rp2, cex = 1.2, font = 2)
trellis.unfocus()
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()



       
       ##########################################################################################################################################################
       
       ##########################################################################################################################################################
       
       ##########################################################################################################################################################
       
       ##########################################################################################################################################################
       




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

cs1 <- corAR1(0.2, form = ~1 | Subject)

cs1AR1 <- corAR1(0.8, form = ~1 | Subject)
cs1AR1. <- Initialize(cs1AR1, data=pro)
corMatrix(cs1AR1.)
##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
#pro$Time.Index[1] <- 1
#pro$Time.Index[2] <- 2
#pro$Time.Index[3] <- 3
#pro$Time.Index[4] <- 4
#pro$Time.Index[5] <- 5
#pro$Time.Index[6] <- 6

main= "rMSSD by Subject and Time"
xyplot(RMSSD ~ Time.Index | Subject, data = pro, type = "l", par.settings = theEconomist.theme(), ylab = "rMSSD (ms)", xlab = "Time Index (Number of 5-Minute Windows)", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


main= "SDNN by Subject and Time"
xyplot(SDANN ~ Time.Index | Subject, data = pro, type = "l", par.settings = theEconomist.theme(), ylab = "SDNN (ms)", xlab = "Time Index (Number of 5-Minute Windows)", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main= "pNN50 by Subject and Time"
xyplot(pNN50 ~ Time.Index | Subject, data = pro, type = "l", par.settings = theEconomist.theme(), ylab = "pNN50 (%)", xlab = "Time Index (Number of 5-Minute Windows)", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main= "LF-HF by Subject and Time"
xyplot(LF.HF ~ Time.Index | Subject, data = pro, type = "l", par.settings = theEconomist.theme(), ylab = "LF-HF", xlab = "Time Index (Number of 5-Minute Windows)", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main= "LFnu by Subject and Time"
xyplot(Lfnu ~ Time.Index | Subject, data = pro, type = "l", par.settings = theEconomist.theme(), ylab = "LFnu", xlab = "Time Index (Number of 5-Minute Windows)", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main= "HFnu by Subject and Time"
xyplot(Hfnu ~ Time.Index | Subject, data = pro, type = "l", par.settings = theEconomist.theme(), ylab = "HFnu", xlab = "Time Index (Number of 5-Minute Windows)", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()



main= "rMSSD by SDNN"
xyplot(RMSSD ~ SDANN, data = pro, type = "p", par.settings = theEconomist.theme(), ylab = "rMSSD (ms)", xlab = "SDNN (ms)", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main= "rMSSD by pNN50"
xyplot(RMSSD ~ pNN50, data = pro, type = "p", par.settings = theEconomist.theme(), ylab = "rMSSD (ms)", xlab = "pNN50 (%)", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main= "rMSSD by LF-HF"
xyplot(RMSSD ~ LF-HF, data = pro, type = "p", par.settings = theEconomist.theme(), ylab = "rMSSD (ms)", xlab = "LF-HF", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()


main= "rMSSD by LFnu"
xyplot(RMSSD ~ Lfnu, data = pro, type = "p", par.settings = theEconomist.theme(), ylab = "rMSSD (ms)", xlab = "LFnu", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main= "rMSSD by HFnu"
xyplot(RMSSD ~ Hfnu, data = pro, type = "p", par.settings = theEconomist.theme(), ylab = "rMSSD (ms)", xlab = "HFnu", main=main)
dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()

main="Correlations Among HRV Measures"
splom(~pro[71:76], data = pro, par.settings = theEconomist.theme(), panel.axis(ticks=FALSE, outside=TRUE), 
      varnames = c("sdnn", "pNN50", "rmssd", "LFnu", "HFnu", "LF/HF"), main=main, xlab="")

dev.copy(png, paste(main, "png", sep="."))
dev.off()
dev.off()



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

#Proposal Time Find Here#

##########################################################################################################################################################

##########################################################################################################################################################
attach(pro)
setwd("C:\\users\\mhughes\\desktop")
##########################################################################################################################################################
abbr <- na.omit(pro)

a1<-log(abbr$Air)
a2 <- log(abbr$AirLag)
a3 <- log(abbr$AirCum)

n1 <- abbr$Sound
n2 <- abbr$SoundLag
n3 <- abbr$SoundCum

xa1 <- abbr$XsAeq
xa2 <- abbr$XsAeqLag
xa3 <- abbr$XsAeqCum

ya1 <- abbr$YsAeq
ya2 <- abbr$YsAeqLag
ya3 <- abbr$YsAeqCum

za1 <- abbr$ZsAeq
za2 <- abbr$ZsAeqLag
za3 <- abbr$ZsAeqCum

va1 <- abbr$VSAeq
va2 <- abbr$VSAeqLag
va3 <- abbr$VSAeqCum


xv1 <- abbr$XsVDV8
xv2 <- abbr$XsVDV8Lag
xv3 <- abbr$XsVDV8Cum

yv1 <- abbr$YsVDV8
yv2 <- abbr$YsVDV8Lag
yv3 <- abbr$YsVDV8Cum

zv1 <- abbr$ZsVDV8
zv2 <- abbr$ZsVDV8Lag
zv3 <- abbr$ZsVDV8Cum

vv1 <- abbr$VSVDV8
vv2 <- abbr$VSVDV8Lag
vv3 <- abbr$VSVDV8Cum

xam <- abbr$FreqAmpXs
zam <- abbr$FreqAmpZs

y = abbr$Hfnu
x = xam
xlag = xam
xcum = xam

mod.dat <- groupedData(y ~ x | Subject, data=abbr) 
mod1 <- lme(y ~ a1 + n1, data = mod.dat, random = ~ 1 | Subject, method="ML")
mod2 <- lme(y ~ a2 + n2, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="ML")
mod3 <- lme(y ~ a3 + n3, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="ML")

mod4 <- lme(y ~ a1 + n1 + a1*n1, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="ML")
mod5 <- lme(y ~ a2 + n2 + a2*n2, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="ML")
mod6 <- lme(y ~ a3 + n3 + a1*n1, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="ML")
mod7 <- lme(y ~ a3 + n3 + a3*n3, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="ML")

names <- c("mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


mod.dat <- groupedData(y ~ x | Subject, data=abbr) 
mod1 <- lme(y ~ a1 + n1, data = mod.dat, random = ~ 1 | Subject, method="ML")
mod2 <- lme(y ~ a2 + n2, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="REML")
mod3 <- lme(y ~ a3 + n3, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="REML")

mod4 <- lme(y ~ a1 + n1 + a1*n1, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="REML")
mod5 <- lme(y ~ a2 + n2 + a2*n2, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="REML")
mod6 <- lme(y ~ a3 + n3 + a1*n1, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="REML")
mod7 <- lme(y ~ a3 + n3 + a3*n3, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit, corr = corAR1(), method="REML")
sqrt(mean(residuals(mod1)^2))
sqrt(mean(residuals(mod2)^2))
sqrt(mean(residuals(mod3)^2))
sqrt(mean(residuals(mod4)^2))
sqrt(mean(residuals(mod5)^2))
sqrt(mean(residuals(mod6)^2))
sqrt(mean(residuals(mod7)^2))
 
bestmod <- mod6
summary(bestmod)

bestmod$coefficients$fixed

plot(fitted(bestmod), resid(bestmod))
abline(h=0)
plot(fitted(bestmod), abs(resid(bestmod)))


plot(resid(bestmod) ~ fitted(bestmod), pch=unclass(pro$Subject))

plot(jitter(fitted(bestmod)), resid(bestmod), xlab="Fitted", ylab="Residuals")

qqnorm(resid(bestmod), ylab="Residuals")
qqline(resid(bestmod))


gi = influence(bestmod)
qqnorml(gi$coef[,2])



dev.copy(png, paste(main, "png", sep="."))

intervals(bestmod)

anova(bestmod)

mod4 <- lm(VSVDV8 ~ Seat, data = pro, na.action = na.omit)
summary(mod4)


plot(bestmod)
qqnorm(resid(bestmod))
qqnorm(bestmod, ~ranef(., level=0))
qqnorm(bestmod, ~ranef(., level=1))
qqnorm(bestmod, ~ resid(., type = "p") | Subject)
qqnorm(bestmod, ~ resid(., type = "p"))
qqnorm(bestmod, ~ranef(.))




##########################################################################################################################################################
#####HERE#######


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

m <- lm(SDANN ~ Sound, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(Sound), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(SDANN ~ Sound, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(SDANN ~ ., subset=c(1, 2))

mod.dat <- groupedData(SDANN ~ Sound | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ Sound, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ Sound + SoundLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ Sound + SoundLag + SoundCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ Sound, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ Sound + SoundLag + SoundCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ Sound + SoundLag + SoundCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(SDANN ~ Sound + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest, "BR")

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod2, mod3)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(SDANN ~ XsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(XsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(SDANN ~ XsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(SDANN ~ ., subset=c(1, 2))

mod.dat <- groupedData(SDANN ~ XsAeq | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ XsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ XsAeq + XsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ XsAeq + XsAeqLag + XsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ XsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ XsAeq + XsAeqLag + XsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ XsAeq + XsAeqLag + XsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(SDANN ~ XsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(SDANN ~ YsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(YsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(SDANN ~ YsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(SDANN ~ ., subset=c(1, 2))

mod.dat <- groupedData(SDANN ~ YsAeq | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ YsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ YsAeq + YsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ YsAeq + YsAeqLag + YsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ YsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ YsAeq + YsAeqLag + YsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ YsAeq + YsAeqLag + YsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(SDANN ~ YsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(SDANN ~ ZsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(ZsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(SDANN ~ ZsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

opar <- par(mfrow = c(2,2), oma = c(0,0,1.1,0))
plot(m, las = 1)

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

mod.dat <- groupedData(SDANN ~ ZsAeq | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ ZsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ ZsAeq + ZsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ ZsAeq + ZsAeqLag + ZsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ ZsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ ZsAeq + ZsAeqLag + ZsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ ZsAeq + ZsAeqLag + ZsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(SDANN ~ ZsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ VSAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(VSAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(SDANN ~ VSAeq, data = pro)
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

mod.dat <- groupedData(SDANN ~ VSAeq | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ VSAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ VSAeq + VSAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ VSAeq + VSAeqLag + VSAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ VSAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ VSAeq + VSAeqLag + VSAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ VSAeq + VSAeqLag + VSAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(SDANN ~ VSAeq + Seat, data = mod.dat, random = ~ 1 | Subject)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(SDANN ~ VSAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ ZfAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(ZfAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(SDANN ~ ZfAeq, data = pro)
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

mod.dat <- groupedData(SDANN ~ ZfAeq | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ ZfAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ ZfAeq + ZfAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ ZfAeq + ZfAeqLag + ZfAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ ZfAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ ZfAeq + ZfAeqLag + ZfAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ ZfAeq + ZfAeqLag + ZfAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(SDANN ~ ZfAeq + Seat, data = mod.dat, random = ~ 1 | Subject)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(SDANN ~ ZfAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ XsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(XsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(SDANN ~ XsVDV8, data = pro)
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

mod.dat <- groupedData(SDANN ~ XsVDV8 | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ XsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ XsVDV8 + XsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ XsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(SDANN ~ XsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(SDANN ~ XsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ YsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(YsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ YsVDV8, data = pro)
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

mod.dat <- groupedData(SDANN ~ YsVDV8 | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ YsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ YsVDV8 + YsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ YsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(SDANN ~ YsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(SDANN ~ YsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ ZsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(ZsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ ZsVDV8, data = pro)
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

mod.dat <- groupedData(SDANN ~ ZsVDV8 | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ ZsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ ZsVDV8 + ZsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ ZsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(SDANN ~ ZsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(SDANN ~ ZsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ ZfVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(ZfVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ ZfVDV8, data = pro)
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

mod.dat <- groupedData(SDANN ~ ZfVDV8 | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ ZfVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ ZfVDV8 + ZfVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ ZfVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(SDANN ~ ZfVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(SDANN ~ ZfVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ ZsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(ZsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ ZsVDV8, data = pro)
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

mod.dat <- groupedData(SDANN ~ ZsVDV8 | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ ZsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ ZsVDV8 + ZsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ ZsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(SDANN ~ ZsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

intervals(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(SDANN ~ ZsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ VSVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(VSVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ VSVDV8, data = pro)
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

mod.dat <- groupedData(SDANN ~ VSVDV8 | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ VSVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ VSVDV8 + VSVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ VSVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(SDANN ~ VSVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(SDANN ~ VSVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ ZfVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(ZfVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ ZfVDV8, data = pro)
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

mod.dat <- groupedData(SDANN ~ ZfVDV8 | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ ZfVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ ZfVDV8 + ZfVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(SDANN ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(SDANN ~ ZfVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(SDANN ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(SDANN ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(SDANN ~ ZfVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod5)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(SDANN ~ ZfVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ MeanPowerXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(MeanPowerXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ MeanPowerXs, data = pro)
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

mod.dat <- groupedData(SDANN ~ MeanPowerXs | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ MeanPowerXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ MeanPowerXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ MeanPowerXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ MeanPowerYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(MeanPowerYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ MeanPowerYs, data = pro)
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

mod.dat <- groupedData(SDANN ~ MeanPowerYs | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ MeanPowerYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ MeanPowerYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ MeanPowerYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ MeanPowerZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(MeanPowerZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ MeanPowerZs, data = pro)
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

mod.dat <- groupedData(SDANN ~ MeanPowerZs | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ MeanPowerZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ MeanPowerZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ MeanPowerZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ MeanPowerZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(MeanPowerZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ MeanPowerZf, data = pro)
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

mod.dat <- groupedData(SDANN ~ MeanPowerZf | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ MeanPowerZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ MeanPowerZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ MeanPowerZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ FreqPeakXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(FreqPeakXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ FreqPeakXs, data = pro)
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

mod.dat <- groupedData(SDANN ~ FreqPeakXs | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ FreqPeakXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ FreqPeakXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ FreqPeakXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ FreqPeakYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(FreqPeakYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ FreqPeakYs, data = pro)
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

mod.dat <- groupedData(SDANN ~ FreqPeakYs | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ FreqPeakYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ FreqPeakYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ FreqPeakYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ FreqPeakZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(FreqPeakZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ FreqPeakZs, data = pro)
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

mod.dat <- groupedData(SDANN ~ FreqPeakZs | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ FreqPeakZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ FreqPeakZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ FreqPeakZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ FreqPeakZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(FreqPeakZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ FreqPeakZf, data = pro)
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

mod.dat <- groupedData(SDANN ~ FreqPeakZf | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ FreqPeakZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ FreqPeakZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ FreqPeakZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ FreqAmpXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(FreqAmpXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ FreqAmpXs, data = pro)
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

mod.dat <- groupedData(SDANN ~ FreqAmpXs | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ FreqAmpXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ FreqAmpXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)


anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ FreqAmpXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ FreqAmpYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(FreqAmpYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ FreqAmpYs, data = pro)
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

mod.dat <- groupedData(SDANN ~ FreqAmpYs | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ FreqAmpYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ FreqAmpYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ FreqAmpYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ FreqAmpZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(FreqAmpZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ FreqAmpZs, data = pro)
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

mod.dat <- groupedData(SDANN ~ FreqAmpZs | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ FreqAmpZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ FreqAmpZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ FreqAmpZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(SDANN ~ FreqAmpZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(SDANN ~ log(FreqAmpZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(SDANN ~ FreqAmpZf, data = pro)
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

mod.dat <- groupedData(SDANN ~ FreqAmpZf | Subject, data=pro) 
mod0 <- lme(SDANN ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(SDANN ~ FreqAmpZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(SDANN ~ FreqAmpZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(SDANN ~ FreqAmpZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(pNN50 ~ Air, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(Air), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(n) ~ fitted(n), pch=unclass(pro$Seat))

plot(jitter(fitted(n)), resid(n), xlab="Fitted", ylab="Residuals")

qqnorm(resid(n), ylab="Residuals")
qqline(resid(n))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(n)))
rug(resid(n))

n <- lm(pNN50 ~ log(Air), data = pro)
gi = influence(n)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(n))
pro[cooks.distance(n) > 0.2, ]
n <- lm(pNN50 ~ ., subset=c(1, 2))

mod.dat <- groupedData(pNN50 ~ log(Air) | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ log(Air), data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ log(Air) + log(AirLag), random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ log(Air) + log(AirLag) + log(AirCum), random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ log(Air), data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ log(Air) + log(AirLag) + log(AirCum) + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ log(Air) + log(AirLag) + log(AirCum) + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(pNN50 ~ Air + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

summary(mod5)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod2, mod3)

anova(mod5)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
plot(mod3, resid(., type = "p") ~ fitted(.) | log(Air), id = 0.05, adj = -0.3 )
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(pNN50 ~ Sound, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(Sound), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(pNN50 ~ Sound, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(pNN50 ~ ., subset=c(1, 2))

mod.dat <- groupedData(pNN50 ~ Sound | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ Sound, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ Sound + SoundLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ Sound + SoundLag + SoundCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ Sound, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ Sound + SoundLag + SoundCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ Sound + SoundLag + SoundCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(pNN50 ~ Sound + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest, "BR")

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod2, mod3)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(pNN50 ~ XsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(XsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(pNN50 ~ XsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(pNN50 ~ ., subset=c(1, 2))

mod.dat <- groupedData(pNN50 ~ XsAeq | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ XsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ XsAeq + XsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ XsAeq + XsAeqLag + XsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ XsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ XsAeq + XsAeqLag + XsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ XsAeq + XsAeqLag + XsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(pNN50 ~ XsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(pNN50 ~ YsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(YsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(pNN50 ~ YsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(pNN50 ~ ., subset=c(1, 2))

mod.dat <- groupedData(pNN50 ~ YsAeq | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ YsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ YsAeq + YsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ YsAeq + YsAeqLag + YsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ YsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ YsAeq + YsAeqLag + YsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ YsAeq + YsAeqLag + YsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(pNN50 ~ YsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(pNN50 ~ ZsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(ZsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(pNN50 ~ ZsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

opar <- par(mfrow = c(2,2), oma = c(0,0,1.1,0))
plot(m, las = 1)

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

mod.dat <- groupedData(pNN50 ~ ZsAeq | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ ZsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ ZsAeq + ZsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ ZsAeq + ZsAeqLag + ZsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ ZsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ ZsAeq + ZsAeqLag + ZsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ ZsAeq + ZsAeqLag + ZsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(pNN50 ~ ZsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ VSAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(VSAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(pNN50 ~ VSAeq, data = pro)
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

mod.dat <- groupedData(pNN50 ~ VSAeq | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ VSAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ VSAeq + VSAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ VSAeq + VSAeqLag + VSAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ VSAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ VSAeq + VSAeqLag + VSAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ VSAeq + VSAeqLag + VSAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(pNN50 ~ VSAeq + Seat, data = mod.dat, random = ~ 1 | Subject)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(pNN50 ~ VSAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ ZfAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(ZfAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(pNN50 ~ ZfAeq, data = pro)
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

mod.dat <- groupedData(pNN50 ~ ZfAeq | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ ZfAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ ZfAeq + ZfAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ ZfAeq + ZfAeqLag + ZfAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ ZfAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ ZfAeq + ZfAeqLag + ZfAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ ZfAeq + ZfAeqLag + ZfAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(pNN50 ~ ZfAeq + Seat, data = mod.dat, random = ~ 1 | Subject)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(pNN50 ~ ZfAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ XsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(XsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(pNN50 ~ XsVDV8, data = pro)
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

mod.dat <- groupedData(pNN50 ~ XsVDV8 | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ XsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ XsVDV8 + XsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ XsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(pNN50 ~ XsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(pNN50 ~ XsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ YsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(YsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ YsVDV8, data = pro)
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

mod.dat <- groupedData(pNN50 ~ YsVDV8 | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ YsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ YsVDV8 + YsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ YsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(pNN50 ~ YsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(pNN50 ~ YsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ ZsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(ZsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ ZsVDV8, data = pro)
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

mod.dat <- groupedData(pNN50 ~ ZsVDV8 | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ ZsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ ZsVDV8 + ZsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ ZsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(pNN50 ~ ZsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(pNN50 ~ ZsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ ZfVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(ZfVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ ZfVDV8, data = pro)
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

mod.dat <- groupedData(pNN50 ~ ZfVDV8 | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ ZfVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ ZfVDV8 + ZfVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ ZfVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(pNN50 ~ ZfVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(pNN50 ~ ZfVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ ZsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(ZsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ ZsVDV8, data = pro)
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

mod.dat <- groupedData(pNN50 ~ ZsVDV8 | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ ZsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ ZsVDV8 + ZsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ ZsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(pNN50 ~ ZsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

intervals(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(pNN50 ~ ZsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ VSVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(VSVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ VSVDV8, data = pro)
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

mod.dat <- groupedData(pNN50 ~ VSVDV8 | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ VSVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ VSVDV8 + VSVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ VSVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(pNN50 ~ VSVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(pNN50 ~ VSVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ ZfVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(ZfVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ ZfVDV8, data = pro)
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

mod.dat <- groupedData(pNN50 ~ ZfVDV8 | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ ZfVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ ZfVDV8 + ZfVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(pNN50 ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(pNN50 ~ ZfVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(pNN50 ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(pNN50 ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(pNN50 ~ ZfVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod5)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(pNN50 ~ ZfVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ MeanPowerXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(MeanPowerXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ MeanPowerXs, data = pro)
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

mod.dat <- groupedData(pNN50 ~ MeanPowerXs | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ MeanPowerXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ MeanPowerXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ MeanPowerXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ MeanPowerYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(MeanPowerYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ MeanPowerYs, data = pro)
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

mod.dat <- groupedData(pNN50 ~ MeanPowerYs | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ MeanPowerYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ MeanPowerYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ MeanPowerYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ MeanPowerZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(MeanPowerZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ MeanPowerZs, data = pro)
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

mod.dat <- groupedData(pNN50 ~ MeanPowerZs | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ MeanPowerZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ MeanPowerZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ MeanPowerZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ MeanPowerZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(MeanPowerZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ MeanPowerZf, data = pro)
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

mod.dat <- groupedData(pNN50 ~ MeanPowerZf | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ MeanPowerZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ MeanPowerZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ MeanPowerZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ FreqPeakXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(FreqPeakXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ FreqPeakXs, data = pro)
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

mod.dat <- groupedData(pNN50 ~ FreqPeakXs | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ FreqPeakXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ FreqPeakXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ FreqPeakXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ FreqPeakYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(FreqPeakYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ FreqPeakYs, data = pro)
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

mod.dat <- groupedData(pNN50 ~ FreqPeakYs | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ FreqPeakYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ FreqPeakYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ FreqPeakYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ FreqPeakZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(FreqPeakZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ FreqPeakZs, data = pro)
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

mod.dat <- groupedData(pNN50 ~ FreqPeakZs | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ FreqPeakZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ FreqPeakZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ FreqPeakZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ FreqPeakZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(FreqPeakZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ FreqPeakZf, data = pro)
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

mod.dat <- groupedData(pNN50 ~ FreqPeakZf | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ FreqPeakZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ FreqPeakZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ FreqPeakZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ FreqAmpXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(FreqAmpXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ FreqAmpXs, data = pro)
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

mod.dat <- groupedData(pNN50 ~ FreqAmpXs | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ FreqAmpXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ FreqAmpXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)


anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ FreqAmpXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ FreqAmpYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(FreqAmpYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ FreqAmpYs, data = pro)
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

mod.dat <- groupedData(pNN50 ~ FreqAmpYs | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ FreqAmpYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ FreqAmpYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ FreqAmpYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ FreqAmpZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(FreqAmpZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ FreqAmpZs, data = pro)
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

mod.dat <- groupedData(pNN50 ~ FreqAmpZs | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ FreqAmpZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ FreqAmpZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ FreqAmpZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(pNN50 ~ FreqAmpZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(pNN50 ~ log(FreqAmpZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(pNN50 ~ FreqAmpZf, data = pro)
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

mod.dat <- groupedData(pNN50 ~ FreqAmpZf | Subject, data=pro) 
mod0 <- lme(pNN50 ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(pNN50 ~ FreqAmpZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(pNN50 ~ FreqAmpZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(pNN50 ~ FreqAmpZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

cs1 <- corAR1(0.2, form = ~1 | Subject)

cs1AR1 <- corAR1(0.8, form = ~1 | Subject)
cs1AR1. <- Initialize(cs1AR1, data=pro)
corMatrix(cs1AR1.)


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(Lfnu ~ Air, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(Air), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(n) ~ fitted(n), pch=unclass(pro$Seat))

plot(jitter(fitted(n)), resid(n), xlab="Fitted", ylab="Residuals")

qqnorm(resid(n), ylab="Residuals")
qqline(resid(n))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(n)))
rug(resid(n))

n <- lm(Lfnu ~ log(Air), data = pro)
gi = influence(n)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(n))
pro[cooks.distance(n) > 0.2, ]
n <- lm(Lfnu ~ ., subset=c(1, 2))

mod.dat <- groupedData(Lfnu ~ log(Air) | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ log(Air), data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ log(Air) + log(AirLag), random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ log(Air) + log(AirLag) + log(AirCum), random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ log(Air), data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ log(Air) + log(AirLag) + log(AirCum) + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ log(Air) + log(AirLag) + log(AirCum) + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod2b <- update(mod2, corr = corARMA(p = 1))

options(warn=-1)

set$x <- zoo(mod2b$residuals)
set$y  <- lag(x, -1, na.pad = TRUE)
set$Sub <- mod2b$groups$Subject
#set <- set[2:67, ]

xyplot(
  set$y ~ set$x,
  data = set, par.settings = theEconomist.theme(),
  group = set$Sub, auto.key=TRUE, na.rm=TRUE, type=c("p", "r"), ylim = c(-100, 200), xlim=c(-20, 90)
)

gatt<-NULL
setun <- NULL
gatt <- mod2$residuals[ ,1]
setun <- data.frame(gatt)
setun$x <- zoo(gatt)
setun$y  <- lag(setun$x, -1, na.pad = TRUE)
setun$Sub <- mod2$groups$Subject
#setun <- setun[2:67, ]

xyplot(
  setun$y ~ setun$x,
  data = setun, par.settings = theEconomist.theme(),
  group = setun$Sub, auto.key=TRUE, na.rm=TRUE, type=c("p", "r"), ylim = c(-100, 200), xlim=c(-20, 90)
)



summary(mod1)
summary(mod2)
summary(mod2b)

intervals(mod2)
intervals(mod2b)

display(mod1)
display(mod2)
display(mod3)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod2b)

seatTest <- aov(Lfnu ~ Air + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod2", "mod2b") 
aictab(cand.set=list(mod2, mod2b), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod2, mod3)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



anova(mod2b)
plot(mod2b)
qqnorm(resid(mod2b))
qqnorm(mod2b, ~ranef(., level=1))
qqnorm(mod2b, ~ resid(., type = "p") | Subject)
qqnorm(mod2b, ~ resid(., type = "p"))
qqnorm(mod2b, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(Lfnu ~ Sound, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(Sound), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(Lfnu ~ Sound, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(Lfnu ~ ., subset=c(1, 2))

mod.dat <- groupedData(Lfnu ~ Sound | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ Sound, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ Sound + SoundLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ Sound + SoundLag + SoundCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ Sound, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ Sound + SoundLag + SoundCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ Sound + SoundLag + SoundCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(Lfnu ~ Sound + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest, "BR")

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod2, mod3)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(Lfnu ~ XsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(XsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(Lfnu ~ XsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(Lfnu ~ ., subset=c(1, 2))

mod.dat <- groupedData(Lfnu ~ XsAeq | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ XsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ XsAeq + XsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ XsAeq + XsAeqLag + XsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ XsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ XsAeq + XsAeqLag + XsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ XsAeq + XsAeqLag + XsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(Lfnu ~ XsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(Lfnu ~ YsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(YsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(Lfnu ~ YsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(Lfnu ~ ., subset=c(1, 2))

mod.dat <- groupedData(Lfnu ~ YsAeq | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ YsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ YsAeq + YsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ YsAeq + YsAeqLag + YsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ YsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ YsAeq + YsAeqLag + YsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ YsAeq + YsAeqLag + YsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(Lfnu ~ YsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(Lfnu ~ ZsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(ZsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(Lfnu ~ ZsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

opar <- par(mfrow = c(2,2), oma = c(0,0,1.1,0))
plot(m, las = 1)

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

mod.dat <- groupedData(Lfnu ~ ZsAeq | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ ZsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ ZsAeq + ZsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ ZsAeq + ZsAeqLag + ZsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ ZsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ ZsAeq + ZsAeqLag + ZsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ ZsAeq + ZsAeqLag + ZsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(Lfnu ~ ZsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ VSAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(VSAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(Lfnu ~ VSAeq, data = pro)
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

mod.dat <- groupedData(Lfnu ~ VSAeq | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ VSAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ VSAeq + VSAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ VSAeq + VSAeqLag + VSAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ VSAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ VSAeq + VSAeqLag + VSAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ VSAeq + VSAeqLag + VSAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(Lfnu ~ VSAeq + Seat, data = mod.dat, random = ~ 1 | Subject)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(Lfnu ~ VSAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ ZfAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(ZfAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(Lfnu ~ ZfAeq, data = pro)
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

mod.dat <- groupedData(Lfnu ~ ZfAeq | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ ZfAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ ZfAeq + ZfAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ ZfAeq + ZfAeqLag + ZfAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ ZfAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ ZfAeq + ZfAeqLag + ZfAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ ZfAeq + ZfAeqLag + ZfAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(Lfnu ~ ZfAeq + Seat, data = mod.dat, random = ~ 1 | Subject)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(Lfnu ~ ZfAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ XsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(XsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(Lfnu ~ XsVDV8, data = pro)
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

mod.dat <- groupedData(Lfnu ~ XsVDV8 | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ XsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ XsVDV8 + XsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ XsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(Lfnu ~ XsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(Lfnu ~ XsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ YsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(YsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ YsVDV8, data = pro)
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

mod.dat <- groupedData(Lfnu ~ YsVDV8 | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ YsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ YsVDV8 + YsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ YsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(Lfnu ~ YsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(Lfnu ~ YsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ ZsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(ZsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ ZsVDV8, data = pro)
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

mod.dat <- groupedData(Lfnu ~ ZsVDV8 | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ ZsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ ZsVDV8 + ZsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ ZsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(Lfnu ~ ZsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(Lfnu ~ ZsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ ZfVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(ZfVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ ZfVDV8, data = pro)
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

mod.dat <- groupedData(Lfnu ~ ZfVDV8 | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ ZfVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ ZfVDV8 + ZfVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ ZfVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(Lfnu ~ ZfVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(Lfnu ~ ZfVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ ZsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(ZsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ ZsVDV8, data = pro)
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

mod.dat <- groupedData(Lfnu ~ ZsVDV8 | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ ZsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ ZsVDV8 + ZsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ ZsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(Lfnu ~ ZsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

intervals(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(Lfnu ~ ZsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ VSVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(VSVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ VSVDV8, data = pro)
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

mod.dat <- groupedData(Lfnu ~ VSVDV8 | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ VSVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ VSVDV8 + VSVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ VSVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(Lfnu ~ VSVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(Lfnu ~ VSVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ ZfVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(ZfVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ ZfVDV8, data = pro)
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

mod.dat <- groupedData(Lfnu ~ ZfVDV8 | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ ZfVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ ZfVDV8 + ZfVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(Lfnu ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(Lfnu ~ ZfVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(Lfnu ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(Lfnu ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(Lfnu ~ ZfVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod5)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(Lfnu ~ ZfVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ MeanPowerXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(MeanPowerXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ MeanPowerXs, data = pro)
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

mod.dat <- groupedData(Lfnu ~ MeanPowerXs | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ MeanPowerXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ MeanPowerXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ MeanPowerXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ MeanPowerYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(MeanPowerYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ MeanPowerYs, data = pro)
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

mod.dat <- groupedData(Lfnu ~ MeanPowerYs | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ MeanPowerYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ MeanPowerYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ MeanPowerYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ MeanPowerZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(MeanPowerZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ MeanPowerZs, data = pro)
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

mod.dat <- groupedData(Lfnu ~ MeanPowerZs | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ MeanPowerZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ MeanPowerZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ MeanPowerZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ MeanPowerZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(MeanPowerZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ MeanPowerZf, data = pro)
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

mod.dat <- groupedData(Lfnu ~ MeanPowerZf | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ MeanPowerZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ MeanPowerZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ MeanPowerZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ FreqPeakXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(FreqPeakXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ FreqPeakXs, data = pro)
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

mod.dat <- groupedData(Lfnu ~ FreqPeakXs | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ FreqPeakXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ FreqPeakXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ FreqPeakXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ FreqPeakYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(FreqPeakYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ FreqPeakYs, data = pro)
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

mod.dat <- groupedData(Lfnu ~ FreqPeakYs | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ FreqPeakYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ FreqPeakYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ FreqPeakYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ FreqPeakZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(FreqPeakZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ FreqPeakZs, data = pro)
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

mod.dat <- groupedData(Lfnu ~ FreqPeakZs | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ FreqPeakZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ FreqPeakZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ FreqPeakZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ FreqPeakZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(FreqPeakZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ FreqPeakZf, data = pro)
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

mod.dat <- groupedData(Lfnu ~ FreqPeakZf | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ FreqPeakZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ FreqPeakZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ FreqPeakZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ FreqAmpXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(FreqAmpXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ FreqAmpXs, data = pro)
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

mod.dat <- groupedData(Lfnu ~ FreqAmpXs | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ FreqAmpXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ FreqAmpXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)


anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ FreqAmpXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ FreqAmpYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(FreqAmpYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ FreqAmpYs, data = pro)
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

mod.dat <- groupedData(Lfnu ~ FreqAmpYs | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ FreqAmpYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ FreqAmpYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ FreqAmpYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ FreqAmpZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(FreqAmpZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ FreqAmpZs, data = pro)
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

mod.dat <- groupedData(Lfnu ~ FreqAmpZs | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ FreqAmpZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ FreqAmpZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ FreqAmpZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(Lfnu ~ FreqAmpZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(Lfnu ~ log(FreqAmpZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(Lfnu ~ FreqAmpZf, data = pro)
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

mod.dat <- groupedData(Lfnu ~ FreqAmpZf | Subject, data=pro) 
mod0 <- lme(Lfnu ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(Lfnu ~ FreqAmpZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(Lfnu ~ FreqAmpZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(Lfnu ~ FreqAmpZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################



##########################################################################################################################################################

##########################################################################################################################################################



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(LF.HF ~ Air, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(Air), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(n) ~ fitted(n), pch=unclass(pro$Seat))

plot(jitter(fitted(n)), resid(n), xlab="Fitted", ylab="Residuals")

qqnorm(resid(n), ylab="Residuals")
qqline(resid(n))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(n)))
rug(resid(n))

n <- lm(LF.HF ~ log(Air), data = pro)
gi = influence(n)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(n))
pro[cooks.distance(n) > 0.2, ]
n <- lm(LF.HF ~ ., subset=c(1, 2))

mod.dat <- groupedData(LF.HF ~ log(Air) | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ log(Air), data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ log(Air) + log(AirLag), random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ log(Air) + log(AirLag) + log(AirCum), random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ log(Air), data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ log(Air) + log(AirLag) + log(AirCum) + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ log(Air) + log(AirLag) + log(AirCum) + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(LF.HF ~ Air + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

summary(mod5)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod2, mod3)

anova(mod5)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
plot(mod3, resid(., type = "p") ~ fitted(.) | log(Air), id = 0.05, adj = -0.3 )
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(LF.HF ~ Sound, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(Sound), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(LF.HF ~ Sound, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(LF.HF ~ ., subset=c(1, 2))

mod.dat <- groupedData(LF.HF ~ Sound | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ Sound, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ Sound + SoundLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ Sound + SoundLag + SoundCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ Sound, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ Sound + SoundLag + SoundCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ Sound + SoundLag + SoundCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(LF.HF ~ Sound + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest, "BR")

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod2, mod3)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(LF.HF ~ XsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(XsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(LF.HF ~ XsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(LF.HF ~ ., subset=c(1, 2))

mod.dat <- groupedData(LF.HF ~ XsAeq | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ XsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ XsAeq + XsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ XsAeq + XsAeqLag + XsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ XsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ XsAeq + XsAeqLag + XsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ XsAeq + XsAeqLag + XsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(LF.HF ~ XsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(LF.HF ~ YsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(YsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(LF.HF ~ YsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

library(faraway)
halfnorm(cooks.distance(m))
pro[cooks.distance(m) > 0.2, ]
m <- lm(LF.HF ~ ., subset=c(1, 2))

mod.dat <- groupedData(LF.HF ~ YsAeq | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ YsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ YsAeq + YsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ YsAeq + YsAeqLag + YsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ YsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ YsAeq + YsAeqLag + YsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ YsAeq + YsAeqLag + YsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(LF.HF ~ YsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

m <- lm(LF.HF ~ ZsAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(ZsAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(LF.HF ~ ZsAeq, data = pro)
gi = influence(m)
qqnorml(gi$coef[,2])

opar <- par(mfrow = c(2,2), oma = c(0,0,1.1,0))
plot(m, las = 1)

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

mod.dat <- groupedData(LF.HF ~ ZsAeq | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ ZsAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ ZsAeq + ZsAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ ZsAeq + ZsAeqLag + ZsAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ ZsAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ ZsAeq + ZsAeqLag + ZsAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ ZsAeq + ZsAeqLag + ZsAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)


summary(mod1)
summary(mod2)
summary(mod3)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(LF.HF ~ ZsAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ VSAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(VSAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(LF.HF ~ VSAeq, data = pro)
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

mod.dat <- groupedData(LF.HF ~ VSAeq | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ VSAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ VSAeq + VSAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ VSAeq + VSAeqLag + VSAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ VSAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ VSAeq + VSAeqLag + VSAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ VSAeq + VSAeqLag + VSAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(LF.HF ~ VSAeq + Seat, data = mod.dat, random = ~ 1 | Subject)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(LF.HF ~ VSAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ ZfAeq, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(ZfAeq), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(LF.HF ~ ZfAeq, data = pro)
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

mod.dat <- groupedData(LF.HF ~ ZfAeq | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ ZfAeq, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ ZfAeq + ZfAeqLag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ ZfAeq + ZfAeqLag + ZfAeqCum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ ZfAeq, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ ZfAeq + ZfAeqLag + ZfAeqCum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ ZfAeq + ZfAeqLag + ZfAeqCum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(LF.HF ~ ZfAeq + Seat, data = mod.dat, random = ~ 1 | Subject)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

display(mod1)
display(mod2)
display(mod3)
display(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)

seatTest <- aov(LF.HF ~ ZfAeq + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

KRmodcomp(mod1, mod0)

AICc(mod3)
AICc(mod2)
AICc(mod1)
AICc(mod0)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

anova(mod3, mod5)
anova(mod2)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ XsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(XsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))

plot(density(resid(m)))
rug(resid(m))

m <- lm(LF.HF ~ XsVDV8, data = pro)
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

mod.dat <- groupedData(LF.HF ~ XsVDV8 | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ XsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ XsVDV8 + XsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ XsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ XsVDV8 + XsVDV8Lag + XsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(LF.HF ~ XsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(LF.HF ~ XsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ YsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(YsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ YsVDV8, data = pro)
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

mod.dat <- groupedData(LF.HF ~ YsVDV8 | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ YsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ YsVDV8 + YsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ YsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ YsVDV8 + YsVDV8Lag + YsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(LF.HF ~ YsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(LF.HF ~ YsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ ZsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(ZsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ ZsVDV8, data = pro)
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

mod.dat <- groupedData(LF.HF ~ ZsVDV8 | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ ZsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ ZsVDV8 + ZsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ ZsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(LF.HF ~ ZsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(LF.HF ~ ZsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ ZfVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(ZfVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ ZfVDV8, data = pro)
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

mod.dat <- groupedData(LF.HF ~ ZfVDV8 | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ ZfVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ ZfVDV8 + ZfVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ ZfVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(LF.HF ~ ZfVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(LF.HF ~ ZfVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))




##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ ZsVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(ZsVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ ZsVDV8, data = pro)
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

mod.dat <- groupedData(LF.HF ~ ZsVDV8 | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ ZsVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ ZsVDV8 + ZsVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ ZsVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ ZsVDV8 + ZsVDV8Lag + ZsVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(LF.HF ~ ZsVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)

intervals(mod5)

anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(LF.HF ~ ZsVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ VSVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(VSVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ VSVDV8, data = pro)
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

mod.dat <- groupedData(LF.HF ~ VSVDV8 | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ VSVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ VSVDV8 + VSVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ VSVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ VSVDV8 + VSVDV8Lag + VSVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(LF.HF ~ VSVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)


summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(LF.HF ~ VSVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ ZfVDV8, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(ZfVDV8), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ ZfVDV8, data = pro)
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

mod.dat <- groupedData(LF.HF ~ ZfVDV8 | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ ZfVDV8, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ ZfVDV8 + ZfVDV8Lag, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod3 <- lme(LF.HF ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)

mod4 <- lme(LF.HF ~ ZfVDV8, data = mod.dat, random = list(Subject = pdDiag(~Time.Index)))
mod5 <- lme(LF.HF ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
mod6 <- lme(LF.HF ~ ZfVDV8 + ZfVDV8Lag + ZfVDV8Cum + Seat, random = ~ 0 + Seat | Subject, data = mod.dat, na.action = na.omit)

mod7 <- lme(LF.HF ~ ZfVDV8 + Seat, data = mod.dat, random = ~ 1 | Subject)

names <- c("mod0", "mod1", "mod2", "mod3", "mod4", "mod5", "mod6", "mod7") 
aictab(cand.set=list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod5)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


anova(mod0)
anova(mod1)
anova(mod2)
anova(mod3)
anova(mod4)
anova(mod5)
anova(mod6)
anova(mod7)

seatTest <- aov(LF.HF ~ ZfVDV8 + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqnorm(mod3, ~ranef(., level=1))
qqnorm(mod3, ~ resid(., type = "p") | Subject)
qqnorm(mod3, ~ resid(., type = "p"))
qqnorm(mod3, ~ranef(.))

anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqnorm(mod5, ~ranef(., level=1))
qqnorm(mod5, ~ resid(., type = "p") | Subject)
qqnorm(mod5, ~ resid(., type = "p") | Seat)
qqnorm(mod5, ~ resid(., type = "p"))
qqnorm(mod5, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
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


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ MeanPowerYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(MeanPowerYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ MeanPowerYs, data = pro)
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

mod.dat <- groupedData(LF.HF ~ MeanPowerYs | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ MeanPowerYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ MeanPowerYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ MeanPowerYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ MeanPowerZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(MeanPowerZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ MeanPowerZs, data = pro)
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

mod.dat <- groupedData(LF.HF ~ MeanPowerZs | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ MeanPowerZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ MeanPowerZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ MeanPowerZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ MeanPowerZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(MeanPowerZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ MeanPowerZf, data = pro)
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

mod.dat <- groupedData(LF.HF ~ MeanPowerZf | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ MeanPowerZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ MeanPowerZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ MeanPowerZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ FreqPeakXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(FreqPeakXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ FreqPeakXs, data = pro)
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

mod.dat <- groupedData(LF.HF ~ FreqPeakXs | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ FreqPeakXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ FreqPeakXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ FreqPeakXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ FreqPeakYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(FreqPeakYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ FreqPeakYs, data = pro)
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

mod.dat <- groupedData(LF.HF ~ FreqPeakYs | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ FreqPeakYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ FreqPeakYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ FreqPeakYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ FreqPeakZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(FreqPeakZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ FreqPeakZs, data = pro)
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

mod.dat <- groupedData(LF.HF ~ FreqPeakZs | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ FreqPeakZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ FreqPeakZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ FreqPeakZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ FreqPeakZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(FreqPeakZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ FreqPeakZf, data = pro)
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

mod.dat <- groupedData(LF.HF ~ FreqPeakZf | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ FreqPeakZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ FreqPeakZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ FreqPeakZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ FreqAmpXs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(FreqAmpXs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ FreqAmpXs, data = pro)
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

mod.dat <- groupedData(LF.HF ~ FreqAmpXs | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ FreqAmpXs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ FreqAmpXs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)


anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ FreqAmpXs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ FreqAmpYs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(FreqAmpYs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ FreqAmpYs, data = pro)
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

mod.dat <- groupedData(LF.HF ~ FreqAmpYs | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ FreqAmpYs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ FreqAmpYs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ FreqAmpYs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))


##########################################################################################################################################################


##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ FreqAmpZs, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(FreqAmpZs), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ FreqAmpZs, data = pro)
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

mod.dat <- groupedData(LF.HF ~ FreqAmpZs | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ FreqAmpZs, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ FreqAmpZs + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ FreqAmpZs + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))



##########################################################################################################################################################

##########################################################################################################################################################

##########################################################################################################################################################
par(mfrow = c(1,1))

m <- lm(LF.HF ~ FreqAmpZf, data = pro)
plot(fitted(m), resid(m))
abline(h=0)

n <- lm(LF.HF ~ log(FreqAmpZf), data = pro)
plot(fitted(n), resid(n))
abline(h=0)
plot(fitted(n), abs(resid(n)))
summary(lm(abs(resid(n)) ~ fitted(n)))

plot(resid(m) ~ fitted(m), pch=unclass(pro$Seat))

qqnorm(resid(m), ylab="Residuals")
qqline(resid(m))


m <- lm(LF.HF ~ FreqAmpZf, data = pro)
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

mod.dat <- groupedData(LF.HF ~ FreqAmpZf | Subject, data=pro) 
mod0 <- lme(LF.HF ~ 1, random = ~ 1 | Subject, data = mod.dat)
mod1 <- lme(LF.HF ~ FreqAmpZf, data = mod.dat, random = ~ 1 | Subject)
mod2 <- lme(LF.HF ~ FreqAmpZf + Seat, random = ~ 1 | Subject, data = mod.dat, na.action = na.omit)
names <- c("mod0", "mod1", "mod2") 
aictab(cand.set=list(mod0, mod1, mod2), modnames=names, second.ord=TRUE, sort=TRUE)

intervals(mod2)

summary(mod1)
summary(mod2)

anova(mod0)
anova(mod1)
anova(mod2)

seatTest <- aov(LF.HF ~ FreqAmpZf + Seat, data=pro)
summary(seatTest)
TukeyHSD(seatTest)

anova(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqnorm(mod2, ~ranef(., level=1))
qqnorm(mod2, ~ resid(., type = "p") | Subject)
qqnorm(mod2, ~ resid(., type = "p"))
qqnorm(mod2, ~ranef(.))