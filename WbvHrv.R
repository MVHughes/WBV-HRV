if (!require("pacman")) install.packages("pacman")
pacman::p_load(stats, AICcmodavg, arm, faraway, lme4, astsa, RHRV, pbkrtest, MASS, ggplot2, ggthemes, dplyr, timeSeries, timeDate, TSA, memisc, lubridate, forecast, zoo, stats, Hmisc, reshape2, RColorBrewer, caTools, xts, scatterplot3d, lattice, latticeExtra, gridExtra, car)

options(max.print = 222)
FileLocation = "G:/HrvValuesOnly.csv"
RecordLocation = "G:/"

par(mfrow=c(1,1))
HrvValuesOnly <- read.csv("G:/HrvValuesOnly.csv", stringsAsFactors=FALSE)
attach(HrvValuesOnly)

HrvValuesOnly$Initial.DateTime <- ymd_hms(HrvValuesOnly$Initial.DateTime)
#rm(HrvValuesOnly)
#Create table of subject measurement durations, 

HrvValuesOnly$Subject <- factor(HrvValuesOnly$Subject)
HrvValuesOnly$Period <- factor(HrvValuesOnly$Period)
HrvValuesOnly$Condition <- factor(HrvValuesOnly$Condition)
HrvValuesOnly$HRV.Segment <- factor(HrvValuesOnly$HRV.Segment)
HrvValuesOnly$id <- factor(HrvValuesOnly$id)

levels(HrvValuesOnly$id)
subsets <- split(HrvValuesOnly, HrvValuesOnly$id, drop=TRUE)

str(subsets)
head(subsets$S10NULL0M1$Actual.time.hr)
str(subsets$S10NULL0M1$DateTime)
subsets$S10NULL0M1$DateTime <- as.POSIXct(subsets$S10NULL0M1$DateTime, format = "%m/%d/%Y %H:%M")
subsets$S10NULL0M1$DateTime <- format(subsets$S10NULL0M1$DateTime, "%d/%m/%Y %H:%M:%S")


export <- data.frame(subsets$S10NULL0M1$RR)
FileName = paste(subsets$S10NULL0M1$id[1], "csv", sep = ".")
write.table(export, FileName, sep=",", col.names=FALSE, row.names = FALSE)


BegTime <- subsets$S10NULL0M1$DateTime[1]

md <- CreateHRVData(Verbose = TRUE)
md <- LoadBeatRR(md, FileName, RecordPath = ".", scale = 0.001, datetime = BegTime)

md <- BuildNIHR(md)
md <- InterpolateNIHR(md, freqhr = 4, method = "linear")
md <- CreateTimeAnalysis(md, size = 300, interval = 7.8125)
md <- CreateFreqAnalysis(md)



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
