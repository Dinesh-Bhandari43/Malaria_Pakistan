

# LOAD THE PACKAGES
library(dlnm) ; library(splines) ; library(MASS)

library(dlnm)
library(mvmeta)
library(splines)
library(tsModel)
library(lubridate)
library(ggplot2)


Lakki <- read.csv("Malaria_Lakki_Final.csv", stringsAsFactors = FALSE)


Lakki$date <- dmy(Lakki$Date)

Lakki$mm <- month(Lakki$date)

Lakki$mm <- as.factor(Lakki$mm)

save(Lakki, file = "data.RData")
load("data.Rdata")


cb <- crossbasis(Lakki$Tmean, lag=5, argvar=list(fun="thr",thr=22.4), 
                 arglag=list(fun="integer"))

model <- glm(total_cases ~ cb + mm + ns(Lakki$date, df=4*10) 
             + offset(log(Population)) + ns(Lakki$Relative_Humidity, df=3), 
             Lakki, family=quasipoisson, maxit = 50, na.action="na.exclude")


pred <- crosspred(cb, model, at = 22.4:35.3, by=0.1) 


allRRfit <- (pred$allRRfit-1)*100
allRRfitlow <- (pred$allRRlow-1)*100
allRRfithigh <- (pred$allRRhigh-1)*100




df <- data.frame(matrix(nrow=length(pred$predvar), ncol=4))
colnames(df) <- c("temp", "RRfit", "RRlow", "RRhigh")
df$temp <- pred$predvar
df$RRfit <- (pred$allRRfit-1)*100
df$RRlow <- (pred$allRRlow-1)*100
df$RRhigh <- (pred$allRRhigh-1)*100


predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))

predvar <- quantile(Lakki$Tmean,1:99/100,na.rm=T)

# Exposure-lag response
xlab <- expression(paste("Temperature (",degree,"C)"))

#pdf("Figure 1.pdf",height=4.5,width=10)
layout(matrix(c(1,2),ncol=2,byrow=T))

# PLOT - 3D
par(mar=c(2,3,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"3d",ltheta=150,xlab="Temperature (C)",ylab="Lag",zlab="RR", 
     col=gray(0.9), main="Exposure-lag-response")


par(mar=c(5,4,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"overall",col="red",ylim=c(0.9,1.40),axes=T,lab=c(6,5,7),xlab=xlab,
     ylab="RR",main="Overall")

##########################################3
# PLOT - FIGURE 1

xlab <- expression(paste("Temperature (",degree,"C)"))

#pdf("Figure 1.pdf",height=4.5,width=10)
layout(matrix(c(1,2),ncol=2,byrow=T))

# PLOT - 3D
par(mar=c(2,3,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"3d",ltheta=150,xlab="Temperature (C)",ylab="Lag",zlab="RR", 
     col=gray(0.9), main="Exposure-lag-response")

# OVERALL
# The plots show the cumulative exposure-response association, in terms of 
#    relative risks (RR) and centered in the MMT, across the 21 days of lag.
par(mar=c(5,4,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"overall",col="red",ylim=c(0.5,2.5),axes=T,lab=c(6,5,7),xlab=xlab,
     ylab="RR",main="Overall")


layout(1)

summary(model)
#dev.off()

