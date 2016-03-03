gnp_t <- read.table("~/Desktop/linear/coursedata/gnp96.dat")
plot(gnp_t)
#time series
gnp_t.ts <-ts(gnp_t$V2)
plot(gnp_t.ts)
acf(gnp_t.ts)
pacf(gnp_t.ts)
#model trend
time.gnp<-c(1:length(gnp_t.ts))
gnp_t.trend <-lm(gnp_t.ts~time.gnp)
e.ts.gnp <- gnp_t.trend$residuals
acf(e.ts.gnp)
pacf(e.ts.gnp)
# first order difference
par(mfrow=c(1,2))
acf(diff(e.ts.gnp), main="ACF of Residuals from spam.trend")
pacf(diff(e.ts.gnp),main="PACF of Residuals from spam.trend")
par(mfrow=c(1,1))
#take a log for time series
plot(log(gnp_t.ts))
acf(log(gnp_t.ts))
#first order diff
#diff after log is better
gnp.log<-diff(log(gnp_t.ts))
plot(diff(gnp_t.ts))
plot(gnp.log)
#acf and pacf for gnp.log
par(mfrow=c(2,1))
acf(gnp.log, 24)
pacf(gnp.log, 24)
#fit arima models from the acf and pacf:
gnp.log.ar = arima(gnp.log, order = c(1,0,0))

library(forecast)
gnp.auto <- auto.arima(gnp.log)
gnp.auto

#diagnostics
tsdiag(gnp.auto, gof.lag = 20)
tsdiag(gnp.log.ar, gof.lag = 20)
