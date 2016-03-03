
#************************************
#
#	Data & Code
#
#************************************
# Get xplant and xplantC from the S24Poisson.R and S23lme.R

r11xplant <- read.table("~/Desktop/linear/transplant/R11xplant.csv", sep = ",", header = T)

r11donor<-read.table("~/Desktop/linear/transplant/R11donor.csv", sep = ",", header = T)

uva <- read.table("~/Desktop/linear/transplant/UVAxplant.csv", sep = ",", header = T)

duke <- read.table("~/Desktop/linear/transplant/Dukexplant.csv", sep = ",", header = T)

mcv <- read.table("~/Desktop/linear/transplant/MCVxplant.csv", sep = ",", header = T)

unc <- read.table("~/Desktop/linear/transplant/UNCxplant.csv", sep = ",", header = T)

uvaeth <- read.table("~/Desktop/linear/transplant/UVAethnic.csv", sep = ",", header =T)

mcveth <- read.table("~/Desktop/linear/transplant/MCVethnic.csv", sep = ",", header =T)

dukeeth <- read.table("~/Desktop/linear/transplant/DukeEthnic.csv", sep = ",", header =T)

# Note - get the unc data

unceth <- read.table("~/Desktop/linear/transplant/UNCethnic.csv", sep = ",", header =T)

r11donor<-read.table("~/Desktop/linear/transplant/R11donor.csv", sep = ",", header = T)

r11xplant<-read.table("~/Desktop/linear/transplant/R11xplant.csv", sep = ",", header = T)


#***************************************************************
#
#  Source files and libraries
#
#***************************************************************

# Source 
library(boot) #If you don't have this library, install it by: install.packages('boot')
source("~/Desktop/linear/code/TSbootfunctions.R")
source("~/Desktop/linear/code/Transplant.plots.R")
source("~/Desktop/linear/code/DemographicFunctions.R")
source("~/Desktop/linear/code/LMEfunctions.R")


library(boot)
library(forecast)
library(MASS)
library(lattice)
library(lme4)


summary(xplant)

# xplantC, The version without 2014

xplantC <- xplant[-c(seq(28,112,28)),]

dim(xplant)
dim(xplantC)

# Source files and packages

# Source 

source("TSbootfunctions.R") 
source("Transplant.plots.R")
source("DemographicFunctions.R") 
source("LMEfunctions.R") 

library(boot)
library(forecast)
library(MASS)
library(lattice)
library(lme4) 


#************************************
#
#	Lagged 1 Data Frame
#
#************************************

# Need a new DF with lag 1
# for both Kidney and Liver

xplantL1 <- data.frame(xplantC[-seq(1,4*27, by =27),], LiverL1 = xplantC$Liver[-seq(27, 4*27, by = 27)], KidneyL1 = xplantC$Kidney[-seq(27, 4*27, by = 27)])

head(xplantL1[,c("Kidney", "KidneyL1","Liver", "LiverL1")])


#************************************
#
#	Mixture Model
#	for Kidney Transplants
#
#************************************

# Plot of centers

center.plot(cbind(xplantC[which(xplantC$School == "UVA"), "Kidney"],xplantC[which(xplantC$School == "MCV"), "Kidney"], xplantC[which(xplantC$School == "Duke"), "Kidney"], xplantC[which(xplantC$School == "UNC"), "Kidney"]), Year = seq(1988,max(xplantC$Year )), title = "Kidney Transplants")

# What we should do to improve UVA?


# The model that fits best

xplant.kid.lme <- lmer(Kidney~ nYears + KidneyL1+ (nYears|School), data = xplantL1, REML = F)

# Get summary

# plot of random effects
ranef(xplant.kid.lme)


dotplot(ranef(xplant.kid.lme,  condVar = T), scales = list(relation = "free"))



#*****************************
#
#	Prediction
#
#*****************************
#extract new data

newkdata <- data.frame(Kidney =0, nYears = 28, KidneyL1 = xplantC$Kidney[27], School = "UVA")

################

# prediction formulas

mm <- model.matrix(terms(xplant.kid.lme), newkdata)

xplant.kid.pred <- mm %*% t(coef(xplant.kid.lme)[[1]][which(row.names(ranef(xplant.kid.lme)$School) == newkdata$School),])


# prediction function

mypredict <- function(.)
{
	mm <- model.matrix(terms(.), newkdata)
	mm %*% t(coef(.)[[1]][which(row.names(ranef(.)$School) == "UVA"),])
}


kid.boot <- bootMer(xplant.kid.lme, mypredict, nsim =200)

# The confidence intervals

kid.boot.pred <- apply(kid.boot$t,2, quantile, c(0.025, 0.5, 0.975))

#kid.boot.pred <- quantile(kid.boot$t, c(0.025, 0.5, 0.975))[1]

# plot the prediction 
# with the TS
#how uva done in the past several years 
plot(xplantC$Year[1:27], xplantC$Kidney[1:27], type = "b", pch = 19, xlim = c(1988, 2015), ylab = "Number of Transplants", xlab = "Year", main = "Kidney Transplants at UVA")


# bootstrap  prediction

points(2015, kid.boot.pred[2] , col = "green", pch = 19)

# CI

segments(2015, kid.boot.pred[2], 2015, kid.boot.pred[1], col = "green")
segments(2015, kid.boot.pred[2], 2015, kid.boot.pred[3], col = "green")


# Twice current

points(2015, 2*xplant$Kidney[28], col = "orange", pch = 19)


#************************************
#
#	Estimating improved results
#	Simulating MCV like performance
#
#************************************

# Can we make UvA more like MCV?
#Roanoke
# Coefficients


# Get 50% of MCV advantage?

plusup <- 0.5 * (coef(xplant.kid.lme)[[1]][which(row.names(ranef(xplant.kid.lme)$School) == "MCV"),2] - coef(xplant.kid.lme)[[1]][which(row.names(ranef(xplant.kid.lme)$School) == "UVA"),2]) 

(coef(xplant.kid.lme)[[1]][which(row.names(ranef(xplant.kid.lme)$School) == "MCV"),2])
(coef(xplant.kid.lme)[[1]][which(row.names(ranef(xplant.kid.lme)$School) == "UVA"),2])
# Change coefficient change uva like mcv, use the rginal data to see the result

newcoef <- t(coef(xplant.kid.lme)[[1]][which(row.names(ranef(xplant.kid.lme)$School) == newkdata$School),]) + c(0, plusup, 0)
newcoef
# prediction formulas

# model matrix 
#extract the coefficient 
mm <- model.matrix(terms(xplant.kid.lme), newkdata)

# Simulation of new performance

uva.kid.sim <- rnorm(2000, mm %*% newcoef, sd(kid.boot$t))

# Plotting the simulation results
 
points(2015.5, median(uva.kid.sim), col = "magenta", pch=19)
length(uva.kid.sim)
#actionable suggestion mimic the mcv, cause mcv is the best in the past ten years
# CI
segments(2015.5, median(uva.kid.sim), 2015.5, quantile(uva.kid.sim, 0.975), col = "magenta")
segments(2015.5, median(uva.kid.sim), 2015.5, quantile(uva.kid.sim, 0.025), col = "magenta")

# Add a legend 

legend(1989,100, legend = c("No Change", "MCV-Like"), lwd = 2, col = c("green","magenta"))

#************************************
#
#	Poisson Regression - GLM
#	for Liver Transplants
#
#************************************

# No-pooled model
#if the centern was build, what the  performance of uva, ronoc centern was build

uva.liv.np <- glm(Liver~ LiverL1 , data = xplantL1[which(xplantC$School == "UVA"), ], family = quasipoisson)

# summary?

summary(uva.liv.np)

uva.liv.npp <- glm(Liver~ LiverL1 , data = xplantL1[which(xplantC$School == "UVA"), ], family = poisson)

summary(uva.liv.npp)

sum((residuals(uva.liv.np, type = "pearson"))^2)/uva.liv.np$df.residual

# Pooled Model

uva.liv.pool <- glm(Liver~ LiverL1 + School , data = xplantL1, family = quasipoisson)

# summary?
summary(uva.liv.pool)

##################
# predictions 

# No-pooled prediction

# New data

newLdata <- data.frame(LiverL1 = xplantC[27, "Liver"], School = "UVA")

# Prediction with SE

liv.pred.NP <- predict(uva.liv.np, type = "response", newdata = newLdata, se.fit = T)

liv.pred.NPP <- predict(uva.liv.npp, type = "response", newdata = newLdata, se.fit = T)

# prediction plot

# Base data plot 

plot(xplantC$Year[which(xplantC$School == "UVA")], xplantC$Liver[which(xplantC$School == "UVA")], type = "b", xlim = c(1988, 2015), xlab = "Year", ylab = "Number of Transplants", main = "Liver Transplants")

# Prediction plots

# No-pooled plot 

points(2015, liv.pred.NP$fit, pch = 19, col = "blue")

segments(2015, liv.pred.NP$fit, 2015, liv.pred.NP$fit + 1.96*liv.pred.NP$se.fit, col = "blue")

segments(2015, liv.pred.NP$fit, 2015, liv.pred.NP$fit - 1.96*liv.pred.NP$se.fit, col = "blue")

# Pooled model prediction

# new data

newLdata.pool <- data.frame(nYears = 27, LiverL1 = xplantC[27, "Liver"], School = "UVA")

liv.pred.pool <- predict(uva.liv.pool, type = "response", newdata = newLdata.pool, se.fit = T)

# Plot

points(2015.5, liv.pred.pool$fit, pch = 19, col = "green")

segments(2015.5, liv.pred.pool$fit, 2015.5, liv.pred.pool$fit + 1.96*liv.pred.pool$se.fit, col = "green")

segments(2015.5, liv.pred.pool$fit, 2015.5, liv.pred.pool$fit - 1.96*liv.pred.pool$se.fit, col = "green")

# Twice current 

points(2015, 2*xplant$Liver[28], col = "orange", pch = 19)


#*********************************************
#	
#		Center Model
#
#*********************************************

# Look at the plot of UVA liver transplants

# How can we improve?

#In 2005, they constucted a new center to market to transplant patients... Thus, the increase
#in transplant patients

# New UVA center opens in 2006


# Add the Center variable

xplantL1$Center <- c(rep(0, 2006-1989), rep(1, max(xplantL1$Year) - 2006 + 1), rep(0, (max(xplantL1$Year) - 1989 + 1)*3))

# No-pooled
# Model with Center 

#pooled
uva.liv.np.cent <- glm(Liver~ LiverL1 + School + Center , data = xplantL1, family = quasipoisson)

#no pooled
uva.liv.np.cent <- glm(Liver~ LiverL1 + Center , data = xplantL1[which(xplantC$School == "UVA"), ], family = quasipoisson)

# summary?
summary(uva.liv.np.cent)

# diagnostics?


# No-pooled Center prediction

# New data adding center

newLdata.center <- data.frame(LiverL1 = xplantC[28, "Liver"], School = "UVA", Center = 1)

# Prediction with SE

liv.pred.NP.cent <- predict(uva.liv.np.cent, newdata = newLdata.center, type = "response", se.fit= T)

# Plot the prediction

points(2015, liv.pred.NP.cent$fit, pch = 19, col = "tomato")

segments(2015, liv.pred.NP.cent$fit, 2015, liv.pred.NP.cent$fit + 1.96*liv.pred.NP.cent$se.fit, col = "tomato")

segments(2015, liv.pred.NP.cent$fit, 2015, liv.pred.NP.cent$fit - 1.96*liv.pred.NP.cent$se.fit, col = "tomato")

# What happens if we open another center

# New data adding a second center

newLdata.center2 <- data.frame(LiverL1 = xplantC[28, "Liver"], School = "UVA", Center = 2)

# Prediction with SE

liv.pred.NP.2cent <- predict(uva.liv.np.cent, newdata = newLdata.center2, type = "response", se.fit =  T)

# new plot

# Base plot 
# Change the ylim max of 150

plot(xplantC$Year[which(xplantC$School == "UVA")], xplantC$Liver[which(xplantC$School == "UVA")], type = "b", xlim = c(1988, 2015), xlab = "Year", ylab = "Number of Transplants", main = "Liver Transplants", ylim = c(0, 150))


# Prediction plot with new center

points(2015.5, liv.pred.NP.2cent$fit, pch = 19, col = "royalblue")

segments(2015.5, liv.pred.NP.2cent$fit, 2015.5, liv.pred.NP.2cent$fit + 1.96*liv.pred.NP.2cent$se.fit, col = "royalblue")

segments(2015.5, liv.pred.NP.2cent$fit, 2015.5, liv.pred.NP.2cent$fit - 1.96*liv.pred.NP.2cent$se.fit, col = "royalblue")

# Add the plot of the current prediction

# Add a legend




