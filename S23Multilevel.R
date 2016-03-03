
#***************************************************************
#***************************************************************
#
#  Designs to improve numbers of Treated Kidney Transplant Patients
#
#***************************************************************
#***************************************************************


#***************************************************************
#
#  Read in the data
#
#***************************************************************
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
library(lme4) # you will need to install if you don't have it.


#********************************
#
#	Combining Data
#
#********************************

# Years

currentYear <- 2015

nYear <- length(uva$Year)

# Look at the differences in the DF

ncol(uva)
ncol(mcv)
ncol(duke)
ncol(unc)


names(uva)
names(mcv)
names(duke)
names(unc)

# Find commmon variables
comvar <- intersect((intersect(names(uva), names(mcv))),(intersect(names(uva), names(unc))))

# just take the xplants
comvar <- comvar[1:7]

# Create a unified DF
xplant <- rbind(uva[,comvar], mcv[,comvar], duke[,comvar], unc[,comvar])

summary(xplant)

dim(xplant)

# Create a school variable
SchoolNames <- c("UVA", "MCV", "Duke", "UNC")


# Add school to the DF
xplant$School <- c(rep("UVA", 28), rep("MCV", 28), rep("Duke", 28), rep("UNC", 28))

# make it a factor

xplant$School <- factor(xplant$School)

# Add number of years to the DF

xplant$nYears <- xplant$Year - 1987

# Check the new DF

summary(xplant)

# A cleaner DF without the partial year
# So, remove the partial year

xplantC <- xplant[-c(seq(28,112,28)),]

# Check this DF

summary(xplantC)

str(xplantC)

# Ethnic data for kidney transplants

# Get the data for each school

uvaketh <- subdata("Kidney", uvaeth)

mcvketh <- subdata("Kidney", mcveth)

dukeketh <- subdata("Kidney", dukeeth)

uncketh <- subdata("Kidney", unceth)


# Remove year current year and combine all ethnic groups other than white 
# into one category

# UVA
uvake <- data.frame(Kidney.W = uvaketh[-nYear, "Kidney.W"], Kidney.O = apply(uvaketh[-nYear,which(colnames(uvaketh) != "Kidney.W")], 1, sum))

# MCV
mcvke <- data.frame(Kidney.W = mcvketh[-nYear, "Kidney.W"], Kidney.O = apply(mcvketh[-nYear,which(colnames(mcvketh) != "Kidney.W")], 1, sum))

# Duke
dukeke <- data.frame(Kidney.W = dukeeth[-nYear, "Kidney.W"], Kidney.O = apply(dukeketh[-nYear,which(colnames(dukeketh) != "Kidney.W")], 1, sum))

#UNC
uncke <- data.frame(Kidney.W = unceth[-nYear, "Kidney.W"], Kidney.O = apply(uncketh[-nYear,which(colnames(uncketh) != "Kidney.W")], 1, sum))

# Create a diversity score

ethnic <- rbind(uvake, mcvke, dukeke, uncke)
ethnic
ethnic.ratio <- ethnic$Kidney.O/(ethnic$Kidney.O+ethnic$Kidney.W)
ethnic.ratio
# Add to DF 

xplantC <- cbind(xplantC, Ratio = ethnic.ratio)
xplantC
#***************************************************************
#
#  LME models for Kidney Transplants
#
#***************************************************************

# graphics

# What are your conclusions from these plots?

bwplot(Kidney~School, data = xplantC)

bwplot(Ratio~School, data = xplantC)

# Center plot

center.plot(cbind( uva$Kidney[-nYear], unc$Kidney[-nYear], mcv$Kidney[-nYear], duke$Kidney[-nYear]), Year = seq(1988,(1988+nYear - 2)), title = "Kidney")


############################
# Simple model of varying intercepts by school
############################

xplant.lme1 <- lmer(Kidney~ 1 + (1|School) , data = xplantC, REML = T)

# What is the relationship between coefficients, fixed effects and random effects?

coef(xplant.lme1)

fixef(xplant.lme1)

ranef(xplant.lme1)

# plot of random effects

# What do you conclude?

dotplot(ranef(xplant.lme1,  condVar = T))

# Plot the intercepts

# Is this a good model?

center.plot(cbind( uva$Kidney[-nYear], unc$Kidney[-nYear], mcv$Kidney[-nYear], duke$Kidney[-nYear]), Year = seq(1988,(1988+nYear - 2)), title = "Kidney")

sapply(1:4, function(i){
	abline(h = coef(xplant.lme1)$School[i,1])
	text(2015, coef(xplant.lme1)$School[i,1], levels(xplantC$School)[i])
})


# Diagnostics

# is this a good model?

plot(xplant.lme1) # residuals vs. fits

qqnorm(residuals(xplant.lme1))
qqline(residuals(xplant.lme1))

# Any serial correlation?

# ACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	acf(residuals(xplant.lme1)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# PACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	pacf(residuals(xplant.lme1)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# Get the model with maximum likelihood
# Why?

xplant.lme1 <- lmer(Kidney~ 1 + (1|School) , data = xplantC, REML = F)

############################
# A model with number of years
# Varying intercept
############################

xplant.lme2 <- lmer(Kidney~ nYears + (1|School), data = xplantC, REML = T)

# get the coefficients, fixed & random effects

coef(xplant.lme2)

fixef(xplant.lme2)

ranef(xplant.lme2)

# plot of random effects
# what do you observe?

dotplot(ranef(xplant.lme2,  condVar = T))

# Plot the center intercepts with constant slope
# Good model?

center.plot(cbind( uva$Kidney[-nYear], unc$Kidney[-nYear], mcv$Kidney[-nYear], duke$Kidney[-nYear]), Year = seq(1988,(1988+nYear - 2)), title = "Kidney")


sapply(1:4, function(i){
	abline(a = coef(xplant.lme2)$School[i,1] - 1988*coef(xplant.lme2)$School[i,2], b = coef(xplant.lme2)$School[i,2], col = c( "blue3" , "purple" , "lightblue3","orange")[i])
})

	

# Diagnostics
# Good model?

plot(xplant.lme2) # residuals vs. fits
f <- fitted(xplant.lme2)
r <- residuals(xplant.lme2)
plot(f,r)

qqnorm(residuals(xplant.lme2))
qqline(residuals(xplant.lme2))

# Serial correlation?

# ACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	acf(residuals(xplant.lme2)[seq(1,26)+((x-1)*26)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# PACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	pacf(residuals(xplant.lme2)[seq(1,26)+((x-1)*26)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# Get the model with maximum likelihood

xplant.lme2 <- lmer(Kidney~ nYears + (1|School), data = xplantC, REML = F)

AIC(xplant.lme2)
# Compare lme1 to lme2

anova(xplant.lme1, xplant.lme2) #notice no test = "Chi"

#########################
# Model with years
# varying intercepts and slopes
#########################

# xy plots

xyplot(Kidney~nYears|School, data = xplantC,type = c("p", "r"))


xyplot(Kidney~Ratio|School, data = xplantC,type = c("p", "r"))


# the model with varying slopes and intercepts

xplant.lme3 <- lmer(Kidney~ nYears + (nYears|School), data = xplantC, REML = T)
AIC(xplant.lme3)
# get the coefficients, fixed & random effects

coef(xplant.lme3)

fixef(xplant.lme3)

ranef(xplant.lme3)

# plot of random effects
# what do you observe?

dotplot(ranef(xplant.lme3,  condVar = T), scales = list(relation = "free"))

# Diagnostics
# Good model?
AIC(uvamcv.kidney.lm)


plot(xplant.lme3) # residuals vs. fits

qqnorm(residuals(xplant.lme3))
qqline(residuals(xplant.lme3))

# Plot the center intercepts with constant slope

center.plot(cbind( uva$Kidney[-nYear], unc$Kidney[-nYear], mcv$Kidney[-nYear], duke$Kidney[-nYear]), Year = seq(1988,(1988+nYear - 2)), title = "Kidney")


sapply(1:4, function(i){
	abline(a = coef(xplant.lme3)$School[i,1] - 1988*coef(xplant.lme3)$School[i,2], b = coef(xplant.lme3)$School[i,2], col = c( "blue3" , "purple" , "lightblue3","orange")[i])
})


# Serial correlation?
# ACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	acf(residuals(xplant.lme3)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# PACF


par(mfrow = c(2,2))
sapply(1:4, function(x){
	pacf(residuals(xplant.lme3)[seq(1,26)+((x-1)*26)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# Get the model with maximum likelihood

xplant.lme3 <- lmer(Kidney~ nYears + (nYears|School), data = xplantC, REML = F)


# Compare lme2 to lme3

anova(xplant.lme2, xplant.lme3) #notice no test = "Chi"
summary(xplant.lme3)
#########################
# Model with diversity ratio
# varying intercepts and slopes
# for years lme4
#########################
#ethnic!
xplant.lme4 <- lmer(Kidney~nYears + Ratio + (1+ nYears|School), data = xplantC, REML = F)

coef(xplant.lme4)

fixef(xplant.lme4)

ranef(xplant.lme4)

dotplot(ranef(xplant.lme4,  condVar = T), scales = list(relation = "free"))

# Diagnostics
# Good model?

plot(xplant.lme4) # residuals vs. fits

qqnorm(residuals(xplant.lme4))
qqline(residuals(xplant.lme4))

# Plot the center intercepts with constant slope

center.plot(cbind( uva$Kidney[-nYear], unc$Kidney[-nYear], mcv$Kidney[-nYear], duke$Kidney[-nYear]), Year = seq(1988,(1988+nYear - 2)), title = "Kidney")


sapply(1:4, function(i){
  abline(a = coef(xplant.lme4)$School[i,1] - 1988*coef(xplant.lme3)$School[i,2], b = coef(xplant.lme3)$School[i,2], col = c( "blue3" , "purple" , "lightblue3","orange")[i])
})


# Serial correlation?
# ACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
  acf(residuals(xplant.lme4)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# PACF


par(mfrow = c(2,2))
sapply(1:4, function(x){
  pacf(residuals(xplant.lme4)[seq(1,26)+((x-1)*26)], main = SchoolNames[x])
})
par(mfrow = c(1,1))


# Compare with xplant.lme3

anova(xplant.lme3, xplant.lme4)


#########################
# Model with AR1
# varying intercepts and slopes
# for years
#########################


# Need a new DF with lag 1

xplantL1 <- data.frame(xplantC[-seq(1,4*27, by =27),], KidneyL1 = xplantC$Kidney[-seq(27, 4*27, by = 27)])

head(xplantL1[,c("Kidney", "KidneyL1")])

# Model with lag 1

xplant.lme5 <- lmer(Kidney~ nYears + KidneyL1+ (nYears|School), data = xplantL1, REML = F)
xplant.mcvuva.lme5 <- lm(uva$Kidney[2:27]-mcv$Kidney[2:27] ~KidneyL1[2:27], data = xplantL1)
# get the coefficients, fixed & random effects

coef(xplant.lme5) 

# compare with xplant.lme3

fixef(xplant.lme5)

ranef(xplant.lme5)

# plot of random effects

dotplot(ranef(xplant.lme5,  condVar = T), scales = list(relation = "free"))

# Diagnostics
# Good model?

plot(xplant.lme5) # residuals vs. fits

qqnorm(residuals(xplant.lme5))
qqline(residuals(xplant.lme5))


# Serial correlation?
# ACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	acf(residuals(xplant.lme5)[seq(1,26)+((x-1)*26)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# PACF


par(mfrow = c(2,2))
sapply(1:4, function(x){
	pacf(residuals(xplant.lme5)[seq(1,25)+((x-1)*25)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# Get the model with maximum likelihood

xplant.lme5 <- lmer(Kidney~ nYears + KidneyL1 +(nYears|School), data = xplantL1, REML = F)


# Compare lme3 to lme5

# need xplant.lme3 to use the same DF

xplant.lme3 <- lmer(Kidney~ nYears + (nYears|School), data = xplantL1, REML = F)


anova(xplant.lme3, xplant.lme5) #notice no test = "Chi"


#*****************************
#
#	Prediction
#
#*****************************


newkdata <- data.frame(Kidney =0, nYears = 28, KidneyL1 = xplantC$Kidney[27], School = "UVA")

################
# xplant.lme5

# prediction formulas

mm <- model.matrix(terms(xplant.lme5), newkdata)

xplant.lme5.pred <- mm %*% t(coef(xplant.lme5)[[1]][which(row.names(ranef(xplant.lme5)$School) == newkdata$School),])
xplant.lme5.pred <- mm %*% t(coef(xplant.lme5)[[1]][which(row.names(ranef(xplant.lme5)$School) == newkdata$School),])

xplant.lme5.pred 
#######
# Bootstrapping prediction and CI
# ignore warnings of convergence problems
#after bootstrip more prediction

kidboot.pred <- lmeboot("Kidney", xplantL1, formula(xplant.lme5), predict(xplant.lme5), residuals(xplant.lme5), newkdata, "School", "UVA", 200)

quantile(kidboot.pred, c(.025, .5, .975))

hist(kidboot.pred)
abline(v = median(kidboot.pred), col = "red", lwd = 2)


# Second method for bootstrapping
# the prediction and CI


# prediction function

mypredict <- function(.)
{
	mm <- model.matrix(terms(.), newkdata)
	mm %*% t(coef(.)[[1]][which(row.names(ranef(.)$School) == "UVA"),])
}

lme5.boot <- bootMer(xplant.lme5, mypredict, nsim =200)

lme5.boot.pred <- apply(lme5.boot$t,2, quantile, c(0.025, 0.5, 0.975))


# plot the prediction 
# with the TS

plot(xplantC$Year[1:27], xplantC$Kidney[1:27], type = "b", pch = 19, xlim = c(1988, 2015), ylab = "Number of Transplants", xlab = "Year", main = "Kidney Transplants at UVA")

# bootstrap 1 prediction

points(2015,quantile(kidboot.pred, c(.5)) , col = "red", pch = 19)

# CI

segments(2015, quantile(kidboot.pred, c(.5)), 2015, quantile(kidboot.pred, c(.975)), col = "red")
segments(2015, quantile(kidboot.pred, c(.5)), 2015, quantile(kidboot.pred, c(.025)), col = "red")


# bootstrap 2 prediction

points(2015.5, lme5.boot.pred[2] , col = "green", pch = 19)

# CI

segments(2015.5, lme5.boot.pred[2], 2015.5, lme5.boot.pred[1], col = "green")
segments(2015.5, lme5.boot.pred[2], 2015.5, lme5.boot.pred[3], col = "green")


# Current

points(2015, xplant$Kidney[28], col = "blue", pch = 19)
xplant$Kidney[28]


# Twice current

points(2015, 2*xplant$Kidney[28], col = "cyan", pch = 19)

# legend

legend(1990, 100, legend = c("Current", "2 X Current", "Bootstrap 1", "Bootstrap 2"), col = c("blue", "cyan", "red", "green"), pch = 19)

#***************************************************************
#
#  LME models for Liver Transplants
#
#***************************************************************

# Ethnic data for kidney transplants

# Get the data for each school

uvaleth <- subdata("Liver", uvaeth)

mcvleth <- subdata("Liver", mcveth)

dukeleth <- subdata("Liver", dukeeth)

uncleth <- subdata("Liver", unceth)

# Remove year current year and combine all ethnic groups other than white 
# into one category

# UVA
uvale <- data.frame(Liver.W = uvaleth[-nYear, "Liver.W"], Liver.O = apply(uvaleth[-nYear,which(colnames(uvaleth) != "Liver.W")], 1, sum))

# MCV
mcvle <- data.frame(Liver.W = mcvleth[-nYear, "Liver.W"], Liver.O = apply(mcvleth[-nYear,which(colnames(mcvleth) != "Liver.W")], 1, sum))

# Duke
#dukele <- data.frame(Liver.W = dukeeth[-nYear, "Liver.W"], Liver.O = apply(dukeleth[-nYear,which(colnames(dukeleth) != "Liver.W")], 1, sum))
dukele <- data.frame(Liver.W = dukeleth[-nYear, "Liver.W"], Liver.O = apply(dukeleth[-nYear,which(colnames(dukeleth) != "Liver.W")], 1, sum))

#UNC
#uncle <- data.frame(Liver.W = uncleth[-nYear, "Liver.W"], Liver.O = apply(uncleth[-nYear,which(colnames(uncleth) != "Liver.W")], 1, sum))
uncle <- data.frame(Liver.W = uncleth[-nYear, "Liver.W"], Liver.O = apply(uncleth[-nYear,which(colnames(uncleth) != "Liver.W")], 1, sum))
# Create a diversity score

ethnic <- rbind(uvale, mcvle, dukele, uncle)

Lethnic.ratio <- ethnic$Liver.O/(ethnic$Liver.O+ethnic$Liver.W)

# Add to DF 

xplantC <- cbind(xplantC, LRatio = Lethnic.ratio)


# graphics

bwplot(Liver~School, data = xplantC)

bwplot(LRatio~School, data = xplantC)

# Plot of centers


center.plot(cbind( uva$Liver[-nYear], unc$Liver[-nYear], mcv$Liver[-nYear], duke$Liver[-nYear]), Year = seq(1988,(1988+nYear - 2)), title = "Liver")



############################
# Simple model of varying intercepts by school
############################


Lxplant.lme1 <- lmer(Liver~ 1 + (1|School) , data = xplantC, REML = T)

# Get coefficients, fixed and random effects

coef(Lxplant.lme1)

fixef(Lxplant.lme1)

ranef(Lxplant.lme1)

# plot of random effects

dotplot(ranef(Lxplant.lme1,  condVar = T))

# Diagnostics
# Good model?

plot(Lxplant.lme1) # residuals vs. fits

qqnorm(residuals(Lxplant.lme1))
qqline(residuals(Lxplant.lme1))

# Serial correlation?
# ACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	acf(residuals(Lxplant.lme1)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# PACF


par(mfrow = c(2,2))
sapply(1:4, function(x){
	pacf(residuals(Lxplant.lme1)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# Get the model with maximum likelihood

Lxplant.lme1 <- lmer(Liver~ 1 + (1|School) , data = xplantC, REML = F)

############################
# A model with number of years
# Varying intercept
############################


Lxplant.lme2 <- lmer(Liver~ nYears + (1|School), data = xplantC, REML = T)

# Get coefficients, fixed and random effects

coef(Lxplant.lme2)

fixef(Lxplant.lme2)

ranef(Lxplant.lme2)

# plot of random effects

dotplot(ranef(Lxplant.lme2,  condVar = T))

# Diagnostics
# Good model?

plot(Lxplant.lme2) # residuals vs. fits

qqnorm(residuals(Lxplant.lme2))
qqline(residuals(Lxplant.lme2))

# Serial correlation?
# ACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	acf(residuals(Lxplant.lme2)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# PACF


par(mfrow = c(2,2))
sapply(1:4, function(x){
	pacf(residuals(Lxplant.lme2)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# Get the model with maximum likelihood

Lxplant.lme2 <- lmer(Liver~ nYears + (1|School), data = xplantC, REML = F)


# Compare lme1 to lme2

anova(Lxplant.lme1, Lxplant.lme2) #notice no test = "Chi"

#########################
# Model with years
# varying intercepts and slopes
#########################

# xy plots

xyplot(Liver~nYears|School, data = xplantC,type = c("p", "r"))


xyplot(Liver~Ratio|School, data = xplantC,type = c("p", "r"))


# the model with varying slopes and intercepts

Lxplant.lme3 <- lmer(Liver~ nYears + (nYears|School), data = xplantC, REML = T)

# Get coefficients, fixed and random effects

coef(Lxplant.lme3)

fixef(Lxplant.lme3)

ranef(Lxplant.lme3)

# plot of random effects

dotplot(ranef(Lxplant.lme3,  condVar = T))

# Diagnostics
# Good model?

plot(Lxplant.lme3) # residuals vs. fits

qqnorm(residuals(Lxplant.lme3))
qqline(residuals(Lxplant.lme3))

# Serial correlation?
# ACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	acf(residuals(Lxplant.lme3)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# PACF


par(mfrow = c(2,2))
sapply(1:4, function(x){
	pacf(residuals(Lxplant.lme3)[seq(1,27)+((x-1)*27)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# Get the model with maximum likelihood

Lxplant.lme3 <- lmer(Liver~ nYears + (nYears|School), data = xplantC, REML = F)


# Compare lme2 to lme3

anova(Lxplant.lme2, Lxplant.lme3) #notice no test = "Chi"

#########################
# Model with diversity ratio
# varying intercepts and slopes
# for years
#########################

Lxplant.lme4 <- lmer(Liver~nYears + Ratio + (1|School), data = xplantC, REML = F)

# Compare to Lxplant.lme2

anova(Lxplant.lme2, Lxplant.lme4)

