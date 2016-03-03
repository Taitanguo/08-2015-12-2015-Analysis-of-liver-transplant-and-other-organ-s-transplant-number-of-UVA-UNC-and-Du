

#***************************************************************
#***************************************************************
#
#  Predicting the Number of Liver Transplants
#	Poisson Regression
#
#***************************************************************
#***************************************************************

# Ref
# http://glmm.wikidot.com/faq
# http://permalink.gmane.org/gmane.comp.lang.r.lme4.devel/12284
#***************************************************************
#
#  Read in the data
#
#***************************************************************


# Note - get the unc data

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

#********************************
#
#	Combining Data
#
#********************************

# Years

currentYear <- 2015

nYear <- length(uva$Year)


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


# Add school to the DF data frame

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

#**************************************
#
#		Mixture Models with Poisson
#			Regression
#
#
#**************************************

# Mixture model with varying intercepts

Lxplant.lme1 <- glmer(Liver~ 1 + (1|School), data = xplantC, family = poisson)

# Get coefficients, fixed and random effects, difference between them multilevel
#want to find a middle effect model use, is the goal.
coef(Lxplant.lme1) # this is the average of each centern of four school
ranef(Lxplant.lme1)#just the diviation of every centern from the fixed
fixef(Lxplant.lme1)#fixed:  combination from four into one, this is fixed effect

# plot the random effects with SE

dotplot(ranef(Lxplant.lme1, condVar = T))

# Residual vs. fitted plot
#performed poorly we dont want to see any pattern or trend
plot(Lxplant.lme1)

## How does this model perform with regards to model diagnostics?

# Adding number of years as a predictor

Lxplant.lme2 <- glmer(Liver~ nYears + (1|School), data = xplantC, family = poisson)

# Get summary

summary(Lxplant.lme2)

# Get coefficients, fixed and random effects
coef(Lxplant.lme2)#intercept
ranef(Lxplant.lme2)#slope
fixef(Lxplant.lme2)#

# Dotplot of ranef

dotplot(ranef(Lxplant.lme2, condVar= T)) # uve and mcv performs better

# Residual vs. fitted plot
#better
plot(Lxplant.lme2)

# Serial correlation 

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
#determine whether we should add lag
## Is there any serial correlation we need to model?

# Compare to Lxplant.lme1
anova(Lxplant.lme1, Lxplant.lme2) 
#second one better
##  Which model do you choose between lme1 and lme2 and why?


##########################################
# Model with varying intercepts and slopes

Lxplant.lme3 <- glmer(Liver~ nYears + (nYears|School), data = xplantC, family = poisson)


# Compare to varying intercepts only

anova(Lxplant.lme2, Lxplant.lme3)

##  Which model do you choose between lme2 and lme3 and why?

#########################
# Model with AR1
# varying intercepts and slopes
# for years
#########################

# Need a new DF with lag 1

xplantL1 <- data.frame(xplantC[-seq(1,4*27, by =27),], LiverL1 = xplantC$Liver[-seq(27, 4*27, by = 27)])

head(xplantL1[,c("Liver", "LiverL1")])

# Model with lag 1 no nYears

Lxplant.lme5 <- glmer(Liver~ LiverL1+ (LiverL1|School), data = xplantL1, family = poisson)

# coefficients, random and fixed effects
coef(Lxplant.lme5)
ranef(Lxplant.lme5)
fixef(Lxplant.lme5)

# plot of random effects

dotplot(ranef(Lxplant.lme5,  condVar = T),scales = list(relation = "free"))
#intercept higher the liver lower
# Diagnostics

plot(Lxplant.lme5) # residuals vs. fits

qqnorm(residuals(Lxplant.lme5))
qqline(residuals(Lxplant.lme5))

## How does this model perform with regards to model diagnostics?

# Serial correlation?
# ACF

par(mfrow = c(2,2))
sapply(1:4, function(x){
	acf(residuals(Lxplant.lme5)[seq(1,26)+((x-1)*26)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

# PACF
#no lag here great 
par(mfrow = c(2,2))
sapply(1:4, function(x){
	pacf(residuals(Lxplant.lme5)[seq(1,26)+((x-1)*26)], main = SchoolNames[x])
})
par(mfrow = c(1,1))

## Is there any serial correlation we need to model after including a lag of 1 in our model?

# Get the model with maximum likelihood

Lxplant.lme5 <- glmer(Liver~ LiverL1+ (LiverL1|School), data = xplantL1, family = poisson, REML = F)

# Test the model against just a varying intercepts model

Lxplant.lme4 <- glmer(Liver~ LiverL1+ (1|School), data = xplantL1, family = poisson)

anova(Lxplant.lme4,Lxplant.lme5)

## Which model do we choose between lme4 and lme5 and why?
#the lme5 is better 
#*****************************
#
#	Prediction
#
#*****************************

newLdata <- data.frame(Liver = rep(0, 4), LiverL1 = xplantC$Liver[seq(27,108,27)], School = c("UVA","MCV", "Duke", "UNC"))
#for four different schools 
predict(Lxplant.lme5, newdata = newLdata, type = "response")

#************************

# Overdispersion

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(Lxplant.lme5)
#dispertion ratio should be close to 1

#************************************
#
#	Poisson Regression - GLM
#	Overdispersion with Quasipoisson
#
#************************************

# A basic poisson model

uva.liv.glm1  <- glm(Liver~. , data = xplantL1[which(xplantL1$School == "UVA"), c("Liver" ,"LiverL1")], family = poisson)

# Get a summary

summary(uva.liv.glm1)

# Model utility test

uva.liv.glm1.null <- glm(Liver~1 , data = xplantL1[which(xplantL1$School == "UVA"), c("Liver" ,"LiverL1")], family = poisson)
anova(uva.liv.glm1.null,uva.liv.glm1 ,test="Chi")


# Dispersion test

sum(resid(uva.liv.glm1, type = "pearson")^2/uva.liv.glm1$df.residual)


#################################
# Correcting for overdispersion

# Quasi-poisson model
# unpooled model

uva.liv.glm2  <- glm(Liver~. , data = xplantL1[which(xplantL1$School == "UVA"), c("Liver" ,"LiverL1")], family = quasipoisson)

# Get the summary

summary(uva.liv.glm2)

# Model utility test

uva.liv.glm.null  <- glm(Liver~1 , data = xplantL1[which(xplantL1$School == "UVA"), c("Liver" ,"LiverL1")], family = quasipoisson)

anova(uva.liv.glm.null, uva.liv.glm2, test = "Chi")

# Diagnostics

par(mfrow = c(2,2))
plot(uva.liv.glm2)
par(mfrow = c(1,1))

# Serial Corrlation?

uva.liv.res <- residuals(uva.liv.glm2, type = "pearson")

par(mfrow = c(1,2))
acf(uva.liv.res)
pacf(uva.liv.res)
par(mfrow = c(1,1))

# prediction 

newLdata <- data.frame(LiverL1 = xplantC[27, "Liver"], School = "UVA")


liv.pred <- predict(uva.liv.glm2, newdata = newLdata, type = "response", se.fit = T)

# prediction plot

plot(xplantC$Year[which(xplantC$School == "UVA")], xplantC$Liver[which(xplantC$School == "UVA")], type = "b", xlim = c(1988, 2015))

lines(xplantC$Year[which(xplantC$School == "MCV")], xplantC$Liver[which(xplantC$School == "MCV")], type = "b", col = "green")


lines(xplantC$Year[which(xplantC$School == "Duke")], xplantC$Liver[which(xplantC$School == "Duke")], type = "b", col = "red2")


lines(xplantC$Year[which(xplantC$School == "UNC")], xplantC$Liver[which(xplantC$School == "UNC")], type = "b", col = "blue2")


points(2015, liv.pred$fit, pch = 19)

segments(2015, liv.pred$fit, 2015, liv.pred$fit + 1.96*liv.pred$se.fit)

segments(2015, liv.pred$fit, 2015, liv.pred$fit - 1.96*liv.pred$se.fit)


# Pooled Model
uva.liv.glm4 <- glm(Liver~ LiverL1 + School , data = xplantL1, family = quasipoisson)


# Diagnostics
par(mfrow = c(2,2))
plot(uva.liv.glm4)
par(mfrow = c(1,1))

# Serial correlation?

uva.liv.res <- residuals(uva.liv.glm4, type = "pearson")

par(mfrow = c(1,2))
acf(uva.liv.res)
pacf(uva.liv.res)
par(mfrow = c(1,1))



#*******************

#prediction

newLdata <- data.frame(nYears = 28, LiverL1 = xplantC[27, "Liver"], School = "UVA")

liv.pred <- predict(uva.liv.glm4, newdata = newLdata, type = "response", se.fit = T)

# Plot

plot(xplantC$Year[1:27], xplantC$Liver[1:27], type = "b", xlim = c(1988, 2015))

points(2015, liv.pred$fit, pch = 19)

segments(2015, liv.pred$fit, 2015, liv.pred$fit + 1.96*liv.pred$se.fit)

segments(2015, liv.pred$fit, 2015, liv.pred$fit - 1.96*liv.pred$se.fit)


#********************

# Compare pooled, unpooled and mixture model predictions

