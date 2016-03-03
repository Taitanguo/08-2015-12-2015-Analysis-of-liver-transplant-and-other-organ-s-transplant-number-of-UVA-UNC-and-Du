#***************************************************************
#
# Bootstrap Regression and Time Series
#
#***************************************************************

#***************************************************************
#
#  Read in the data
#
#***************************************************************

# Read data
r11xplant <- read.table("~/Desktop/linear/transplant/R11xplant.csv", sep = ",", header = T)

r11donor<-read.table("~/Desktop/linear/transplant/R11donor.csv", sep = ",", header = T)

uva <- read.table("~/Desktop/linear/transplant/UVAxplant.csv", sep = ",", header = T)

duke <- read.table("~/Desktop/linear/transplant/Dukexplant.csv", sep = ",", header = T)

mcv <- read.table("~/Desktop/linear/transplant/MCVxplant.csv", sep = ",", header = T)

unc <- read.table("~/Desktop/linear/transplant/UNCxplant.csv", sep = ",", header = T)

#setwd(sourcedir)

#  Source the bootstrapping functions
library(boot) #If you don't have this library, install it by: install.packages('boot')
source("~/Desktop/linear/code/TSbootfunctions.R")

#Kidney Model:

##Build a linear model, uva.kid.lm that predicts uva kidney transplants by region 11 kidney donors from 1988-2014
uva.kid.lm <- lm(uva$Kidney[-28]~r11donor$Kidney[-28])

##  Is your model significant?
summary(uva.kid.lm)

# Plot diagnostics uva.kid.lm.  What do you examine?
par(mfrow=c(2,2))
plot(uva.kid.lm)
par(mfrow=c(1,1))

#  Bootstrapping the linear model

# Get the fitted values from the regression model uva.kid.lm and store in uva.kfit
# Hint: use the fitted() function with the linear model as the option
uva.kfit <- fitted(uva.kid.lm)

# Get the residuals from the regression model uva.kid.lm and store in uva.ke
# Hint: use the residuals() function with the linear model as the option
uva.ke <- residuals(uva.kid.lm)

#  Get the regression model for uva.kid.lm- model.matrix()
uva.mod <- model.matrix(uva.kid.lm)


# Use the RTSB function (from TSbootfunctions.R) to obtain the bootstrap
# distribution for the coefficients. This may take a few minutes.
# Be patient.	

uva.kid.boot <- RTSB(uva$Kidney[-28], r11donor$Kidney[-28], uva.kfit, uva.ke, uva.mod,2000)


# What are the estimates in uva.kid.boot?
uva.kid.boot
#the first one is a intersect and second is the another one
summary(uva.kid.lm)
summary(uva.kid.boot)
summary(uva.kid.boot$t[,2])
# Get the 99% CI for uva.kid.boot

boot.ci(uva.kid.boot, .99)

# Get the 95% CI for uva.kid.boot 
#BCA is adjusted and the mean is not skewed

boot.ci(uva.kid.boot, .95)

#	Plot the results for the coeffiecient for region 11 donors
plot(uva.kid.boot, index = 2) 

#	A set of configurable plots
par(mfrow = c(1,2))
hist(uva.kid.boot$t[,2], main = "Region 11 Donors",xlab ="Coefficient Values",   col = "steelblue", breaks = 50)
qqnorm(uva.kid.boot$t[,2])
qqline(uva.kid.boot$t[,2])
par(mfrow = c(1,1))

#  Bootstrapping TS

# Evaluating residual correlation from the model uva.kid.lm
# Hint: use the acf() and pcf()

par(mfrow = c(1,2))
acf(uva.kid.lm$residuals)
pacf(uva.kid.lm$residuals)
par(mfrow = c(1,1))

#	Fit an ar model to the residuals using the yule-walker method
diff.ar.kid <- ar(uva.kid.lm$residuals, method = "yule-walker") 

# How many autoregressive terms are needed?
diff.ar.kid

# If we use diff.ar.kid why we choose the 3:27 and 2:26 and 1:25
uva.kid.lm2<- lm(uva$Kidney[3:27]~r11donor$Kidney[3:27]+ uva.kid.lm$residuals[2:26] + uva.kid.lm$residuals[1:25])

summary(uva.kid.lm2)

# The problem here is we have only a few observations (26)

par(mfrow=c(2,2))
plot(uva.kid.lm2)
par(mfrow=c(1,1))

# Let's also try an AR(1) model of the residuals to regression, linear model.
uva.kid.lm2<- lm(uva$Kidney[2:27]~r11donor$Kidney[2:27]+ uva.kid.lm$residuals[1:26])

summary(uva.kid.lm2)
par(mfrow=c(2,2))
plot(uva.kid.lm2)
par(mfrow=c(1,1))

#ar(2) is the same as the ar(1)

#try model ar(2)

#get the fitted values from the regression model
uva.kfit2<-fitted(uva.kid.lm2)

#get the residuals from the regressoin model

uva.ke2 <- residuals(uva.kid.lm2)

#get the regression model
uva.mod2 <- model.matrix(uva.kid.lm2)

#Use the RTSB function as before on the new model to obtain the bootstrap
uva.kid.boot2<-RTSB(uva$Kidney[3:27], r11donor$Kidney[3:27],uva.kfit2,uva.ke2,uva.mod2,2000)

#if we use diff.ar?? what is diff.ar 

#plot the results for the coeffiecient for the region 11 donors
plot(uva.kid.boot2,iindex = 2)

plot(uva.kid.boot2,iindex = 3)
plot(uva.kid.boot2,iindex = 4)
plot(uva.kid.boot2,iindex = 5)

# get the 95% CI for uva.kid.boot2  what this index means???
boot.ci(uva.kid.boot2,.95,index=1)
boot.ci(uva.kid.boot2,.95)

######## Repeat for liver transplants

# Build a linear model to predict uva liver transplants in terms of region 11 donors




# Model significance?




# Generate the diagnostic plots. Do you see any problems?





# Estimate the liver model with bootstrapping (by residuals). Is b1 significant?



#   Bootstrapping LM




# What is the 95% CI of r11donor?




# Plot the distribution of B1



# Time series models for residuals of liver model




# Generate the ACF and PACF plots of the residuals from uva.liver.lm. What's your conclusion?  



# Fit an ar model to the residuals, what order do you select?




# Bootstrap your time series model. Are the coefficients significant?





# Plot the results for the coeffiecient for region 11 donors and time series components




# What are the confidence intervals for each of the parameters?



