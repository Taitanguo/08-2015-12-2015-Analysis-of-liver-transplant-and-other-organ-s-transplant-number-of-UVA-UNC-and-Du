#***************************************************************
#
#      Session 22: Simulation & Bootstrapping
#     	 A Case Study: Transplant Center
#   		   
#***************************************************************

#***************************************************************
#
#  Read in the data
#
#***************************************************************

#Set working directory

#Read data
r11xplant <- read.table("~/Desktop/linear/transplant/R11xplant.csv", sep = ",", header = T)

r11donor<-read.table("~/Desktop/linear/transplant/R11donor.csv", sep = ",", header = T)

uva <- read.table("~/Desktop/linear/transplant/UVAxplant.csv", sep = ",", header = T)

duke <- read.table("~/Desktop/linear/transplant/Dukexplant.csv", sep = ",", header = T)

mcv <- read.table("~/Desktop/linear/transplant/MCVxplant.csv", sep = ",", header = T)

unc <- read.table("~/Desktop/linear/transplant/UNCxplant.csv", sep = ",", header = T)

#setwd(sourcedir)
#    Source the bootstrapping functions
library(boot) #If you don't have this library, install it by: install.packages('boot')
source("~/Desktop/linear/code/TSbootfunctions.R")

#***************************************************************
#
# Part 1: Basic Statistics (10 mins)
#
#***************************************************************

#Step 1.1 Compare the performance of UVa with MCV kidney transplants
uva.kidney<-uva$Kidney
mcv.kidney<-mcv$Kidney

#   Get the distribution of uva$Kidney, mcv$Kidney, r11donor$Kidney. What do you observe? 
kidney <-data.frame(uva$Kidney,mcv$Kidney, r11donor$Kidney)
uva.pairs(as.matrix(kidney))

#   On average, how many kidney transplants are performed at UVa per year? MCV?
mean(uva$Kidney)
#   Perform a paired t-test between uva and mcv kidney transplants
t.test(uva$Kidney, mcv$Kidney, paired = T)
#Step 1.2 Compare the performance of UVa with Duke kidney transplants 

#   Get the distribution of uva$Kidney, Duke$Kidney, r11xplant$Kidney. What do you observe?
kidney2 <-data.frame(uva$Kidney,duke$Kidney, r11xplant$Kidney)
uva.pairs(as.matrix(kidney2))
#   On average, how many kidney transplants are performed at UVa per year? Duke?
mean(duke$Kidney)
#   Perform a paired t-test between uva and mcv kidney transplants
t.test(uva$Kidney, duke$Kidney, paired = T)
#Step 1.3 Use bootstrapping to test the hypothesis: there is not a significant difference between UVa and MCV kidney transplants.
kidney.diff<-ts(uva$Kidney-mcv$Kidney,1988,2014)

bs.mean<-function(x,i)
{
  return(mean(x[i]))
}

bs.kidney.diff<-boot(kidney.diff,bs.mean,R=1000)


#   What are the standard errors of the mean? 
bs.kidney.diff
#   What're the 95% confidence intervals? 
boot.ci(bs.kidney.diff,0.95,type=c('bca','perc'))
#   Do you accept or reject the null hypothesis?
#Reject because 0 is not in interval
#Step 1.4 Use bootstrapping to test the hypothesis: There is not a significant difference between UVa and Duke kidney transplants.
kidney.diff1<-ts(uva$Kidney-duke$Kidney,1988,2014)

bs.mean<-function(x,i)
{
  return(mean(x[i]))
}

bs.kidney.diff1<-boot(kidney.diff1,bs.mean,R=1000)
#   What are the standard errors of the mean? 
bs.kidney.diff1
#   What are the 95% confidence intervals? 
boot.ci(bs.kidney.diff1,0.95,type=c('bca','perc'))
#   Do you accept or reject the null hypothesis?
#Reject
#Step 1.5 Get the scatter plot matrix with the above 6 variables. Describe what you observe. 
#You can use either uva.pairs() {source("SPM_Panel.R")} or pairs().
Kidney<-data.frame(uva$Kidney,duke$Kidney,mcv$Kidney,unc$Kidney, r11donor$Kidney)
uva.pairs(as.matrix(Kidney))


#***************************************************************
#
# Part 2: Linear Regression Models (15 mins)
#
#***************************************************************

# Test the hypothesis: There is no difference between the forecasted numbers of kidney 
# transplants that will be performed at UVA and at MCV in 2014.
uvamcv.kidney.lm<- lm((uva$Kidney-mcv$Kidney)~r11donor$Kidney)
summary(uvamcv.kidney.lm)
### Low adj- R^2, Model significant at 99% Confidence level, coeeficient is significant
# Step 2.1 Build a linear model: 
# uva$Kidney-mcv$Kidney = b0+b1*r11donor$Kidney+e. Call it uvamcv.kidney.lm.
# Analyze the result: R^2, model utility test, t-tests, etc.
uvamcv.kidney.lm<- lm((uva$Kidney-mcv$Kidney)~r11donor$Kidney)
summary(uvamcv.kidney.lm)
### Low adj- R^2, Model significant at 99% Confidence level, coeeficient is significant


#Step 2.2 Generate the diagnostic plots. Do you see any problems?

par(mfrow=c(2,2))
plot(uvamcv.kidney.lm)
par(mfrow=c(1,1))
## Due to few points, the resiual vs fitted plot is not symmetrical around 0. 
## QQ plot is not normal near the tails; Further, there are influential points
# Step 2.3 Estimate the model with bootstrapping (by residuals). Is b1 significant?
#    Get the regression model
# Get the fitted values from the regression model
uvamcv.fit <- fitted(uvamcv.kidney.lm)

#    Get the residuals from the regression model
uvamcv.re <- residuals(uvamcv.kidney.lm)

#    Get the regression model
uvamcv.mod <- model.matrix(uvamcv.kidney.lm)

#   Bootstrapping LM
uvamcv.kidney.boot <- RTSB(uva$Kidney[-28]-mcv$Kidney[-28], r11donor$Kidney[-28], uvamcv.fit, uvamcv.re, uvamcv.mod,5000)   #5000 samples
uvamcv.kidney.boot$t
sqrt(abs(var(uvamcv.kidney.boot$t)))

#    95% CI of r11donor
boot.ci(uvamcv.kidney.boot, .95, index=2)

#    Distribution of b1
par(mfrow = c(1,2))
hist(uvamcv.kidney.boot$t[,2], main = "Region 11 Donors",xlab ="Coefficient Values",   col = "steelblue", breaks = 50)
qqnorm(uvamcv.kidney.boot $t[,2])
qqline(uvamcv.kidney.boot $t[,2])
par(mfrow = c(1,1))





#Step 2.4* What about Duke? Repeat the above steps and compare the results. 

# Test the hypothesis: There is no difference between the forecasted numbers of kidney 
# transplants that will be performed at UVA and at Duke in 2014. use this
# test that last year's data can predict next year's data. 
par(mfrow = c(1,2))
acf(uvamcv.le)
pacf(uvamcv.le)
par(mfrow = c(1,1))
kidney <-data.frame(uva$Kidney,duke$Kidney,mcv$Kidney, r11donor$Kidney, r11xplant$Kidney)
uva.pairs(as.matrix(kidney))
#ACF - significant lags stop after lag 1, sinusoidal decay
#PACF - no significant lags



#***************************************************************
#
# Part 3: Time Series Models (20 mins)
#
#***************************************************************
#Step 3.1 Generate the ACF and PACF plots of the residuals from your part 2 linear model 
#for UVA and MCV kidney transplants. What's your conclusion?  


#Step 3.2 Based on the above ACF and PACF plots, what time series do you suggest to model the residuals?

#    Fit an ar model to the residuals

diff.ar.kid <- ar(uvamcv.kidney.lm$residuals, method = "yule-walker")
diff.ar.kid
#Use AR(1)

#    Add the AR model of the residuals to regression linear model. 
#    Call this model uvamcv.kidney.lm2. Analyze the regression results.

uvamcv.kidney.lm2<- lm(uva$Kidney[2:27]-mcv$Kidney[2:27]~r11donor$Kidney[2:27]+ uvamcv.kidney.lm$residuals[1:26])

summary(uvamcv.kidney.lm2)
# build new data frame 
xplantL1 <- data.frame(xplantC[-seq(1,4*27, by =27),], KidneyL1 = xplantC$Kidney[-seq(27, 4*27, by = 27)])
head(xplantL1[,c("Kidney", "KidneyL1")])
#build a prediction model

xplant.mcvuva.lme5 <- lm(uva$Kidney[2:27]-mcv$Kidney[2:27] ~KidneyL1[2:27], data = xplantL1)

par(mfrow=c(2,2))
plot(xplant.mcvuva.lme5)
par(mfrow=c(1,1))
# Then build a bootstrap model
# Get the fitted values from the regression model
uvamcv.fit <- fitted(xplant.mcvuva.lme5)

#    Get the residuals from the regression model
uvamcv.re <- residuals(xplant.mcvuva.lme5)

#    Get the regression model
uvamcv.mod <- model.matrix(xplant.mcvuva.lme5)
uvamcv.kidneylag1.boot <- RTSB(uva$Kidney[-28]-mcv$Kidney[-28],KidneyL1[2:27], uvamcv.fit, uvamcv.re, uvamcv.mod,5000)

#   Generate the diagnostic plots. Do you see any problems?

par(mfrow=c(2,2))
plot(uvamcv.kidney.lm2)
par(mfrow=c(1,1))

#Step 3.3 Bootstrap the above time series model. Are the coefficients significant?

#    Get the fitted values from the regression model


#    Get the residuals from the regression model


#    Get the regression model


#     Use the RTSB function to obtain the bootstrap


#     The estimates


#    Plot the results for the coeffiecient for region 11 donors


#    Plot the results for the coeffiecient for time series components



#Step 3.5* What about Duke? Repeat the above steps and compare the results. 


#***************************************************************
#
# Part 4: Predicting Differences in Kidney Transplants Part 1 (15 mins)
#
#***************************************************************

#Step 4.1 Build an AR model to predict the difference in 2014


#Step 4.2 Use the predict function with ar model to forecast 2014 differences between UVA and MCV


#Step 4.3 Bootstrapping the difference of UVa and MCV in 2014
#    To obtain a bootstrap estimate of the prediction for 2014
#    use the TSB function in the source file.
#    It takes three arguments:
#    tsint - the time series
#    oth.arg - the data for the new estimate
#    boot.number- number of replications (default=1000)


#Step 4.4 What about Duke? Repeat the above steps and compare for Duke.


#***************************************************************
#
# Part 5*: Predicting Differences in Kidney Transplants Part 2 (15 mins)
#
#***************************************************************

#Step 5.1 Develop an AR model of region 11 kidney donors


#Step 5.2 Forecast the R11 donors and standard errors for 2014 using your ar model from step 5.1


#Step 5.3 Use the linear model from part 3.3 combined with the forecast of region 11 kidney donors
#to forecast the differences in number of kidney transplants between UVa and MCV for 2014.
#Use the predict() function

#   Creating the new data frame


#   Predict the linear model with the time series


#Step 5.4* Bootstrap the Forecast from the linear model combined with the forecast of region 11 kidney donors to forecast the differences
#in number of kidney transplants between UVa and MCV for 2014.


#   Bootstrap prediction using RFB function


#   Bootstrap plot 


#   Bootstrap confidence intervals


#Step 5.5 Plot the current and predictions for each value along with the confidence intervals


#Step 5.6 What about Duke? Repeat the above steps for Duke.