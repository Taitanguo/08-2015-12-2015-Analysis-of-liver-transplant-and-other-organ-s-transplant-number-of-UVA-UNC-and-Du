#************************************************************
#
#				SPAM Filter 1 
#			Graphical Analysis of SPAM 

#cancerdata filter
cancer <- read.table("~/Desktop/linear/cancerdata/wdbc.data.txt", sep = " ", header = F)

#**************************************************************
source("~/Desktop/linear/code/SPM_Panel.r")
source("~/Desktop/linear/code/PCAplots.R")
source("~/Desktop/linear/code/FactorPlots.r")
source("~/Desktop/linear/code/pc.glm.R")
summary(cancer)

uva.pairs(spam[,c(1:10,58)])

uva.pairs(log(spam[,c(1:10,58)]+.00001))

uva.pairs(log(spam[,c(48:57,58)]+.00001))

par(mfrow = c(3,3))
for(i in 1:9)
{
  boxplot(spam[,i]~spam[,58], xlab = "Type of Email", ylab = "Count", main = paste("Variable", i))  
}
par(mfrow = c(1,1))

par(mfrow = c(3,3))
for(i in 49:57)
{
  boxplot(spam[,i]~spam[,58], xlab = "Type of Email", ylab = "Count", main = paste("Variable", i))  
}
par(mfrow = c(1,1))
#****************************************************
#
#		Principal Components
#
#****************************************************

spam.pca = princomp(spam[,1:57], cor = T)
biplot(spam.pca)

#1754
summary(spam[,which(spam.pca$loadings[,2]>0.2)])
spam[1754,which(spam.pca$loadings[,2]>0.2)]
# Obtain the biplot.fact of your principal components???? what is a factplot?
biplot.fact(spam.pca,spam[,58])
legend(20,15, legend = c("Spam", "Ham"), pch = c(18,19), col = c("red", "blue"))
#***************************************************************
#
#		Log of Predictors
#
#***************************************************************
lspam <- log(spam[,-58] + .1)
lspam[,58] <- spam[,58]

#1. Obtain box plots for log transforms of variables 1-9 and then 49-57.

par(mfrow = c(3,3))
for(i in 1:9)
{
  boxplot(lspam[,i]~lspam[,58], xlab = "Type of Email", ylab = "Count", main = paste("Variable", i))  
}
par(mfrow = c(1,1))

par(mfrow = c(3,3))
for(i in 49:57)
{
  boxplot(lspam[,i]~lspam[,58], xlab = "Type of Email", ylab = "Count", main = paste("Variable", i))  
}
par(mfrow = c(1,1))



#2. Obtain scatter plot matices for log transforms of variables 1-9 with variable 58 and then log transforms of variables 49-57 with variable 58.
par(mfrow = c(3,3))
for(i in 1:9)
{
  scatter.smooth(lspam[,i]~lspam[,58], xlab = "Type of Email", ylab = "Count", main = paste("Variable", i))  
}
par(mfrow = c(1,1))
#3. Obtain the principal components for the log transform of variables 1-57. 
#   Look at the biplot and explain what you see.

lspam.pc = princomp(lspam[,1:57], cor = T)
biplot(lspam.pc)

#4. Obtain the biplot.fact of your principal components for the log transformed variables.
#   Explain what you see.
biplot.fact(lspam.pc,lspam[,58])
legend(-7,10, legend = c("Spam", "Ham"), pch = c(18,19), col = c("red", "blue"))
#***************************************************************
#get a glm model 
spam.glm <- glm(V58~., data = spam, family = binomial)
spam.null <- glm(V58~1, data = spam, family = binomial)
anova(spam.null, spam.glm, test = "Chi")
# what's this mean 58 means that 58 is a response variable. and data points out the variable that we should choose from the 55-58,we must contain the 58
spam.cap<-glm(V58~.,data = spam[,55:58], family = binomial )
summary(spam.cap)
#remove one variable in order to test wether a model without that is useful
spam.no57 <- update(spam.glm, .~.-V57, data = spam)
anova(spam.no57, spam.glm, test = "Chi")
# predict model logodds
predict(spam.glm, newdata = spam[1,])
# odds 1
exp(predict(spam.glm, newdata = spam[1,]))

# compare the drop1 chi square test to the approximate Gaussian test in summary.
library(MASS)
# drop one by one
drop1(spam.glm, response~., test = "Chi", data = spam)
#  Compare the step model with capital letter predictors to the capital letter model
# both delete
step.cap <- step(spam.cap, data = spam, family = binomial)


# GLM with Interactions
spam.glm2 <- glm(V58~. + (V5+V6+V7)^2, data = spam, family = binomial)
AIC(spam.glm)
BIC(spam.glm)

#  Principal Components Regression     ?????????????????????????????????

# obtain the principal components for the predictors with a correlation matrix

spam.pca <- princomp(spam[,-58], cor = T)
summary(spam.pca)
# Scree plot

screeplot(spam.pca)

# to see how many components are needed for 90% of the variance we use

var.comp(spam.pca, 90)
# Use pc.glm() to get the principal component regression results.

spampca.glm98 <- pc.glm(spam.pca, 98, spam[,58])

# Do a model utility test starting with pc.null()

spampc.null <- pc.null(spam.pca, 98, spam[,58])

anova(spampc.null, spampca.glm98, test = "Chi")


#	Evaluation of Generalized Lnear Models 

library(lattice)
#xyplot ??? interaction plot??????????????? ??????????????????interaction???

xyplot(V58~V1 | cut(V57, breaks = 3), data = Spam$train )


interaction.plot(cut(spam$V1, breaks = 2), cut(spam$V52, breaks = 3), spam$V58, ylim = c(0,1))


# Here is an example of how to          ??????????????????????????????
# obtain a GLM with principal components
# accounting for 90% of the variance
# using the training data.


spam.pca <- princomp(Spam$train[,-58], cor = T)

spampca.glm90 <- pc.glm(spam.pca, 90, Spam$train[,58])

summary(spampca.glm90)
