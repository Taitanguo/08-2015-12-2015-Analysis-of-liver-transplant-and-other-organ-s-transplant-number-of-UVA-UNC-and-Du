# Coded By Prof. Donald Brown, SIE, UVA

#******************************************
#
#		PC Variance Functions
#
#******************************************

# I have written these functions to find the 
# principal components that capture the desired 
# amount of variance in the data.

seq.sum <- function(x)
{
	s <- c(x[1],rep(0, length(x)-1))
	for(i in 2:length(x))
	{s[i] <- x[i] + s[i-1]}
	s
	}

var.comp <- function(pc.obj, goal.percent)
{
	var.pc <- pc.obj$sdev^2
	var.prop <- var.pc/sum(var.pc)
	max(which(round(100*(seq.sum(var.prop))) <= goal.percent))
	}

#******************************************
#
#		PC Regression
#
#******************************************


# I have written this function to make it easy to generate
# the regression from the principal components that 
# capture the desired amount of variance.

pc.reg <- function(pc.obj, goal.percent, response)
{
	k <- var.comp(pc.obj, goal.percent)
	lm(r ~ ., data = data.frame(pc.obj$scores[,1:k], r = response))
	}
