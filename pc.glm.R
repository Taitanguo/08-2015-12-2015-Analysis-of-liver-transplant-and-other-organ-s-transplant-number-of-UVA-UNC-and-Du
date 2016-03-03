# I have written these functions to find the 
# principal components that capture the desired 
# amount of variance in the data.

# Run these functions first (one time).

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
	max(which(round(100*(seq.sum(var.prop)),2) <= goal.percent))
	}


# I have written this function to make it easy to generate
# the logistic regression from the principal components that 
# capture the desired amount of variance.

pc.glm <- function(pc.obj, goal.percent, response)
{
	k <- var.comp(pc.obj, goal.percent)
	glm(r ~ ., data = data.frame(pc.obj$scores[,1:k], r = response), family = binomial)
	}
	
	
pc.null <- function(pc.obj, goal.percent, response)
{
	k <- var.comp(pc.obj, goal.percent)
	glm(r ~ 1, data = data.frame(pc.obj$scores[,1:k], r = response), family = binomial)
	}



# This function predicts with the pc.glm object


predict.pc.glm <- function(pcglm.obj, pc.obj, ndata, type = "response")
{
	
	M <- scale(ndata, center = pc.obj$center, scale = pc.obj$scale)
	PC <- data.frame(M %*% pc.obj$loadings)
	predict(pcglm.obj, newdata = PC, type = type)
	}

#
