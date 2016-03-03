

#****************************************
#
#	Functions for Bootstrap Predictions 
#	with LME objects
#
#****************************************



lmepred <- function(ResponseName, DataF, modelform, fit, resid, predata, RE, REV)
{
	i <- which(DataF[,RE] == REV)
	resid[i] <- sample(resid[i], replace = T)
	DataF[i,ResponseName] <- fit[i] + resid[i]
	M1 <- lmer(modelform, data = DataF, REML = T)
	mm <- model.matrix(terms(M1), predata)
	mm %*% t(coef(M1)[[1]][which(row.names(ranef(M1)[[1]]) == REV),])
	
}


lmeboot <- function(ResponseName, DataF, modelform, fit, resid, predata, RE, REV, num)
{
	sapply(1:num, function(x)
{
myboot <- rep(NA, length(x))
myboot[x] <- lmepred(ResponseName, DataF, modelform, fit, resid, predata, RE, REV)	
}  )
}




