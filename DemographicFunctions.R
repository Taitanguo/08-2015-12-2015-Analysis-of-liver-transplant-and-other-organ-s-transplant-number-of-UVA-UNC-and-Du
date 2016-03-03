

#**************************************
#
#  Functions to calculate entropy
#
#**************************************

#	M the data matrix


standard <- function(X, S)
{
	X/S
	}

entropy <- function(M)
{
	S <- apply(M, 1, sum)
	V <- standard(M, S= S)
	ent <- -V*log(V)
	ent[which(V == 0)] <- 0	
	apply(ent, 1, sum)
	}

#**************************************
#
#  Function to format the data
#  for use in entropy calculations
#
#**************************************

#	RN the list of all possible names
# 	N the names in the select data set
#	M the data set

formdata <- function(RN, N, M)
{
	Z <- matrix(nrow = nrow(M), ncol = length(N))
	for(i in 1:length(N))
	{
		Z[,i] <- M[,paste(N[i])]
		}
	L<- setdiff(RN, N)
	if(length(L) > 0)
	{
		Z <- cbind(Z, matrix(rep(0, length(L)*nrow(M)), nrow = nrow(M)))
		}
	Z
	}
	
	
#**************************************
#
#  Function to format the data by organ name
#  for use in modeling
#
#**************************************

#	oname is the organ name to be subsetted
# 	data is a data frame with organ transplants by 
#	category, e.g., ethnic group or age

subdata <- function(oname, M)
{
	N <- nchar(oname)
	Mnames <- colnames(M)[which(substr(colnames(M), 1, N) == oname & nchar(colnames(M)) < 11)]
	M[,Mnames]
	
	}



#**************************************
#
#  Function to calculate the % white 
#  in the population
#
#**************************************

#	M the matrix of ethnic groups 
#   assumed in the order specified in OPTN



pw <- function(M)
{
	S <- apply(M, 1, sum)
	M[,1]/S
	}
