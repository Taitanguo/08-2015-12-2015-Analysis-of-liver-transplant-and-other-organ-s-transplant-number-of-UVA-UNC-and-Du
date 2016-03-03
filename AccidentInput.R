

#******************************************************************
#
#			Input Functions for Accident Data
#
#*******************************************************************


# function to input files into a list



file.inputl <- function(my.path)
{
	my.dir <- getwd()
	setwd(my.path)
	my.files <- list.files(my.path)
	sapply(1:length(my.files), function(i)
	{
		read.csv(my.files[i])
	})
	
}

# Function to create a data frame as the combination of multiple data frames in a list.


combine.data <- function(Data.List, Vars)
	{
		DF <- rbind(Data.List[[1]][, Vars])
		for(i in 2:length(Data.List))
	{
		DF <- rbind(DF, acts[[i]][,Vars])
	}
	DF
	}


# Code to input files into the global environment

# First get the path for the files, here it is the working directory

# path <- getwd()

# Get the vector of data files

# my.files <- list.files(path)

# Read in all files with names "acts"

#  for(i in my.files) 
# { 
#   	nam <- paste("acts", substr(i, 14, 15), sep ='')
# 	 	assign(nam, read.csv(i))
# }
 



