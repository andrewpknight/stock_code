####################################################################
## This function just aggregates the sensor data very quickly ## 

# aggvars = a vector containing the variables you want to aggregate (e.g., e, t, zyx)
# groupvars = a vector containing the grouping variables that you want to use to aggregate
# d = the dataset
# stats = a vector containing which stats you want (e.g., mean, sd)
# suffs = a vector containing the suffix that you want appended for each metric (e.g., "_x", "_sd")

sensor.agg <- function(aggvars,groupvars,d, stats, suffs) {
	require(multilevel)
	group_list <- vector("list", length(groupvars))
	for(i in 1:length(groupvars)) {
		group_list[[i]] <- d[,groupvars[i]]
	}
	for(s in 1:length(stats)) {
		agg <- aggregate(d[,aggvars], by=group_list, as.character(stats[s]), na.rm=TRUE)
		colnames(agg) <- c(groupvars, paste(aggvars, as.character(suffs[s]), sep=""))
		if(s == 1) {
			ret <- agg
		} else {
			ret <- merge(ret, agg, by=groupvars)
		}
	}
	return(ret)	
}

####################################################################
## This function adds a variable called "button" to a sensor dataset. 

## The function assumes that the infile is in the structure that is created by the parse_eda_data.rb script. The function finds any row with missing data across records and marking as a button push. NOTE: This presumes that there is not other missing data in the dataset. So, I'd run this first, before merging with other survey data ##

sensor.add.button <- function(infile) {
	infile$button <- ifelse((is.na(infile$b) & is.na(infile$t)), 1, 0)
	return(infile)
}

####################################################################
## This function standardizes data for each individual in the dataset
## The function assumes that you have a file that contains a unique identifier
## for each individual in the dataset.
## stand_vars = a vector of variables that you want to standardize
## ind_id = the unique id variable for individuals
## d = the datafile
## suff = the suffix that you want to append to the standardized variable name

sensor.standardize <- function(stand_vars, ind_id, d, suff) {
	require(multilevel)
	# First, get the individual mean and standard deviation
	agg.x <- aggregate(d[,stand_vars], by=list(d[,ind_id]), mean, na.rm=TRUE)
	agg.sd <- aggregate(d[,stand_vars], by=list(d[,ind_id]), sd, na.rm=TRUE)	
	colnames(agg.x) <- c(ind_id, paste(stand_vars, "_x", sep=""))
	colnames(agg.sd) <- c(ind_id, paste(stand_vars, "_sd", sep=""))	
	agg <- merge(agg.x, agg.sd, by=ind_id)
	
	# Second, scale the original data, subtracting from the mean and dividing by the SD
	d <- merge(d, agg, by=ind_id, all.x=TRUE)
	
	for(i in 1:length(stand_vars)) {
		vname <- paste(stand_vars[i], suff, sep="")
		d[,vname] <- (d[,stand_vars[i]] - d[,paste(stand_vars[i], "_x", sep="")]) / d[,paste(stand_vars[i], "_sd", sep="")]
	}
	return(d)
}

####################################################################
## This function creates a motion variable based on the z, y, and x coordinates provided
## by the Q sensors. The motion variable is computed as the Euclidean distance between
## successive rows in the dataset for these three positions in space. 
## This function is used in the data import script. 
## d = data file that you're supplying
## The script returns the datafile with a new variable called "a" (for accelerometer)

sensor.create.motion <- function(d) {	
	vrs <- c("z", "y", "x")	
	# Go by individual; dataset must be properly sorted by individual and time before coming in here.
	d$a <- NA
	ids <- unique(d$indiv_b_id)
	count <- 1
	for(i in ids) {
		wrk <- d[d$indiv_b_id == i, ]
		wrk[2:length(wrk$indiv_b_id), c("a")] <- sqrt(rowSums(diff(as.matrix(wrk[,vrs]))^2, na.rm=TRUE))			
		if(count == 1) {
			dat <- wrk
		} else {
			dat <- rbind(dat, wrk)
		}
		count <- count + 1
	}
	return(dat)	
}

####################################################################
## This function imports a raw sensor dataset, which has been output from the 
## Ruby script that parses the eda files. This file will create a unique identifier for
## individuals based on the date the data were collected and the log file from the
## Q sensor. The script will output a datafile containing formatted time/date variables
## (i.e., bracelet_time, bracelet_start_time), as well as a motion variable calculated 
## using the function above.
## fname = the path to the csv file that you're importing (e.g., "./sensor_data.csv"
## second = Boolean indicating whether you want the script to aggregate to the second
## 			level. This is usually a good idea, unless you're at a fine resolution
## standardize = Boolean indicating whether you want the script to standardize across
## 				 individuals. If your indiv id is clean, this is a good idea. 

sensor.import.raw <- function(fname, second, standardize) {

	# Which variables are returned by default			
	ret.vars <- c("indiv_b_id", "bracelet_start_time", "bracelet_time", "second", "z", "y", "x", "b", "t", "e", "a")
	
	d <- read.delim(fname, sep=",", header=TRUE, as.is=TRUE, colClasses=c("character", "numeric", "character", rep("numeric",9)))

	# Create a unique individual id variable (based on the log file number)	
	d$bracelet_start_date <- substr(as.character(d$bracelet_start_time), 0, 10)
	d$indiv_b_id <- paste(d$bracelet_start_date, d$bracelet_id, d$file_number, sep="_")

	# Create the time variables
	d$second <- floor(d$sample_number/d$sampling_rate) + 1	
	d$bracelet_start_time <- as.POSIXct(d$bracelet_start_time, origin="1970-01-01 00:00:00")
	d$bracelet_time <- d$bracelet_start_time + d$second
	
	d <- sensor.create.motion(d)
	
	# If the user wants to bring this up to the second level of analysis, do so
	if(second) {
		d.agg <- sensor.agg(aggvars=c("z", "y", "x", "b", "t", "e", "a"),groupvars=c("indiv_b_id", "bracelet_time"),d=d, stats=c("mean"), suffs=c(""))	
		
		d.rec <- unique(d[,c("indiv_b_id", "bracelet_start_time", "second", "bracelet_time")])
		
		d <- merge(d.rec, d.agg, by=c("indiv_b_id", "bracelet_time"))
		d <- d[,ret.vars]
	}
	
	if(standardize) {
		d <- sensor.standardize(stand_vars=c("t", "e", "a"), ind_id=c("indiv_b_id"), d=d, suff="_z")	
	}	
	return(d)
}

####################################################################
## This function creates a decent looking plot of sensor data. 
## x = the vector of x values (e.g., bracelet_time or second)
## y = the metric you want to graph (e.g., e_z, t, a)
## ylim = a vector containing the min and the max for the y axis

sensor.graph.line <- function(x, y, fname, ylim) {
	quartz(width=10, height=5, type="pdf", file=fname)
	par(family="Tahoma", bg="white", mar=c(4,4,1,1), mgp=c(3,1,0))
	
	plot(x,y,type="n", ylim=ylim, xlim=c(min(x), max(x)), axes=TRUE, ylab="Electrodermal Activity", xlab="Time", col="dark green", bty="n", cex.axis=.75)
	
	lines(x,y, lwd=2.5, col="dark green")
	xx <- c(x, rev(x))
	yy <- c(rep(min(y), length(y)), rev(y))
	polygon(xx, yy, col="light green", border=FALSE) 
	
	dev.off()

}

