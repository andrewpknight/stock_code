## Setup ##
get.cross.ts <- function(i,j,d,y) {
	# This function inputs two id numbers, a dataset, and a metric and outputs a matched time series between the two individuals.
	# i = ego id
	# j = alter id
	# d = dataset
	# y = metric / variable
	
	# ego time series
	d1 <- d[d$indiv_id == i & d$bracelet_time >= d$preso_start & d$bracelet_time <= d$qa_end, c("bracelet_time", y)]
	
	# alter time series
#	d2 <- d[d$indiv_id == j & d$bracelet_time >= d$preso_start & d$bracelet_time <= d$qa_end, c("bracelet_time", y)]	
	d2 <- d[d$indiv_id == j, c("bracelet_time", y)]
	
	# Merge together ego and alter into a single dataset
	d <- merge(d1,d2,by=c("bracelet_time"), sort=T)
	colnames(d) <- c("bracelet_time", "ts1", "ts2")
		
	# Return the merged time series
	return(d[,c("ts1", "ts2")])	
}

get.ts <- function(i,d,y) {
	# This function inputs one id number, a dataset, and a metric and outputs the time series for the individual
	# i = ego id
	# d = dataset
	# y = metric / variable
	
	# ego time series
	d1 <- d[d$indiv_id == i & d$bracelet_time >= d$preso_start & d$bracelet_time <= d$qa_end, c("bracelet_time", y)]
	colnames(d1) <- c("bracelet_time", "ts1")		
	# Return the merged time series
	return(d1[,c("ts1")])	
}

get.delay.parm <- function(ts, max.delay) {
	# This function inputs one timeseries and, using the average mutual information method, outputs the delay that should be used for embedding the time series in phase space through the time-delayed embedding method.
	require(tseriesChaos)
	# ts = time series	
	# Get the delay using the average mutual information approach
	o <- mutual(ts, partitions=floor(length(ts)*.01), lag.max=100, plot=F)
	delay <- min(which(diff(sign(diff(as.numeric(o))))==-2))
	if(delay > max.delay ) delay <- max.delay
	return(delay)	
}

get.embed.parm <- function(ts, del, max.embed) {
	# This function inputs one timeseries and an appropriate delay. Using the false nearest neighbor method, the function outputs the number of embedding dimensions that are needed to adequately represent the dynamics in phase space
	require(tseriesChaos)
	# Get the false nearest neighbor for 15 dimensions
    fnn = false.nearest(ts, m = max.embed, d = del, t = 0, 
        eps = sd(ts)/10)
    
    # Get the number of dimensions that yields the first minimum
    em.dims = as.numeric(which(fnn[1, ] == min(fnn[1, 
        ], na.rm = TRUE)))
    
    # if there are multiple dimensions with the same value, return only one
    if (length(em.dims) > 1) {
        em.dims = em.dims[1]
    }    
    return(em.dims)
}

time.delay.embed <- function(x,m,d) {
	n = length(x) - (m-1)*d
	res = matrix(0,n,m)
	for (i in 1:m) res[,i] = x[((i-1)*d+1):(n+(i-1)*d)]
	res
}

rp.write <- function(i,j, inp,v) {
	fname1 <- "~/Dropbox/startup_weekend/analyses/recur/mem1.dat"
	fname2 <- "~/Dropbox/startup_weekend/analyses/recur/mem2.dat"
	d <-	get.cross.ts(i,j,inp,v)	
	if(length(d$ts1) >= 120) {
		write.table(d$ts1, fname1, sep="", col.names=F, row.names=F) 
		write.table(d$ts2, fname2, sep="", col.names=F, row.names=F)  
	} else {
		write.table("", fname1, sep="", col.names=F, row.names=F) 
		write.table("", fname2, sep="", col.names=F, row.names=F)  		
	}		
}

get.threshold.parm <- function(i,j,d,y, parms.i, parms.j, include.diag=TRUE) {
	require(fields)
#	parms.i <- parms[parms$indiv_id == i, c("delay", "embed")]
#	parms.j <- parms[parms$indiv_id == j, c("delay", "embed")]	
	embed.parm <- max(parms.i$embed, parms.j$embed)
	delay.parm <- mean(c(parms.i$delay, parms.j$delay), na.rm=TRUE)			
	d <- get.cross.ts(i,j,d,y)
	n <- length(d$ts1)
	if(n > (embed.parm-1)*delay.parm) {
		# embed the cross recurrence time series
		ts1.xyz <- time.delay.embed(d$ts1, embed.parm, delay.parm)
		ts2.xyz <- time.delay.embed(d$ts2, embed.parm, delay.parm)
		D <- rdist(ts1.xyz, ts2.xyz)
		D.rescale <- 100*D/max(D)
		recurrence.rate <- NULL
		for(parm in 1:100) {
			ind = which(D.rescale < parm, arr.ind=T)			
			# exclude the diagonal if requested
			if(include.diag == FALSE) {
				# remove the diagonals from the recurrence count
				ind <- data.frame(ind)
				ind <- ind[ind$row != ind$col, ]
				numer <- length(ind$row)
				denom <- length(D)-sqrt(D)
				recurrence.rate[parm] <- 100*numer/denom
			} else {
				recurrence.rate[parm] = 100*length(ind[,1])/length(D)				
			}
			if(recurrence.rate[parm] >= 1.7) break
		}
		threshold.parm <- parm
		threshold.parm.unscaled <- parm*max(D)/100
		rr <- recurrence.rate[parm]
	} else {
		threshold.parm <- NA
		threshold.parm.unscaled <- NA		
		rr <- NA
	}			
	res <- c(i,j,n,delay.parm, embed.parm, threshold.parm, threshold.parm.unscaled, rr)	
	res
}

cross.recurr.matrix <- function(d, parms, distance=FALSE, ret.matrix=FALSE) {
	require(fields)
	n <- length(d$ts1)
	# embed the time series
	ts1.xyz <- time.delay.embed(d$ts1, parms$embed, parms$delay)
	ts2.xyz <- time.delay.embed(d$ts2, parms$embed, parms$delay)	
	D <- rdist(ts1.xyz, ts2.xyz)
	if(distance == TRUE) {
		return(D)
	} else {
		if(ret.matrix==TRUE) {
			D <- ifelse(D < parms$radius_unscaled, 1, 0)		
			return(D)
		} else {
			ind <- data.frame(which(D < parms$radius_unscaled, arr.ind=T))		
			return(ind)
		}
	}
}


cross.recurr.plot <- function(d,parms, fname, saveplot, addline) {
	require(fields)	
	make.plot <- function(d, ind, axis.length, addline) {
		if(addline) {
#			nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(3,1), c(1,3), TRUE)
			nf <- layout(matrix(c(2,1), 2, 1, byrow=TRUE), widths=c(5), heights=c(5,1), TRUE)
#			layout.show(nf)
			par(mar = c(0,3,0,1))
			plot(d[1:axis.length, c("ts1")], type="l", col="blue", axes=FALSE, ylim=c(min(d[,c("ts1", "ts2")]), max(d[,c("ts1", "ts2")])))
			lines(d[1:axis.length, c("ts2")], lty=1, col="red")						
			
			par(mar = c(3,3,1,1))		
		}
		
		plot(0, xlim=c(0, axis.length), ylim=c(0, axis.length), type="n", main="Cross Recurrence Plot", xlab="i", ylab="j")
		points(ind$row, ind$col, pch=".", col="#333333")		
		
	}

	n <- length(d$ts1)
	# embed the time series
	ts1.xyz <- time.delay.embed(d$ts1, parms$embed, parms$delay)
	ts2.xyz <- time.delay.embed(d$ts2, parms$embed, parms$delay)	
	D <- rdist(ts1.xyz, ts2.xyz)
	ind <- data.frame(which(D < parms$radius_unscaled, arr.ind=T))
	axis.length <- length(D[1,])
	
	if(addline) {
		width=6
		height=6
	} else {
		width=6
		height=6
	}
	
	if(saveplot == T) {
		quartz(width=width,height=height,type="pdf",file=fname)	
			make.plot(d,ind, axis.length, addline)
		dev.off()
	} else {
		quartz(width=width,height=height)			
			make.plot(d,ind, axis.length, addline)
	}	
}	

recurr.plot <- function(d,parms, fname, saveplot) {
	require(fields)
	n <- length(d$ts1)
	# embed the time series
	ts1.xyz <- time.delay.embed(d$ts1, parms$embed, parms$delay)
	D <- rdist(ts1.xyz, ts1.xyz)
	ind <- data.frame(which(D < parms$radius_unscaled, arr.ind=T))
	axis.length <- length(D[1,])
	if(saveplot == T) {
		quartz(width=8,height=8,type="pdf",file=fname)	
		plot(0, xlim=c(0, axis.length), ylim=c(0, axis.length), type="n", main="Recurrence Plot", xlab="i", ylab="i")
		points(ind$row, ind$col, pch=15, col="#333333", cex=.1)	
		dev.off()
	} else {
		plot(0, xlim=c(0, axis.length), ylim=c(0, axis.length), type="n", main="Recurrence Plot", xlab="i", ylab="i")
		points(ind$row, ind$col, pch=15, col="#333333", cex=.1)	
	}	
}	

# This function returns a diagonal at an offset from the main
superdiag <- function(A,i) {
	# how many rows are there in the matrix
	n<-nrow(A);
	
	if(i >= n) stop("the index cannot be greater than the matrix width")
  if(i < 0) {
  	i <- abs(i);
		len<-n-i;
		c <- 1:len; 
		r <- (i+1):n;   
  } else {
		len<-n-i;
		r <- 1:len; 
		c <- (i+1):n; 
	}
  indices<-(c-1)*n+r; 
  A[indices]
}

ak_crqa <- function (ts1, ts2, delay, embed, rescale, radius, normalize, 
    mindiagline, minvertline, tw = 0, whiteline = F, recpt = F, 
    side = "both", checkl = list(do = F)) 
{
	require(crqa)
    v11 = v21 = NULL
    if (recpt == FALSE) {
        ts1 = as.vector(as.matrix(ts1))
        ts2 = as.vector(as.matrix(ts2))
        if (is.matrix(ts1)) {
            stop("Your data must consist of a single column of data.")
        }
        if (is.matrix(ts2)) {
            stop("Your data must consist of a single column of data.")
        }
        if (checkl$do == TRUE) {
            tsnorm = checkts(ts2, ts1, checkl$datatype, checkl$thrshd, 
                checkl$pad)
            if (tsnorm[[2]] == FALSE) {
                stop("Time-series difference longer than threshold. Increase threshold, or set checkl$do = FALSE avoiding normalization of ts")
            }
            else {
                ts1 = tsnorm[[1]][, 1]
                ts2 = tsnorm[[1]][, 2]
            }
        }
        if (length(ts1) < embed * delay) {
            stop("Phase-space (embed*delay) longer than ts1")
        }
        if (length(ts2) < embed * delay) {
            stop("Phase-space (embed*delay) longer than ts2")
        }
        if (normalize > 0) {
            switch(normalize, {
                1
                ts1 = (ts1 - min(ts1))
                ts1 = ts1/max(ts1)
                ts2 = (ts2 - min(ts2))
                ts2 = ts2/max(ts2)
            }, {
                2
                ts1 = (ts1 - mean(ts1))/sd(ts1)
                ts2 = (ts2 - mean(ts2))/sd(ts2)
            })
        }
        for (loop in 1:embed) {
            vectorstart = (loop - 1) * delay + 1
            vectorend = length(ts1) - ((embed - loop) * delay)
            assign(paste("v1", loop, sep = ""), ts1[vectorstart:vectorend])
        }
        for (loop in 1:embed) {
            vectorstart = (loop - 1) * delay + 1
            vectorend = length(ts2) - ((embed - loop) * delay)
            assign(paste("v2", loop, sep = ""), ts2[vectorstart:vectorend])
        }
        dimts1 = dimts2 = vector()
        for (loop in 1:embed) {
            if (loop == 1) {
                dimts1 = v11
            }
            else {
                eval(parse(text = paste("dimts1 = cbind(dimts1,", 
                  paste("v1", loop, sep = ""), ", deparse.level = 0)", 
                  sep = "")))
            }
        }
        for (loop in 1:embed) {
            if (loop == 1) {
                dimts2 = v21
            }
            else {
                eval(parse(text = paste("dimts2 = cbind(dimts2,", 
                  paste("v2", loop, sep = ""), ", deparse.level = 0)", 
                  sep = "")))
            }
        }
        v1l = length(v11)
        v2l = length(v21)
        dm = rdist(dimts1, dimts2)
        if (rescale > 0) {
            switch(rescale, {
                1
                rescaledist = mean(dm)
                dmrescale = (dm/rescaledist) * 100
            }, {
                2
                rescaledist = max(dm)
                dmrescale = (dm/rescaledist) * 100
            })
        }
        else {
            dmrescale = dm
        }
        ind = which(dmrescale <= radius, arr.ind = TRUE)
        r = ind[, 1]
        c = ind[, 2]
    }
    else {
        ts1 = matrix(as.logical(ts1), ncol = ncol(ts1))
        v1l = nrow(ts1)
        v2l = ncol(ts1)
        ind = which(ts1 > 0, arr.ind = TRUE)
        r = ind[, 1]
        c = ind[, 2]
    }
    if (length(r) != 0 & length(c) != 0) {
        S = sparseMatrix(r, c, dims = c(v1l, v2l))
        S = t(S)
        S = theiler(S, tw)
        if (side == "upper") {
            S = as.matrix(S)
            S[lower.tri(S, diag = TRUE)] = 0
            S = Matrix(S, sparse = TRUE)
        }
        if (side == "lower") {
            S = as.matrix(S)
            S[upper.tri(S, diag = TRUE)] = 0
            S = Matrix(S, sparse = TRUE)
        }
        if (side == "both") {
            S = S
        }
        spdiagonalize = spdiags(S)
        B = spdiagonalize$B
        numrecurs = length(which(B == TRUE))
        percentrecurs = (numrecurs/((v1l * v2l)-v1l)) * 100
        if (is.vector(B)) {
            false = rep(FALSE, length(B))
            B = rbind(false, B, false, deparse.level = 0)
        }
        else {
            false = rep(FALSE, ncol(B))
            B = as.matrix(B)
            B = rbind(false, B, false, deparse.level = 0)
        }
        diaglines = sort(diff(which(B == FALSE)) - 1, decreasing = TRUE)
        diaglines = diaglines[-which(diaglines < mindiagline)]
        if (length(diaglines) != 0) {
            numdiaglines = length(diaglines)
            maxline = max(diaglines)
            meanline = mean(diaglines)
            tabled = as.data.frame(table(diaglines))
            total = sum(tabled$Freq)
            p = tabled$Freq/total
            del = which(p == 0)
            if (length(del) > 0) {
                p = p[-del]
            }
            entropy = -sum(p * log(p))
            relEntropy = entropy/(-1 * log(1/nrow(tabled)))
            pdeter = sum(diaglines)/numrecurs * 100
            restt = tt(S, minvertline, whiteline)
            lam = restt$lam
            TT = restt$TT
        }
        else {
            numdiaglines = 0
            maxline = 0
            pdeter = NA
            entropy = NA
            relEntropy = NA
            meanline = 0
            lam = 0
            TT = 0
            RP = NA
        }
        results = list(RR = percentrecurs, DET = pdeter, NRLINE = numdiaglines, 
            maxL = maxline, L = meanline, ENTR = entropy, rENTR = relEntropy, 
            LAM = lam, TT = TT, RP = S)
    }
    else {
        results = list(RR = 0, DET = NA, NRLINE = 0, maxL = 0, 
            L = 0, ENTR = NA, rENTR = NA, LAM = NA, TT = NA, 
            RP = NA)
    }
    return(results)
}