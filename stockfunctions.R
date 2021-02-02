# Setup the environment
apk.r.setup <- function() {
#	pckgs <- c("lme4", "psych", "quantreg", "xtable", "reshape2", "arm", "metafor", "MCMCpack", "multilevel", "sqldf", "RMySQL", "lavaan", "Hmisc", "sna", "tseriesChaos", "fNonlinear", "timeSeries", "tm", "pscl", "xlsx", "foreign", "nFactors", "pwr", "amen", "TripleR", "ggplot2", "gplots", "lattice", "crqa", "matrixStats", "data.table", "AER", "sas7bdat", "clusterSEs", "RCurl", "sandwich", "RSA", "lavaan.survey", "survey", "mclust", "tidyLPA", "stargazer", "nnet")
	pckgs = c("lme4", "reshame2", "multilevel", "lavaan", "Hmisc", "openxlsx", "lattice", "matrixStats", "data.table", "RSA", "lavaan.survey", "survey", "paws")
	install.packages(pckgs, dependencies=TRUE, repos="http://cran.wustl.edu")
}

# FUNCTION FOR OUTPUTTING A CORRELATION MATRIX
correlation.matrix <- function(x, diagval=NA){
 require(Hmisc)
 require(matrixStats)
 x <- as.matrix(x)
 #j <- as.matrix(j)
 R <- rcorr(x)$r
 p <- rcorr(x)$P
 M <- formatC(round(colMeans(x, na.rm=T),2), digits=2, format="f")

 SD <- formatC(round(colSds(x, na.rm=T),2),digits=2, format="f")

 ## define notations for significance levels; spacing is important.
 mystars <- ifelse(p < .01, "", ifelse(p < .05, "", ifelse(p < .10, "", "")))


 ## trunctuate the matrix that holds the correlations to two decimal
 R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]

 ## build a new matrix that includes the correlations with their apropriate stars
 Rnew <- R
 Rnew<-matrix(paste(Rnew, mystars, sep=""), ncol=ncol(x))
 diag(Rnew) <- diagval
 rownames(Rnew) <- colnames(x)
 colnames(Rnew) <- paste(colnames(x), "", sep="")

 ## remove upper triangle
 Rnew <- as.matrix(Rnew)
 Rnew[upper.tri(Rnew, diag = FALSE)] <- ""

 ## Add the mean and the SD
 Rnew0<-cbind(M, SD, Rnew)

 Rnew1 <- as.data.frame(Rnew0)

 ## remove last column and return the matrix (which is now a data frame)
 #Rnew <- cbind(Rnew[1:length(Rnew)-1])
 return(Rnew1)
}  

plot.interaction <- function(dv, v1, v2, controls,dfile){
eq<-paste(dv,v1,sep="~")
eq<-paste(eq,v2,sep="*")
if(length(controls) > 0)
{
for(i in 1:length(controls))
{
	eq<-paste(eq,controls[i],sep="+")	
}
}
arg<-model.frame(as.formula(eq), dfile)
out<-glm(arg)

v1.x<-mean(dfile[,c(v1)], na.rm=TRUE)
v1.sd<-sd(dfile[,c(v1)], na.rm=TRUE)
v1.lo<-v1.x-v1.sd
v1.hi<-v1.x+v1.sd

v2.x<-mean(dfile[,c(v2)], na.rm=TRUE)
v2.sd<-sd(dfile[,c(v2)], na.rm=TRUE)
v2.lo<-v2.x-v2.sd
v2.hi<-v2.x+v2.sd

lolo<-out$coefficients[1]+out$coefficients[2]*v1.lo+out$coefficients[3]*v2.lo+out$coefficients[length(out$coefficients)]*v1.lo*v2.lo
lohi<-out$coefficients[1]+out$coefficients[2]*v1.lo+out$coefficients[3]*v2.hi+out$coefficients[length(out$coefficients)]*v1.lo*v2.hi
hilo<-out$coefficients[1]+out$coefficients[2]*v1.hi+out$coefficients[3]*v2.lo+out$coefficients[length(out$coefficients)]*v1.hi*v2.lo
hihi<-out$coefficients[1]+out$coefficients[2]*v1.hi+out$coefficients[3]*v2.hi+out$coefficients[length(out$coefficients)]*v1.hi*v2.hi

loline<-c(lolo, lohi)
hiline<-c(hilo, hihi)

#quartz(type="pdf", file=filename)
plot(hiline, type="o",col="red", ylim=c(1,2.5), axes=FALSE, main="", cex.main=.75, xlab="v1", pch=16, ylab="yvar", lty=1, lwd=2)
axis(1, at=1:2, lab=c("Minus 1 SD", "Plus 1 SD"), tck=0, cex.axis=.6)
axis(2, at=c("2.25", "2.75", "3.25"), lab=c("2.25", "2.75", "3.25"), tck=0, cex.axis=.7)
lines(loline, type="o", col="blue", pch=18, lwd=2, lty=2)
legend("bottomright", legend=c("+1 SD v2", "-1 SD v2"), cex=.8,col=c("red", "blue"), pch=c(16,18), lty=c(1,2), bty="n")
#dev.off()
return(out)
}

blau<-function(catvar, groupid, dfile){
require(sqldf)
q0<-paste("SELECT ",groupid,", COUNT(",groupid,") as groupsize FROM dfile WHERE ",catvar," IS NOT NULL GROUP BY ",groupid, sep=" ")
out0<-data.frame(sqldf(q0))
new<-merge(dfile, out0, by=c(groupid))
q1<-paste("SELECT ",groupid,",groupsize, COUNT(", catvar,") AS catsize FROM new GROUP BY ",groupid,", ",catvar, sep=" ")
out1<-data.frame(sqldf(q1))
out1$prop2<-(out1$catsize/out1$groupsize)*(out1$catsize/out1$groupsize)
q2<-paste("SELECT ",groupid,", (1-SUM(prop2)) AS ",catvar,"_blau FROM out1 GROUP BY ",groupid,sep="")
d<-data.frame(sqldf(q2))
return<-d[,c(groupid,paste(catvar,"_blau",sep=""))]
}

regtable<-function(models, rnames)
{
	numrows<-length(rnames)	
	rnames<-data.frame(rnames)
	colnames(rnames)<-c("Var")
	df<-c("DF",rep("",length(models)*3))
    rsq<-c("R2",rep("",length(models)*3))
    fstat<-c("F",rep("",length(models)*3))	
	for(i in 1:length(models))
	{	
		# Get the coefficient information - doing as much info as would be needed in a table. 
		s<-mapply(summary, models[i], USE.NAMES = TRUE)
		coef<-s[4,1]$coefficients
		m.b <- round(coef[,1],2)
		m.se <- round(coef[,2],2)		
		m.p <- round(coef[,4],3)
		coefs<-data.frame(cbind(row.names(coef), m.b, m.se, m.p))
		colnames(coefs)<-c("Var",paste("b",i,sep=""), paste("se",i,sep=""), paste("p", i, sep=""))
		
		# Get model summary information. 
		index <- (i-1)*3 + 2
		df[index]<-paste(s[10,1]$fstatistic[2], s[10,1]$fstatistic[3], sep=",")
		rsq[index]<-formatC(round(s[8,1]$r.squared[1],2), digits=2, format="f")
		fstat[index]<-formatC(round(s[10,1]$fstatistic[1],2), digits=2, format="f")
		fstat[index+2]<-round(pf(s[10,1]$fstatistic[1], s[10,1]$fstatistic[2], s[10,1]$fstatistic[3], lower.tail=FALSE), 3)		
		if(i == 1)
		{
			tab<-merge(rnames, coefs, by=c("Var"), all.x=TRUE, sort=FALSE)
		} else {
			tab<-merge(tab, coefs, by=c("Var"), all.x=TRUE, sort=FALSE)
		}
	}
	
	# Add the summary statistics for each model
	
	data.frame(tab)
	# Now make sure that the order is as desired
	tab<-merge(rnames, tab, by=c("Var"), all.x=TRUE, sort=FALSE)
	
	sumstats<-data.frame(rbind(df, rsq, fstat))
	colnames(sumstats)<-colnames(tab)
	tab<-rbind(tab,sumstats)
	return(tab)
}

hlmtable<-function(models, rnames)
{
	numrows<-length(rnames)	
	rnames<-data.frame(rnames)
	colnames(rnames)<-c("Var")
	df<-c("DF",rep("",length(models)))
    rsq<-c("R2",rep("",length(models)))
    fstat<-c("F",rep("",length(models)))	
	for(i in 1:length(models))
	{	
		s<-mapply(summary, models[i], USE.NAMES = TRUE)
		coef<-s[4,1]$coefficients
		column<-rep("",length(coef[,1]))
		name<-rep("",length(coef[,1]))
		for(z in 1:length(coef[,1]))
		{
			name[z]<-rownames(coef)[z]		
			star <- ifelse(coef[z,4] < .01, "**", ifelse(coef[z,4] < .05, "*", ifelse(coef[z,4] < 0.10, "+","")))
   			val<-formatC(round(coef[z,1],2),digits=2, format="f")
			column[z]<-paste(val,star,sep="")
		}
		coefs<-data.frame(cbind(name, column))
		colnames(coefs)<-c("Var",paste("M",i,sep=""))
		
		df[(i+1)]<-paste(s[10,1]$fstatistic[2], s[10,1]$fstatistic[3], sep=",")
		rsq[(i+1)]<-formatC(round(s[8,1]$r.squared[1],2), digits=2, format="f")
		fval<-formatC(round(s[10,1]$fstatistic[1],2), digits=2, format="f")
		pval<-pf(s[10,1]$fstatistic[1], s[10,1]$fstatistic[2], s[10,1]$fstatistic[3], lower.tail=FALSE)
		star <- ifelse(pval < .001, "***", ifelse(pval < .01, "** ", ifelse(pval < .05, "* ", " ")))
		fstat[(i+1)]<-paste(fval, star, sep="")		
		if(i == 1)
		{
			tab<-merge(coefs, rnames, by=c("Var"), all.y=TRUE, sort=FALSE)
		} else {
			tab<-merge(tab, coefs, by=c("Var"), all.y=TRUE, sort=FALSE)
		}
	}
	
	# Add the summary statistics for each model
	
	data.frame(tab)
	# Now make sure that the order is as desired
	tab<-merge(rnames, tab, by=c("Var"), all.x=TRUE, sort=FALSE)
	
	sumstats<-data.frame(rbind(df, rsq, fstat))
#	sumstats<-sumstats[,c(1,rev(2:length(colnames(tab))))]
	colnames(sumstats)<-colnames(tab)
	tab<-rbind(tab,sumstats)
	return(tab)
}

specify_decimal <- function(x, k) format(round(x, k), nsmall=k)

lmer.table <- function(models, rnames) {
	require(arm)
	rnames <- c(rnames, "Random Intercept", "Random Residual", "AIC", "Deviance")
	rnames<-data.frame(rnames, stringsAsFactors=F)
	colnames(rnames)<-c("Var")
   
	for(i in 1:length(models)) {
		vals <- display(models[[i]])
		rbind(summary(models[[i]])$varcor,NA)

		int <- round(attr(summary(models[[i]])$varcor[[1]], "stddev"),4)^2
		res <- round(attr(summary(models[[i]])$varcor, "sc"),4)^2
		aic <- round(vals$AIC,4)
		dev <- round(vals$deviance,4)

		Var <- names(vals$coef)
		b <- round(vals$coef,4)
		se <- round(vals$se, 4)
		t <- round(b/se, 4)
		new.col <- data.frame(Var,b,se, t,stringsAsFactors=F)
		new.col <- rbind(new.col, c("Random Intercept", int, NA, NA))
		new.col <- rbind(new.col, c("Random Residual", res, NA, NA))
		new.col <- rbind(new.col, c("AIC", aic, NA, NA))				
		new.col <- rbind(new.col, c("Deviance", dev, NA, NA))						
		colnames(new.col)[2:4] <- c(paste("b",i,sep=""), paste("se",i,sep=""), paste("t",i,sep=""))
		if(i == 1 ) {
			tab <- merge(rnames, new.col, by=c("Var"), all.x=TRUE, sort=FALSE)
		} else {
			tab <- merge(tab, new.col, by=c("Var"), all.x=TRUE, sort=FALSE)
		}
	}
	
	# Order everything
	tab<-merge(rnames, tab, by=c("Var"), all.x=TRUE, sort=FALSE)	
	return(tab)	
}

# Output a nice table of results for bayesian models

bayestable <- function(vars) {
	p<-matrix(rep(NA,(6*ncol(vars))), nrow=ncol(vars), ncol=6)
	for(i in 1:ncol(vars)){
		p[i,1] <- mean(vars[,i], na.rm=TRUE)
		p[i,2] <- sd(vars[,i], na.rm=TRUE)
		p[i,4] <- quantile(vars[,i], .025)
		p[i,5] <- quantile(vars[,i], .975)
		if(p[i,1] > 0)
		{
			num <- length(which(vars[,i] < 0))
		}
		else
		{
			num <- length(which(vars[,i] > 0))
		}
		p[i,6] <- num/length(vars[,i])
	}	
	p[,3] <- p[,1] / p[,2]
	rownames(p) <- colnames(vars)
	colnames(p) <- c("X","SD","t",".025",".975","p")
	return(p)
}

# A function for calculating simple slopes given an object output from lme, z is a vector of three points at which you want slopes


# A function for calculating simple slopes given coefficients, standard errors, and the positions in the array of the relevant variables

ss.lin.ols <- function(coefs, ses, x.n, z.n, xz.n,z, zvals) {
	b.0 <- coefs[1]
	b.1 <- coefs[x.n]
	b.2 <- coefs[z.n]
	b.3 <- coefs[xz.n]
	
#	z.hi <- mean(z, na.rm=TRUE) + sd(z, na.rm=TRUE)
#	z.av <- mean(z, na.rm=TRUE)
#	z.lo <- mean(z, na.rm=TRUE) - sd(z, na.rm=TRUE)	

	z.hi <- zvals[3]
	z.av <- zvals[2]
	z.lo <- zvals[1]		
	
	int.hi <- (b.0 + b.2*z.hi)
	int.av <- (b.0 + b.2*z.av)
	int.lo <- (b.0 + b.2*z.lo)		
	
	ss.hi <- (b.1+b.3*z.hi)
	ss.av <- (b.1+b.3*z.av)
	ss.lo <- (b.1+b.3*z.lo)
	
	se.hi <- sqrt(ses[x.n,x.n] + 2*z.hi*ses[x.n,xz.n] + (z.hi^2)*ses[xz.n,xz.n])
	se.av <- sqrt(ses[x.n,x.n] + 2*z.av*ses[x.n,xz.n] + (z.av^2)*ses[xz.n,xz.n])
	se.lo <- sqrt(ses[x.n,x.n] + 2*z.lo*ses[x.n,xz.n] + (z.lo^2)*ses[xz.n,xz.n])
	
	se.int.hi <- sqrt(ses[1,1] + 2*z.hi*ses[1,z.n] + (z.hi^2)*ses[z.n, z.n])
	se.int.av <- sqrt(ses[1,1] + 2*z.av*ses[1,z.n] + (z.av^2)*ses[z.n, z.n])	
	se.int.lo <- sqrt(ses[1,1] + 2*z.lo*ses[1,z.n] + (z.lo^2)*ses[z.n, z.n])			
	
	tval.hi <- ss.hi/se.hi
	tval.av <- ss.av/se.av 
	tval.lo <- ss.lo/se.lo	
	
	tval.int.hi <- int.hi/se.int.hi
	tval.int.av <- int.av/se.int.av	
	tval.int.lo <- int.lo/se.int.lo			
	
	dfs <- length(z) - (length(coefs)-1) - 1
	
	pval.hi <- 2*pt(-abs(tval.hi),df=dfs)
	pval.av <- 2*pt(-abs(tval.av),df=dfs)
	pval.lo <- 2*pt(-abs(tval.lo),df=dfs)	
	
	pval.int.hi <- 2*pt(-abs(tval.int.hi),df=dfs)
	pval.int.av <- 2*pt(-abs(tval.int.av),df=dfs)	
	pval.int.lo <- 2*pt(-abs(tval.int.lo),df=dfs)	
		
	ss <- rbind(ss.hi, ss.av, ss.lo)
	se <- rbind(se.hi, se.av, se.lo)
	tval <- rbind(tval.hi, tval.av, tval.lo)
	dfs <- rbind(dfs, dfs, dfs)
	pval <- rbind(pval.hi, pval.av, pval.lo)
	
	int <- rbind(int.hi, int.av, int.lo)
	se.int <- rbind(se.int.hi, se.int.av, se.int.lo)
	tval.int <- rbind(tval.int.hi, tval.int.av, tval.int.lo)
	pval.int <- rbind(pval.int.hi, pval.int.av, pval.int.lo)		
	
	ret <- cbind(ss, se, tval, dfs, pval, int, se.int, tval.int, pval.int)
	row.names(ret) <- c("+1SD Z","Mean Z","-1SD Z")
	colnames(ret) <- c("SS","SE","t","df","p", "INT", "INT.SE", "INT.t", "INT.p")
	return(ret)
}

ss.lin.rcm <- function(obj, x.n, z.n, xz.n, z) {
	coefs <- fixef(obj)
	ses <- vcov(obj)
	dfs <- length(obj@frame[,1])-length(coefs)
	
	b.1 <- coefs[x.n]
	b.2 <- coefs[z.n]
	b.3 <- coefs[xz.n]
	
	z.hi <- mean(z, na.rm=TRUE) + sd(z, na.rm=TRUE)
	z.av <- mean(z, na.rm=TRUE)
	z.lo <- mean(z, na.rm=TRUE) - sd(z, na.rm=TRUE)		
	
	ss.hi <- (b.1+b.3*z.hi)
	ss.av <- (b.1+b.3*z.av)
	ss.lo <- (b.1+b.3*z.lo)
	
	se.hi <- sqrt(ses[x.n,x.n] + 2*z.hi*ses[x.n,xz.n] + (z.hi^2)*ses[xz.n,xz.n])
	se.av <- sqrt(ses[x.n,x.n] + 2*z.av*ses[x.n,xz.n] + (z.av^2)*ses[xz.n,xz.n])
	se.lo <- sqrt(ses[x.n,x.n] + 2*z.lo*ses[x.n,xz.n] + (z.lo^2)*ses[xz.n,xz.n])
	
	tval.hi <- ss.hi/se.hi
	tval.av <- ss.av/se.av 
	tval.lo <- ss.lo/se.lo	
	
	dfs <- dfs
	
	pval.hi <- 2*pt(-abs(tval.hi),df=dfs)
	pval.av <- 2*pt(-abs(tval.av),df=dfs)
	pval.lo <- 2*pt(-abs(tval.lo),df=dfs)	
	
	ss <- rbind(ss.hi, ss.av, ss.lo)
	se <- rbind(se.hi, se.av, se.lo)
	tval <- rbind(tval.hi, tval.av, tval.lo)
	dfs <- rbind(dfs, dfs, dfs)
	pval <- rbind(pval.hi, pval.av, pval.lo)
	
	ret <- cbind(ss, se, tval, dfs, pval)
	row.names(ret) <- c("+1SD Z","Mean Z","-1SD Z")
	colnames(ret) <- c("SS","SE","t","df","p")
	return(ret)
}


ss.lin.rcm.bin <- function(obj, x.n, z.n, xz.n) {
	coefs <- fixef(obj)
	ses <- vcov(obj)
	dfs <- length(obj@frame[,1])-length(coefs)
	
	b.0 <- coefs[1]	
	b.1 <- coefs[x.n]
	b.2 <- coefs[z.n]
	b.3 <- coefs[xz.n]
	
	z.hi <- 1
	z.lo <- 0
	
	int.hi <- (b.0 + b.2*z.hi)
	int.lo <- (b.0 + b.2*z.lo)		
	
	
	ss.hi <- (b.1+b.3*z.hi)
	ss.lo <- (b.1+b.3*z.lo)
	
	se.hi <- sqrt(ses[x.n,x.n] + 2*z.hi*ses[x.n,xz.n] + (z.hi^2)*ses[xz.n,xz.n])
	se.lo <- sqrt(ses[x.n,x.n] + 2*z.lo*ses[x.n,xz.n] + (z.lo^2)*ses[xz.n,xz.n])
		
	se.int.hi <- sqrt(ses[1,1] + 2*z.hi*ses[1,z.n] + (z.hi^2)*ses[z.n, z.n])
	se.int.lo <- sqrt(ses[1,1] + 2*z.lo*ses[1,z.n] + (z.lo^2)*ses[z.n, z.n])			
	
	tval.hi <- ss.hi/se.hi
	tval.lo <- ss.lo/se.lo	
	
	tval.int.hi <- int.hi/se.int.hi
	tval.int.lo <- int.lo/se.int.lo			
		
	dfs <- dfs
	
	pval.hi <- 2*pt(-abs(tval.hi),df=dfs)
	pval.lo <- 2*pt(-abs(tval.lo),df=dfs)	
	
	pval.int.hi <- 2*pt(-abs(tval.int.hi),df=dfs)
	pval.int.lo <- 2*pt(-abs(tval.int.lo),df=dfs)	
	
	
	ss <- rbind(ss.hi, ss.lo)
	se <- rbind(se.hi, se.lo)
	tval <- rbind(tval.hi, tval.lo)
	dfs <- rbind(dfs, dfs)
	pval <- rbind(pval.hi,pval.lo)
	
	int <- rbind(int.hi, int.lo)
	se.int <- rbind(se.int.hi, se.int.lo)
	tval.int <- rbind(tval.int.hi,  tval.int.lo)
	pval.int <- rbind(pval.int.hi,  pval.int.lo)		
	
	
	ret <- cbind(ss, se, tval, dfs, pval, int, se.int, tval.int, pval.int)
	row.names(ret) <- c("Z=1","Z=0")
	colnames(ret) <- c("SS","SE","t","df","p", "INT", "INT.SE", "INT.t", "INT.p")
	return(ret)
}



# Function to return typical scale statistics for group-level constructs
scale.scores <- function(vars, dat) {
	count <- 1	
	for(i in vars) {		
		if(length(i[[3]]) > 1) {
		
			val = data.frame(rowMeans(dat[,c(i[[3]])], na.rm=T))
			colnames(val) = i[[2]]
			if(count == 1) {
				return.dat <- val	
			} else {
				return.dat <- cbind(return.dat, val)
			}
			count = count + 1					
		}	
	}
	dat = cbind(dat, return.dat)
	return(dat)	
}

scale.stats <- function (vars, groupvar, dat) {
	require(psych)
	require(multilevel)
	
	# initialize a matrix to hold the results of this


	
	icc1 <- rep(NA, length(vars))
	icc1p <- rep(NA, length(vars))
	icc2 <- rep(NA, length(vars))
	alpha <- rep(NA, length(vars))
	alpha_group <- rep(NA, length(vars))	
	rwgjmed <- rep(NA, length(vars))
	rwgjx <- rep(NA, length(vars))
	rwgjsd <- rep(NA, length(vars))
	rwgjmin <- rep(NA, length(vars))
	rwgjmax <- rep(NA, length(vars))
#	indiv_n <- rep(NA, length(vars))	
#	group_n <- rep(NA, length(vars))		
#	group_size <- rep(NA, length(vars))			
	name <- rep(NA, length(vars))
	count <- 1
	for(i in vars) {
			name[count] <- i[[1]]
			dv <- dat[, c(i[[2]])]
			gv <- dat[, c(groupvar)]    
			items <- dat[,c(i[[3]])]
			out <- aov(dv ~ as.factor(gv), data=dat)
			icc1[count] <-round(ICC1(out),3)
			icc1p[count] <- round(summary(out)[[1]][["Pr(>F)"]], 4)[1]
			icc2[count] <-round(ICC2(out),3)
			
			## Get the individual-level n
#			indiv_n[count] <- sum(!is.na(dv))
			
			## Get the group-level n			
#			group_n[count] <- length(unique(gv))
						
			if(length(i[[3]]) <= 1) {
				alpha[count] <- NA
				alpha_group[count] <- NA				
				out<-rwg(items,gv, i[[4]])
				rwgjmed[count]<-round(median(out$rwg, na.rm=TRUE),3)			
				rwgjx[count]<-round(mean(out$rwg, na.rm=TRUE),3)
				rwgjsd[count]<-round(sd(out$rwg, na.rm=TRUE),3)
				rwgjmin[count]<-round(min(out$rwg, na.rm=TRUE),3)
				rwgjmax[count]<-round(max(out$rwg, na.rm=TRUE),3)        								
			} else {
				# aggregate the items for the group-level alpha
				agg <- aggregate(dat[,c(i[[3]])], by=list(dat[,groupvar]), mean, na.rm=TRUE)
				colnames(agg) <- c(groupvar,c(i[[3]]))
				out <- alpha(agg[,c(i[[3]])], check.keys=FALSE)				
				alpha_group[count] <- round(out$total['raw_alpha'],3)				
				# Produce the individual-level alpha
				out <- alpha(items, na.rm=TRUE, check.keys=FALSE)
				alpha[count] <- round(out$total['raw_alpha'],3)
				out<-rwg.j(items,gv, i[[4]])
				rwgjmed[count]<-round(median(out$rwg.j, na.rm=TRUE),3)			
				rwgjx[count]<-round(mean(out$rwg.j, na.rm=TRUE),3)
				rwgjsd[count]<-round(sd(out$rwg.j, na.rm=TRUE),3)
				rwgjmin[count]<-round(min(out$rwg.j, na.rm=TRUE),3)
				rwgjmax[count]<-round(max(out$rwg.j, na.rm=TRUE),3)        				
			}			
			count <- count + 1
	}
	v <- data.frame(name,  alpha=unlist(alpha), alpha_group=unlist(alpha_group), icc1, icc1p, icc2, rwgjx, rwgjsd, rwgjmed, rwgjmin, rwgjmax, stringsAsFactors=F)
	# Clean up the formatting of this return table
	return(v)
}

# This function takes the construct listing and generates two tables - one containing fit indices, the second containing standardized loadings
scale.cfa <- function(vars,dat) {
	require(lavaan)
	count <- 1
	for(construct in vars) {
		itemlist <- paste0(construct$items, collapse=" + ")
		model <- paste('f1 =~ ', itemlist, sep="")
		mod <- cfa(model, data=dat)
		mod.fit <- fitMeasures(mod)
		mod.est <- standardizedSolution(mod)
		fit.ind <- c(round(mod.fit[c(3,4,5,9,10,19,20,23,29)],3))
		ind.names <- mod.est[1:length(construct$items),3]
		std.loads <- round(mod.est[1:length(construct$items),4],3)
		var.name <- c(construct$name)
		lds <- cbind(rep(var.name), length(ind.names), ind.names, std.loads)
		if(count == 1) {
			fit.out <- c(var.name, fit.ind)
			load.out <- lds
		} else {
			fit.out <- rbind(fit.out, c(var.name, fit.ind))
			load.out <- rbind(load.out, lds)
		}
		count <- count + 1
	}
	row.names(fit.out) <- 1:length(fit.out[,1])
	fit.out <- data.frame(fit.out)
	load.out <- data.frame(load.out)
	return(list(fit.out, load.out))
}


# Aggregation Function

ak_aggregate <- function(dat, constructs, groupvars, stdev=FALSE, items=FALSE) {
	require(multilevel)
	require(sqldf)
	aggvars <- rep(NA, length(constructs))
	for(i in 1:length(constructs)) {
		aggvars[i] <- constructs[[i]]$var	
	}	
	
	if(items) {
		for(i in 1:length(constructs)) {
			aggvars <- c(aggvars, constructs[[i]]$items)
		}
	}
		
	group_list <- vector("list", length(groupvars))
	for(i in 1:length(groupvars)) {
		group_list[[i]] <- dat[,groupvars[i]]
	}
	
	q <- paste("SELECT ",paste(groupvars, collapse=", "), ", COUNT(", groupvars[length(groupvars)],") AS groupsize FROM dat GROUP BY ", paste(groupvars,collapse=","))

	n <- sqldf(q)
	
	agg.x <- aggregate(dat[,aggvars], group_list, mean, na.rm=TRUE)
	colnames(agg.x) <- c(groupvars, paste(aggvars, "_x", sep=""))
	
	if(stdev) {
		agg.sd <- aggregate(dat[,aggvars], group_list, sd, na.rm=TRUE)
		colnames(agg.sd) <- c(groupvars, paste(aggvars, "_sd", sep=""))	
		agg.dat <- merge(agg.x, agg.sd, by=groupvars)			
	} else {
		agg.dat <- agg.x
	}
	agg.dat <- merge(agg.dat, n, by=groupvars)
	return(agg.dat)
}

# Function to ge ta correlation lookup table

cor.table <- function(d) {
	varlist <- rep(NA, length(colnames(d)))
	count <- 1
	for(i in 1:length(colnames(d))) {
		if(is.numeric(d[,i]) || is.integer(d[,i])) {
			varlist[count] <- colnames(d)[i]
			count <- count + 1
		}
	}
	cvars <- varlist[!is.na(varlist)]	
	
	size <- length(cvars)*length(cvars)
	xvar <- rep(NA, size)
	yvar <- rep(NA, size)
	r <- rep(NA, size)
	p <- rep(NA, size)
	nobs <- rep(NA, size)	
	count <- 1
	for(x in 1:length(cvars)) {
		xname <- cvars[x]
		for(y in 1:length(cvars)) {
			yname <- cvars[y]	
			xvar[count] <- xname
			yvar[count] <- yname						
			nobs[count] = nrow(d[!is.na(d[,xname]) & !is.na(d[,yname]), ]) 
			## See if there are enough observations ##
			if(nobs[count] >= 2) {						
				cout <- cor.test(d[,xname], d[,yname], na.rm=TRUE)
				r[count] <- cout$estimate
				p[count] <- cout$p.value					
			} else {
				r[count] <- NA
				p[count] <- NA			
			}
			count <- count + 1
		}
	}
	ctab <- data.frame(cbind(xvar, yvar, r, p, nobs), stringsAsFactors=FALSE)
	ctab[,3] <- round(as.numeric(ctab[,3]),3)
	ctab[,4] <- round(as.numeric(ctab[,4]),3)
	ctab[,5] <- as.numeric(ctab[,5])
	return(ctab[!is.na(ctab$p) & ctab$xvar != ctab$yvar, ])
}

# function for r-squared for lmer

# R-squared information
r2.corr.lmer <- function(lmer.object) {
	summary(lm(attr(lmer.object, "y") ~ fitted (lmer.object) ))$r.squared
}

# R-squared information
Rsq <- function(reml.mod) { 
  ## Based on 
   ## N. J. D. Nagelkerke. A note on a general definition 
   ## of the coefficient of determination. Biometrika, 78:691–692, 1991. 
   ml.mod <- update(reml.mod, REML=FALSE) 
   l.B <- logLik(ml.mod) 
   l.0 <- logLik( lm(ml.mod@y ~ 1) ) 
   Rsq <- 1 - exp( - ( 2/length(ml.mod@y) ) * (l.B - l.0) ) 
Rsq[1] 
} 

##################################################################
# This function conducts the Dawson & Richter probing of a three-way interaction #
# Input: An output object either generated via lm or lmer, indication of the positioning in the variable order of the different components #
# Output: A table of slope differences and test statistics 
##################################################################

three.way <- function(obj, x,w,z,xw,xz,wz,xwz,w.vals, z.vals) {
	# Figure out if this is a multilevel model or a single level model
	# to extract the parameters estimates and the standard errors
	if( class(obj) == "lm" ) {
		b <- coef(obj)
		se <- vcov(obj)
		df <- obj$df.residual
	} else {
		b <- fixef(obj)
		se <- vcov(obj)
		df <- length(obj@frame[,1])-length(b)		
	}
	# Compute the simple slopes
	b1 = b[x] + b[xz]*z.vals[2] + b[xw]*w.vals[2] + b[xwz]*z.vals[2]*w.vals[2]
	b2 = b[x] + b[xz]*z.vals[2] + b[xw]*w.vals[1] + b[xwz]*z.vals[2]*w.vals[1]
	b3 = b[x] + b[xz]*z.vals[1] + b[xw]*w.vals[2] + b[xwz]*z.vals[1]*w.vals[2]
	b4 = b[x] + b[xz]*z.vals[1] + b[xw]*w.vals[1] + b[xwz]*z.vals[1]*w.vals[1]	
	
	# Compute the standard errors for the simple slopes
	s1 = sqrt( se[x,x] + se[xz,xz]*z.vals[2]^2 + se[xw,xw]*w.vals[2]^2 + se[xwz,xwz]*(w.vals[2]^2)*(z.vals[2]^2) + se[x,xz]*z.vals[2]*2 + se[x,xw]*w.vals[2]*2 + se[x,xwz]*z.vals[2]*w.vals[2]*2 + se[xz,xw]*z.vals[2]*w.vals[2]*2 + se[xz,xwz]*2*w.vals[2]*z.vals[2]^2 + se[xw,xwz]*z.vals[2]*2*w.vals[2]^2 )
	
	s2 = sqrt( se[x,x] + se[xz,xz]*z.vals[2]^2 + se[xw,xw]*w.vals[1]^2 + se[xwz,xwz]*(w.vals[1]^2)*(z.vals[2]^2) + se[x,xz]*z.vals[2]*2 + se[x,xw]*w.vals[1]*2 + se[x,xwz]*z.vals[2]*w.vals[1]*2 + se[xz,xw]*z.vals[2]*w.vals[1]*2 + se[xz,xwz]*2*w.vals[1]*z.vals[2]^2 + se[xw,xwz]*z.vals[2]*2*w.vals[1]^2 )	

	s3 = sqrt( se[x,x] + se[xz,xz]*z.vals[1]^2 + se[xw,xw]*w.vals[2]^2 + se[xwz,xwz]*(w.vals[2]^2)*(z.vals[1]^2) + se[x,xz]*z.vals[1]*2 + se[x,xw]*w.vals[2]*2 + se[x,xwz]*z.vals[1]*w.vals[2]*2 + se[xz,xw]*z.vals[1]*w.vals[2]*2 + se[xz,xwz]*2*w.vals[2]*z.vals[1]^2 + se[xw,xwz]*z.vals[1]*2*w.vals[2]^2 )

	s4 = sqrt( se[x,x] + se[xz,xz]*z.vals[1]^2 + se[xw,xw]*w.vals[1]^2 + se[xwz,xwz]*(w.vals[1]^2)*(z.vals[1]^2) + se[x,xz]*z.vals[1]*2 + se[x,xw]*w.vals[1]*2 + se[x,xwz]*z.vals[1]*w.vals[1]*2 + se[xz,xw]*z.vals[1]*w.vals[1]*2 + se[xz,xwz]*2*w.vals[1]*z.vals[1]^2 + se[xw,xwz]*z.vals[1]*2*w.vals[1]^2 )
	
	# Compute the t values for these simple slopes
	t1 = b1/s1
	t2 = b2/s2
	t3 = b3/s3
	t4 = b4/s4
	
	# Get the p values for these simple slopes
	p1 = 2*pt(-abs(t1),df=df)
	p2 = 2*pt(-abs(t2),df=df)
	p3 = 2*pt(-abs(t3),df=df)
	p4 = 2*pt(-abs(t4),df=df)				
	
	# Compute the simple slope differences
	b12 = b[xw]*(w.vals[2]-w.vals[1]) + b[xwz]*z.vals[2]*(w.vals[2]-w.vals[1])
	
	b13 = b[xz]*(z.vals[2]-z.vals[1]) + b[xwz]*w.vals[2]*(z.vals[2]-z.vals[1])
	
	b14 = b[xz]*(z.vals[2]-z.vals[1]) + b[xw]*(w.vals[2]-w.vals[1]) + b[xwz]*(z.vals[2]*w.vals[2]-z.vals[1]*w.vals[1])
	
	b23 = b[xz]*(z.vals[2]-z.vals[1]) + b[xw]*(w.vals[1]-w.vals[2]) + b[xwz]*(z.vals[2]*w.vals[1]-z.vals[1]*w.vals[2])
	
	b24 = b[xz]*(z.vals[2]-z.vals[1]) + b[xwz]*w.vals[1]*(z.vals[2]-z.vals[1])
	
	b34 = b[xw]*(w.vals[2]-w.vals[1]) + b[xwz]*z.vals[1]*(w.vals[2]-w.vals[1])
		
	# Get the standard errors for the differences
	s12 = (w.vals[2]-w.vals[1])*sqrt( se[xw,xw] + se[xwz,xwz]*z.vals[2]^2 + 2*z.vals[2]*se[xw,xwz] )
	
	s13 = (z.vals[2]-z.vals[1])*sqrt( se[xz,xz] + se[xwz,xwz]*w.vals[2]^2 + 2*w.vals[2]*se[xz, xwz] )
	
	s14 = sqrt( se[xz,xz]*(z.vals[2]-z.vals[1])^2 + se[xw,xw]*(w.vals[2]-w.vals[1])^2 + se[xwz,xwz]*(z.vals[2]*w.vals[2]-z.vals[1]*w.vals[1])^2 + 2*( se[xz,xw]*(z.vals[2]-z.vals[1])*(w.vals[2]-w.vals[1]) + se[xz,xwz]*(z.vals[2]-z.vals[1])*(z.vals[2]*w.vals[2]-z.vals[1]*w.vals[1]) + se[xw,xwz]*(w.vals[2]-w.vals[1])*(z.vals[2]*w.vals[2]-z.vals[1]*w.vals[1]) ) )
	
	s23 = sqrt( se[xz,xz]*(z.vals[2]-z.vals[1])^2 + se[xw,xw]*(w.vals[1]-w.vals[2])^2 + se[xwz,xwz]*(z.vals[2]*w.vals[1]-z.vals[1]*w.vals[2])^2 + 2*( se[xz,xw]*(z.vals[2]-z.vals[1])*(w.vals[1]-w.vals[2]) + se[xz,xwz]*(z.vals[2]-z.vals[1])*(z.vals[2]*w.vals[1]-z.vals[1]*w.vals[2]) + se[xw,xwz]*(w.vals[1]-w.vals[2])*(z.vals[2]*w.vals[1]-z.vals[1]*w.vals[2]) ) ) 
	
	s24 = (z.vals[2]-z.vals[1])*sqrt( se[xz,xz] + se[xwz,xwz]*w.vals[1]^2 + 2*w.vals[1]*se[xz, xwz] )
	
	s34 = (w.vals[2]-w.vals[1])*sqrt( se[xw,xw] + se[xwz,xwz]*z.vals[1]^2 + 2*z.vals[1]*se[xw,xwz] )
	
	# Get the test statistics
	t12 = b12/s12
	t13 = b13/s13
	t14 = b14/s14
	t23 = b23/s23
	t24 = b24/s24
	t34 = b34/s34	
	
	# Get the p values
	p12 = 2*pt(-abs(t12),df=df)	
	p13 = 2*pt(-abs(t13),df=df)	
	p14 = 2*pt(-abs(t14),df=df)	
	p23 = 2*pt(-abs(t23),df=df)	
	p24 = 2*pt(-abs(t24),df=df)	
	p34 = 2*pt(-abs(t34),df=df)						
	
	# Create the output tables
	#t1 = info for the four simple slopes
	zlab = c("High", "High", "Low", "Low")
	wlab = c("High", "Low", "High", "Low")	
	slab = c(1:4)
	bval = c(b1,b2,b3,b4)
	seval = c(s1,s2,s3,s4)
	tval = c(t1,t2,t3,t4)
	dfval= c(df,df,df,df)
	pval = c(p1,p2,p3,p4)	
	
	table1 = data.frame(zlab,wlab, slab,bval,seval, tval,dfval,pval, stringsAsFactors=F)
	colnames(table1)<- c("Z", "W", "slope_label", "b", "se", "t", "df", "p")
	
	#t2 = info for slope difference tests
	slope1 = c(1,1,1,2,2,3)
	slope2 = c(2,3,4,3,4,4)	
	slab = c("a", "b", "c", "d", "e", "f")
	bval = c(b12,b13,b14,b23,b24,b34)
	seval = c(s12,s13,s14,s23,s24,s34)
	tval = c(t12,t13,t14,t23,t24,t34)
	dfval= c(df,df,df,df,df,df)
	pval = c(p12,p13,p14,p23,p24,p34)	
	table2 = data.frame(slope1,slope2, slab,bval,seval, tval,dfval,pval, stringsAsFactors=F)
	colnames(table2)<- c("slope_1", "slope_2", "slope_label", "b", "se", "t", "df", "p")
		
	# Return the values
	return(list(table1,table2))
}

indirect.effects <- function(m1, m2, coef1, coef2, nsims, d) {	
	# Add an id variable to the file
	d$boot_id <- 1:length(d[,1])
	i <- 1
	ie <- rep(NA,nsims)
	for(i in 1:nsims) {
		# pull the bootstrap sample
		unt <- sample(1:length(d[,1]), length(d[,1]), replace=TRUE)
		unt.num <- d$boot_id
		pulls <- data.frame(cbind(unt.num, unt))		
		bset <- merge(pulls, d, by.x=c("unt"), by.y=c("boot_id"), sort=TRUE)	
		
		# Run the regressions
		o1 <- lm(m1, data=bset)
		o2 <- lm(m2, data=bset)
		
		# Pull the coefficients and calculate the indirect effect
		ie[i] <- o1$coefficients[coef1]*o2$coefficients[coef2]				
	}
	return(ie)	
}

indirect.effects.lmer <- function(m1, m2, coef1, coef2, nsims, d) {	
	# Add an id variable to the file
	d$boot_id <- 1:length(d[,1])
	i <- 1
	ie <- rep(NA,nsims)
	for(i in 1:nsims) {
		# pull the bootstrap sample
		unt <- sample(1:length(d[,1]), length(d[,1]), replace=TRUE)
		unt.num <- d$boot_id
		pulls <- data.frame(cbind(unt.num, unt))		
		bset <- merge(pulls, d, by.x=c("unt"), by.y=c("boot_id"), sort=TRUE)	
		
		# Run the regressions
		o1 <- lmer(m1, data=bset)
		o2 <- lmer(m2, data=bset)
		
		# Pull the coefficients and calculate the indirect effect
		ie[i] <- o1$coefficients[coef1]*o2$coefficients[coef2]				
	}
	return(ie)	
}

boot.table <- function(eff) {
	est <- mean(eff,na.rm=TRUE)
	se <- sd(eff,na.rm=TRUE)	
	lcl <- quantile(eff, 0.025, na.rm=TRUE)
	ucl <- quantile(eff, 0.975, na.rm=TRUE)	
	cltab <- data.frame(cbind(est, se, lcl, ucl))
	cltab <- round(cltab, 4)
	return(data.frame(cltab))
}

indirect.three.way <- function(m1, m2, x,w,z,xw,xz,wz,xwz, m, w.vals, z.vals, nsims, d) {

	# Add an id variable to the file
	d$boot_id <- 1:length(d[,1])
	ie1 <- rep(NA,nsims)
	ie2 <- rep(NA,nsims)
	ie3 <- rep(NA,nsims)
	ie4 <- rep(NA,nsims)			
	for(i in 1:nsims) {
		# pull the bootstrap sample
		unt <- sample(1:length(d[,1]), length(d[,1]), replace=TRUE)
		unt.num <- d$boot_id
		pulls <- data.frame(cbind(unt.num, unt))		
		bset <- merge(pulls, d, by.x=c("unt"), by.y=c("boot_id"), sort=TRUE)	
		
		# Run the regressions
		o1 <- lm(m1, data=bset)
		o2 <- lm(m2, data=bset)
		
		# Figure out if this is a multilevel model or a single level model
		# to extract the parameters estimates and the standard errors
		if( class(o1) == "lm" ) {
			b <- coef(o1)
			b.m <- coef(o2)
			df <- o1$df.residual
		} else {
			b <- fixef(o1)
			df <- length(o1@frame[,1])-length(b)		
		}		
		
		# Compute the simple slopes
		b1 = b[x] + b[xz]*z.vals[2] + b[xw]*w.vals[2] + b[xwz]*z.vals[2]*w.vals[2]
		b2 = b[x] + b[xz]*z.vals[2] + b[xw]*w.vals[1] + b[xwz]*z.vals[2]*w.vals[1]
		b3 = b[x] + b[xz]*z.vals[1] + b[xw]*w.vals[2] + b[xwz]*z.vals[1]*w.vals[2]
		b4 = b[x] + b[xz]*z.vals[1] + b[xw]*w.vals[1] + b[xwz]*z.vals[1]*w.vals[1]				
		
		# Pull the coefficients and calculate the indirect effect
		ie1[i] <- b1*b.m[m]				
		ie2[i] <- b2*b.m[m]				
		ie3[i] <- b3*b.m[m]				
		ie4[i] <- b4*b.m[m]														
	}

	ie <- cbind(ie1, ie2, ie3, ie4)
	est <- apply(ie, 2, mean, na.rm=TRUE)
	se <- apply(ie, 2, sd, na.rm=TRUE)	
	lcl <- apply(ie, 2, function(x) quantile(x, 0.025, na.rm=TRUE))
	ucl <- apply(ie, 2, function(x) quantile(x, 0.975, na.rm=TRUE))	

	# Format a nicer table of results	
	tab <- data.frame(round(t(rbind(est, se, lcl, ucl)),3)	)
	tab$w <- c("High", "Low", "High", "Low")	
	tab$z <- c("High", "High", "Low", "Low")
	return(tab[,c(6,5,1:4)])	
}

## This little blip is also for the mod.med stuff below - - this is old code that I need to update to save face for the embarrassing way that I coded it. 

run.model.stage1 <- function(d, form1, form2, val) {
	# first model
	m1 <- lm(form1, data=d[!is.na(d$y), ])
	a1 <- m1$coefficients['x']
	a3 <- m1$coefficients['xw']

	# second model
	m2 <- lm(form2, data=d)
	b1 <- m2$coefficients['m']

	ind <- (a1 + a3*val)*b1
	return(ind)
}


## This function is for use the below - this just bootstraps the modmed

modmed.boot.stage1 <- function(d, form1, form2, n.samp, val) {
	b.res <- rep(NA, n.samp)
	for(i in 1:n.samp) {
		# Get dataset, with replacement
		d.samp <- d[sample(nrow(d), length(d$y), replace=TRUE), ]
		b.res[i] <- run.model.stage1(d.samp, form1, form2, val)
	}
	b <- mean(b.res, na.rm=TRUE)
	se <- sd(b.res, na.rm=TRUE)
	z <- b/se
	p <- 2*pnorm(-abs(z))
	ci.ub <- quantile(b.res, 0.975)
	ci.lb <- quantile(b.res, 0.025)

	ln <- round(cbind(b,se,z,p,ci.lb, ci.ub),3)
	rownames(ln) <- round(val,2)
	return(ln)
}

## This function computes conditional indirect effects 

modmed.stage1 <- function(data, x,w,m,y, controls, n.samp, wvals) {
	form1 <- formula(paste("m ~ ", paste(controls, collapse=" + "), " + ", "x + w + xw", sep=""))
		
	form2 <- formula(paste("y ~ ", paste(controls, collapse=" + "), " + ", "x + w + xw + m", sep=""))

	# rename the variables in the dataset that matter into what they should be
	colnames(data)[colnames(data) == x] <- "x"
	colnames(data)[colnames(data) == w] <- "w"
	colnames(data)[colnames(data) == m] <- "m"
	colnames(data)[colnames(data) == y] <- "y"			
	
	# Create the interaction variable for easy usage
	data$xw <- data$x*data$w
	
	# subset the data
	data <- data[,c(controls, "x", "w", "xw", "m", "y")]

	count <- 1
	for(val in wvals) {
		ln <- modmed.boot.stage1(data, form1, form2, n.samp, val)
		if(count == 1) {
			tab <- ln
		} else {
			tab <- rbind(tab, ln)
		}	
		count <- count + 1
	}
	return(tab)
}



### Doing a dirty re-write of the functions to accoutn for second stage moderation in the psych ownership => ownership conflict effect through marking*helpseeking

run.model.stage2 <- function(d, form1, form2, val) {
	# first model
	m1 <- lm(form1, data=d[!is.na(d$y), ])
	a1 <- m1$coefficients['x']
#	a3 <- m1$coefficients['xw']

	# second model
	m2 <- lm(form2, data=d)
	b1 <- m2$coefficients['m']
	b3 <- m2$coefficients['mw']	

	ind <- (b1 + b3*val)*a1
	return(ind)
}


## This function is for use the below - this just bootstraps the modmed

modmed.boot.stage2 <- function(d, form1, form2, n.samp, val) {
	b.res <- rep(NA, n.samp)
	for(i in 1:n.samp) {
		# Get dataset, with replacement
		d.samp <- d[sample(nrow(d), length(d$y), replace=TRUE), ]
		b.res[i] <- run.model.stage2(d.samp, form1, form2, val)
	}
	b <- mean(b.res, na.rm=TRUE)
	se <- sd(b.res, na.rm=TRUE)
	z <- b/se
	p <- 2*pnorm(-abs(z))
	ci.ub <- quantile(b.res, 0.975)
	ci.lb <- quantile(b.res, 0.025)

	ln <- round(cbind(b,se,z,p,ci.lb, ci.ub),3)
	rownames(ln) <- round(val,2)
	return(ln)
}

## This function computes conditional indirect effects 

modmed.stage2 <- function(data, x,w,m,y, controls, n.samp, wvals) {
	form1 <- formula(paste("m ~ ", paste(controls, collapse=" + "), " + ", "x", sep=""))
		
	form2 <- formula(paste("y ~ ", paste(controls, collapse=" + "), " + ", "x + m + w + mw", sep=""))

	# rename the variables in the dataset that matter into what they should be
	colnames(data)[colnames(data) == x] <- "x"
	colnames(data)[colnames(data) == w] <- "w"
	colnames(data)[colnames(data) == m] <- "m"
	colnames(data)[colnames(data) == y] <- "y"			
	
	# Create the interaction variable for easy usage
	data$mw <- data$m*data$w
	
	# subset the data
	data <- data[,c(controls, "x", "w", "mw", "m", "y")]

	count <- 1
	for(val in wvals) {
		ln <- modmed.boot.stage2(data, form1, form2, n.samp, val)
		if(count == 1) {
			tab <- ln
		} else {
			tab <- rbind(tab, ln)
		}	
		count <- count + 1
	}
	return(tab)
}

######################################################
## This function calculates ICC values for a LMER results object
######################################################

lmer.icc <- function(obj) {
	residual <- attr(VarCorr(obj), "sc")^2
	int <- sapply(VarCorr(obj), function(x) attr(x, "stddev")["(Intercept)"]^2)
	icc <- sapply(int, function(x) x/(sum(int) + residual)) 
	return(icc)
}

indirect.effects.lmer <- function(m1, m2, coef1, coef2, nsims, d) {	
	# Add an id variable to the file
	d$boot_id <- 1:length(d[,1])
	i <- 1
	ie <- rep(NA,nsims)
	for(i in 1:nsims) {
		# pull the bootstrap sample
		unt <- sample(1:length(d[,1]), length(d[,1]), replace=TRUE)
		unt.num <- d$boot_id
		pulls <- data.frame(cbind(unt.num, unt))		
		bset <- merge(pulls, d, by.x=c("unt"), by.y=c("boot_id"), sort=TRUE)	
		
		# Run the regressions
		o1 <- lmer(m1, data=bset)
		o2 <- lmer(m2, data=bset)
		b1 <- fixef(o1)[coef1]
		b2 <- fixef(o2)[coef2]		
		# Pull the coefficients and calculate the indirect effect
		ie[i] <- b1*b2 
	}
	return(ie)	
}




######################################################
## The code below is used to estimate conditional indirect effects in lmer models
######################################################

## This little blip is also for the mod.med stuff below - - this is old code that I need to update to save face for the embarrassing way that I coded it. 

run.model.lmer.stage1 <- function(d, form1, form2, val) {
	# first model
	m1 <- lmer(form1, data=d[!is.na(d$y), ])
	a1 <- fixef(m1)['x']
	a3 <- fixef(m1)['xw']

	# second model
	m2 <- lmer(form2, data=d)
	b1 <- fixef(m2)['m']

	ind <- (a1 + a3*val)*b1
	return(ind)
}


## This function is for use the below - this just bootstraps the modmed

modmed.lmer.boot.stage1 <- function(d, form1, form2, n.samp, val) {
	b.res <- rep(NA, n.samp)
	for(i in 1:n.samp) {
		# Get dataset, with replacement
		d.samp <- d[sample(nrow(d), length(d$y), replace=TRUE), ]
		b.res[i] <- run.model.lmer.stage1(d.samp, form1, form2, val)
	}
	b <- mean(b.res, na.rm=TRUE)
	se <- sd(b.res, na.rm=TRUE)
	z <- b/se
	p <- 2*pnorm(-abs(z))
	ci.ub <- quantile(b.res, 0.975)
	ci.lb <- quantile(b.res, 0.025)

	ln <- round(cbind(b,se,z,p,ci.lb, ci.ub),3)
	rownames(ln) <- round(val,2)
	return(ln)
}

## This function computes conditional indirect effects 

modmed.lmer.stage1 <- function(data, x,w,m,y, controls, groupvar, n.samp, wvals) {
	form1 <- formula(paste("m ~ ", paste(controls, collapse=" + "), " + ", "x + w + xw", " + (1|",groupvar,")", sep=""))
		
	form2 <- formula(paste("y ~ ", paste(controls, collapse=" + "), " + ", "x + w + xw + m"," + (1|",groupvar,")", sep=""))

	# rename the variables in the dataset that matter into what they should be
	colnames(data)[colnames(data) == x] <- "x"
	colnames(data)[colnames(data) == w] <- "w"
	colnames(data)[colnames(data) == m] <- "m"
	colnames(data)[colnames(data) == y] <- "y"			
	
	# Create the interaction variable for easy usage
	data$xw <- data$x*data$w
	
	# subset the data
	data <- data[,c(controls, groupvar, "x", "w", "xw", "m", "y")]

	count <- 1
	for(val in wvals) {
		ln <- modmed.lmer.boot.stage1(data, form1, form2, n.samp, val)
		if(count == 1) {
			tab <- ln
		} else {
			tab <- rbind(tab, ln)
		}	
		count <- count + 1
	}
	return(tab)
}



### Doing a dirty re-write of the functions to accoutn for second stage moderation in the psych ownership => ownership conflict effect through marking*helpseeking

run.model.lmer.stage2 <- function(d, form1, form2, val) {
	# first model
	m1 <- lmer(form1, data=d[!is.na(d$y), ])
	a1 <- fixef(m1)['x']
#	a3 <- fixef(m1)['xw']

	# second model
	m2 <- lmer(form2, data=d)
	b1 <- fixef(m2)['m']
	b3 <- fixef(m2)['mw']	

	ind <- (b1 + b3*val)*a1
	return(ind)
}


## This function is for use the below - this just bootstraps the modmed

modmed.lmer.boot.stage2 <- function(d, form1, form2, n.samp, val) {
	b.res <- rep(NA, n.samp)
	for(i in 1:n.samp) {
		# Get dataset, with replacement
		d.samp <- d[sample(nrow(d), length(d$y), replace=TRUE), ]
		b.res[i] <- run.model.lmer.stage2(d.samp, form1, form2, val)
	}
	b <- mean(b.res, na.rm=TRUE)
	se <- sd(b.res, na.rm=TRUE)
	z <- b/se
	p <- 2*pnorm(-abs(z))
	ci.ub <- quantile(b.res, 0.975)
	ci.lb <- quantile(b.res, 0.025)

	ln <- round(cbind(b,se,z,p,ci.lb, ci.ub),3)
	rownames(ln) <- round(val,2)
	return(ln)
}

## This function computes conditional indirect effects 

modmed.lmer.stage2 <- function(data, x,w,m,y, controls, groupvar,n.samp, wvals) {
	form1 <- formula(paste("m ~ ", paste(controls, collapse=" + "), " + ", "x"," + (1|",groupvar,")", sep=""))
		
	form2 <- formula(paste("y ~ ", paste(controls, collapse=" + "), " + ", "x + m + w + mw"," + (1|",groupvar,")", sep=""))

	# rename the variables in the dataset that matter into what they should be
	colnames(data)[colnames(data) == x] <- "x"
	colnames(data)[colnames(data) == w] <- "w"
	colnames(data)[colnames(data) == m] <- "m"
	colnames(data)[colnames(data) == y] <- "y"			
	
	# Create the interaction variable for easy usage
	data$mw <- data$m*data$w
	
	# subset the data
	data <- data[,c(controls,groupvar, "x", "w", "mw", "m", "y")]

	count <- 1
	for(val in wvals) {
		ln <- modmed.lmer.boot.stage2(data, form1, form2, n.samp, val)
		if(count == 1) {
			tab <- ln
		} else {
			tab <- rbind(tab, ln)
		}	
		count <- count + 1
	}
	return(tab)
}

# Reverse score any negatively worded items
rev.items <- function(x, item_max){
	y<-item_max-x
	return(y)
}

# ------------------------------ #
# This is a simple function for returning the group size of a clustered dataset
# ------------------------------ #

get.groupsize = function(d, g_id) {
	require(data.table)	
	d$new_group_id = d[,g_id]
	d.dt = data.table(d)	
	d.ag = data.frame(d.dt[,list(group_size = .N), by=list(new_group_id)])
	names(d.ag)[1] = g_id
	return(d.ag)
}



# ------------------------------ #
# This stupid function does the scale stats for a single-level construct
# ------------------------------ #
scale.stats.single = function(d, constructs) {
	require(psych)
	require(lavaan)
	
	# initialize a matrix to hold the results of this
	o = matrix(, nrow=length(constructs), ncol=22)
	for(r in 1:length(constructs)) {
		v_info = constructs[[r]]
		num_items = length(v_info$items)
		
		if(num_items > 1) {		
			# Get the alpha for this variable		
			alpha_out = alpha(d[,v_info$items], na.rm=T, check.keys=F)		
			alpha_val = as.numeric(round(alpha_out$total['raw_alpha'],3))
		} else {
			alpha_val = NA
		}
		
		# If there are at least three items, do a single factor confirmatory factor analysis
		if(num_items >= 3) {
			itemlist <- paste0(v_info$items, collapse=" + ")		
			model <- paste('f1 =~ ', itemlist, sep="")
			mod <- cfa(model, data=d)
			mod.fit <- fitMeasures(mod)
			fit.ind <- as.numeric(round(mod.fit[c(3,4,5,9,10,19,20,23,29)],3))
		} else {
		
			fit.ind = rep(NA, 9)
		
		}
		o[r, 1] = v_info$name
		o[r, 2] = v_info$var
		o[r, 3] = num_items
		o[r, 4] = sum(!is.na(d[,v_info$var]))
		o[r, 5] = sum(is.na(d[,v_info$var]))		
		o[r, 6] = round(min(d[,v_info$var], na.rm=T),3)
		o[r, 7] = round(max(d[,v_info$var], na.rm=T),3)
		o[r, 8] = round(mean(d[,v_info$var], na.rm=T),3)
		o[r, 9] = round(sd(d[,v_info$var], na.rm=T),3)
		o[r, 10] = round(quantile(d[,v_info$var], .25, na.rm=T),3)
		o[r, 11] = round(quantile(d[,v_info$var], .5, na.rm=T),3)
		o[r, 12] = round(quantile(d[,v_info$var], .75, na.rm=T),3)				
		o[r, 13] = alpha_val
		o[r,14:22] = fit.ind	
	}
	
	o.d = data.frame(o, stringsAsFactors=F)
	names(o.d) = c("var_desc", "var_name", "num_items", "num_responses", "num_missing", "min", "max", "mean", "sd", "pct_25", "pct50", "pct_75", "alpha", "chisq", "df", "pvalue", "cfi", "tli", "aic", "bic", "rmsea", "srmr")
	o.d[,3:22] = lapply(o.d[,3:22], as.numeric)
	return(o.d)
}


sample.stats = function(vars, data, out.file=NA) {
	n = unlist(lapply(data[, vars],function(x) sum(!is.na(x))))	
	x = unlist(lapply(data[, vars],mean, na.rm=T))
	sd = unlist(lapply(data[, vars],sd, na.rm=T))
	pct25 = unlist(lapply(data[, vars],function(x) quantile(x, na.rm=T)[2]))	
	med = unlist(lapply(data[, vars],median, na.rm=T))
	pct75 = unlist(lapply(data[, vars],function(x) quantile(x, na.rm=T)[4]))		
	min = unlist(lapply(data[, vars],min, na.rm=T))
	max = unlist(lapply(data[, vars],max, na.rm=T))
	o = data.frame(rbind(t(n), t(x), t(sd), t(pct25), t(med), t(pct75), t(min), t(max)))
	o$param = c("n","mean", "sd", "pct25", "median", "pct75", "min", "max")
	o = o[,c("param", vars)]
	if(!is.na(out.file)) {
		write.table(o, out.file, sep=",", row.names=F)
	}
	return(o)
}





