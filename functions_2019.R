
## This function takes in a survey structure file and outputs
## a dataframe containing the constructs and their items, appropriately reverse
## scored as needed

getConstructItems = function(ck) {

	# create lists containing the item names and the construct names 
	construct_names = unique(ck[!is.na(ck$scale_group) & ck$scale_group != "", "scale_group"])

	# loop through the construct names, pulling the items, changing the name if it is reverse score to create a list of items for that particular construct name
	all_items = ck[ck$scale_group %in% construct_names, "new_colname"]
	all_constructs = matrix(nrow=length(all_items), ncol=2)
	count = 1
	for(i in construct_names) {
		items = NULL
		neg_items = NULL	
		items = ck[ck$scale_group == i & !is.na(ck$scale_group) & ck$col_rev == 0, "new_colname"]
		neg_items = ck[ck$scale_group == i & !is.na(ck$scale_group) & ck$col_rev == 1, "new_colname"]
		if(length(neg_items) > 0) {
			items = c(items, paste(neg_items, "_r", sep=""))
		}
		all_constructs[count:(count+length(items)-1), 1] = i
		all_constructs[count:(count+length(items)-1), 2] = items	
		count = count+length(items)
	}
	# This is the construct file
	all_constructs = data.frame(all_constructs, stringsAsFactors=F)
	names(all_constructs) = c("construct_name", "item_name")
	return(all_constructs)
}

### ### ### ### ### ### ### ### ### ### ### 
# This function reverse scores negatively worded items
# in a column key file
### ### ### ### ### ### ### ### ### ### ### 

reverseNegativeItems = function(ck, d) {
	neg.items = ck[ck$col_rev == 1, "new_colname"]	
	for(i in neg.items) {	
		rev_var_name = paste(i, "r", sep="_")	
		d[,rev_var_name] = (ck[ck$new_colname == i, "scale_max"]+1) - d[,i]			
	}
	return(d)
}

### ### ### ### ### ### ### ### ### ### ### 
# This function checks the item range based on the scale setup file
### ### ### ### ### ### ### ### ### ### ### 
checkItemRange = function(ck, d) {

	# Loop through the column key file (this is the cumbersome approach; for larger files and smaller actions use lapply
	outOfRange = c()
	for(i in 1:nrow(ck)){

		# get the current column name
		var_name = ck[i,"new_colname"]

		# Check to make sure that items are not out of bounds if they have a scale_max and scale_min
		if(!is.na(ck[i,"scale_max"])) {
	
			# are any scores outside of the range
			if(min(d[,var_name], na.rm=T) < ck[i,"scale_min"] | max(d[,var_name], na.rm=T) > ck[i,"scale_max"]) {		
				outOfRange = c(outOfRange, var_name)
			}	
		}
	}
	return(outOfRange)
}

### ### ### ### ### ### ### ### ### ### ### 
# This function calculates a scale score based on a column key file.
# It uses the get constructs function
### ### ### ### ### ### ### ### ### ### ### 

getScaleScores = function(ck, d) {
	all_constructs = getConstructItems(ck) 
	for(i in unique(all_constructs$construct_name)) {
		items = all_constructs[all_constructs$construct_name == i, "item_name"]
		if(length(items) > 1 ) {
			d[,i] = rowMeans(d[,c(items)], na.rm=T)
		} else {
			d[,i] = d[, items]
		}
	}
	return(d)
}


# Reverse score any negatively worded items
rev.items <- function(x, item_max){
	y<-(item_max+1)-x
	return(y)
}


# ------------------------------ #
# This old (and stupidly written, clunky) function does the scale stats for a single-level construct
# ------------------------------ #
getScaleStats = function(ck, d) {
	require(psych)
	require(lavaan)
	
	all_constructs = getConstructItems(ck) 	
	
	# initialize a matrix to hold the results of this
	constructs = unique(all_constructs$construct_name)
	o = matrix(, nrow=length(constructs), ncol=22)

	for(r in 1:length(constructs)) {
		alpha_out = NA
		alpha_val = NA
		num_items = nrow(all_constructs[all_constructs$construct_name == constructs[r], ])

		# If number of items is greater than one, get an alpha
		if(num_items > 1) {
			alpha_out = alpha(d[,all_constructs[all_constructs$construct_name == constructs[r], "item_name"]], na.rm=TRUE, check.keys=F)		
			alpha_val = as.numeric(round(alpha_out$total['raw_alpha'],3))		
		}
		
		# if there are at least three items, do a single factor confirmatory factor analysis
		if(num_items >= 3) {
			item_rh <- paste0(all_constructs[all_constructs$construct_name == constructs[r], "item_name"], collapse=" + ")				
			model <- paste('f1 =~ ', item_rh, sep="")
			mod <- cfa(model, data=d)
			mod.fit <- fitMeasures(mod)
			fit.ind <- as.numeric(round(mod.fit[c(3,4,5,9,10,19,20,23,29)],3))
		} else {
			fit.ind = rep(NA, 9)		
		}
		# clunky way of outputting the information I want 
		o[r, 1] =constructs[r]
		o[r, 2] = constructs[r]
		o[r, 3] = num_items
		o[r, 4] = sum(!is.na(d[,constructs[r]]))
		o[r, 5] = sum(is.na(d[,constructs[r]]))		
		o[r, 6] = round(min(d[,constructs[r]], na.rm=T),3)
		o[r, 7] = round(max(d[,constructs[r]], na.rm=T),3)
		o[r, 8] = round(mean(d[,constructs[r]], na.rm=T),3)
		o[r, 9] = round(sd(d[,constructs[r]], na.rm=T),3)
		o[r, 10] = round(quantile(d[,constructs[r]], .25, na.rm=T),3)
		o[r, 11] = round(quantile(d[,constructs[r]], .5, na.rm=T),3)
		o[r, 12] = round(quantile(d[,constructs[r]], .75, na.rm=T),3)				
		o[r, 13] = alpha_val
		o[r,14:22] = fit.ind				
	}	
	o.d = data.frame(o, stringsAsFactors=F)
	names(o.d) = c("var_desc", "var_name", "num_items", "num_responses", "num_missing", "min", "max", "mean", "sd", "pct_25", "pct50", "pct_75", "alpha", "chisq", "df", "pvalue", "cfi", "tli", "aic", "bic", "rmsea", "srmr")
	o.d[,3:22] = lapply(o.d[,3:22], as.numeric)
	return(o.d)
}

getGroupScaleStats = function(ck, d, g) {
	# ck is the colkey file
	# d is the dataset
	# g is the grouping variable
	require(psych)
	require(lavaan)
	require(multilevel)
	
	all_constructs = getConstructItems(ck) 	
	conItems = c(unique(all_constructs$construct_name), all_constructs$item_name)
	constructs = unique(all_constructs$construct_name)	

	# Need to adjust the number of columns to fit with the additional stats we're going to drop
	o = matrix(, nrow=length(constructs), ncol=32)	

	# Create a group-level dataset containing just the focal items and construct scale scores
	d.agg = aggregate(d[,conItems], by=list(d[,g]), mean, na.rm=T)

	


	for(r in 1:length(constructs)) {

		num_items = nrow(all_constructs[all_constructs$construct_name == constructs[r], ])

		# If number of items is greater than one, get an alpha
		alpha_out = NA
		alpha_val = NA
		alpha_out_grp = NA
		alpha_val_grp = NA		

		if(num_items > 1) {
			alpha_out = alpha(d[,all_constructs[all_constructs$construct_name == constructs[r], "item_name"]], na.rm=TRUE, check.keys=F)		
			alpha_val = as.numeric(round(alpha_out$total['raw_alpha'],3))		

			alpha_out_grp = alpha(d.agg[,all_constructs[all_constructs$construct_name == constructs[r], "item_name"]], na.rm=TRUE, check.keys=F)		
			alpha_val_grp = as.numeric(round(alpha_out_grp$total['raw_alpha'],3))		
		}
		
		# if there are at least three items, do a single factor confirmatory factor analysis
		if(num_items >= 3) {
			item_rh <- paste0(all_constructs[all_constructs$construct_name == constructs[r], "item_name"], collapse=" + ")				
			model <- paste('f1 =~ ', item_rh, sep="")
			mod <- cfa(model, data=d)	
			if(lavInspect(mod, what="converged")) {
				mod.fit <- fitMeasures(mod)
				fit.ind <- as.numeric(round(mod.fit[c(3,4,5,9,10,19,20,23,29)],3))				
			} else {
				fit.ind = rep(NA, 9)			
			}
		} else {
			fit.ind = rep(NA, 9)		
		}


		# Calculate the group statistics

		icc.formula = as.formula(paste(constructs[r], "~ as.factor(",g,")",sep=""))
		icc.out <- aov(icc.formula, data=d)
		icc1 <-round(ICC1(icc.out),3)
		icc1p <- round(summary(icc.out)[[1]][["Pr(>F)"]], 4)[1]
		icc2 <-round(ICC2(icc.out),3)		

		num.options = unique(ck[ck$scale_group == constructs[r], "scale_max"])-unique(ck[ck$scale_group == constructs[r], "scale_min"]) + 1
		uniform.variance = (num.options^2-1)/12
		if(num_items <=1) {
			rwg.out<-rwg(d[,all_constructs[all_constructs$construct_name == constructs[r], "item_name"]],d[,g], uniform.variance)
			rwgjmed<-round(median(rwg.out$rwg, na.rm=TRUE),3)			
			rwgjx<-round(mean(rwg.out$rwg, na.rm=TRUE),3)
			rwgjsd<-round(sd(rwg.out$rwg, na.rm=TRUE),3)
			rwgjmin<-round(min(rwg.out$rwg, na.rm=TRUE),3)
			rwgjmax<-round(max(rwg.out$rwg, na.rm=TRUE),3)        								
		} else {
			rwg.out<-rwg.j(d[,all_constructs[all_constructs$construct_name == constructs[r], "item_name"]],d[,g], uniform.variance)
			rwgjmed<-round(median(rwg.out$rwg, na.rm=TRUE),3)			
			rwgjx<-round(mean(rwg.out$rwg, na.rm=TRUE),3)
			rwgjsd<-round(sd(rwg.out$rwg, na.rm=TRUE),3)
			rwgjmin<-round(min(rwg.out$rwg, na.rm=TRUE),3)
			rwgjmax<-round(max(rwg.out$rwg, na.rm=TRUE),3)        								
		}

		# clunky way of outputting the information I want 
		o[r, 1] =constructs[r]
		o[r, 2] = constructs[r]
		o[r, 3] = num_items
		o[r, 4] = sum(!is.na(d[,constructs[r]]))
		o[r, 5] = sum(is.na(d[,constructs[r]]))		
		o[r, 6] = round(min(d[,constructs[r]], na.rm=T),3)
		o[r, 7] = round(max(d[,constructs[r]], na.rm=T),3)
		o[r, 8] = round(mean(d[,constructs[r]], na.rm=T),3)
		o[r, 9] = round(sd(d[,constructs[r]], na.rm=T),3)
		o[r, 10] = round(quantile(d[,constructs[r]], .25, na.rm=T),3)
		o[r, 11] = round(quantile(d[,constructs[r]], .5, na.rm=T),3)
		o[r, 12] = round(quantile(d[,constructs[r]], .75, na.rm=T),3)				
		o[r, 13] = alpha_val
		o[r,14:22] = fit.ind	

		o[r, 23] = sum(!is.na(d.agg[,constructs[r]]))
		o[r, 24] = alpha_val_grp
		o[r, 25] = icc1
		o[r, 26] = icc1p		
		o[r, 27] = icc2		
		o[r, 28] = rwgjmed
		o[r, 29] = rwgjx
		o[r, 30] = rwgjsd
		o[r, 31] = rwgjmin
		o[r, 32] = rwgjmax
	}	
	o.d = data.frame(o, stringsAsFactors=F)
	names(o.d) = c("var_desc", "var_name", "num_items", "num_responses", "num_missing", "min", "max", "mean", "sd", "pct_25", "pct50", "pct_75", "alpha", "chisq", "df", "pvalue", "cfi", "tli", "aic", "bic", "rmsea", "srmr", "num_groups", "alpha_group", "icc1", "icc1p", "icc2", "rwg_med", "rwg_mean", "rwg_sd", "rwg_min", "rwg_max")
	o.d[,3:32] = lapply(o.d[,3:32], as.numeric)
	return(o.d)
}




# ------------------------------ #
# This is an update of the correlation matrix 
# ------------------------------ #
getCorrelationMatrix = function(ck, d) {
	require(Hmisc)
	require(matrixStats)	
	all_constructs = getConstructItems(ck) 	
	vars = unique(all_constructs$construct_name)
	x = d[,vars]
	x <- as.matrix(x)
	R <- rcorr(x)$r
	p <- rcorr(x)$P
	M <- formatC(round(colMeans(x, na.rm=T),2), digits=2, format="f")
	SD <- formatC(round(colSds(x, na.rm=T),2),digits=2, format="f")

	## trunctuate the matrix that holds the correlations to two decimal
	R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]

	## build a new matrix that includes the correlations with their apropriate stars
	Rnew <- R
	diag(Rnew) <- ""
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


# ------------------------------ #
# This old (and stupidly written, clunky) function does a correlation matrix
# ------------------------------ #

correlation.matrix <- function(x, diagval=NA){
	require(Hmisc)
	require(matrixStats)
	x <- as.matrix(x)
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
