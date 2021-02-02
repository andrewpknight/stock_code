############################################################
# This file contains general functions for processing sociometric badge data
############################################################
#fname <- "./b_data/test_dat_1min.xlsx"
require(xlsx)


read.b.data <- function(fname, report) {
	# Some computers will require the following option to be set. So far, I know my MBP requires this, but my Mac Mini does not. These options set the memory of the java virtual machine a bit higher than it might be otherwise. This will be important for very large datasets. Once we know more about differences across environments, we can standardize the call to do this or not.
#	 options( java.parameters = "-Xmx4g" )
	# options( java.parameters = "-Xmx2g" )
	require(xlsx)
	# First, read in the summary sheet to get the information that is loaded there needed to parse the other columns
	summary_sheet <- read.xlsx2(fname, as.data.frame=TRUE, header=TRUE, sheetName="summary", stringsAsFactors=FALSE, colClasses=c("character", "character", rep("numeric",20)))

	# create new column names and store the individual names in a different vector
	people_n <- length(summary_sheet)-2
	ssi_ids <- colnames(summary_sheet)[3:length(summary_sheet)]
#	colnames(summary_sheet)[3:length(summary_sheet)] <- paste("p",1:people_n, sep="_")
#	people_column_names <- paste("p",1:people_n, sep="")
	people_column_names <- ssi_ids


	#### Specify the sheet information to pull ####
	
	# Here's a little function to create column names when things are in wide-form
	b_get_col_names <- function(y, people_column_names, sval) {
		v <- merge(y,people_column_names)
		nms <- paste(v$y,v$x, sep=sval)
		return(nms)
	}

	# This is the data structure for the different sheets
	
	t_BM_activity1 <- list(
		name="t_BM_activity1", 
		var_prefix=c("BM_activity"), 
		desc="BM Activity", 
		level="individual", 
		temporal="dynamic",
		cols=c("time_stamp", paste(people_column_names, "BM_activity",sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)+1)
	)
	
	t_BM_rate1 <- list(
		name="t_BM_rate1", 
		var_prefix=c("BM_rate"), 
		desc="BM Rate", 
		level="individual", 
		temporal="dynamic",
		cols=c("time_stamp", paste(people_column_names, "BM_rate",sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)+1)
	)
	
	t_BM_consistency1 <- list(
		name="t_BM_consistency1", 
		var_prefix=c("BM_consistency"), 
		desc="BM Consistency", 
		level="individual", 
		temporal="dynamic",
		cols=c("time_stamp", paste(people_column_names, "BM_consistency",sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)+1)
	)	

	t_BM_mirroring1 <- list(
		name="t_BM_mirroring1", 
		var_prefix=c("bm_mirroring"),
		desc="Body Motion Mirroring", 
		level="dyad",
		temporal="dynamic",
		cols=c("time_stamp", paste(b_get_col_names(people_column_names, people_column_names, sval="_"), "bm_mirroring", sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)^2+1)
	)

	t_posture_activity1 <- list(
		name="t_posture_activity1", 
		var_prefix=c("posture_activity"), 
		desc="Posture Activity", 
		level="individual", 
		temporal="dynamic",
		cols=c("time_stamp", paste(people_column_names, "posture_activity",sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)+1)
	)

	t_posture_rate1 <- list(
		name="t_posture_rate1", 
		var_prefix=c("posture_rate"), 
		desc="Posture Rate", 
		level="individual", 
		temporal="dynamic",
		cols=c("time_stamp", paste(people_column_names, "posture_rate",sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)+1)
	)

	t_posture_mirroring1 <- list(
		name="t_posture_mirroring1", 
		var_prefix=c("posture_mirroring"), 
		desc="Body Posture Mirroring", 
		level="dyad", 
		temporal="dynamic",
		cols=c("time_stamp", paste(b_get_col_names(people_column_names, people_column_names, sval="_"), "posture_mirroring", sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)^2+1)	
	)

	n_tt_turntaking <- list(
		name="n_tt_turntaking", 
		var_prefix=c("turntaking"),
		desc="Turntaking Matrix", 
		level="dyad", 
		temporal="static", 
		cols=c("ssi_id", paste(people_column_names, "turntaking", sep="%")),
		data_row_start=2, 
		column_in_types=c("character", rep("numeric", length(people_column_names)))
	)

	r_tt_turntaking1 <- list(
		name="r_tt_turntaking1", 
		var_prefix=c("turntaking"),
		desc="Turntaking Detail",
		level="individual", 
		temporal="static",
		cols=c("ssi_id", "n_turns_total", "n_turns_by_badge", "n_turns_after", "n_turns_self", "n_segments_speaking", "n_segments_silent", "n_turns_per_second", "avg_segment_speaking", "avg_segment_silent", "n_interruptions_successful", "n_interruptions_unsuccessful"),
		data_row_start=2, 
		column_in_types=c("character", rep("numeric", 11))	
	)
	# Create the speech_profile variables
	t_speech_profile1 <- list(
		name="t_speech_profile1", 
		var_prefix= c("speaking", "overlap", "listening", "silent", "total_speaking", "total_silent"),
		desc="Speech Profile Detail", 
		level="individual", 
		temporal="dynamic", 
		cols=c("time_stamp", b_get_col_names(c("speaking", "overlap", "listening", "silent", "total_speaking", "total_silent"), people_column_names, sval="%")),
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)*6+1)
	)

	t_audio_front_vol_mirroring1 <- list(
		name="t_audio_front_vol_mirroring1", 
		var_prefix=c("audio_front_vol_mirroring"), 
		desc="Audio Front Volume Mirroring", 
		level="dyad", 
		temporal="dynamic",
		cols=c("time_stamp", paste(b_get_col_names(people_column_names, people_column_names, sval="_"), "audio_front_vol_mirroring", sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)^2+1)	
	)

	t_audio_front_f0_mirroring1 <- list(
		name="t_audio_front_f0_mirroring1", 
		var_prefix=c("audio_front_f0_mirroring"), 
		desc="Audio Front Mirroring", 
		level="dyad", 
		temporal="dynamic",
		cols=c("time_stamp", paste(b_get_col_names(people_column_names, people_column_names, sval="_"), "audio_front_f0_mirroring", sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)^2+1)	
	)	

	t_audio_back_vol_mirroring1 <- list(
		name="t_audio_back_vol_mirroring1", 
		var_prefix=c("audio_back_vol_mirroring"), 
		desc="Audio back Volume Mirroring", 
		level="dyad", 
		temporal="dynamic",
		cols=c("time_stamp", paste(b_get_col_names(people_column_names, people_column_names, sval="_"), "audio_back_vol_mirroring", sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)^2+1)	
	)

	t_audio_back_f0_mirroring1 <- list(
		name="t_audio_back_f0_mirroring1", 
		var_prefix=c("audio_back_f0_mirroring"), 
		desc="Audio back Mirroring", 
		level="dyad", 
		temporal="dynamic",
		cols=c("time_stamp", paste(b_get_col_names(people_column_names, people_column_names, sval="_"), "audio_back_f0_mirroring", sep="%")), 
		data_row_start=3, 
		column_in_types=rep("numeric", length(people_column_names)^2+1)	
	)		

	t_audio_front_pitch1 <- list(
		name="t_audio_front_pitch1", 
		var_prefix=c("pitch", "speech_volume"),
		desc="Audio Front Pitch", 
		level="individual", 
		temporal="dynamic", 
		cols=c("time_stamp", b_get_col_names(c("pitch", "speech_volume"), people_column_names, sval="%")),
		data_row_start=3,
		column_in_types=rep("numeric", length(people_column_names)*2+1)	
	)

	t_audio_back_pitch1 <- list(
		name="t_audio_back_pitch1", 
		var_prefix=c("pitch", "speech_volume"),
		desc="Audio Back Pitch", 
		level="individual", 
		temporal="dynamic",
		cols=c("time_stamp", b_get_col_names(c("pitch", "speech_volume"), people_column_names, sval="%")),
		data_row_start=3,
		column_in_types=rep("numeric", length(people_column_names)*2+1)		
	)
	if(report == "group") {
		sheet_names <- list(n_tt_turntaking, r_tt_turntaking1, t_speech_profile1, t_BM_mirroring1, t_posture_mirroring1, n_tt_turntaking, r_tt_turntaking1, t_speech_profile1, t_audio_front_vol_mirroring1, t_audio_front_f0_mirroring1, t_audio_back_vol_mirroring1, t_audio_back_f0_mirroring1)
	} else 
	if(report == "experimental"){
		sheet_names <- list(t_BM_mirroring1, t_BM_activity1, t_BM_rate1, t_BM_consistency1, t_posture_rate1, t_posture_activity1, t_posture_mirroring1, n_tt_turntaking, r_tt_turntaking1, t_speech_profile1, t_audio_front_pitch1, t_audio_back_pitch1)	
	} else
	if(report == "dyad") {
		sheet_names <- list(t_BM_mirroring1, t_posture_mirroring1, n_tt_turntaking, r_tt_turntaking1, t_speech_profile1, t_audio_front_vol_mirroring1, t_audio_front_f0_mirroring1, t_audio_back_vol_mirroring1, t_audio_back_f0_mirroring1)			
	}

	dat_lst <- list("summary_sheet"=summary_sheet)

	for(i in sheet_names) {
		d <- read.xlsx2(fname, as.data.frame=TRUE, header=FALSE, sheetName=i$name, stringsAsFactors=FALSE, startRow=i$data_row_start, colIndex=1:length(i$cols), colClasses=i$column_in_types)
		colnames(d) <- i$cols
	
		### Convert the timestamp if this is a dynamic dataset ###
		if(i$temporal == "dynamic") {
			d$time_stamp <- as.POSIXct((as.numeric(d$time_stamp)-25569)*86400, tz="CDT", origin="1970-01-01")
		}	
	
		assign(i$name, d)		
		## Go through a few specific formats to produce better-structured datasets
		
		if(i$level == "dyad" && i$temporal == "dynamic") {
			count <- 1
			for( ego in people_column_names ) {
				for( alter in people_column_names ) {
					var_names <- paste(ego, "_", alter, "%", i$var_prefix, sep="")
					tmp <- d[, c("time_stamp", var_names)]
					colnames(tmp) <- c("time_stamp", i$var_prefix)
					tmp$ego <- ego
					tmp$alter <- alter
					if( count == 1 ) {
						new.d <- tmp
					} else {
						new.d <- rbind(new.d, tmp)
					}
					count <- count + 1
				}
			}
			assign(i$name, new.d)
		}
	
		if(i$level == "individual" && i$temporal == "dynamic") {
			count <- 1
			for( ego in people_column_names ) {
					var_names <- paste(ego, "%", i$var_prefix, sep="")
				tmp <- d[, c("time_stamp", var_names)]
				colnames(tmp) <- c("time_stamp", i$var_prefix)
				tmp$ssi_id <- ego
				if( count == 1 ) {
					new.d <- tmp
				} else {
					new.d <- rbind(new.d, tmp)
				}
				count <- count + 1
			}
			assign(i$name, new.d)	
		}
		dat_lst[[i$name]] <- eval(as.name(i$name))
	}
	return(dat_lst)
}


###### SIMPLE BAR CHART FUNCTION FOR PRODUCING HORIZONTAL BAR CHARTS OF BADGE METRICS ######
b.barchart <- function(metric, d, out_name, xlabel, ylabel, singlepage) {
	if(singlepage) {
#		quartz(width=8, height=4)	
		par(mar=c(5,10,10,1))
	} else {
		quartz(width=9,height=3,type="pdf",file=out_name,dpi=300)
		par(mar= c(5,10,1,1))
		par(lwd = 1)
		par(xpd=TRUE)
		par(mfcol=c(1,1), family="Arial", mgp=c(3,1,0), bg="white")		
	}

	out<-barplot(d[,metric], xlim=c(0,max(d[,metric])*1.2), main="", ylab="", beside=TRUE, las=1, border=NA, cex.names=1.5, horiz=TRUE, axes=TRUE, col="#777777", xlab=ylabel)

	# These are the row labels
	text(0, out, d[,ylabel], xpd=TRUE, pos=2, col="#222222", cex=.75)
	if(!singlepage) {
		dev.off()			
	}
}

b.balancechart <- function(metric, d, balance_point, out_name, xlabel, point_label, singlepage) {
	if(singlepage) {
#		quartz(width=9, height=3)	
		par(mar=c(5,10,10,1))
	} else {
		quartz(width=9,height=3,type="pdf",file=out_name,dpi=300)
		par(mar= c(1,0,0,0))
		par(lwd = 1)
		par(xpd=TRUE)
		par(mfcol=c(1,1), family="Arial", mgp=c(3,1,0), bg="white")		
	}
	people_n <- length(d[,metric])
	xmin <- min(d[,metric], na.rm=TRUE)-balance_point/people_n
	xmax <- max(d[,metric], na.rm=TRUE)+balance_point/people_n
	plot(d[,metric], 1:people_n, xlim=c(xmin, xmax), axes=F, xlab=xlabel, ylab="", ylim=c(-1,people_n+1))
	grad.length <- 100
	# do the gradient from the far left to the border of the goal zone
	goal_left <- balance_point-balance_point*.30
	xvals <- seq(xmin, goal_left, length.out=grad.length)
	colfunc <- colorRampPalette(c("red", "yellow"))
	xcols <- colfunc(grad.length)
	for(i in 1:length(xvals)) {
		rect(xleft=xvals[i], ybottom=0, xright=xvals[i+1], ytop=people_n+1, border=NA, col=xcols[i])
	}
	
	goal_right <- balance_point+balance_point*.30
	xvals <- seq(goal_right, xmax, length.out=grad.length)
	colfunc <- colorRampPalette(c("yellow", "red"))
	xcols <- colfunc(grad.length)
	for(i in 1:length(xvals)) {
		rect(xleft=xvals[i], ybottom=0, xright=xvals[i+1], ytop=people_n+1, border=NA, col=xcols[i])
	}
	
	xvals <- seq(goal_left, balance_point, length.out=grad.length)
	colfunc <- colorRampPalette(c("yellow", "green"))
	xcols <- colfunc(grad.length)
	for(i in 1:length(xvals)) {
		rect(xleft=xvals[i], ybottom=0, xright=xvals[i+1], ytop=people_n+1, border=NA, col=xcols[i])
	}	

	xvals <- seq(balance_point,goal_right, length.out=grad.length)
	colfunc <- colorRampPalette(c("green", "yellow"))
	xcols <- colfunc(grad.length)
	for(i in 1:length(xvals)) {
		rect(xleft=xvals[i], ybottom=0, xright=xvals[i+1], ytop=people_n+1, border=NA, col=xcols[i])
	}	
	
	text(d[,metric], 1:people_n, d[,point_label], cex=.75)	
	arrows(goal_right, -.5, xmax, -.5)
	lines(c(goal_right, goal_right), c(-.25,-.75))
	arrows(xmin, -.5, goal_left, -.5, code=1)	
	lines(c(goal_left, goal_left), c(-.25,-.75))	
	text(balance_point, -.5, "Goal Zone")
	text(xmin+(goal_left-xmin)/2, -.85, "Not Enough Speaking", cex=.8)	
	text(goal_right+(xmax-goal_right)/2, -.85, "Too Much Speaking", cex=.8)		
	if(!singlepage) {
		dev.off()			
	}	
}

b.balancechart2 <- function(metric, d, balance_point, out_name, xlabel, point_label, singlepage) {
	if(singlepage) {
#		quartz(width=9, height=3)	
		par(mar=c(5,10,10,1))
	} else {
		quartz(width=9,height=3,type="pdf",file=out_name,dpi=300)
		par(mar= c(1,0,0,0))
		par(lwd = 1)
		par(xpd=TRUE)
		par(mfcol=c(1,1), family="Arial", mgp=c(3,1,0), bg="white")		
	}
	people_n <- length(d[,metric])
	xmin <- min(d[,metric], na.rm=TRUE)-balance_point/people_n
	xmax <- max(d[,metric], na.rm=TRUE)+balance_point/people_n
	
	# Put people in different categories
	inc <- .20
	r1 <- balance_point + balance_point*inc
	r2 <- balance_point + balance_point*inc*2	
	l1 <- balance_point - balance_point*inc
	l2 <- balance_point - balance_point*inc*2	
	
	# This is inside the balance point box
	d$xpt <- ifelse(d[,metric] >= l1 & d[,metric] <= r1, 0, NA)
	# This is in the l2 box
	d$xpt <- ifelse(d[,metric] < l2, -2, d$xpt)
	# This is in the l1 box
	d$xpt <- ifelse(d[,metric] < l1 & d[,metric] >=l2, -1, d$xpt)	
	# This is in the r1 box
	d$xpt <- ifelse(d[,metric] > r1 & d[,metric] <= r2, 1, d$xpt)
	# This is in the r2 box
	d$xpt <- ifelse(d[,metric] > r2, 2, d$xpt)				
	
	
	plot(d$xpt, 1:people_n, xlim=c(-2.5, 2.5), axes=F, xlab=xlabel, ylab="", ylim=c(-1,people_n+1))
	
	red <- "#F59759"
	yellow <- "#F5CA59"
	green <- "#93BC69"
	
	# Draw rectangles
	rect(xleft=-2.5, ybottom=0, xright=-1.5, ytop=people_n+1, border=NA, col=red)
	rect(xleft=-1.5, ybottom=0, xright=-0.5, ytop=people_n+1, border=NA, col=yellow)	
	rect(xleft=-0.5, ybottom=0, xright=0.5, ytop=people_n+1, border=NA, col=green)		
	rect(xleft=.5, ybottom=0, xright=1.5, ytop=people_n+1, border=NA, col=yellow)			
	rect(xleft=1.5, ybottom=0, xright=2.5, ytop=people_n+1, border=NA, col=red)				
	

	text(d$xpt, 1:people_n, d[,point_label], col="#FFFFFF", cex=.75)	
	arrows(0.5, -.5, 2.5, -.5)
	lines(c(0.5, 0.5), c(-.25,-.75))
	arrows(-2.5, -.5, -0.5, -.5, code=1)	
	lines(c(-0.5, -0.5), c(-.25,-.75))	
	text(0, -.5, "Goal Zone")
	text(-1.5, -.85, "Not Enough Speaking", cex=.8)	
	text(1.5, -.85, "Too Much Speaking", cex=.8)		
	if(!singlepage) {
		dev.off()			
	}	
}


b.pacechart <- function(x,y, d, out_name, ylabel, singlepage) {
	indiv_ids <- unique(d$indiv_id)
	people_n <- length(indiv_ids)
	if(singlepage) {
#		quartz(width=9, height=3)	
		par(mar=c(5,10,10,1))
	} else {
		quartz(width=9,height=3,type="pdf",file=out_name,dpi=300)
		par(mar= c(0,5,0,0))
		par(lwd = 1)
		par(xpd=TRUE)
		par(mfcol=c(1,1), family="Arial", mgp=c(3,1,0), bg="white")		
	}
		# construct a blank plot 
	plot(d[,x], d[,y], xlim=c(min(d[,x], na.rm=TRUE), max(d[,x], na.rm=TRUE)), ylim=c(0, people_n+1), type="n", axes=FALSE, xlab="Timeline of Group Meeting", ylab="")
	# Loop through team members and add their lines to the graph
	count <- 1
	d$bin <- ifelse(d[,y] > 15 | d[,y], 1, 0)
	d$bin <- ifelse(is.na(d[,y]), 0, d$bin)	
	for(indiv in indiv_ids) {
		d.sub <- d[d$indiv_id == indiv, ]
		symbols(d.sub[,x], rep(count, length(d.sub[,x])), circles=d.sub[,y], inches=.10, fg="white", bg="dark blue", add=TRUE)   	
		# These are the row labels
		count <- count + 1
	}
	labs <- unique(d[,c("indiv_id", ylabel)])
	text(min(d.sub[,x], na.rm=TRUE), 1:length(indiv_ids), labs[,ylabel], xpd=TRUE, pos=2, col="#222222", cex=0.75)			
	if(!singlepage) {
		dev.off()			
	}
}

create.group.report <- function(tm_id, indiv_labels, focal_indiv, out_prefix, out_suffix, singlepage=FALSE) {

	### SETUP THE OUTPUT DEPENDING ON THE USER INPUT ###
	if(singlepage) {
		out_name <- paste(out_prefix, "report_", out_suffix, ".pdf", sep="")

		quartz(width=11,height=8.5,type="pdf",file=out_name,dpi=300)
		# 1/2 inch margins all around
		par(mar= c(0,0,1,1), xpd=TRUE, family="Arial", bg="white")
		m<-matrix(c(1,1,2,3,4,5), 3,2,byrow=TRUE)
		layout(m, widths=c(.8,rep(.8,4)), heights=c(.1, rep(.35,4)))
		
		# Get the name for this person to write the header information
		name <- roster[roster$indiv_id == focal_indiv, indiv_labels]
		plot.new()
		text(.5, .5, name)							
	}

	# Read the data for the selected team
	fname <- paste("./team_",tm_id, ".xlsx", sep="")
	o <- read.b.data(fname, report="group")
	### GET SUMMARY INFORMATION, WHICH CONTAINS PARTICIPANT NAMES ###
	d <- o$summary_sheet
	people_n <- length(d)-2
	ssi_ids <- colnames(d)[3:length(colnames(d))]
	people_column_names <- paste("p",1:people_n, sep="")
	peeps <- data.frame(ssi_ids, people_column_names, stringsAsFactors=F)
	colnames(peeps) <- c("ssi_id", "ssi_temp_id")
	
	# Add the names from the roster to the peeps set
	peeps <- merge(roster, peeps, by=c("ssi_id"), all.x=TRUE)
	peeps$label <- peeps[,indiv_labels]	
	if(!is.na(focal_indiv)) {
		peeps$label <- ifelse(peeps$indiv_id == focal_indiv, peeps$label, " ")
	} 
	
	### NEED TO FIGURE OUT THE RIGHT TIMEFRAME FOR THIS SPECIFIC TEAM ###
	# First, roll-up the speech_profile sheet, looking for the latest time when speaking is zero across all participants (treat the first moment of speech as the start of the data)
	d <- o$t_speech_profile1
	d <- merge(d, peeps, by=c("ssi_temp_id"))	
	spch_agg <- aggregate(d[, c("total_speaking")], by=list(d$time_stamp), sum, na.rm=TRUE)
	colnames(spch_agg) <- c("time_stamp", "total_speaking_sum")
	start_time <- NA
	end_time <- NA
	for(i in 1:length(spch_agg$time_stamp)) {
		if(spch_agg[i, c("total_speaking_sum")] > 0 && is.na(start_time)) {
			start_time <- spch_agg[i, c("time_stamp")]
		}
		if(spch_agg[i, c("total_speaking_sum")] == 0 && !is.na(start_time) && is.na(end_time)) {
			end_time <- spch_agg[i, c("time_stamp")]					
		}	
		
		if(spch_agg[i, c("total_speaking_sum")] > 0 && !is.na(start_time) && !is.na(end_time)) {
			end_time <- spch_agg[i, c("time_stamp")]					
		}			
	}
#	spch_agg[,c("time_stamp", "total_speaking_sum")]
	# Replace the longer speech profile dataset with this one
	d <- d[d$time_stamp > start_time & d$time_stamp < end_time, ]
	
	### FIRST GRAPH: WHAT IS YOUR AVERAGE SPEAKING SEGMENT LENGTH ###
	d <- o$r_tt_turntaking1
	colnames(d)[1] <- "ssi_id"
	d <- merge(d, peeps, by=c("ssi_id"))
	out_name1 <- paste(out_prefix, "avg_segment_speaking_", out_suffix, ".pdf", sep="")	
	b.barchart("avg_segment_speaking", d, out_name1, "Average Speaking Segment Length", "label", singlepage)
	
	### SECOND GRAPH: WHAT IS YOUR BALANCE OF SPEAKING AND LISTENING
	##KW: "total speaking" here incorporates the speech overlap segments which we know are grossly overestimated --> I propose using "speaking" instead##
	d <- data.frame(t(o$summary_sheet[,3:(length(people_column_names)+2)]), stringsAsFactors=F)
	colnames(d) <- o$summary_sheet[,1]
	d$ssi_id <- rownames(d)
	d <- merge(d, peeps, by=c("ssi_id"))

# 	First, calculate each individual's percentage of the team's total speaking time
	d$tm_sum <- round(100*d$speech_profile_speaking/sum(d$speech_profile_speaking, na.rm=TRUE), 0)

#	Second, calculate the balance window, which is just the proportion of the team
	balance_point <- 100/length(d$tm_sum)	

	metric <- "tm_sum"
	xlabel <- ""	
	out_name2 <- paste(out_prefix, "speaking_balance_", out_suffix, ".pdf", sep="")			
	b.balancechart2("tm_sum", d, balance_point, out_name2, xlabel, "label", singlepage)
	
	### THIRD GRAPH: WHAT IS YOUR PACING OF CONTRIBUTIONS OVER TIME?
	d <- o$t_speech_profile1

	# cut out any null points of data using the info above
	d <- merge(d, peeps, by=c("ssi_temp_id"))
	d <- d[d$time_stamp > start_time & d$time_stamp < end_time, ]	
	x <- "time_stamp"
	y <- "speaking"
	out_name3 <- paste(out_prefix, "pacing_", out_suffix, ".pdf", sep="")	
	b.pacechart(x,y,d,out_name3, "label", singlepage)
	
	### FOURTH GRAPH: who do you engage with during a group discussion
	# Do this as a heatmap.
	d <- o$n_tt_turntaking
	
	# First get the ssi_temp_ids
	ssi_temp_ids <- unlist(lapply(strsplit(colnames(d)[2:(length(d))], "%"), function(x) x[1]))	
	colnames(d) <- c("ssi_id", ssi_temp_ids)
	# Check for any missing values and remove if there are any
	d <- merge(d, peeps, by.x=c("ssi_id"), all.x=TRUE)	
	b_rem <- d[is.na(d$team_id), c("ssi_id")]
	# Get the position in the temp of these
	pos <- match(b_rem, ssi_ids)
	b_temps <- people_column_names[pos]
	d <- d[,!(names(d) %in% b_temps)]
	d <- d[d$ssi_id %in% peeps[peeps$team_id == tm_id, c("ssi_id")], ]
	# Rename the columns using the values of the rows
	colnames(d)[2:(length(d$ssi_id)+1)] <- d$ssi_id

	mat <- data.matrix(d[,2:(length(d$ssi_id)+1)])
	max_val <- max(mat, na.rm=TRUE)
	cols <- rev(gray((0:max_val/max_val)))	
	out_name4 <- paste(out_prefix, "turntaking_", out_suffix, ".pdf", sep="")							
	
	if(singlepage) {
#		quartz(width=6,height=6	)
		par(mar= c(0,0,0,0))
	} else {
		quartz(width=9,height=3,type="pdf",file=out_name4,dpi=300)	
		par(mar= c(0,0,0,0))
		par(lwd = 1)
		par(xpd=TRUE)
		par(mfcol=c(1,1), family="Arial", mgp=c(3,1,0), bg="white")		
	}	

	plot.new()
	# draw the outside box
#	rect(0,0,1,1)	
	# draw the boxes
	boxsize <- .8/length(d$ssi_id)
	start.x <- .2
	start.y <- .2
	ego <- d$ssi_id[2]
	alter <- d$ssi_id[1]
	for(ego in d$ssi_id) {
		text(.1, start.y+boxsize/2, d[d$ssi_id == ego, c("label")], pos=1, cex=0.75)
		for(alter in d$ssi_id) {
			col.val <- d[d$ssi_id == ego, alter]
			rect(start.x, start.y, start.x+boxsize, start.y+boxsize, col=cols[col.val])
			if(ego == d$ssi_id[1]) {
				text(start.x+boxsize/2, .1, d[d$ssi_id == alter, c("label")], cex=0.75)				
			}
			start.x <- start.x + boxsize			
		}
		start.y <- start.y + boxsize	
		start.x <- .2	
	}	
	if(!singlepage) {
		dev.off()			
	}
	
	if(singlepage) {
		dev.off()
	}
}
#spch_profile <- d
b.temporal.boundaries <- function(spch_profile) {
	d <- spch_profile
	spch_agg <- aggregate(d[, c("total_speaking")], by=list(d$time_stamp), sum, na.rm=TRUE)
	colnames(spch_agg) <- c("time_stamp", "total_speaking_sum")
	start_time <- NA
	end_time <- NA
	for(i in 1:length(spch_agg$time_stamp)) {
		if(spch_agg[i, c("total_speaking_sum")] > 0 && is.na(start_time)) {
			start_time <- spch_agg[i, c("time_stamp")]
		}
		if(spch_agg[i, c("total_speaking_sum")] == 0 && !is.na(start_time) && is.na(end_time)) {
			end_time <- spch_agg[i, c("time_stamp")]					
		}	
	
		if(spch_agg[i, c("total_speaking_sum")] > 0 && !is.na(start_time) && !is.na(end_time)) {
			end_time <- spch_agg[i, c("time_stamp")]					
		}			
	}
	tmp.bound <- c(start_time, end_time)
	attr(tmp.bound, "tzone") <- "CDT"	
	return(tmp.bound)
}


##################################################################
# These functions will parse mock interview data and prepare for upload to the server
##################################################################

read.roster <- function(fileroot, report_name) {
	fname <- paste(fileroot, "roster.xlsx", sep="_")
	r <- read.xlsx2(fname, as.data.frame=TRUE, header=TRUE, sheetName="Sheet1", stringsAsFactors=FALSE, colClasses=c("numeric", "character", "character", "character", "character", "character", "numeric", "character", "numeric"))
	r$username <- r$indiv_id
	r$ssi_id <- toupper(r$ssi_id)
	r <- r[r$username != "" & r$report_name == report_name, ]
	return(r)
}

## This function takes four mirroring variables and creates aggregate and combined forms
create.mirroring <- function(o) {
	d <- o$t_audio_front_vol_mirroring1
	agg1 <- aggregate(d[,c("audio_front_vol_mirroring")], by=list(d$ego, d$alter), median, na.rm=TRUE)
	colnames(agg1) <- c("ego", "alter", "audio_front_vol_mirroring")

	d <- o$t_audio_back_vol_mirroring1
	agg2 <- aggregate(d[,c("audio_back_vol_mirroring")], by=list(d$ego, d$alter), median, na.rm=TRUE)
	colnames(agg2) <- c("ego", "alter", "audio_back_vol_mirroring")

	d <- o$t_audio_front_f0_mirroring1
	agg3 <- aggregate(d[,c("audio_front_f0_mirroring")], by=list(d$ego, d$alter), median, na.rm=TRUE)
	colnames(agg3) <- c("ego", "alter", "audio_front_f0_mirroring")

	d <- o$t_audio_back_f0_mirroring1
	agg4 <- aggregate(d[,c("audio_back_f0_mirroring")], by=list(d$ego, d$alter), median, na.rm=TRUE)
	colnames(agg4) <- c("ego", "alter", "audio_back_f0_mirroring")

	d <- o$t_BM_mirroring1
	agg5 <- aggregate(d[,c("bm_mirroring")], by=list(d$ego, d$alter), median, na.rm=TRUE)
	colnames(agg5) <- c("ego", "alter", "bm_mirroring")

	d <- o$t_posture_mirroring1
	agg6 <- aggregate(d[,c("posture_mirroring")], by=list(d$ego, d$alter), median, na.rm=TRUE)
	colnames(agg6) <- c("ego", "alter", "posture_mirroring")

	# merge these together
	mrg <- merge_recurse(list(agg1,agg2,agg3,agg4,agg5, agg6), by=c("ego", "alter"))
	mirroring <- mrg[mrg$ego != mrg$alter, ]
	mirroring$physical_mirroring <- rowMeans(mirroring[,c("bm_mirroring", "posture_mirroring")], na.rm=TRUE)
	mirroring$vocal_mirroring <- rowMeans(mirroring[,c("audio_front_vol_mirroring", "audio_back_vol_mirroring", "audio_front_f0_mirroring", "audio_back_f0_mirroring")], na.rm=TRUE)
	
	return(mirroring[,c("ego", "alter", "physical_mirroring", "vocal_mirroring")])
}

create.indiv.static <- function(o, roster, report, fileroot) {	
	if(report == "dyad") {
		# Get the turntaking information
		d <- data.frame(t(o$summary_sheet[,3:length(o$summary_sheet)]), stringsAsFactors=F)
		colnames(d) <- o$summary_sheet[,1]
		d$ssi_id <- colnames(o$summary_sheet[,3:length(o$summary_sheet)])
		# Add the turntaking to the indiv_static
		d1 <- merge(d, o$r_tt_turntaking1, by=c("ssi_id"))			
		# merge with the roster information
		indiv_static <- merge(d1, roster, by=c("ssi_id"))	
		
		indiv_static_vars <- c("speech_profile_speaking", "speech_profile_total_speaking", "avg_segment_speaking", "n_interruptions_successful", "n_interruptions_unsuccessful")
		count <- 1
		for(v in indiv_static_vars) {
			for(r in 1:length(indiv_static$username)) {
				c1 <- indiv_static[r, c("username", "report_name")]
				c2 <- v
				c3 <- round(indiv_static[r,v],2)
				d <- data.frame(c1,c2,c3, stringsAsFactors=FALSE)
				colnames(d) <- c("username", "reportdesc", "metric", "value")
				if(count == 1 ) {
					indiv_static_out <- d
				} else {
					indiv_static_out <- rbind(indiv_static_out, d)
				}
				count <- count + 1		
			}
		}
		write.table(indiv_static_out, file=paste("./web/",fileroot, "_indiv_static.txt", sep=""), sep="\t",  row.names=FALSE, quote=FALSE)
	}
}

create.indiv.dynamic <- function(o, roster, report, fileroot) {
	if(report == "dyad") {	
		# Now do the indiv_dynamic datafile
		d <- o$t_speech_profile1	
		tmp.bound <- b.temporal.boundaries(d)

		# Replace the longer speech profile dataset with this one
		d <- d[d$time_stamp > tmp.bound[1] & d$time_stamp < tmp.bound[2], ]

		# merge with the roster information
		indiv_dynamic <- merge(d, roster, by=c("ssi_id"))
		
		indiv_dynamic_vars <- c("speaking", "total_speaking")
		count <- 1
		for(v in indiv_dynamic_vars) {
			for(r in 1:length(indiv_dynamic$time_stamp)) {
				c1 <- indiv_dynamic[r, c("username", "report_name")]
				c2 <- indiv_dynamic[r, c("time_stamp")]
				c3 <- v
				c4 <- round(indiv_dynamic[r,v], 2)
				d <- data.frame(c1,c2,c3,c4, stringsAsFactors=FALSE)
				colnames(d) <- c("username", "reportdesc", "time_stamp", "metric", "value")
				if(count == 1) {
					indiv_dynamic_out <- d
				} else {
					indiv_dynamic_out <- rbind(indiv_dynamic_out, d)
				}
				count <- count + 1	
			}
		}
		write.table(indiv_dynamic_out, file=paste("./web/",fileroot, "_indiv_dynamic.txt", sep=""), sep="\t",  row.names=FALSE, quote=FALSE)		
	}
}
create.dyad.static <- function(o, roster, report, fileroot) {
	if(report == "dyad") {
		d <- o$n_tt_turntaking
		# First get the alter ssi_ids
		alter_ssi_id <- unlist(lapply(strsplit(colnames(d)[2:(length(d))], "%"), function(x) x[1]))	

		# Loop through egos and alters to build a dyadic dataset
		count <- 1
		for(ego in d$ssi_id) {
			for(alter in alter_ssi_id) {
				# Get the value
				turntaking <- d[d$ssi_id == ego, paste(alter, "turntaking", sep="%")]
				add.row <- data.frame(ego,alter,turntaking, stringsAsFactors=FALSE)
				if(count == 1) {
					o.d <- add.row
				} else {
					o.d <- rbind(o.d, add.row)
				}
				count <- count + 1
			}	
		}

		# Add the mirroring variables
		mirroring <- create.mirroring(o)
		o.d <- merge(o.d, mirroring[,c("ego", "alter", "physical_mirroring", "vocal_mirroring")], by=c("ego", "alter"))

		# merge with the roster information
		sub.roster <- roster[,c("ssi_id", "username", "report_name")]
		o.d.ego <- merge(o.d, sub.roster, by.x=c("ego"), by.y=c("ssi_id"))
		colnames(o.d.ego) <- c("ego_ssi_id", "alter_ssi_id", "turntaking", "physical_mirroring", "vocal_mirroring", "ego_username", "report_name")
		o.d.both <- merge(o.d.ego, sub.roster[,c("ssi_id", "username")], by.x=c("alter_ssi_id"), by.y=c("ssi_id"))
		colnames(o.d.both)[8] <- "alter_username"

		dyad_static <- o.d.both[,c("ego_ssi_id", "alter_ssi_id", "ego_username", "alter_username", "report_name", "turntaking", "physical_mirroring", "vocal_mirroring")]

		## Create the out files
		dyad_static_vars <- c("turntaking", "physical_mirroring", "vocal_mirroring")
		count <- 1
		for(v in dyad_static_vars) {
			for(r in 1:length(dyad_static$ego_username)) {
				c1 <- dyad_static[r,c("report_name", "ego_username")]
				c2 <- dyad_static[r,c("alter_username")]
				c3 <- v
				c4 <- dyad_static[r,v]
				d <- data.frame(c1,c2,c3,c4, stringsAsFactors=FALSE)
				colnames(d) <- c("reportdesc", "ego_username", "alter_username", "metric", "value")		
				if(count == 1) {
					dyad_static_out <- d
				} else {
					dyad_static_out <- rbind(dyad_static_out, d)
				}
				count <- count + 1			
			}
		}
		write.table(dyad_static_out, file=paste("./web/",fileroot, "_dyad_static.txt", sep=""), sep="\t",  row.names=FALSE, quote=FALSE)
	}
}

## This function is for use in the Rube Goldberg project. 
## This is not yet a generalized function; it is still in progress. 

read.rube.data <- function(fname, data_structure) {
	# Some computers will require the following option to be set. So far, I know my MBP requires this, but my Mac Mini does not. These options set the memory of the java virtual machine a bit higher than it might be otherwise. This will be important for very large datasets. Once we know more about differences across environments, we can standardize the call to do this or not.
	options( java.parameters = "-Xmx4g" )
	# options( java.parameters = "-Xmx2g" )
	require(xlsx)
	# First, read in the summary sheet to get the information that is loaded there needed to parse the other columns
	summary_sheet <- read.xlsx2(fname, as.data.frame=TRUE, header=TRUE, sheetName="summary", stringsAsFactors=FALSE, colClasses=c("character", "character", rep("numeric",20)))

	# create new column names and store the individual names in a different vector
	people_n <- length(summary_sheet)-2
	ssi_ids <- colnames(summary_sheet)[3:length(summary_sheet)]


	dat_lst <- list("summary_sheet"=summary_sheet)

	for(i in data_structure) {
		d <- read.xlsx2(fname, as.data.frame=TRUE, header=FALSE, sheetName=i$name, stringsAsFactors=FALSE, startRow=i$data_row_start, colIndex=1:length(i$cols), colClasses=i$column_in_types)
		colnames(d) <- i$cols

		### Convert the timestamp if this is a dynamic dataset ###
		if(i$temporal == "dynamic") {
			d$time_stamp <- as.POSIXct((as.numeric(d$time_stamp)-25569)*86400, tz="America/New_York", origin="1970-01-01") + 60*60*4
		}	

		assign(i$name, d)		

		## Go through a few specific formats to produce better-structured datasets

		if(i$level == "individual" && i$temporal == "dynamic") {
			count <- 1
			for( ego in ssi_ids ) {
					var_names <- paste(ego, "%", i$var_prefix, sep="")
				tmp <- d[, c("time_stamp", var_names)]
				colnames(tmp) <- c("time_stamp", i$var_prefix)
				tmp$ssi_id <- ego
				if( count == 1 ) {
					new.d <- tmp
				} else {
					new.d <- rbind(new.d, tmp)
				}
				count <- count + 1
			}
			assign(i$name, new.d)	
		}
		dat_lst[[i$name]] <- eval(as.name(i$name))
	}
	return(dat_lst)
}



