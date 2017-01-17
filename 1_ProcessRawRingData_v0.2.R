GetName <- function(){
	# get the filename from the current working directory
	directory <- basename(getwd())

	# directory naming from MRR: "CHIPNAME_gaskGASETTYPE_DATE"
	# extracts and returns GASKETTYPE from directory name
	name <- unlist(strsplit(directory, split = "_"))
	name <- name[2]

	# define name as global variable for use in other functions
	name <<- gsub('gask','',name) # removes "gask" from name
}

AggData <- function(loc = 'plots') {
	# load relevant libraries
	library(readr)
    
        # get working directory to reset at end of function
	directory <- getwd()
	
	filename <- "../groupNames_XPP.csv"

	# get information of chip layout from github repository
	if (!file.exists("../groupNames_XPP.csv")){
	        url <- "https://raw.githubusercontent.com/jwade1221/XenograftProteinProfiling/master/groupNames_XPP.csv"
	        filename <- paste("../", basename(url), sep="")
	        download.file(url, filename)
	}
	
	# define recipe as global variable for use in other functions
	recipe <<- read_csv(filename, col_types = cols())
	targets <- recipe$Target

	# generate list of rings to analyze (gets all *.csv files)
	rings <- list.files(directory, pattern = ".csv", recursive = FALSE)
	removeFiles <- c("comments.csv")
	rings <- rings[!rings %in% removeFiles]

	# create empty data frame to store data
	df <- data.frame()
	
	# add data to data frame corresponding for each ring in rings
	for (i in rings) {
		ring <- as.vector(i)
		dat <- read_csv(ring, col_types = cols(), col_names = FALSE)
		time_shift <- dat[ ,1]
		shift <- dat[ ,2]
		ringStr <- strsplit(i, "\\.")[[1]]
		ringNum <- as.numeric(ringStr[1])
		groupNum <- (ringNum - 1) %/% 4 + 1
		ring <- rep(ringNum, nrow(dat))
		group <- rep(groupNum, nrow(dat))
		groupName <- as.character(recipe$Target[[groupNum]])
		groupName <- rep(groupName, nrow(dat))
		channel <- recipe$Channel[[groupNum]]
		channel <- rep(channel, nrow(dat))
		run <- rep(name, nrow(dat))
		tmp <- data.frame(ring, group, time_shift, shift, groupName, 
		        channel, run)
		df <- rbind(df, tmp)
	}
    
	# renames columns in df
	names(df) <- c("Ring", "Group", "Time", "Shift", "Target", "Channel",
	        "Experiment")
    
	# creates "plots" directory if one does not exist
	if (!file.exists(loc)){dir.create(loc)}
    
	# saves aggregated data with name_allRings.csv
	write_csv(df, paste(loc, '/', name, "_allRings.csv", sep=""))
	
	# returns working directory to top level
	setwd(directory)
}

SubtractControl <- function(loc = 'plots', ch, cntl){
	#load relevant libraries
	library(readr)
	library(dplyr)
	
        # get working directory to reset at end of function
	directory = getwd()
	
	# get ring data and filter by channel
	dat <- read_csv(paste(loc, "/", name, "_", "allRings.csv", sep = ''), 
		col_types = cols())
	dat <- filter(dat, Channel == ch)

	# get thermal control averages
	controls <- filter(dat, Target == cntl)
	ringList <- unique(controls$Ring)
	
	# gets times from first thermal control
	times <- filter(controls, Ring == ringList[1]) %>% select(Time)
	df.controls <- data.frame(times)
	
	# create dataframe with all controls
	for (i in ringList){
		ringShift <- filter(controls, Ring == i) %>% select(Shift)
		names(ringShift) <- paste('Ring', i, sep='')
		df.controls <- cbind(df.controls, ringShift)
	}
	
	# averages thermal controls
	cols <- ncol(df.controls)
	if (length(unique(controls$Ring)) != 1) {
	        df.controls$avgControls <- rowMeans(df.controls[,c(2:cols)])
	} else {
	        df.controls$avgControls <- df.controls[,c(2:cols)]
	}
	avgControls <- as.vector(df.controls$avgControls)

	#subtracts thermal controls from each ring
	ringNames <- unique(dat$Ring)
	for(i in ringNames){
		ringDat <- filter(dat, Ring == i) %>% select(Shift)
		ringTC <- ringDat - avgControls
		dat[dat$Ring == i, 4] <- ringTC
	}

	write_csv(dat, paste(loc,"/", name, "_", cntl, "Control", "_ch", ch, 
	        ".csv", sep = ''))   
}

PlotRingData <- function(cntl, ch, loc = 'plots'){
	# loads relevant libraries
	library(ggplot2)
	library(readr)
	library(RColorBrewer)

	# get working directory to reset at end of function
	directory <- getwd()
	
	# use thermally controlled data if desired
	if (cntl != "raw"){
		dat <- read_csv(paste(loc, "/", name, "_", cntl, "Control", 
		        "_ch", ch,".csv", sep=''), col_types = cols())
	} else if (cntl == "raw") {
		dat <- read_csv(paste(loc, "/", name, "_allRings.csv", sep=''), 
			col_types = cols())
	}
	
	#set colors for plot
	colorCount <- length(unique(dat$Target))
	getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)

	#configure plot and legend
	plots <- ggplot(dat, aes(Time, Shift, colour = factor(Target))) + 
		xlab("Time (min)") + 
		ylab(expression(paste("Relative Shift (",Delta,"pm)"))) +
		scale_colour_manual(values = getPalette, name = 'Target') +
		theme_bw() + theme(panel.grid = element_blank(), 
			axis.title.x=element_blank()) + 
		theme(legend.key = element_rect(colour = 'white',
			fill = 'white'), legend.key.size = unit(0.4, "cm"))
	
	if (cntl == "raw"){
	        plots <- plots + geom_point(size = 1) + facet_grid(.~ Channel)
	} else {plots <- plots + geom_point(size = 1)}

	#plot figure, uncomment to plot
	# plots

	#save plot, uncomment to save
	filename <- paste(name, "_", cntl, "Control", "_ch", ch, ".png", sep="")
	setwd(loc)
	ggsave(plots, file = filename, width = 8, height = 6)
	setwd(directory)
}

GetNetShifts <- function(cntl, ch, loc = 'plots', time1, time2){
	# load relevant libraries
	library(readr)
	library(dplyr)

	# get working directory to reset at end of function
	directory <- getwd()
	
	# use thermally controlled data if desired
	if (cntl != "raw"){
	        dat <- read_csv(paste(loc, "/", name, "_", cntl, "Control", 
	                "_ch", ch, ".csv", sep=""), col_types = cols())
	} else {
	        dat <- read_csv(paste(loc, "/", name, "_", "allRings.csv", 
	                sep=''), col_types = cols())
	}
	
	# generate list of rings and empty dataframe to store net shift data
	ringList <- unique(dat$Ring)
	dat.rings <- data.frame()
	
	# locations for each time is determined using which, min, and abs func
	for (i in ringList){
		dat.ring <- filter(dat, Ring == i)
		time1.loc <- which.min(abs(dat.ring$Time - time1))
		time1.val <- dat.ring$Shift[time1.loc]
		time2.loc <- which.min(abs(dat.ring$Time - time2))
		time2.val <- dat.ring$Shift[time2.loc]
		ring <- i
		group <- unique(dat.ring$Group)
		target <- unique(dat.ring$Target)
		experiment <- unique(dat.ring$Experiment)
		channel <- unique(dat.ring$Channel)
		tmp <- data.frame(i, group, target, time1.val, time2.val,
		        experiment, channel)
		dat.rings <- rbind(dat.rings, tmp)
	}
	
	# renames dat.rings columns
	names(dat.rings) <- c("Ring", "Group", "Target", "Shift.1", "Shift.2", 
                "Experiment", "Channel")
	
	# calculate nat shift and create new column in dataframe
	dat.rings$`Net Shift` <- dat.rings$Shift.1 - dat.rings$Shift.2

	# save net shift data
	setwd(loc)
	write_csv(dat.rings, paste(name, "_netShifts_ch", ch, ".csv", sep=""))
	setwd(directory)
}

PlotNetShifts <- function(cntl, ch, loc = 'plots'){
	# load relevant libraries
	library(ggplot2)
	library(readr)
	library(RColorBrewer)

	# get working directory to reset at end of function
	directory <- getwd()
	
	# get net shift data
	if (cntl != "raw"){
	        dat <- read_csv(paste(loc, "/", name, "_netShifts_ch", ch,
	                ".csv", sep=""), col_types = cols())
	} else {
	        dat <- read_csv(paste(loc,"/", name, "_netShifts_chU.csv", 
	                sep=""), col_types = cols())

	}

	# configure colors
	colorCount <- nrow(dat)
	getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)

	# configure plot and legend
	plots <- ggplot(dat, aes(Target, `Net Shift`, fill = factor(Ring))) +
		geom_bar(stat = 'identity', position = 'dodge') +
		scale_fill_manual(values = getPalette, name = 'Ring') +
		ylab(expression(paste("Net Shift ( ", Delta,"pm)"))) +
		theme_bw() + theme(panel.grid = element_blank(), 
		       axis.title.x=element_blank()) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		theme(legend.key = element_rect(colour = 'white',
			fill = 'white'), legend.key.size = unit(0.3, "cm"))

	# save plot, uncomment to save
	filename <- paste(name, "_NetShift_ch", ch, ".png", sep="")
	setwd(loc)
	ggsave(plots, file = filename, width = 10, height = 6)
	setwd(directory)
}

GetAvgShifts <- function(cntl, ch, loc = 'plots'){
	# load relevant libraries
	library(readr)

	# get working directory to reset at end of function
	directory <- getwd()
	
	# get net shift data
	if (cntl != "raw"){
	        dat <- read_csv(paste(loc, "/", name, "_netShifts_ch", ch,
	                ".csv", sep=""), col_types = cols())
	} else {
	        dat <- read_csv(paste(loc,"/", name, "_netShifts_chU.csv", 
	                sep=""), col_types = cols())

	}
	
	# generate list of targets and empty dataframe
	targets <- unique(dat$Target)
	df <- data.frame()

	# average for each target and calculate sd, se, and cv
	for(i in targets){
		dat.group <- dat[dat$Target == i,]
		dat.shifts <- dat.group$`Net Shift`
		avgShift <- mean(dat.shifts)
		sdShift <- sd(dat.shifts)
		seShift <- sd(dat.shifts) / sqrt(length(dat.shifts))
		cvShift <- sd(dat.shifts) / avgShift * 100
		target <- unique(dat.group$Target)
		experiment <- unique(dat.group$Experiment)
		channel <- unique(dat.group$Channel)
		tmp <- data.frame(target, avgShift, sdShift, seShift, cvShift,
		        experiment, channel)
		df <- rbind(df, tmp)
	}
	
	# renamte dataframe columns
	names(df) <- c("Target", "Average Shift", "SD", "SE", "CV", 
                "Experiment", "Channel")

	# save data
	setwd(loc)
        write_csv(df, paste(name, "_avgShifts_ch", ch, ".csv", sep=""))
	setwd(directory)
}

PlotAvgShifts <- function(cntl, ch, loc = 'plots'){
	# load relevant libraries
	library(ggplot2)
	library(readr)
	library(RColorBrewer)

	# get working directory to reset at end of function
	directory <- getwd()
	
	# get avg shifts data
	if (cntl != "raw"){
	        dat <- read_csv(paste(loc, "/", name, "_avgShifts_ch", ch,
	                ".csv", sep=""), col_types = cols())
	} else {
	        dat <- read_csv(paste(loc,"/", name, "_avgShifts_chU.csv", 
	                sep=""), col_types = cols())
	}
	# configure colors
	colorCount <- nrow(dat)
	getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)

	# set error bars
	limits <- aes(ymax = `Average Shift` + SD, ymin = `Average Shift` - SD)

	# configure plot and legend
	plots <- ggplot(dat, 
		aes(x = Target, y = `Average Shift`, fill = Target)) +
		geom_errorbar(limits, width = 0.3) +
		scale_fill_manual(values = getPalette) +
		ylab(expression(paste("Net Shift ( ",Delta,"pm)"))) +
		theme_bw() + theme(panel.grid = element_blank(), 
			axis.title.x=element_blank()) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		guides(fill = FALSE)
	
	if (cntl == "raw"){
	        plots <- plots + geom_bar(stat = "identity") + 
	                facet_grid(.~ Channel)
	} else {plots <- plots + geom_bar(stat = "identity")}

	#save figure, uncomment to save
	filename <- paste(name, "_AvgShift_ch", ch, ".png", sep="")
	setwd(loc)
	ggsave(plots, file = filename, width = 8, height = 6)
	setwd(directory)
}

FindBadClusters <- function(loc = 'plots', cvCutoff = 20){
	# load relevant libraries
	library(readr)
	library(dplyr)
	
	# get working directory to reset at end of function
	directory <- getwd()
	
	# get avg shift data
	dat <- read_csv(paste(loc, "/", name, "_", "avgShifts.csv", sep=''), 
		    col_types = cols())
	
	# generate group list and vector to hold bad groups
	groupList <- dat$Target
	badGroups <- vector("numeric")
	counter <- 0
	
	# if group is outside of acceptable range & if true, add to badGroups
	for (i in 2:length(groupList)){
		groupCV <- filter(dat, Target == groupList[i]) %>% 
			select(CV) %>% unlist()
		if (groupCV > cvCutoff){
			counter <- counter + 1
			badGroups[counter] <- groupList[i]
		}
	}
	
	# convert badGroups to dataframe and save data if bad groups exist
	badGroups <- as.data.frame(badGroups)
	if(length(badGroups) > 0){
		write_csv(badGroups, paste(loc, "/", name, "_", 
		"badGroups.csv", sep=""))
	}
	setwd(directory)
}

PlotBadClusters <- function(loc = 'plots'){
	# load relevant libraries
	library(readr)
	library(dplyr)
	library(RColorBrewer)
	library(ggplot2)

	# get working directory to reset at end of function
	directory <- getwd()
	
	# get bad ring data
	badRings <- read_csv(paste(loc, "/", name, "_", "badGroups.csv", 
		sep = ""), col_types = cols())
	badRings <- unlist(badRings)
	
	# only execute if bad rings exist
	if (length(badRings) > 0){
		dat <- read_csv(paste(loc, "/", name, "_", "netShifts.csv", 
		sep=""), col_types = cols())
		dat.bad <- filter(dat, groupName %in% badRings)
		
		#configure colors
		colorCount <- nrow(dat.bad)
		getPalette <- colorRampPalette(brewer.pal(8, 
			"Paired"))(colorCount)

		#configure plot
		plots <- ggplot(dat.bad) +
			geom_bar(aes(groupName, netShifts, fill = factor(ring)), 
				stat = 'identity', position = 'dodge') +
			scale_fill_manual(values = getPalette, name = 'Rings') +
			ylab(expression(paste("Net Shift ( ", Delta,"pm)"))) +
			theme_bw() + theme(panel.grid = element_blank(), 
				axis.title.x=element_blank()) + 
			theme(axis.text.x = element_text(angle = 45, 
				hjust = 1)) +
			theme(legend.key = element_rect(colour = 'white',
				fill = 'white'), 
				legend.key.size = unit(0.3, "cm"))

		#save plot, uncomment to save
		filename <- paste(name, "BadRings.png", sep="_")
		setwd(loc)
		ggsave(plots, file = filename, width = 10, height = 6)
		write_csv(dat.bad, paste(name, 'badRings.csv', sep="_"))
	}
	setwd(directory)
}

Grubbs <- function (x, type = 10, opposite = FALSE, two.sided = FALSE){
    # this is slightly modified from the "outliers" library
    # added value output
    library(outliers)
    
    if (sum(c(10, 11, 20) == type) == 0) 
        stop("Incorrect type")
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    if (type == 11) {
        g <- (x[n] - x[1])/sd(x)
        u <- var(x[2:(n - 1)])/var(x) * (n - 3)/(n - 1)
        pval = 1 - pgrubbs(g, n, type = 11)
        method <- "Grubbs test for two opposite outliers"
        alt = paste(x[1], "and", x[n], "are outliers")
    }
    else if (type == 10) {
        if (xor(((x[n] - mean(x)) < (mean(x) - x[1])), opposite)) {
            alt = paste("lowest value", x[1], "is an outlier")
            o <- x[1]
            d <- x[2:n]
        }
        else {
            alt = paste("highest value", x[n], "is an outlier")
            o <- x[n]
            d <- x[1:(n - 1)]
        }
        g <- abs(o - mean(x))/sd(x)
        u <- var(d)/var(x) * (n - 2)/(n - 1)
        pval <- 1 - pgrubbs(g, n, type = 10)
        method <- "Grubbs test for one outlier"
    }
    else {
        if (xor(((x[n] - mean(x)) < (mean(x) - x[1])), opposite)) {
            alt = paste("lowest values", x[1], ",", x[2], "are outliers")
            u <- var(x[3:n])/var(x) * (n - 3)/(n - 1)
        }
        else {
            alt = paste("highest values", x[n - 1], ",", x[n], 
                        "are outliers")
            u <- var(x[1:(n - 2)])/var(x) * (n - 3)/(n - 1)
        }
        g <- NULL
        pval <- pgrubbs(u, n, type = 20)
        method <- "Grubbs test for two outliers"
    }
    if (two.sided) {
        pval <- 2 * pval
        if (pval > 1) 
            pval <- 2 - pval
    }
    # added "val = o"
    RVAL <- list(statistic = c(G = g, U = u), alternative = alt, 
                 p.value = pval, method = method, data.name = DNAME, val = o)
    class(RVAL) <- "htest"
    return(RVAL)
    return(o)
}

IdentifyOutliers <- function(loc = 'plots'){
	# load relevant libraries
	library(readr)
	library(dplyr)

	# get working directory to reset at end of function
	directory <- getwd()
	
	# get bad groups list
	badRings <- read_csv(paste(loc, "/", name, "_", "badGroups.csv", sep=""),
			 col_types = cols())
	badRings <- unlist(badRings)

	# identify outliers with Grubbs test and save list
	if (length(badRings) > 0){
		dat <- read_csv(paste(loc, "/", name, "_", 
			"badRings.csv", sep=""), col_types = cols())
		badClusters <- unique(dat$group)
		counter <- 1
		tossedRings <- vector("numeric")

		for(i in 1:length(badClusters)){
			dat.cluster <- filter(dat, group == badClusters[i]) %>% 
			select(netShifts) %>% unlist()
			outliers <- Grubbs(dat.cluster)
			ring.outlier <- filter(dat, netShifts == 
				       outliers$val) %>% select(ring)
			tossedRings[counter] <- ring.outlier
			counter <- counter + 1
		}
		tossedRings <- as.data.frame(tossedRings)
		write_csv(tossedRings, "plots/tossedRings.csv")
	}
    
	setwd(directory)
}

RemoveRings <- function(loc = 'corrected'){
	# load relevant libraries
	library(readr)

	# ensure that this function does not fun in plots directory
	if(loc == 'plots'){break}

	# get working directory to reset at end of function
	directory <- getwd()
	
	# get tossed rings
	RemoveRings <- as.numeric(read_csv("plots/tossedRings.csv"), 
	        col_types = cols())
	
	# repeat AggData without outlier rings
	groupNames <- recipe$groupNames
	holder <- recipe[,c(1,2)]
	rings <- list.files(directory)
	RemoveRings <- paste(RemoveRings, '.csv', sep='')
	remove <- c(RemoveRings, 'plots', 'comments.csv')
	rings <- rings[!rings %in% remove]

	# create empty data frame
	df <- data.frame()
	
	# add data to data frame corresponding for each ring in rings
	for (i in rings) {
		ring <- as.vector(i)
		dat <- read_csv(ring, col_types = cols(), col_names = FALSE)
		time_shift <- dat[ ,1]
		shift <- dat[ ,2]
		ringStr <- strsplit(i, "\\.")[[1]]
		ringNum <- as.numeric(ringStr[1])
		groupNum <- (ringNum - 1) %/% 4 + 1
		ring <- rep(ringNum, nrow(dat))
		group <- rep(groupNum, nrow(dat))
		groupName <- as.character(holder$groupNames[[groupNum]])
		groupName <- rep(groupName, nrow(dat))
		tmp <- data.frame(ring, group, time_shift, shift, groupName)
		df <- rbind(df, tmp)
	}
    
	# renames columns in df
	names(df) <- c("ring", "group", "time", "shift", "groupName")
	
	# save data
	dir.create('corrected')
	setwd('corrected')
	write_csv(df, paste(name, "allRings.csv", sep="_"))
	setwd(directory)
}

RawGo <- function(location = 'plots', getName = TRUE, control, channel){
	# load relevant libraries
	library(ggplot2)
	library(readr)
	library(RColorBrewer)
	library(dplyr)
	library(outliers)

	# get working directory to reset at end of function
	directory <- getwd()

	name <- GetName()

	# only run RemoveRings if location is not plots
	if(location != 'plots'){
		RemoveRings(loc = location)
	} else{
		AggData()
	}
	
	if (control != "raw"){
	        SubtractControl(cntl = control, ch = channel, loc = location)
	}
	PlotRingData(cntl = control, ch = channel, loc = location)
	PlotRingData(cntl = control, ch = channel, loc = location)
	GetNetShifts(cntl = control, ch = channel, 
	        loc = location, time1 = 53, time2 = 39)
	PlotNetShifts(cntl = control, ch = channel, loc = location)
	GetAvgShifts(cntl = control, ch = channel, loc = location)
	PlotAvgShifts(cntl = control, ch = channel, loc = location)
	#FindBadClusters(cntl = control, ch = channel, loc = location)
	#PlotBadClusters(cntl = control, ch = channel, loc = location)
	#IdentifyOutliers(cntl = control, ch = channel, loc = location)

	setwd(directory)
}

AnalyzeData <- function(){
        RawGo(control = "raw", channel = "U")
        RawGo(control = "thermal", channel = 1)
        RawGo(control = "thermal", channel = 2)
}


AnalyzeFullDataSet <- function(){
	directory <- getwd()
	listData <- list.dirs(recursive = FALSE)
	for (i in 1:length(listData)){
		tempDir <- listData[i]
		setwd(tempDir)
		RawGo(control = "raw", channel = "U")
		RawGo(control = "thermal", channel = 1)
		RawGo(control = "thermal", channel = 2)
		#if (file.exists('plots/tossedRings.csv')){
		#    RawGo(location = 'corrected')
		#}
		setwd(directory)
	}
}

AnalyzeFullDataSetAlt <- function(){
        directory <- getwd()
        listData <- list.dirs(recursive = FALSE)
        for (i in 1:length(listData)){
                tempDir <- listData[i]
                setwd(tempDir)
                RawGo()
                if (file.exists('plots/tossedRings.csv')){
                        RawGo(location = 'corrected')
                }
                if (file.exists('corrected/tossedRings.csv')){
                        RawGo(location = 'corrected2')
                }
                setwd(directory)
        }
}