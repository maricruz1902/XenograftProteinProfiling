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

