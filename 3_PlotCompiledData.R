GetData <- function( filename = "Combined_Averaged.csv"){
	library(readr)
	# uncomment to use github repository as source of data
	# url <- "https://raw.githubusercontent.com/jwade1221/XenograftProteinProfiling/master/Combined_Averaged.csv"
	# filename <- basename(url)
	# download.file(url,destfile=filename)
	
        filename <- 'Combined_Averaged.csv'
        dat <<- read_csv(filename, col_types = cols())
}

PlotEachTarget <- function(i, filename = "Combined_Averaged.csv"){
	# load relevant libraries
	library(readr)
	library(ggplot2)
	library(dplyr)
	library(RColorBrewer)

	# import data
	dat <- read_csv(filename, col_types = cols())

	# create lists for data filter and filter data
	targetList <- unique(dat$Target)
	dat.Target <- filter(dat, Target == targetList[i])

	# set error bars
	limits <- aes(ymax = `Average Shift` + SD, ymin = `Average Shift` - SD)

	# configure colors
	colorCount <- length(unique(dat$`Cell Line`)) + 
		length(unique(dat$`Time Point`))
	getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)

	# configure plot and legend
	plots <- ggplot(data = dat.Target,
		aes(x = Treatment, y = `Average Shift`, 
			fill = interaction(factor(`Time Point`), 
			factor(`Cell Line`)), group = 
			interaction(factor(`Time Point`), 
			factor(`Cell Line`)))) +
		geom_bar(stat = "identity", position = "dodge") +
		ylab("Net Shift (pm)") +
		geom_errorbar(limits,position = position_dodge(width = 0.9), 
			width = 0.4) +
		theme_bw() + theme(panel.grid = element_blank(), 
			axis.title.x=element_blank()) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1))+
		scale_fill_manual(values = getPalette, name = 
			"Cell Line and Treatment", labels = 
			c("GBM-6 1 h Treatment", "GBM-6 24 h Treatment",
			"GBM-26 1 h Treatment", "GBM-26 24 h Treatment")) + 
		theme(legend.key = element_rect(colour = 'white',
			fill = 'white'), legend.key.size = unit(0.4, "cm"))+
		ggtitle(paste('Target: ', targetList[i], sep=''))

	# plot figure, uncomment to plot
	# plots
		
	# save figure, uncomment to save
	filename = paste(targetList[i], 'png', sep='.')
	ggsave(plots, file = filename, width = 8, height = 6)
}

PlotEachTimepoint <- function(i, j, filename = "Combined_Averaged.csv"){
	library(readr)
	library(ggplot2)
	library(RColorBrewer)

	#import data
	dat <- read_csv(filename, col_types = cols())

	#create lists for data islation
	timepointList <- unique(dat$`Time Point`)
	cellLineList <- unique(dat$`Cell Line`)
	dat.Timepoint <- dat[dat$`Time Point` == timepointList[i] & 
			     dat$`Cell Line` == cellLineList[j],]

	#set error bars
	limits <- aes(ymax = `Average Shift` + SD, ymin = `Average Shift` - SD)

	#set colors for plot
	colorCount <- length(unique(dat$Target))
	getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)

	#configure plot, colors, and legend
	plots <- ggplot(data = dat.Timepoint,
		aes(x = Treatment, y = `Average Shift`,
			fill = factor(Target))) +
		geom_bar(stat = "identity", 
			position = position_dodge(width = 0.9)) +
		ylab("Net Shift (pm)") +
		scale_fill_manual(values = getPalette, name = 'Target') +
		geom_errorbar(limits, width = 0.4, size = 0.3, colour = 'grey40',
			position = position_dodge(width = 0.9)) +
		theme_bw() + 
		theme(panel.grid = element_blank(), 
			axis.title.x=element_blank()) + 
	        theme(axis.text.x = element_text(angle = 45, hjust = 1),
	              panel.background = element_blank()) +
		theme(legend.key = element_rect(colour = 'white', 
			fill = 'white'), legend.key.size = unit(0.4, "cm")) +
		ggtitle(paste('Cell Line: ', cellLineList[j],
			' Time Point: ', timepointList[i], ' h', sep=''))

	#plot figure, uncomment to plot
	#plots

	#save figure, uncomment to save
	filename = paste(timepointList[i],cellLineList[j], 'png', sep='.')
	ggsave(plots, file = filename, width = 8, height = 6)
}

PlotEachCellLine <- function(i, filename = "Combined_Averaged.csv"){
	library(readr)
	library(ggplot2)
	library(RColorBrewer)

	# import data
	dat <- read_csv(filename, col_types = cols())

	# create lists for data islation
	cellLineList <- unique(dat$`Cell Line`)
	dat.Timepoint <- dat[dat$`Cell Line` == cellLineList[i],]

	# set error bars
	limits <- aes(ymax = `Average Shift` + SD, ymin = `Average Shift` - SD)

	# set colors for plot
	colorCount <- length(unique(dat$Target)) * length(unique(dat$`Time Point`))
	getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)

	# configure plot, colors, and legend
	plots <- ggplot(data = dat.Timepoint,
		aes(x = Treatment, y = `Average Shift`,
			fill = interaction(factor(Target),
			factor(`Time Point`)))) +
		geom_bar(stat = "identity", 
			position = position_dodge(width = 0.9)) +
		ylab("Net Shift (pm)") +
		scale_fill_manual(values = getPalette, name = 'Target') +
		geom_errorbar(limits, width = 0.4, size = 0.3, 
			colour = 'grey40', position = 
			position_dodge(width = 0.9)) +
		theme_bw() + theme(panel.grid = element_blank(),
			axis.title.x=element_blank()) + 
	        theme(axis.text.x = element_text(angle = 45, hjust = 1),
	              panel.background = element_blank()) +
		theme(legend.key = element_rect(colour = 'white',
			fill = 'white'), legend.key.size = unit(0.4, "cm")) +
		ggtitle(paste('Cell Line: ', cellLineList[i], sep=''))

	#plot figure, uncomment to plot
	#plots

	#save figure, uncomment to save
	filename = paste(cellLineList[i], 'png', sep='.')
	ggsave(plots, file = filename, width = 12, height = 6)
}

HeatmapCellTime <- function(i, filename = "Combined_Averaged.csv"){
	library(readr)
	library(ggplot2)
	library(dplyr)

	# get data
	dat <- read_csv(filename, col_types = cols())

	# create necessary lists and matrix
	cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
	timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
	targetList <- select(dat, `Target`) %>% unlist() %>% unique()
	treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
	cLL <- rep(1:length(cellLineList), length(timePointList))
	tPL <- rep(1:length(timePointList), each = length(cellLineList))
	m <- as.matrix(cbind(tPL, cLL))

	cLL.i <- cellLineList[m[i,2]]
	tPL.i <- timePointList[m[i,1]]
	dat.1 <- filter(dat, `Cell Line` == cLL.i & `Time Point` == tPL.i)
	df <- data.frame()

	for (i in 1:length(targetList)){
		dat.2 <- filter(dat.1, `Target` == targetList[i])
		targetShift <-  select(dat.2, `Average Shift`) %>% unlist()
		treatment <- select(dat.2, Treatment) %>% unlist() 
		meanShift <- mean(targetShift)
		target <- rep(targetList[i], length(treatmentList))
		holder <- numeric()
		for (i in 1:length(targetShift)){
			holder[i] <- (targetShift[i] - meanShift) / meanShift
		}
		tmp <- data.frame(holder, target, treatment)
		df <- rbind(df, tmp)
	}
	treatmentList <- select(dat.1, `Treatment`) %>% unique()
	names(df) <- c('Change', 'Target', 'Treatment')

	#configure heatmap
	plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
		geom_tile(stat = "identity") +
		scale_fill_gradient2(low = "blue", mid = "white", 
		high = "red", midpoint = 0, space = "rgb", 
		na.value = "grey50", guide = "colourbar") +
		ggtitle(paste('Cell Line: GBM-', cLL.i, ' Time Point: ', tPL.i,
		' h', sep=''))

	#plot figure, uncomment to plot
	#plots

	#write csv file and save figure, uncomment to save
	filename <- paste("heatmap_GBM", cLL.i, '_', tPL.i, 'h', sep='')
	write_csv(df, paste(filename, '.csv', sep=''))
	ggsave(plots, file = paste(filename, '.png', sep=''), 
	width = 10, height = 6)
}

HeatmapCellLine <- function(i, filename = "Combined_Averaged.csv"){
	library(readr)
	library(ggplot2)
	library(dplyr)

	#get data
	dat <- read_csv(filename, col_types = cols())

	#create necessary lists and matrix
	cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
	timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
	targetList <- select(dat, `Target`) %>% unlist() %>% unique()
	treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
	cLL <- rep(1:length(cellLineList), length(timePointList))
	tPL <- rep(1:length(timePointList), each = length(cellLineList))
	m <- as.matrix(cbind(tPL, cLL))

	cLL.i <- cellLineList[m[i,2]]
	dat.1 <- filter(dat, `Cell Line` == cLL.i)
	df <- data.frame()

	for (i in 1:length(targetList)){
		dat.2 <- filter(dat.1, `Target` == targetList[i])
		targetShift <-  select(dat.2, `Average Shift`) %>% unlist
		treatment.a <- select(dat.2, Treatment) %>% unlist() 
		treatment.b <- select(dat.2, `Time Point`) %>% unlist()
		treatment <- paste(treatment.a, '_', treatment.b, 'h', sep='') 
		meanShift <- mean(targetShift)
		target <- rep(targetList[i], length(treatmentList))
		holder <- numeric()
		for (i in 1:length(targetShift)){
			holder[i] <- (targetShift[i] - meanShift) / meanShift
		}
		tmp <- data.frame(holder, target, treatment)
		df <- rbind(df, tmp)
	}
	treatmentList <- select(dat.1, `Treatment`) %>% unique()
	names(df) <- c('Change', 'Target', 'Treatment')

	#configure heatmap
	plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
		geom_tile(stat = "identity") + 
	        theme(axis.text.x = element_text(angle = 45, hjust = 1),
	              panel.background = element_blank()) +
		scale_fill_gradient2(low = "blue", mid = "white", 
			high = "red", midpoint = 0, space = "rgb", 
			na.value = "grey50", guide = "colourbar") +
		ggtitle(paste('Cell Line: GBM-', cLL.i, sep=''))

	#plot figure, uncomment to plot
	#plots

	#write csv file and save figure, uncomment to save
	filename <- paste("heatmap_GBM", cLL.i, sep='')
	write_csv(df, paste(filename, '.csv', sep=''))
	ggsave(plots, file = paste(filename, '.png', sep=''), 
		width = 10, height = 6)
	}

HeatmapTimePoint <- function(i, filename = "Combined_Averaged.csv"){
	library(readr)
	library(ggplot2)
	library(dplyr)

	#get data
	dat <- read_csv(filename, col_types = cols())

	#create necessary lists and matrix
	cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
	timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
	targetList <- select(dat, `Target`) %>% unlist() %>% unique()
	treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
	cLL <- rep(1:length(cellLineList), length(timePointList))
	tPL <- rep(1:length(timePointList), each = length(cellLineList))
	m <- as.matrix(cbind(tPL, cLL))

	tPL.i <- timePointList[i]
	dat.1 <- filter(dat, `Time Point` == tPL.i)
	df <- data.frame()

	for (i in 1:length(targetList)){
		dat.2 <- filter(dat.1, `Target` == targetList[i])
		targetShift <-  select(dat.2, `Average Shift`) %>% unlist
		treatment.a <- select(dat.2, Treatment) %>% unlist() 
		treatment.b <- select(dat.2, `Cell Line`) %>% unlist()
		treatment <- paste(treatment.a, '_GBM-', treatment.b, sep='') 
		meanShift <- mean(targetShift)
		target <- rep(targetList[i], length(treatmentList))
		holder <- numeric()
		for (i in 1:length(targetShift)){
			holder[i] <- (targetShift[i] - meanShift) / meanShift
		}
		tmp <- data.frame(holder, target, treatment)
		df <- rbind(df, tmp)
	}
	treatmentList <- select(dat.1, `Treatment`) %>% unique()
	names(df) <- c('Change', 'Target', 'Treatment')

	#configure heatmap
	plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
		geom_tile(stat = "identity") + 
	        theme(axis.text.x = element_text(angle = 45, hjust = 1),
	              panel.background = element_blank()) +
		scale_fill_gradient2(low = "blue", mid = "white", 
			high = "red", midpoint = 0, space = "rgb", 
			na.value = "grey50", guide = "colourbar") +
		ggtitle(paste('Time Point: ', tPL.i, ' h', sep=''))

	# plot figure, uncomment to plot
	# plots

	# write csv file and save figure, uncomment to save
	filename <- paste("heatmap_", tPL.i, 'h', sep='')
	write_csv(df, paste(filename, '.csv', sep=''))
	ggsave(plots, file = paste(filename, '.png', sep=''), 
		width = 10, height = 6)
}

HeatmapAll <- function(filename = "Combined_Averaged.csv"){
	library(readr)
	library(ggplot2)
	library(dplyr)

	#get data
	dat <- read_csv(filename, col_types = cols())

	#create necessary lists and matrix
	cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
	timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
	targetList <- select(dat, `Target`) %>% unlist() %>% unique()
	treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
	
	dat.1 <- dat
	df <- data.frame()

	for (i in 1:length(targetList)){
		dat.2 <- filter(dat.1, `Target` == targetList[i])
		targetShift <-  select(dat.2, `Average Shift`) %>% unlist
		treatment.a <- select(dat.2, Treatment) %>% unlist() 
		treatment.b <- select(dat.2, `Time Point`) %>% unlist()
		treatment.c <- select(dat.2, `Cell Line`) %>% unlist()
		treatment <- paste(treatment.a, '_', treatment.b, 'h_GBM', 
			treatment.c, sep='') 
		meanShift <- mean(targetShift)
		target <- rep(targetList[i], length(treatmentList))
		holder <- numeric()
		for (i in 1:length(targetShift)){
			holder[i] <- (targetShift[i] - meanShift) / meanShift
		}
		tmp <- data.frame(holder, target, treatment)
		df <- rbind(df, tmp)
	}
	treatmentList <- select(dat.1, `Treatment`) %>% unique()
	names(df) <- c('Change', 'Target', 'Treatment')

	#configure heatmap
	plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
		geom_tile(stat = "identity") + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1),
		      panel.background = element_blank()) +
		scale_fill_gradient2(low = "blue", mid = "white", 
			high = "red", midpoint = 0, space = "rgb", 
			na.value = "grey50", guide = "colourbar")

	# plot figure, uncomment to plot
	# plots

	# write csv file and save figure, uncomment to save
	filename <- "heatmap_all"
	write_csv(df, paste(filename, '.csv', sep=''))
	ggsave(plots, file = paste(filename, '.png', sep=''),
		width = 10, height = 6)
}

HeatmapCellTimeAlt <- function(i, filename = "Combined_Averaged.csv", 
        removeTxt = c('(+)-Serum', '(-)-Serum')){
        library(readr)
        library(ggplot2)
        library(dplyr)
        
        # get data
        dat <- read_csv(filename, col_types = cols())
        dat <- filter(dat, !Treatment %in% removeTxt)
        
        # create necessary lists and matrix
        cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
        timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
        targetList <- select(dat, `Target`) %>% unlist() %>% unique()
        treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
        cLL <- rep(1:length(cellLineList), length(timePointList))
        tPL <- rep(1:length(timePointList), each = length(cellLineList))
        m <- as.matrix(cbind(tPL, cLL))
        
        cLL.i <- cellLineList[m[i,2]]
        tPL.i <- timePointList[m[i,1]]
        dat.1 <- filter(dat, `Cell Line` == cLL.i & `Time Point` == tPL.i)
        df <- data.frame()
        
        for (i in 1:length(targetList)){
                dat.2 <- filter(dat.1, `Target` == targetList[i])
                targetShift <-  select(dat.2, `Average Shift`) %>% unlist()
                treatment <- select(dat.2, Treatment) %>% unlist() 
                meanShift <- mean(targetShift)
                target <- rep(targetList[i], length(treatmentList))
                holder <- numeric()
                for (i in 1:length(targetShift)){
                        holder[i] <- (targetShift[i] - meanShift) / meanShift
                }
                tmp <- data.frame(holder, target, treatment)
                df <- rbind(df, tmp)
        }
        treatmentList <- select(dat.1, `Treatment`) %>% unique()
        names(df) <- c('Change', 'Target', 'Treatment')
        
        #configure heatmap
        plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
                geom_tile(stat = "identity") +
                scale_fill_gradient2(low = "blue", mid = "white", 
                                     high = "red", midpoint = 0, space = "rgb", 
                                     na.value = "grey50", guide = "colourbar") +
                ggtitle(paste('Cell Line: GBM-', cLL.i, ' Time Point: ', tPL.i,
                              ' h', sep=''))
        
        #plot figure, uncomment to plot
        #plots
        
        #write csv file and save figure, uncomment to save
        filename <- paste("heatmap_GBM", cLL.i, '_', tPL.i, 'h_alt', sep='')
        write_csv(df, paste(filename, '.csv', sep=''))
        ggsave(plots, file = paste(filename, '.png', sep=''), 
               width = 10, height = 6)
}

HeatmapTimePointAlt <- function(i, filename = "Combined_Averaged.csv",
        removeTxt = c('(+)-Serum', '(-)-Serum')){
        library(readr)
        library(ggplot2)
        library(dplyr)
        
        #get data
        dat <- read_csv(filename, col_types = cols())
        dat <- filter(dat, !Treatment %in% removeTxt)
        
        #create necessary lists and matrix
        cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
        timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
        targetList <- select(dat, `Target`) %>% unlist() %>% unique()
        treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
        cLL <- rep(1:length(cellLineList), length(timePointList))
        tPL <- rep(1:length(timePointList), each = length(cellLineList))
        m <- as.matrix(cbind(tPL, cLL))
        
        tPL.i <- timePointList[i]
        dat.1 <- filter(dat, `Time Point` == tPL.i)
        df <- data.frame()
        
        for (i in 1:length(targetList)){
                dat.2 <- filter(dat.1, `Target` == targetList[i])
                targetShift <-  select(dat.2, `Average Shift`) %>% unlist
                treatment.a <- select(dat.2, Treatment) %>% unlist() 
                treatment.b <- select(dat.2, `Cell Line`) %>% unlist()
                treatment <- paste(treatment.a, '_GBM-', treatment.b, sep='') 
                meanShift <- mean(targetShift)
                target <- rep(targetList[i], length(treatmentList))
                holder <- numeric()
                for (i in 1:length(targetShift)){
                        holder[i] <- (targetShift[i] - meanShift) / meanShift
                }
                tmp <- data.frame(holder, target, treatment)
                df <- rbind(df, tmp)
        }
        treatmentList <- select(dat.1, `Treatment`) %>% unique()
        names(df) <- c('Change', 'Target', 'Treatment')
        
        #configure heatmap
        plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
                geom_tile(stat = "identity") + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      panel.background = element_blank()) +
                scale_fill_gradient2(low = "blue", mid = "white", 
                                     high = "red", midpoint = 0, space = "rgb", 
                                     na.value = "grey50", guide = "colourbar") +
                ggtitle(paste('Time Point: ', tPL.i, ' h', sep=''))
        
        #plot figure, uncomment to plot
        #plots
        
        #write csv file and save figure, uncomment to save
        filename <- paste("heatmap_", tPL.i, 'h_alt', sep='')
        write_csv(df, paste(filename, '.csv', sep=''))
        ggsave(plots, file = paste(filename, '.png', sep=''), 
               width = 10, height = 6)
}

HeatmapCellLineAlt <- function(i, filename = "Combined_Averaged.csv",
        removeTxt = c('(+)-Serum', '(-)-Serum')){
        library(readr)
        library(ggplot2)
        library(dplyr)
        
        #get data
        dat <- read_csv(filename, col_types = cols())
        dat <- filter(dat, !Treatment %in% removeTxt)
        
        #create necessary lists and matrix
        cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
        timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
        targetList <- select(dat, `Target`) %>% unlist() %>% unique()
        treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
        cLL <- rep(1:length(cellLineList), length(timePointList))
        tPL <- rep(1:length(timePointList), each = length(cellLineList))
        m <- as.matrix(cbind(tPL, cLL))
        
        cLL.i <- cellLineList[m[i,2]]
        dat.1 <- filter(dat, `Cell Line` == cLL.i)
        df <- data.frame()
        
        for (i in 1:length(targetList)){
                dat.2 <- filter(dat.1, `Target` == targetList[i])
                targetShift <-  select(dat.2, `Average Shift`) %>% unlist
                treatment.a <- select(dat.2, Treatment) %>% unlist() 
                treatment.b <- select(dat.2, `Time Point`) %>% unlist()
                treatment <- paste(treatment.a, '_', treatment.b, 'h', sep='') 
                meanShift <- mean(targetShift)
                target <- rep(targetList[i], length(treatmentList))
                holder <- numeric()
                for (i in 1:length(targetShift)){
                        holder[i] <- (targetShift[i] - meanShift) / meanShift
                }
                tmp <- data.frame(holder, target, treatment)
                df <- rbind(df, tmp)
        }
        treatmentList <- select(dat.1, `Treatment`) %>% unique()
        names(df) <- c('Change', 'Target', 'Treatment')
        
        #configure heatmap
        plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
                geom_tile(stat = "identity") + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                scale_fill_gradient2(low = "blue", mid = "white", 
                                     high = "red", midpoint = 0, space = "rgb", 
                                     na.value = "grey50", guide = "colourbar") +
                ggtitle(paste('Cell Line: GBM-', cLL.i, sep=''))
        
        #plot figure, uncomment to plot
        #plots
        
        #write csv file and save figure, uncomment to save
        filename <- paste("heatmap_GBM", cLL.i, '_alt', sep='')
        write_csv(df, paste(filename, '.csv', sep=''))
        ggsave(plots, file = paste(filename, '.png', sep=''), 
               width = 10, height = 6)
}

HeatmapAllAlt <- function(filename = "Combined_Averaged.csv",
        removeTxt = c('(+)-Serum', '(-)-Serum')){
        library(readr)
        library(ggplot2)
        library(dplyr)
        
        #get data
        dat <- read_csv(filename, col_types = cols())
        dat <- filter(dat, !Treatment %in% removeTxt)
        
        #create necessary lists and matrix
        cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
        timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
        targetList <- select(dat, `Target`) %>% unlist() %>% unique()
        treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
        
        dat.1 <- dat
        df <- data.frame()
        
        for (i in 1:length(targetList)){
                dat.2 <- filter(dat.1, `Target` == targetList[i])
                targetShift <-  select(dat.2, `Average Shift`) %>% unlist
                treatment.a <- select(dat.2, Treatment) %>% unlist() 
                treatment.b <- select(dat.2, `Time Point`) %>% unlist()
                treatment.c <- select(dat.2, `Cell Line`) %>% unlist()
                treatment <- paste(treatment.a, '_', treatment.b, 'h_GBM', 
                                   treatment.c, sep='') 
                meanShift <- mean(targetShift)
                target <- rep(targetList[i], length(treatmentList))
                holder <- numeric()
                for (i in 1:length(targetShift)){
                        holder[i] <- (targetShift[i] - meanShift) / meanShift
                }
                tmp <- data.frame(holder, target, treatment)
                df <- rbind(df, tmp)
        }
        treatmentList <- select(dat.1, `Treatment`) %>% unique()
        names(df) <- c('Change', 'Target', 'Treatment')
        
        #configure heatmap
        plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
                geom_tile(stat = "identity") + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      panel.background = element_blank()) +
                scale_fill_gradient2(low = "blue", mid = "white", 
                                     high = "red", midpoint = 0, space = "rgb", 
                                     na.value = "grey50", guide = "colourbar")
        
        # plot figure, uncomment to plot
        # plots
        
        # write csv file and save figure, uncomment to save
        filename <- "heatmap_all"
        write_csv(df, paste(filename, '.csv', sep=''))
        ggsave(plots, file = paste(filename, '.png', sep=''),
               width = 10, height = 6)
}

GoPlot <- function(){
	dat <- GetData()
	numTargets <- 1:length(unique(dat$Target))
	for(i in numTargets){
		PlotEachTarget(i)
	}
	numTimepoints <- length(unique(dat$`Time Point`))
	numCellLines <- length(unique(dat$`Cell Line`))
	for(i in 1:numTimepoints){
		for (j in 1:numCellLines){
			PlotEachTimepoint(i, j)
		}
	}
	for(i in 1:2){
		PlotEachCellLine(i)
		HeatmapCellLine(i)
		HeatmapCellLineAlt(i)
		HeatmapTimePoint(i)
		HeatmapTimePointAlt(i)
	}
	for(i in 1:4){
		HeatmapCellTime(i)
	        HeatmapCellTimeAlt(i)
	}
	HeatmapAll()
	HeatmapAllAlt()
}