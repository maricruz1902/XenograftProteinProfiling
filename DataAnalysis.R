getData <- function(){
    library(readr)
    url <- "https://raw.githubusercontent.com/jwade1221/XenograftProteinProfiling/master/Combined_Averaged_NoSerumNoMinSerum.csv"
    filename <- basename(url)
    download.file(url,destfile=filename)
    dat <- read_csv(filename, col_types = cols())
    dat
}

plotEachTarget <- function(i, filename = "Combined_Averaged.csv"){
    library(readr)
    library(ggplot2)
    library(RColorBrewer)
    
    #import data
    dat <- read_csv(filename, col_types = cols())
    
    #create lists for data islation
    targetList <- unique(dat$Target)
    dat.Target <- dat[dat$Target == targetList[i],]
    
    #set error bars
    limits <- aes(ymax = `Relative Shift` + `Standard Deviation`, 
                  ymin = `Relative Shift` - `Standard Deviation`)
    
    #get colors
    colorCount <- length(unique(dat$`Cell Line`)) + 
                length(unique(dat$`Time Point`))
    getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)
    
    #plot figure
    plots <- ggplot(data = dat.Target,
				aes(x = Treatment, y = `Relative Shift`,
				fill = interaction(factor(`Time Point`), factor(`Cell Line`)),
				group = interaction(factor(`Time Point`), factor(`Cell Line`)))) +
        geom_bar(stat = "identity", position = "dodge") +
        ylab("Net Shift (pm)") +
        geom_errorbar(limits,position = position_dodge(width = 0.9), 
                      width = 0.4) +
        theme_bw() + theme(panel.grid = element_blank(), 
				axis.title.x=element_blank()) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        scale_fill_manual(values = getPalette, name="Cell Line and Treatment", 
                labels=c("GBM-6 1 h Treatment", "GBM-6 24 h Treatment", 
                "GBM-26 1 h Treatment", "GBM-26 24 h Treatment")) + 
        theme(legend.key = element_rect(colour = 'white',
                fill = 'white'), legend.key.size = unit(0.4, "cm"))+
        ggtitle(paste('Target: ', targetList[i], sep=''))
    
    #plot figure, uncomment to plot
    #plots
                
    #save figure, uncomment to save
    filename = paste(targetList[i], 'png', sep='.')
    ggsave(plots, file = filename, width = 8, height = 6)
}

plotEachTimepoint <- function(i, j, filename = "Combined_Averaged.csv"){
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
    limits <- aes(ymax = `Relative Shift` + `Standard Deviation`, 
                  ymin = `Relative Shift` - `Standard Deviation`)
    
    #set colors for plot
    colorCount <- length(unique(dat$Target))
    getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)
    
    #configure plot, colors, and legend
    plots <- ggplot(data = dat.Timepoint,
        # group lets `ggplot` know we want different errorbars/bars for each day
        aes(x = Treatment, y = `Relative Shift`, fill = factor(Target))) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
        ylab("Net Shift (pm)") +
        scale_fill_manual(values = getPalette, name = 'Target') +
        geom_errorbar(limits, width = 0.4, size = 0.3, colour = 'grey40',
                      position = position_dodge(width = 0.9)) +
        theme_bw() + theme(panel.grid = element_blank(), 
                           axis.title.x=element_blank()) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.key = element_rect(colour = 'white', 
                fill = 'white'), legend.key.size = unit(0.4, "cm")) +
        ggtitle(paste('Cell Line: ', cellLineList[j], 
                      ' Time Point: ', timepointList[i], sep=''))
    
    #plot figure, uncomment to plot
    #plots
    
    #save figure, uncomment to save
    filename = paste(timepointList[i],cellLineList[j], 'png', sep='.')
    ggsave(plots, file = filename, width = 8, height = 6)
}

plotEachCellLine <- function(i, filename = "Combined_Averaged.csv"){
    library(readr)
    library(ggplot2)
    library(RColorBrewer)
    
    #import data
    dat <- read_csv(filename, col_types = cols())
    
    #create lists for data islation
    cellLineList <- unique(dat$`Cell Line`)
    dat.Timepoint <- dat[dat$`Cell Line` == cellLineList[i],]
    
    #set error bars
    limits <- aes(ymax = `Relative Shift` + `Standard Deviation`, 
                  ymin = `Relative Shift` - `Standard Deviation`)
    
    #set colors for plot
    colorCount <- length(unique(dat$Target)) * length(unique(dat$`Time Point`))
    getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)
    
    #configure plot, colors, and legend
    plots <- ggplot(data = dat.Timepoint,
        aes(x = Treatment, y = `Relative Shift`, fill = interaction(factor(Target),
                factor(`Time Point`)))) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
        ylab("Net Shift (pm)") +
        scale_fill_manual(values = getPalette, name = 'Target') +
        geom_errorbar(limits, width = 0.4, size = 0.3, colour = 'grey40',
                      position = position_dodge(width = 0.9)) +
        theme_bw() + theme(panel.grid = element_blank(), 
                           axis.title.x=element_blank()) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.key = element_rect(colour = 'white', 
                                        fill = 'white'), legend.key.size = unit(0.4, "cm")) +
        ggtitle(paste('Cell Line: ', cellLineList[i], sep=''))
    
    #plot figure, uncomment to plot
    #plots
    
    #save figure, uncomment to save
    filename = paste(cellLineList[i], 'png', sep='.')
    ggsave(plots, file = filename, width = 12, height = 6)
}


heatmapCellTime <- function(i, filename = "Combined_Averaged.csv"){
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
    tPL.i <- timePointList[m[i,1]]
    dat.1 <- filter(dat, `Cell Line` == cLL.i & `Time Point` == tPL.i)
    df <- data.frame()
    
    for (i in 1:length(targetList)){
        dat.2 <- filter(dat.1, `Target` == targetList[i])
        targetShift <-  select(dat.2, `Relative Shift`) %>% unlist()
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
    filename <- paste("heatmap_", cLL.i, '_', tPL.i, sep='')
    write_csv(df, paste(filename, '.csv', sep=''))
    ggsave(plots, file = paste(filename, '.png', sep=''), 
          width = 10, height = 6)
}

heatmapCellLine <- function(i, filename = "Combined_Averaged.csv"){
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
        targetShift <-  select(dat.2, `Relative Shift`) %>% unlist
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
    filename <- paste("heatmap_", cLL.i, sep='')
    write_csv(df, paste(filename, '.csv', sep=''))
    ggsave(plots, file = paste(filename, '.png', sep=''), 
           width = 10, height = 6)
}

heatmapTimePoint <- function(i, filename = "Combined_Averaged.csv"){
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
    
    tPL.i <- timePointList[m[i,1]]
    dat.1 <- filter(dat, `Time Point` == tPL.i)
    df <- data.frame()
        
    for (i in 1:length(targetList)){
        dat.2 <- filter(dat.1, `Target` == targetList[i])
        targetShift <-  select(dat.2, `Relative Shift`) %>% unlist
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
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_gradient2(low = "blue", mid = "white", 
            high = "red", midpoint = 0, space = "rgb", 
            na.value = "grey50", guide = "colourbar") +
        ggtitle(paste('Time Point: ', tPL.i, ' h', sep=''))
    
    #plot figure, uncomment to plot
    #plots
    
    #write csv file and save figure, uncomment to save
    filename <- paste("heatmap_", tPL.i, sep='')
    write_csv(df, paste(filename, '.csv', sep=''))
    ggsave(plots, file = paste(filename, '.png', sep=''), 
           width = 10, height = 6)
}

go <- function(){
    dat <- getData()
    numTargets <- length(unique(dat$Target))
    for(i in 1:numTargets){
        plotEachTarget(i, filename = "Combined_Averaged_NoSerumNoMinSerum.csv")
    }
    numTimepoints <- length(unique(dat$`Time Point`))
    numCellLines <- length(unique(dat$`Cell Line`))
    for(i in 1:numTimepoints){
        for (j in 1:numCellLines){
            plotEachTimepoint(i, j, filename = "Combined_Averaged_NoSerumNoMinSerum.csv")
        }
    }
    for(i in 1:2){
        plotEachCellLine(i, filename = "Combined_Averaged_NoSerumNoMinSerum.csv")
        heatmapCellLine(i, filename = "Combined_Averaged_NoSerumNoMinSerum.csv")
        heatmapTimePoint(i, filename = "Combined_Averaged_NoSerumNoMinSerum.csv")
    }
    for(i in 1:4){
        heatmapCellTime(i, filename = "Combined_Averaged_NoSerumNoMinSerum.csv")
    }
}