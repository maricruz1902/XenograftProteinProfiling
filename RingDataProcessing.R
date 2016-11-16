aggData <- function() {
    library(readr)
    
    ## get the working directory
    directory <- getwd()
    url <- "https://raw.githubusercontent.com/jwade1221/XenograftProteinProfiling/master/groupNames_allClusters.csv"
    filename <- basename(url)
    download.file(url,destfile=filename)
    recipe <- read_csv(filename, col_types = cols())
    #recipe <- read_csv("D:/Google Drive/Research/GitRepositories/XenograftProteinProfiling/groupNames_XPP.csv", col_types = cols())
    groupNames <- recipe$groupNames
    holder <- recipe[,c(1,2)]
    ## delete unnecessary files
    unlink("comments.csv", recursive = TRUE)
    unlink("groupNames*", recursive = TRUE)
    unlink("plots", recursive = TRUE)
    rings <- list.files(directory)
    
    ## create empty data frame
    df <- data.frame()
    ## add data to data frame corresponding to id
    for (i in rings) {
        ring <- as.vector(i)
        dat <- read_csv(ring, col_types = cols(), col_names = FALSE)
        time <- dat[ ,1]
        shift <- dat[ ,2]
        ringStr <- strsplit(i, "\\.")[[1]]
        ringNum <- as.numeric(ringStr[1])
        groupNum <- (ringNum - 1) %/% 4 + 1
        ring <- rep(ringNum, nrow(dat))
        group <- rep(groupNum, nrow(dat))
        if (groupNum == 35){groupNum <- 34}
        groupName <- as.character(holder$groupNames[[groupNum]])
        groupName <- rep(groupName, nrow(dat))
        tmp <- data.frame(ring, group, time, shift, groupName)
        df <- rbind(df, tmp)
    }
    names(df) <- c("ring", "group", "time", "shift", "groupName")
    
    if (file.exists("plots")){
        setwd("plots")
    } else {
        dir.create("plots")
        setwd("plots")
    }
    
    write_csv(df, "allRings.csv")
    setwd(directory)
}

thermalControl <- function(){
    library(readr)
    library(dplyr)
    
    directory = getwd()
    dat <- read_csv("plots/allRings.csv", col_types = cols())
    
    #get thermal control averages
    thermals <- filter(dat, groupName == "thermal")
    times <- filter(thermals, ring == 3) %>% select(time)
    ringList <- unique(thermals$ring)
    df.thermals <- data.frame(times)
    for (i in ringList){
        ringShift <- filter(thermals, ring == i) %>% select(shift)
        names(ringShift) <- paste('ring', i, sep='')
        df.thermals <- cbind(df.thermals, ringShift)
    }
    cols <- ncol(df.thermals)
    df.thermals$avgThermals <- rowMeans(df.thermals[,c(2:cols)])
    avgThermals <- as.vector(df.thermals$avgThermals)
    
    #subtract thermal controls
    ringNames <- unique(dat$ring)
    for(i in ringNames){
        ringDat <- filter(dat, ring == i) %>% select(shift)
        ringTC <- ringDat - avgThermals
        dat[dat$ring == i, 4] <- ringTC
    }
    
    write_csv(dat, "plots/allRings_tc.csv")
    
}

plotRingData <- function(thermal = TRUE){
    library(ggplot2)
    library(readr)
    library(RColorBrewer)
    
    directory <- getwd()
    if (thermal == TRUE){
    dat <- read_csv("plots/allRings_tc.csv", col_types = cols())
    } else {
        dat <- read_csv("plots/allRings.csv", col_types = cols())
    }
    #set colors for plot
    colorCount <- length(unique(dat$groupName))
    getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)
    
    #configure plot and legend
    plots <- ggplot(dat, aes(time, shift, colour = factor(groupName))) + 
        geom_point(size = 1) +
        xlab("Time (min)") + 
        ylab(expression(paste("Relative Shift (", Delta,"pm)"))) +
        scale_colour_manual(values = getPalette, name = 'Target') +
        theme_bw() + theme(panel.grid = element_blank(), 
            axis.title.x=element_blank()) + 
        theme(legend.key = element_rect(colour = 'white',
            fill = 'white'), legend.key.size = unit(0.4, "cm"))
        
    #plot figure, uncomment to plot
    #plots
    
    #save plot, uncomment to save
    filename <- "AllRingsTest.png"
    newDir <- "plots"
    if (file.exists("plots")){
        setwd("plots")
    } else {
        dir.create("plots")
        setwd("plots")
    }
    
    ggsave(plots, file = filename, width = 8, height = 6)
    setwd(directory)
}

getNetShifts <- function(thermal = TRUE){
    library(readr)
    library(dplyr)
    
    directory <- getwd()
    if (thermal == TRUE){
        dat <- read_csv("plots/allRings_tc.csv", col_types = cols())
    } else {
        dat <- read_csv("plots/allRings.csv", col_types = cols())
    }
    ringList <- unique(dat$ring)
    time1 <- 72
    time2 <- 53
    dat.rings <- data.frame()
    for (i in ringList){
        dat.ring <- filter(dat, ring == i)
        time1.loc <- which.min(abs(dat.ring$time - time1))
        time1.val <- dat.ring$shift[time1.loc]
        time2.loc <- which.min(abs(dat.ring$time - time2))
        time2.val <- dat.ring$shift[time2.loc]
        ring <- i
        group <- unique(dat.ring$group)
        groupName <- unique(dat.ring$groupName)
        tmp <- data.frame(i, group, groupName, time1.val, time2.val)
        dat.rings <- rbind(dat.rings, tmp)
    }
    names(dat.rings) <- c("ring", "group", "groupName", "shift.1", "shift.2")
    dat.rings$netShifts <- dat.rings$shift.1 - dat.rings$shift.2
    
    #write 'netShifts.csv' in plots/
    if (file.exists("plots")){
        setwd("plots")
    } else {
        dir.create("plots")
        setwd("plots")
    }
    write_csv(dat.rings, "netShifts.csv")
    setwd(directory)
}

plotNetShifts <- function(){
    library(ggplot2)
    library(readr)
    library(RColorBrewer)
    
    directory <- getwd()
    dat <- read_csv("plots/netShifts.csv", col_types = cols())
    
    colorCount <- nrow(dat)
    getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)
    
    plots <- ggplot(dat) +
        geom_bar(aes(groupName, netShifts, fill = factor(ring)), 
            stat = 'identity', position = 'dodge') +
        scale_fill_manual(values = getPalette, name = 'Rings') +
        ylab(expression(paste("Net Shift ( ", Delta,"pm)"))) +
        theme_bw() + theme(panel.grid = element_blank(), 
                       axis.title.x=element_blank()) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.key = element_rect(colour = 'white',
            fill = 'white'), legend.key.size = unit(0.3, "cm"))
    
    #save plot, uncomment to save
    filename <- "NetShiftTest.png"
    newDir <- "plots"
    if (file.exists("plots")){
        setwd("plots")
    } else {
        dir.create("plots")
        setwd("plots")
    }
    ggsave(plots, file = filename, width = 10, height = 6)
    setwd(directory)
}

getAvgShifts <- function(){
    library(readr)
    
    directory <- getwd()
    dat <- read_csv("plots/netShifts.csv", col_types = cols())
    targets <- unique(dat$groupName)
    targets
    df <- data.frame()
    
    for(i in targets){
        dat.group <- dat[dat$groupName == i,]
        dat.shifts <- dat.group$netShifts
        avgShift <- mean(dat.shifts)
        sdShift <- sd(dat.shifts)
        seShift <- sd(dat.shifts) / sqrt(length(dat.shifts))
        cvShift <- sd(dat.shifts) / avgShift * 100
        target <- unique(dat.group$groupName)
        tmp <- data.frame(target, avgShift, sdShift, seShift, cvShift)
        df <- rbind(df, tmp)
    }

    names(df) <- c("Target", "Average Shift", "SD", "SE", "CV")
    
    if (file.exists("plots")){
        setwd("plots")
    } else {
        dir.create("plots")
        setwd("plots")
    }
    write_csv(df, "avgShifts.csv")
    setwd(directory)
}

findBadClusters <- function(){
    library(readr)
    library(dplyr)
    dat <- read_csv("plots/avgShifts.csv", col_types = cols())
    groupList <- dat$Target
    groupList
    cvCutoff <- 20
    badGroups <- vector("numeric")
    counter <- 0
    for (i in 2:length(groupList)){
        groupCV <- filter(dat, Target == groupList[i]) %>% select(CV) %>% unlist()
        if (groupCV > cvCutoff){
            counter <- counter + 1
            badGroups[counter] <- groupList[i]
        }
    }
    badGroups <- as.data.frame(badGroups)
    if(length(badGroups) > 0){
        write_csv(badGroups, 'plots/badGroups.csv')
    }
}

plotBadClusters <- function(){
    library(readr)
    library(dplyr)
    library(RColorBrewer)
    library(ggplot2)
    
    directory <- getwd()
    badRings <- read_csv("plots/badGroups.csv", col_types = cols())
    badRings <- unlist(badRings)
    if (length(badRings) > 0){
        dat <- read_csv("plots/netShifts.csv", col_types = cols())
        dat.bad <- filter(dat, groupName %in% badRings)
        #configure colors
        colorCount <- nrow(dat.bad)
        getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)
        
        #configure plot
        plots <- ggplot(dat.bad) +
            geom_bar(aes(groupName, netShifts, fill = factor(ring)), 
                     stat = 'identity', position = 'dodge') +
            scale_fill_manual(values = getPalette, name = 'Rings') +
            ylab(expression(paste("Net Shift ( ", Delta,"pm)"))) +
            theme_bw() + theme(panel.grid = element_blank(), 
                               axis.title.x=element_blank()) + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            theme(legend.key = element_rect(colour = 'white',
                        fill = 'white'), legend.key.size = unit(0.3, "cm"))
        
        #save plot, uncomment to save
        filename <- "BadRings.png"
        newDir <- "plots"
        if (file.exists("plots")){
            setwd("plots")
        } else {
            dir.create("plots")
            setwd("plots")
        }
        ggsave(plots, file = filename, width = 10, height = 6)
        setwd(directory)
    }
}

plotAvgShifts <- function(){
    library(ggplot2)
    library(readr)
    library(RColorBrewer)
    
    directory <- getwd()
    dat <- read_csv("plots/avgShifts.csv", col_types = cols())
    
    colorCount <- nrow(dat)
    getPalette <- colorRampPalette(brewer.pal(8, "Paired"))(colorCount)
    
    #set error bars
    limits <- aes(ymax = `Average Shift` + SD,
                  ymin = `Average Shift` - SD)
    
    #configure plot and legend
    plots <- ggplot(dat, aes(x = Target, y = `Average Shift`, fill = Target)) +
        geom_bar(stat = 'identity') + geom_errorbar(limits, width = 0.3) +
        scale_fill_manual(values = getPalette) +
        ylab(expression(paste("Net Shift ( ", Delta,"pm)"))) +
        theme_bw() + theme(panel.grid = element_blank(), 
            axis.title.x=element_blank()) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(fill = FALSE)
    
    #plot figure, uncomment to plot
    plots
    
    #save figure, uncomment to save
    filename <- "AvgShiftTest.png"
    newDir <- "plots"
    
    if (file.exists("plots")){
        setwd("plots")
    } else {
        dir.create("plots")
        setwd("plots")
    }
    
    ggsave(plots, file = filename, width = 8, height = 6)
    setwd(directory)
}

go <- function(){
    library(ggplot2)
    library(readr)
    library(RColorBrewer)
    library(dplyr)
    
    directory <- getwd()
    aggData()
    thermalControl()
    plotRingData()
    getNetShifts()
    plotNetShifts()
    getAvgShifts()
    findBadClusters()
    plotBadClusters()
    plotAvgShifts()
    setwd(directory)
}

analyzeFullDataSet <- function(){
    directory <- getwd()
    listData <- list.dirs(recursive = FALSE)
    for (i in 1:length(listData)){
        tempDir <- listData[i]
        setwd(tempDir)
        #print(getwd())
        go()
        setwd(directory)
    }
}