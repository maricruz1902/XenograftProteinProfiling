aggData <- function() {
    ## get the working directory
    directory <- getwd()
    groupList <- seq(1:34)
    groupNames <- read.csv("C:/Users/james_000/Desktop/groupNames.csv")
    holder <- cbind(groupNames, groupList)
    setwd(directory)
    
    ## delete unnecessary files
    unlink("comments.csv", recursive = TRUE)
    unlink("allRings.csv", recursive = TRUE)
    unlink("netShifts.csv", recursive = TRUE)
    unlink("*.txt", recursive = TRUE)
    unlink("*.png", recursive = TRUE)
    files_list <- list.files(directory, full.names=TRUE)
    rings <- list.files(directory)
    x <- length(list.files(directory))
    id <- 1:x
    
    ## create empty data frame
    df <- data.frame()
    
    ## add data to data frame corresponding to id
    library(readr)
    for (i in rings) {
        ring <- as.vector(i)
        dat <- read_csv(ring, col_types = cols(), col_names = FALSE)
        time <- dat[ ,1]
        shift <- dat[ ,2]
        ringStr <- strsplit(ring, "\\.")[[1]]
        ringNum <- as.numeric(ringStr[1])
        groupNum <- (ringNum - 1) %/% 4 + 1
        ring <- rep(ringNum, length(time))
        group <- rep(groupNum, length(time))
        if (groupNum == 35){groupNum <- 34}
        gN <- as.character(holder$groupName[[groupNum]])
        groupName <- rep(gN, length(time))
        tmp <- data.frame(ring, group, time, shift, groupName)
        df <- rbind(df, tmp)
    }
    
    names(df) <- c("ring", "group", "time", "shift", "groupName")
    
    if (file.exists("../plots")){setwd("../plots")
    } else {
        dir.create("../plots")
        setwd("../plots")
    }
    write_csv(df, "../plots/allRings.csv")
    setwd(directory)
}

plotRingData <- function(){
    library(ggplot2)
    library(readr)
    directory <- getwd()
    dat <- read_csv("../plots/allRings.csv", col_types = cols())
    
    #set colors for plot
    colorCount <- length(unique(dat$groupName))
    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))(colorCount)
    
    #configure plot and legend
    plots <- ggplot(dat) + 
        geom_point(aes(time, shift, colour = factor(groupName))) +
        xlab("Time (min)") + 
        ylab(expression(paste("Relative Shift ( ", Delta,"pm)"))) +
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
    if (file.exists("../plots")){setwd("../plots")
    } else {
            dir.create("../plots")
            setwd("../plots")
        }
    
    ggsave(plots, file = filename, width = 8, height = 6)
    setwd(directory)
}

getNetShifts <- function(){
    directory <- getwd()
    dat <- read_csv("../plots/allRings.csv", col_types = cols())
    time1 <- 52
    time2 <- 39
    time1.low <- time1 - 0.07
    time1.high <- time1 + 0.07
    time2.low <- time2 - 0.07
    time2.high <- time2 + 0.07
    dat.tmp <- as.data.frame(c(dat[dat$time < time1.high & dat$time > time1.low, ], 
                               dat[dat$time < time2.high & dat$time > time2.low, ]))
    keepCols <- c('ring', 'group', 'time', 'shift','groupName', 'time.1', 'shift.1')
    dat.rings <- dat.tmp[keepCols]
    colnames(dat.rings)[6] <- c("time.1")
    colnames(dat.rings)[7] <- c("shift.1")
    dat.rings$netShifts <- dat.rings$shift - dat.rings$shift.1
    
    #write 'netShifts.csv' in ../plots/
    if (file.exists("../plots")){setwd("../plots")
    } else {
        dir.create("../plots")
        setwd("../plots")
    }
    write_csv(dat.rings, "../plots/netShifts.csv")
    setwd(directory)
}

plotNetShifts <- function(){
    library(ggplot2)
    
    directory <- getwd()
    dat <- read_csv("../plots/netShifts.csv", col_types = cols())
    
    colorCount <- nrow(dat)
    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))(colorCount)
    
    plots <- ggplot(dat) +
        geom_bar(aes(groupName, netShifts, fill = factor(ring)), 
            stat = 'identity', position = 'dodge') +
        scale_fill_manual(values = getPalette, name = 'Rings') +
        ylab(expression(paste("Net Shift ( ", Delta,"pm)"))) +
        theme_bw() + theme(panel.grid = element_blank(), 
                       axis.title.x=element_blank()) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.key = element_rect(colour = 'white',
            fill = 'white'), legend.key.size = unit(0.4, "cm"))
    
    #save plot, uncomment to save
    filename <- "NetShiftTest.png"
    newDir <- "plots"
    if (file.exists("../plots")){setwd("../plots")
    } else {
            dir.create("../plots")
            setwd("../plots")
    }
    
    ggsave(plots, file = filename, width = 8, height = 6)
    setwd(directory)
}

getAvgShifts <- function(){
    library(ggplot2)
    library(readr)
    
    directory <- getwd()
    dat <- read_csv("../plots/netShifts.csv", col_types = cols())
    targets <- unique(dat$groupName)
    targets
    df <- data.frame()
    
    for(i in targets){
        dat.group <- dat[dat$groupName == i,]
        dat.shifts <- dat.group$netShifts
        avgShift <- mean(dat.shifts)
        sdShift <- sd(dat.shifts)
        seShift <- sd(dat.shifts) / length(dat.shifts)
        target <- unique(dat.group$groupName)
        tmp <- data.frame(target, avgShift, sdShift, seShift)
        df <- rbind(df, tmp)
    }

    names(df) <- c("Target", "Average Shift", "SD", "SE")
    
    head(df)
    if (file.exists("../plots")){setwd("../plots")
    } else {
        dir.create("../plots")
        setwd("../plots")
    }
    write_csv(df, "../plots/avgShifts.csv")
    setwd(directory)
}

plotAvgShifts <- function(){
    library(ggplot2)
    
    directory <- getwd()
    dat <- read_csv("../plots/avgShifts.csv", col_types = cols())
    
    colorCount <- nrow(dat)
    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))(colorCount)
    
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
    if (file.exists("../plots")){setwd("../plots")
    } else {
        dir.create("../plots")
        setwd("../plots")
    }
    
    ggsave(plots, file = filename, width = 8, height = 6)
    setwd(directory)
}

go <- function(){
    aggData()
    plotRingData()
    getNetShifts()
    plotNetShifts()
    getAvgShifts()
    plotAvgShifts()
}
