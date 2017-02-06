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
        
        # change this file name to use alternative ring or group labels
        filename <- "groupNames_allClusters.csv"
        
        # get information of chip layout from github repository
        if (!file.exists("groupNames_allClusters.csv")){
        url <- "https://raw.githubusercontent.com/jwade1221/XenograftProteinProfiling/master/groupNames_allClusters.csv"
        filename <- basename(url)
        download.file(url, filename)
        }
        
        # define recipe as global variable for use in other functions
        recipe <<- read_csv(filename, col_types = cols())
        targets <- recipe$Target
        
        # generate list of rings to analyze (gets all *.csv files)
        rings <- list.files(directory, pattern = ".csv", recursive = FALSE)
        idfile <- grepl("group", rings)
        removeFiles <- c("comments.csv", rings[idfile])
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
        time_point <- seq(1:nrow(dat))
        tmp <- data.frame(ring, group, time_shift, shift, groupName, 
                      channel, run, time_point)
        df <- rbind(df, tmp)
        }
        
        # renames columns in df
        names(df) <- c("Ring", "Group", "Time", "Shift", "Target", "Channel",
                 "Experiment", "Time Point")
        
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
        dat <- read_csv(paste0(loc, "/", name, "_", "allRings.csv"))
        if (ch != "U"){
        dat <- filter(dat, Channel == ch)
        }
        
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
                                    fill = 'white'), legend.key.size = unit(0.4, "cm")) + 
        geom_vline(xintercept = c(5, 10), linetype = "longdash")
        
        if (cntl == "raw"){
        plots <- plots + geom_point(size = 1) + facet_grid(.~ Channel)
        } else {plots <- plots + geom_point(size = 1)}
        
        #plot figure, uncomment to plot
        # plots
        
        # alternative plot
        dat.2 <- dat %>% group_by(Target, `Time Point`) %>% summarise_each(funs(mean, sd), c(Time, Shift))
        head(dat.2)
        g <- ggplot(dat.2, aes(Time_mean, Shift_mean, color = Target))
                
        plot2 <- g + geom_line(size = 1) + 
                geom_ribbon(aes(ymin = Shift_mean - Shift_sd, 
                        ymax = Shift_mean + Shift_sd, linetype = NA), 
                        fill = "slategrey", alpha = 1/3) + 
                theme_bw() + theme(panel.grid = element_blank()) +
                theme(legend.key = element_rect(colour = 'white', fill = 'white'), 
                        legend.key.size = unit(0.4, "cm")) +
                xlab("Time (min)") + 
                ylab(expression(paste("Relative Shift (",Delta,"pm)")))
        
        
        #save plot, uncomment to save
        filename <- paste0(name, "_", cntl, "Control", "_ch", ch)
        filename2 <- paste0(name, "_", cntl, "Control", "_ch", ch, "_avg")
        setwd(loc)
        ggsave(plots, file = paste0(filename, ".png"), width = 8, height = 6)
        ggsave(plots, file = paste0(filename, ".pdf"), width = 8, height = 6)
        ggsave(plot2, file = paste0(filename2, ".png"), width = 8, height = 6)
        ggsave(plot2, file = paste0(filename2, ".pdf"), width = 8, height = 6)
        setwd(directory)
}

GetNetShifts <- function(cntl, ch, loc = 'plots', time1, time2, step = 1){
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
                                  experiment, channel, step)
                dat.rings <- rbind(dat.rings, tmp)
        }
        
        # renames dat.rings columns
        names(dat.rings) <- c("Ring", "Group", "Target", "Shift.1", "Shift.2", 
                              "Experiment", "Channel", "Step")
        
        # calculate nat shift and create new column in dataframe
        dat.rings$`Net Shift` <- dat.rings$Shift.1 - dat.rings$Shift.2
        
        # save net shift data
        setwd(loc)
        write_csv(dat.rings, paste0(name, "_netShifts_ch", ch, "step_", step, ".csv"))
        setwd(directory)
}

PlotNetShifts <- function(cntl, ch, loc = 'plots', step = 1){
        # load relevant libraries
        library(ggplot2)
        library(readr)
        library(dplyr)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # get net shift data
        if (cntl != "raw"){
                dat <- read_csv(paste0(loc, "/", name, "_netShifts_ch", ch, "step_", step,
                                       ".csv"), col_types = cols())
        } else {
                dat <- read_csv(paste(loc,"/", name, "_netShifts_chU.csv", 
                                      sep=""), col_types = cols())
                
        }
        
        # configure plot and legend
        dat.nothermal <- filter(dat, Target != "thermal")
        
        plots <- ggplot(dat.nothermal, aes(Target, `Net Shift`, color = Target)) +
                geom_boxplot() +
                ylab(expression(paste("Net Shift ( ", Delta,"pm)"))) +
                theme_bw() + theme(panel.grid = element_blank(), 
                                   axis.title.x=element_blank()) + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                theme(legend.key = element_rect(colour = 'white',
                                                fill = 'white'), legend.key.size = unit(0.3, "cm"))
        
        # save plot, uncomment to save
        filename <- paste0(name, "_NetShift_ch", ch, "_step", step, ".png", sep="")
        setwd(loc)
        ggsave(plots, file = filename, width = 10, height = 6)
        setwd(directory)
}

AnalyzeData <- function() {
        GetName()
        AggData()
        SubtractControl(ch = 1, cntl = "thermal")
        SubtractControl(ch = 2, cntl = "thermal")
        SubtractControl(ch = "U", cntl = "thermal")
        PlotRingData(cntl = "thermal", ch = 1)
        PlotRingData(cntl = "thermal", ch = 2)
        PlotRingData(cntl = "thermal", ch = "U")
        PlotRingData(cntl = "raw", ch = "U")
        GetNetShifts(cntl = "thermal", ch = 1, time1 = 53, time2 = 39, step = 1)
        GetNetShifts(cntl = "thermal", ch = 1, time1 = 25, time2 = 14, step = 2)
        GetNetShifts(cntl = "thermal", ch = 1, time1 = 13, time2 = 2, step = 3)
        GetNetShifts(cntl = "thermal", ch = 2, time1 = 53, time2 = 39, step = 1)
        GetNetShifts(cntl = "thermal", ch = 2, time1 = 25, time2 = 14, step = 2)
        GetNetShifts(cntl = "thermal", ch = 2, time1 = 13, time2 = 2, step = 3)
        GetNetShifts(cntl = "thermal", ch = "U", time1 = 53, time2 = 39, step = 1)
        GetNetShifts(cntl = "thermal", ch = "U", time1 = 25, time2 = 14, step = 2)
        GetNetShifts(cntl = "thermal", ch = "U", time1 = 13, time2 = 2, step = 3)
        PlotNetShifts(cntl = "thermal", ch = "1", step = 1)
        PlotNetShifts(cntl = "thermal", ch = "1", step = 2)
        PlotNetShifts(cntl = "thermal", ch = "1", step = 3)
        PlotNetShifts(cntl = "thermal", ch = "2", step = 1)
        PlotNetShifts(cntl = "thermal", ch = "2", step = 2)
        PlotNetShifts(cntl = "thermal", ch = "2", step = 3)
        PlotNetShifts(cntl = "thermal", ch = "U", step = 1)
        PlotNetShifts(cntl = "thermal", ch = "U", step = 2)
        PlotNetShifts(cntl = "thermal", ch = "U", step = 3)
}


