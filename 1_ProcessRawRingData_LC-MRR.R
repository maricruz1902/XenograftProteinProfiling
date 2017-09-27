GetName <- function(){
        # get the filename from the current working directory
        directory <- basename(getwd())
        
        # directory naming from MRR: "CHIPNAME_gaskGASK_DATE"
        # extracts and returns GASK from directory name
        name <- unlist(strsplit(directory, split = "_"))[1]
        
        # define name as global variable for use in other functions
        name <<- gsub('gask','',name) # removes "gask" from name
}

AggData <- function(loc = 'plots', filename = 'groupNames_unmodified.csv') {
        # load relevant libraries
        library(tidyverse)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # get information of chip layout from github repository
        if (!file.exists(filename)){
                git <- "https://raw.githubusercontent.com/JamesHWade/XenograftProteinProfiling/master/"
                url <- paste0(git, filename)
                filename <- basename(url)
                download.file(url, filename)
        }
        
        # read in recipe/chip layout
        recipe <- read_csv(filename)
        colnames(recipe)[1] <- "Target" # rename column to remove byte order mark
        targets <- recipe$Target
        
        # generate list of rings to analyze (gets all *.csv files)
        rings <- list.files(directory, pattern = ".csv", recursive = FALSE)
        idfile <- grepl("group", rings)
        removeFiles <- c("comments.csv", rings[idfile])
        rings <- rings[!rings %in% removeFiles]
        
        # add data to data frame corresponding for each ring in rings
        df <- lapply(rings, function(i){
                dat <- read_csv(i, col_names = FALSE)
                ringNum <- as.numeric(strsplit(i, "\\.")[[1]][1])
                recipeCol <- which(recipe$Ring == ringNum)
                tmp <- dat[,c(1,2)] # time and shift from raw data
                tmp$ring <- ringNum
                tmp$group <- recipe$Group[recipeCol]
                tmp$groupName <- as.character(recipe$Target[[recipeCol]])
                tmp$channel <- recipe$Channel[[recipeCol]]
                tmp$run <- name
                tmp$timePoint <- seq(1:nrow(dat))
                tmp
        })
        
        # combine data from list into single data frame
        df <- bind_rows(df)
        
        # renames columns in df
        names(df) <- c("Time", "Shift", "Ring", "Group", "Target", "Channel",
                       "Experiment", "Time Point")
        
        # creates "plots" directory
        dir.create(loc, showWarnings = FALSE)
        
        # saves aggregated data with name_allRings.csv
        write_csv(df, paste(loc, '/', name, "_allRings.csv", sep=""))
}

DataSplitting <- function(loc = "plots", smoothed = TRUE){
        library(tidyverse)
        library(reshape2)
        library(pracma)
        library(baseline)
        
        # use thermally controlled data if desired
        dat <- read_csv(paste0("plots/", name, "_allRings.csv"))
        dat <- filter(dat, Target != "thermal")
        dat$Ring <- factor(dat$Ring)
        
        #separates into subsets corrdinating to each PS injection
        deadTime <- 30
        startTime <- deadTime + 30
        incTime <- 30.9
        df0 <- subset(dat, Time > deadTime & Time < deadTime + 30)
        df1 <- subset(dat, Time > startTime & Time < startTime + 30)
        df2 <- subset(dat, Time > startTime + (incTime * 1) & 
                              Time < startTime + (incTime * 1) + 30)
        df3 <- subset(dat, Time > startTime + (incTime * 2) & 
                              Time < startTime + (incTime * 2) + 30)
        df4 <- subset(dat, Time > startTime + (incTime * 3) & 
                              Time < startTime + (incTime * 3) + 30)
        df5 <- subset(dat, Time > startTime + (incTime * 4) & 
                              Time < startTime + (incTime * 4) + 30)
        df6 <- subset(dat, Time > startTime + (incTime * 5) & 
                              Time < startTime + (incTime * 5) + 30)
        df7 <- subset(dat, Time > startTime + (incTime * 6) & 
                              Time < startTime + (incTime * 6) + 30)
        # df8 <- subset(dat, Time > startTime + (incTime * 7) & 
        #                       Time < startTime + (incTime * 7) + 30)
        # df9 <- subset(dat, Time > startTime + (incTime * 8) & 
        #                       Time < startTime + (incTime * 8) + 30)
        # df10 <- subset(dat, Time > startTime + (incTime * 9) & 
        #                       Time < startTime + (incTime * 9) + 30)
        # df11 <- subset(dat, Time > startTime + (incTime * 9) & 
        #                        Time < startTime + (incTime * 9) + 30)
        
        
        #add column with standard name 
        df0$Standard <- "Ethyl Acetate"
        df1$Standard <- "1.3 kDa"
        df2$Standard <- "3.5 kDa"
        df3$Standard <- "8.7 kDa"
        df4$Standard <- "17.6 kDa"
        df5$Standard <- "35 kDa"
        df6$Standard <- "130 kDa"
        df7$Standard <- "304 kDa"
        # df9$Standard <- "Broad 20 kDa"
        # df10$Standard <- "Broad 47 kDa"
        # df11$Standard <- "Broad 60 kDa"
        
        runsList <- list(df2, df3)#, df0, df1, df4, df5, df6, df7)#, df8, df9, df10, df11)
        
        smoothData <- lapply(runsList, function(i){
                i %>% group_by(Ring) %>% 
                        mutate(Loess = loess(Shift~Time, span = 0.03)$fitted)
        })
        
        smoothData2 <- lapply(smoothData, function(i){
                i %>% group_by(Ring) %>%
                        mutate(Loess = Loess - Loess[1],
                               #Whit = Whit - Whit[1], 
                               Time = Time - Time[1],
                               Shift = Shift - Shift[1],
                               # Better baseline correction
                               Baseline =
                                       baseline.fillPeaks(t(Loess),
                                                          lambda = 3,
                                                          hwi = 5,
                                                          it = 200,
                                                          int = 1000)$corrected,
                               
                               # Fast baseline correction
                               Baseline2 =
                                       as.vector(baseline.medianWindow(t(Loess),
                                                       hwm = 300)$corrected),
                                                        
                               # airPLS baseline correction
                               Baseline3 = airPLS(Shift, lambda = 5,
                                                 differences = 2,itermax = 2000),
                               Volume = Time * 0.3)
                
        })
        
        dat.all <- bind_rows(smoothData2)
        
        g <- ggplot(dat.all, aes(x = Volume,
                            y = Baseline,
                            color = Standard,
                            group = interaction(Ring,Standard))) +
                geom_line() + facet_wrap(~Ring) + ylim(-5, 25)
        
        g2 <- ggplot(dat.all, aes(x = Volume,
                                 y = Baseline2,
                                 color = Standard,
                                 group = interaction(Ring,Standard))) + 
                geom_line() + facet_wrap(~Ring) + ylim(-5, 25)
        
        g3 <- ggplot(dat.all, aes(x = Volume,
                                  y = Loess,
                                  color = Standard,
                                  group = interaction(Ring,Standard))) + 
                geom_line() + facet_wrap(~Ring) + ylim(-5, 25)
        
        ggsave(g, filename = "1.png", width = 8, height = 6)
        ggsave(g2, filename = "2.png", width = 8, height = 6)
        ggsave(g3, filename = "3.png", width = 8, height = 6)
        
        write_csv(dat.all, 
                  path = paste0(loc, "/", name, "_allRings_Smoothed.csv"))
        
        dat.avg <- dat.all %>% group_by(Standard, `Time Point`) %>%
                summarise_at(vars(Time, Shift, Volume, Baseline, Loess), mean)
        
        dat.avg <- filter(dat.avg, Standard != "Ethyl Acetate")

        write_csv(dat.avg, 
                  path = paste0(loc, "/", name, "_inj_combined.csv"))
}

PlotIndyRings <- function(loc = "plots", smoothed = TRUE) {
        library(tidyverse)
        library(ggthemes)
        
        # get working directory to reset at end of function
        
        dat <- read_csv(paste0(loc, "/", name, 
                               "_allRings_Smoothed.csv"))
        
        dat$Ring <- factor(dat$Ring)
        
        if (smoothed) {
                plot <- ggplot(dat, aes(y = Baseline, x = Volume, color = Ring))
                plotName <- paste0(name, "_Smooth")
        } else {
                plot <- ggplot(dat, aes(y = Shift, x = Volume, color = Ring))
                plotName <- name
        }
        
        
        
        plot <- plot + geom_line() + facet_grid(Group~Standard) + 
                xlab("Elution Volume (mL)") +
                ylab(expression(paste("Relative Shift (",Delta,"pm)"))) +
                theme_few()
        
        plot
        ggsave(plot, filename = paste0("plots/", plotName, "_IndyRings.png"), 
               width = 12, height = 6)
}

PlotAvgData <- function(loc = 'plots'){
        # loads relevant libraries and plot theme
        library(tidyverse)
        library(ggthemes)

        # use thermally controlled data if desired
        dat <- read_csv(paste0(loc, "/", name, "_inj_combined.csv"))
        
        plot_theme <- theme_few(base_size = 16)
        
        smooth <- ggplot(dat, 
                         aes(y = Baseline, x = Volume, color = Standard)) +
                geom_line(size = 0.5) + xlab("Elution Volume (mL)") +
                ylab(expression(paste("Relative Shift (",Delta,"pm)"))) + 
                plot_theme
        
        smooth.zoom <- smooth + xlim(2, 4)
        
        smooth.comp <- smooth + 
                geom_line(aes(y = Shift, x = Volume), alpha = 1/2)
        
        smooth.wrap <- smooth + facet_wrap(~Standard)
        
        ggsave(smooth, 
               filename = paste0(loc, "/", name, " Smooth.png"), 
               width = 8, height = 6)
        ggsave(smooth.zoom, 
               filename = paste0(loc, "/", name, " Smooth2.png"), 
               width = 8, height = 6)
        ggsave(smooth.comp, 
               filename = paste0(loc, "/", name, " Smoothed and Averaged.png"), 
               width = 8, height = 6)
        ggsave(smooth.wrap, 
               filename = paste0(loc, "/", name, " Smooth Wrap.png"), 
               width = 8, height = 6)
}

FindPeaksCalcPDI <- function(loc = "plots"){
        library(tidyverse)
        library(ggthemes)
        
        dat <- read_csv(paste0(loc, "/", name, "_inj_combined.csv"))
        
        stdList <- unique(dat$Standard)
        stdList <- stdList[stdList != "Ethyl Acetate"]
        
        peaksData <- list()
        PDI <- list()
        
        m <- -2.2
        b <- 7
        
        for( i in stdList ){
                dat.std <- filter(dat, Standard == i)
                peaks <- findpeaks(x = dat.std$Baseline, 
                                   npeaks = 1, 
                                   minpeakheight = 1)
                startTime <- peaks[3]
                endTime <- peaks[4]
                
                peaksData[[i]] <- dat.std[c(startTime:endTime),]
                peaksData[[i]] <- mutate(peaksData[[i]], 
                                         molWeight = 10^(m * Volume + b))
                peaksData[[i]] <- mutate(peaksData[[i]], 
                                         M_n = sum(molWeight * Baseline) / sum(Baseline))
                peaksData[[i]] <- mutate(peaksData[[i]], 
                                         M_w = sum(molWeight ^ 2 * Baseline) / sum(Baseline * molWeight))
                peaksData[[i]] <- mutate(peaksData[[i]], 
                                         PDI = M_w / M_n)
                
                PDI[i] <- unique(peaksData[[i]]$PDI)
        }
        
        peaksData
        
        capture.output(PDI, file = paste0(loc, "/", "PDI.txt"))
        
        dat.a <- bind_rows(peaksData)
        plot <- ggplot(dat.a, aes(y = Baseline, 
                                  x =Volume, 
                                  color = Standard)) + 
                geom_line() + ggtitle(name)
        
        ggsave(plot, filename = paste0(loc, "/", name, "_PDIPeaks.png"), 
               width = 8, height = 6)
        
}

CheckRingQuality <- function(loc = 'plots', smoothed = TRUE) {
        library(tidyverse)
        
        dat <- read_csv(paste0(loc,"/", name, "_allRings_Smoothed.csv"))
        dat <- subset(dat, Time > 20 & Time < 30)
        
        if (smoothed) {
                smoothVar <- "Baseline"
        } else {
                smoothVar <- "Shift"
        }
        
        ggplot(dat, aes(x = Volume, y = Baseline, color = factor(Ring))) +
                 geom_point() + facet_wrap(~Ring)
        
        dat.var <- dat %>% group_by(Ring) %>% 
                summarize_at(vars(smoothVar), funs(Variance = var))
        topRings <- arrange(dat.var, Variance) %>% select(Ring) %>% head(10)
        
        dat.topRings <- filter(dat, Ring %in% topRings)
        
        
        write_csv(topRings, path = paste0(loc, '/', name, "_ringList.csv"))
        write_csv(dat.topRings, paste0(loc, '/', name, "_topRings.csv"))
}

AnalyzeData <- function(loc = "plots"){
        setwd("D:/Google Drive/Research/Manuscripts/Paper06_IsocraticPolymer/Data/MRR_CalCurve Data_n=1-3/0519/")
        GetName()
        AggData()
        DataSplitting(loc = loc)
        PlotIndyRings(loc = loc)
        PlotIndyRings(loc = loc, smoothed = FALSE)
        PlotAvgData(loc = loc)
        FindPeaksCalcPDI(loc = loc)
        CheckRingQuality(loc = loc, smoothed = T)
}

AnalyzeSelectData <- function(loc = "plots", select = "selectRings"){
        library(tidyverse)
        dir.create(select)
        
        ringList <- read_csv(paste0(loc, '/', name, "_ringList.csv"))
        ringList <- mutate(ringList, Files = 
                                   ifelse(Ring < 10,
                                          paste0("0", Ring, ".csv"), 
                                          paste0(Ring, ".csv")))
        fileList <- select(ringList, Files)
        
        lapply(fileList, function(i){
                file.copy(from = i, to = select)
        })
        directory <- getwd()
        setwd(select)
        AnalyzeData()
        setwd(directory)
}

BatchAnalyze <- function(){
        dataList <- list.dirs(recursive = FALSE)
        directory <- getwd()
        for (i in dataList){
                setwd(i)
                AnalyzeData()
                setwd(directory)
        }
}

WhittakerSmooth <- function(x,w,lambda,differences=1) {
        x=matrix(x,nrow = 1, ncol=length(x))
        L=length(x)
        E=spMatrix(L,L,i=seq(1,L),j=seq(1,L),rep(1,L))
        D=as(diff(E,1,differences),"dgCMatrix")
        W=as(spMatrix(L,L,i=seq(1,L),j=seq(1,L),w),"dgCMatrix")
        background=solve((W+lambda*t(D)%*%D),t((w*x)));
        return(as.vector(background))
}

airPLS <- function(x,lambda=10, differences=1, itermax=2000){
  
        x = as.vector(x)
        m = length(x)
        w = rep(1,m)
        control = 1
        i = 1
        while(control==1){
                z = WhittakerSmooth(x,w,lambda,differences)
                d = x-z
                sum_smaller = abs(sum(d[d<0])) 
                if(sum_smaller<0.001*sum(abs(x))||i==itermax)
                {
                       control = 0
                }
                w[d>=0] = 0
                w[d<0] = exp(i*abs(d[d<0])/sum_smaller)
                w[1] = exp(i*max(d[d<0])/sum_smaller)
                w[m] = exp(i*max(d[d<0])/sum_smaller)
                i=i+1
        }
        return(z)
}

arPLS(test, lambda = 100, ratio = 1)
