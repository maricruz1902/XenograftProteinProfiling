CombineNetShifts <- function(){
        library(tidyverse)
        
        directory <- getwd()
        
        netList <- grep('net', 
                        list.files(pattern = '.csv', recursive = TRUE), 
                        value = TRUE)
        
        removeFiles <- grep("Combined", netList, value = TRUE)
        netList <- netList[!netList %in% removeFiles]
        
        netShifts <- lapply(netList, read_csv)
        
        netComb <- bind_rows(netShifts)
        netComb <- filter(netComb, !grepl("thermal|BL", Target))

        write_csv(netComb, path = "netShiftsCombined.csv")
        
        setwd(directory)
}

PlotCombineNetShifts <- function(){
        library(tidyverse)
        library(ggthemes)
        library(reshape2)
        
        dat <- read_csv("netShiftsCombined.csv")
        
        dat.melt <- melt(dat, id = c("Step", "Target", "Ring", "Channel"), measure = "Net Shift")
        
        ggplot(dat.melt, aes(x = Step, y = value, color = Target)) + geom_point() + facet_wrap(~Target) + theme_few(base_size = 16)
        
        plot_theme <- theme_few(base_size = 16)
        
        plot <- ggplot(dat, aes(x = Step, y = `Net Shift`, 
                                color = Target)) +
                geom_point()+
                plot_theme + ggtitle(name) + facet_wrap(~Target) +
                scale_x_log10() +
                labs(x = "Analyte Concentration (pg/mL)", 
                     y = expression(paste("Relative Shift (",Delta,"pm)")),
                     color = "Target")
                
        
        plot

        ggsave(plot = plot, 
               filename = paste0("netShiftsCombined_ch.png"),
               width = 8, height = 6)
}

Fit <- function(loc = 'plots'){
        dat <- read_csv("netShiftsCombined.csv")
        
        dat$Replicate <- as.factor((dat$Ring - 1) %% 4 + 1)
        
        fit <- list()
        
        targetList <- unique(dat$Target)
        
        plot.1 <- ggplot(dat, aes(x = Step, y = `Net Shift`, color = Target)) + geom_point()
                
                
        plot.1a <- plot.1 + scale_x_log10()
        
        # plot.2 <- ggplot(dat, aes(x = Step, y = `Net Shift`, color = Target)) + 
        #         geom_point() +
        #         geom_smooth(mapping = aes(x = Step, y = `Net Shift`), method = 'nls',
        #                     formula = y ~ A.2 + (A.1-A.2)/(1 + (x/x.0)^p),
        #                     method.args = list(start = c(A.2 =max(y),
        #                                                  A.1 = min(y),
        #                                                  x.0 = mean(y),
        #                                                  p = 1)),
        #                     se = FALSE) +
        #         labs(x = "Analyte Concentration (pg/mL)", 
        #              y = expression(paste("Relative Shift (",Delta,"pm)")),
        #              color = "Target")
        # 
        # plot.2a <- plot.2 + scale_x_log10()
        # 
        # plot.2a + geom_smooth()
        
        for(i in 1:length(targetList)) {
                tar <- targetList[i]
                print(tar)
                dat.tar <- filter(dat, Target == tar)
                y <- dat.tar$`Net Shift`
                x <- dat.tar$Step
                #fit[[i]] <- nls(y ~ SSlogis(x, Asym, xmid, scal))
                fit[[i]] <- summary(nls(formula = y ~ A.2 + (A.1-A.2)/(1 + (x/x.0)^p),
                                        start = list(A.2 =max(y),
                                                     A.1 = min(y),
                                                     x.0 = mean(y),
                                                     p = 1)))
        }
        
        summary(fit[[1]])
        
        capture.output(fit, file = "fitInfo.txt")
        ggsave(plot.1, filename = "plot1.png", width = 8, height = 6)
        ggsave(plot.1a, filename = "plot1a.png", width = 8, height = 6)
        # ggsave(plot.2, filename = "plot2.png", width = 8, height = 6)
        
}

CheckRingQuality <- function(loc = 'plots', 
                             varLevel = 500000) {
        library(tidyverse)
        
        dat <- read_csv(paste0(loc,"/", name, "_allRings.csv"))
        
        dat.avg <- dat %>% group_by(Ring) %>%
                summarise_each(funs(var), c(Shift))
        
        ringWinners <- filter(dat.avg, Shift < varLevel) %>% select(Ring)
        ringLosers <- filter(dat.avg, Shift > varLevel) %>% select(Ring)
        
        write_csv(ringWinners, paste0(loc, '/', name, "_ringWinners.csv"))
        write_csv(ringLosers, paste0(loc, '/', name, "_ringLosers.csv"))
}

AnalyzeData <- function(loc = "plots", cntl = "thermal", filename = 'groupNames_allClusters.csv') {
        GetName()
        runConc <- basename(getwd())
        ch1_MCP1 <- unlist(strsplit(runConc, split = "_"))[4]
        # ch1_IL6 <- unlist(strsplit(runConc, split = "_"))[5]
        ch2_MCP1 <- unlist(strsplit(runConc, split = "_"))[6]
        # ch2_IL6 <- unlist(strsplit(runConc, split = "_"))[8]
        AggData()
        SubtractControl(ch = 1, cntl = cntl)
        SubtractControl(ch = 2, cntl = cntl)
        # PlotRingData(cntl = cntl, ch = 1, splitPlot = TRUE)
        # PlotRingData(cntl = cntl, ch = 2, splitPlot = TRUE)
        GetNetShifts(cntl = cntl, ch = 1, target = "MCP-1",
                     time1 = 52, time2 = 41, step = ch1_MCP1)
        # GetNetShifts(cntl = cntl, ch = 1, target = "IL-6",
        #              time1 = 52, time2 = 41, step = ch1_IL6)
        GetNetShifts(cntl = cntl, ch = 2, target = "MCP-1",
                             time1 = 52, time2 = 41, step = ch2_MCP1)
        # GetNetShifts(cntl = cntl, ch = 2, target = "IL-6",
        #                      time1 = 52, time2 = 41, step = ch2_IL6)
        # PlotNetShifts(cntl = cntl, ch = 1, step = ch1_MCP1, target = "MCP-1")
        # PlotNetShifts(cntl = cntl, ch = 2, step = ch2_MCP1, target = "MCP-1")
}

AnalyzeAllData <- function() {
        foldersList <- list.dirs(recursive = FALSE)
        for (i in foldersList){
                directory <- getwd()
                setwd(i)
                AnalyzeData()
                setwd(directory)
        }
        CombineNetShifts()
        PlotCombineNetShifts()
        Fit()
}
