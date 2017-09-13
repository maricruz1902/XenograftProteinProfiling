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

Fit <- function(loc = 'plots', tarList){
        library(scales)
        dat <- read_csv("netShiftsCombined.csv")
        theme_set(theme_few(base_size = 16))
        
        dat <- filter(dat, Target %in% tarList)
        
        
        dat$Replicate <- as.factor((dat$Ring - 1) %% 4 + 1)
        
        fit <- list()
        
        plot.1 <- ggplot(dat, aes(x = Step, 
                                  y = `Net Shift`, 
                                  color = factor(Group))) + 
                geom_point()
                
                
        plot.1a <- plot.1 + scale_x_log10()
        
        plot.2 <- ggplot(dat, aes(x = Step, y = `Net Shift`)) +#, color = Target)) +
                geom_point() +
                geom_smooth(method = 'nls',
                            formula = y ~ A + (B - A) / (1 + (x / C) ^ D),
                            method.args = list(start = c(A = 10000,
                                                         B = 100,
                                                         C = 4000,
                                                         D = 1)),
                            se = FALSE) +
                labs(x = "Analyte Concentration (pg/mL)",
                     y = expression(paste("Relative Shift (",Delta,"pm)")))

        plot.2a <- plot.2 + coord_trans(x = "log10")
        
        for(i in 1:length(tarList)) {
                tar <- tarList[i]
                dat.tar <- filter(dat, Target == tar)
                y <- dat.tar$`Net Shift`
                x <- dat.tar$Step
                #fit[[i]] <- nls(y ~ SSlogis(x, Asym, xmid, scal))
                fit.info <- nls(formula = y ~ A + (B - A) / (1 + (x / C) ^ D),
                                start = list(A = 10000,
                                             B = 100,
                                             C = 4000,
                                             D = 1),
                                control = list(tol = 5e-08, 
                                               minFactor = 0,
                                               maxiter  = 10000))
                fit[[i]] <- summary(fit.info)
                
                A <- as.numeric(coef(fit.info)[1])
                B <- as.numeric(coef(fit.info)[2])
                C <- as.numeric(coef(fit.info)[3])
                D <- as.numeric(coef(fit.info)[4])
                
                testFun <- function(x) {A + (B - A) / (1 + (x / C) ^ D)}
                
                plot.3 <- ggplot(dat.tar, aes(x = Step, y = `Net Shift`)) +
                        geom_boxplot(aes(group = Step), fill = "red") +
                        stat_function(fun = testFun, color = "blue", size = 1) +
                        scale_x_log10(breaks = 
                                        trans_breaks("log10", function(x) 10^x),
                                      labels = 
                                        trans_format("log10", 
                                                     math_format(10 ^ .x))) +
                        labs(x = "Analyte Concentration (pg/mL)",
                             y = expression(paste("Relative Shift (",
                                                  Delta,
                                                  "pm)"))) +
                        annotation_logticks()
                
                ggsave(plot.3, filename = paste0("plot3_", tar, ".png"),
                       width = 8, height = 6)
        }
        
        capture.output(fit, file = "fitInfo.txt")
        ggsave(plot.1, filename = "plot1.png", width = 8, height = 6)
        ggsave(plot.1a, filename = "plot1a.png", width = 8, height = 6)
        ggsave(plot.2, filename = "plot2.png", width = 8, height = 6)
        ggsave(plot.2a, filename = "plot2a.png", width = 8, height = 6)
        
}

AnalyzeCalData <- function() {
        GetName()
        runConc <- basename(getwd())
        ch1_IL6 <- unlist(strsplit(runConc, split = "_"))[4]
        # ch1_IL6 <- unlist(strsplit(runConc, split = "_"))[5]
        ch2_IL6 <- unlist(strsplit(runConc, split = "_"))[6]
        # ch2_IL6 <- unlist(strsplit(runConc, split = "_"))[8]
        AggData(filename = 'groupNames_allClusters.csv')
        SubtractControl(ch = 1, cntl = "thermal")
        SubtractControl(ch = 2, cntl = "thermal")
        # PlotRingData(cntl = cntl, ch = 1, splitPlot = TRUE)
        # PlotRingData(cntl = cntl, ch = 2, splitPlot = TRUE)
        GetNetShifts(cntl = "thermal", ch = 1, 
                     time1 = 53, time2 = 41, step = ch1_IL6)
        # GetNetShifts(cntl = cntl, ch = 1, target = "IL-6",
        #              time1 = 52, time2 = 41, step = ch1_IL6)
        GetNetShifts(cntl = "thermal", ch = 2,
                     time1 = 53, time2 = 41, step = ch2_IL6)
        # GetNetShifts(cntl = cntl, ch = 2, target = "IL-6",
        #                      time1 = 52, time2 = 41, step = ch2_IL6)
        # PlotNetShifts(cntl = cntl, ch = 1, step = ch1_IL6, target = "IL-6")
        # PlotNetShifts(cntl = cntl, ch = 2, step = ch2_IL6, target = "IL-6")
}

CalibrationStation <- function() {
        foldersList <- list.dirs(recursive = FALSE)
        for (i in foldersList){
                directory <- getwd()
                setwd(i)
                AnalyzeData()
                setwd(directory)
        }
        AnalyzeCalData()
        PlotCombineNetShifts()
        Fit(tarList = "IL-6")
}