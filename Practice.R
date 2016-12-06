HeatmapAllAlt <- function(filename = "combinedAveraged.csv",
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
                treatment <- paste(treatment.a, '_GBM', treatment.c, '_', 
                        treatment.b, 'h', sep='') 
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
        
        xPos <- sort(as.vector(unique(df$Treatment)))
        yPos <- sort(as.vector(unique(df$Target)), decreasing = TRUE)
        #configure heatmap
        plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
                geom_raster() + theme(axis.text.x = element_text(angle = 45, 
                        hjust = 1), panel.background = element_blank()) +
                scale_fill_gradient2(low = "blue", mid = "white", 
                        high = "red", midpoint = 0, space = "rgb", 
                        na.value = "grey50", guide = "colourbar") + 
                scale_x_discrete(limits = xPos) + 
                scale_y_discrete(limits = yPos)
        
        # plot figure, uncomment to plot
        plots
        
        # write csv file and save figure, uncomment to save
        # filename <- "heatmap_all"
        # write_csv(df, paste(filename, '.csv', sep=''))
        # ggsave(plots, file = paste(filename, '.png', sep=''),
        #       width = 10, height = 6)
}