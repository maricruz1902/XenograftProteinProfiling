heatmapCellTime <- function(){
    library(readr)
    library(ggplot2)
    library(dplyr)
    
    #get data
    url <- "https://raw.githubusercontent.com/jwade1221/XenograftProteinProfiling/master/Combined_Averaged.csv"
    filename <- basename(url)
    download.file(url,destfile=filename)
    dat <- read_csv(filename, col_types = cols())
    
    #create necessary lists and matrix
    cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
    timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
    targetList <- select(dat, `Target`) %>% unlist() %>% unique()
    treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
    cLL <- rep(1:length(cellLineList), length(timePointList))
    tPL <- rep(1:length(timePointList), each = length(cellLineList))
    m <- as.matrix(cbind(tPL, cLL))
    
    for (i in 1:nrow(m)){
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
        
        #write 'heatMap_'...'.csv'
        filename <- paste("heatmap_", cLL.i, '_', tPL.i, sep='')
        write_csv(df, paste(filename, '.csv', sep=''))
        plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
            geom_tile(stat = "identity") +
            scale_fill_gradient(low = "blue", high = "yellow", 
                    space = "Lab", na.value = "grey50", guide = "colourbar")
        ggsave(plots, file = paste(filename, '.png', sep=''), 
               width = 10, height = 6)
    }
}

heatmapCellLine <- function(){
    library(readr)
    library(ggplot2)
    library(dplyr)
    
    #get data
    url <- "https://raw.githubusercontent.com/jwade1221/XenograftProteinProfiling/master/Combined_Averaged.csv"
    filename <- basename(url)
    download.file(url,destfile=filename)
    dat <- read_csv(filename, col_types = cols())
    
    #create necessary lists and matrix
    cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
    timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
    targetList <- select(dat, `Target`) %>% unlist() %>% unique()
    treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
    cLL <- rep(1:length(cellLineList), length(timePointList))
    tPL <- rep(1:length(timePointList), each = length(cellLineList))
    m <- as.matrix(cbind(tPL, cLL))
    
    for (i in 1:nrow(m)){
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
        
        #write 'heatMap_'...'.csv'
        filename <- paste("heatmap_", cLL.i, sep='')
        write_csv(df, paste(filename, '.csv', sep=''))
        plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
            geom_tile(stat = "identity") + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_fill_gradient(low = "blue", high = "yellow", 
                    space = "Lab", na.value = "grey50", guide = "colourbar")
        ggsave(plots, file = paste(filename, '.png', sep=''), 
               width = 10, height = 6)
    }
}

heatmapTimePoint <- function(){
    library(readr)
    library(ggplot2)
    library(dplyr)
    
    #get data
    url <- "https://raw.githubusercontent.com/jwade1221/XenograftProteinProfiling/master/Combined_Averaged.csv"
    filename <- basename(url)
    download.file(url,destfile=filename)
    dat <- read_csv(filename, col_types = cols())
    
    #create necessary lists and matrix
    cellLineList <- select(dat, `Cell Line`) %>% unlist() %>% unique()
    timePointList <- select(dat, `Time Point`) %>% unlist() %>% unique()
    targetList <- select(dat, `Target`) %>% unlist() %>% unique()
    treatmentList <- select(dat, `Treatment`) %>% unlist() %>% unique()
    cLL <- rep(1:length(cellLineList), length(timePointList))
    tPL <- rep(1:length(timePointList), each = length(cellLineList))
    m <- as.matrix(cbind(tPL, cLL))
    
    for (i in 1:nrow(m)){
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
        
        #write 'heatMap_'...'.csv'
        filename <- paste("heatmap_", tPL.i, sep='')
        write_csv(df, paste(filename, '.csv', sep=''))
        plots <- ggplot(df, aes(x = Treatment, y = Target, fill = Change)) +
            geom_tile(stat = "identity") + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_fill_gradient(low = "blue", high = "yellow", 
                                space = "Lab", na.value = "grey50", guide = "colourbar")
        ggsave(plots, file = paste(filename, '.png', sep=''), 
               width = 10, height = 6)
    }
}