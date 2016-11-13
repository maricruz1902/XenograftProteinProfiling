getData <- function(){
    library(readr)
    url <- "https://raw.githubusercontent.com/jwade1221/XenograftProteinProfiling/master/Combined_Averaged.csv"
    filename <- basename(url)
    download.file(url,destfile=filename)
    dat <- read_csv(filename, col_types = cols())
    dat
}

plotEachTarget <- function(i, filename = "Combined_Averaged.csv"){
    library(readr)
    library(ggplot2)
    
    #import data
    #dat <- read_csv(filename)
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

plotEachTimepoint <- function(i, j, filename = "CompiledAveraged.csv"){
    library(readr)
    library(ggplot2)
    library(RColorBrewer)
    
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
        geom_errorbar(limits, width = 0.4, size = 0.3, colour = 'grey',
                      position = position_dodge(width = 0.9)) +
        theme_bw() + theme(panel.grid = element_blank(), 
                           axis.title.x=element_blank()) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.key = element_rect(colour = 'white', 
                fill = 'white'), legend.key.size = unit(0.4, "cm")) +
        ggtitle(paste('Cell Line: ', cellLineList[j], 
                      'Time Point: ', timepointList[i], sep=''))
    
    #plot figure, uncomment to plot
    #plots
    
    #save figure, uncomment to save
    filename = paste(timepointList[i],cellLineList[j], 'png', sep='.')
    ggsave(plots, file = filename, width = 8, height = 6)
}

go <- function(){
    dat <- getData()
    numTargets <- length(unique(dat$Target))
    for(i in 1:numTargets){
        plotEachTarget(i)
    }
    numTimepoints <- length(unique(dat$`Time Point`))
    numCellLines <- length(unique(dat$`Cell Line`))
    for(i in 1:numTimepoints){
        for (j in 1:numCellLines){
            plotEachTimepoint(i, j)
        }
    }
}