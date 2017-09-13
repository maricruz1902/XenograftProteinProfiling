# Split data for heatmaps and correlation plots
DataSplitHeat <- function(){
        setwd("D:/Box Sync/Data/")
        library(tidyverse)
        library(reshape2)
        library(gplots) # heatmap.2
        library(Hmisc) # for rcorr function
        library(corrplot)
        
        dat <- read_csv("compiledLabeled.csv")
        
        dat.1 <- filter(dat, !grepl("p53|Abl", Target))
        dat.6 <- filter(dat.1, `Cell Line` == "GBM 6")
        dat.26 <- filter(dat.1, `Cell Line` == "GBM 26")
        dat.1h <- filter(dat.1, `Time Point` == "1 h")
        dat.24h <- filter(dat.1, `Time Point` == "24 h")
        dat.6.1 <- filter(dat.1, `Time Point` == "1 h",
                          `Cell Line` == "GBM 6")
        dat.6.24 <- filter(dat.1, `Time Point` == "24 h",
                           `Cell Line` == "GBM 6")
        dat.26.1 <- filter(dat.1, `Time Point` == "1 h",
                           `Cell Line` == "GBM 26")
        dat.26.24 <- filter(dat.1, `Time Point` == "24 h",
                            `Cell Line` == "GBM 26")
        
        casting.1 <- dcast(dat.1, Treatment + `Cell Line` + `Time Point` ~
                             Target, mean, value.var = "Net Shift")
        casting.6 <- dcast(dat.6, Treatment + `Cell Line` + `Time Point` ~
                             Target, mean, value.var = "Net Shift")
        casting.26 <- dcast(dat.26, Treatment + `Cell Line` + `Time Point` ~
                              Target, mean, value.var = "Net Shift")
        casting.1h <- dcast(dat.1h, Treatment + `Cell Line` + `Time Point` ~
                              Target, mean, value.var = "Net Shift")
        casting.24h <- dcast(dat.24h, Treatment + `Cell Line` + `Time Point` ~
                              Target, mean, value.var = "Net Shift")
        casting.6.1 <- dcast(dat.6.1, Treatment + `Cell Line` + `Time Point` ~
                                  Target, mean, value.var = "Net Shift")
        casting.6.24 <- dcast(dat.6.24, Treatment + `Cell Line` + `Time Point` ~
                               Target, mean, value.var = "Net Shift")
        casting.26.1 <- dcast(dat.26.1, Treatment + `Cell Line` + `Time Point` ~
                                Target, mean, value.var = "Net Shift")
        casting.26.24 <- dcast(dat.26.24, Treatment + `Cell Line` +
                               `Time Point` ~ Target, mean,
                               value.var = "Net Shift")
        
        rownames(casting.1) = paste(casting.1$Treatment,
                                    casting.1$`Cell Line`,
                                    casting.1$`Time Point`,
                                    sep = " ")
        rownames(casting.6) = paste(casting.6$Treatment,
                                    casting.6$`Cell Line`,
                                    casting.6$`Time Point`,
                                    sep = " ")
        rownames(casting.26) = paste(casting.26$Treatment,
                                     casting.26$`Cell Line`, 
                                     casting.26$`Time Point`,
                                     sep = " ")
        rownames(casting.1h) = paste(casting.1h$Treatment,
                                     casting.1h$`Cell Line`,
                                     casting.1h$`Time Point`,
                                     sep = " ")
        rownames(casting.24h) = paste(casting.24h$Treatment,
                                      casting.24h$`Cell Line`,
                                      casting.24h$`Time Point`,
                                      sep = " ")
        rownames(casting.6.1) = paste(casting.6.1$Treatment,
                                      casting.6.1$`Cell Line`, 
                                      casting.6.1$`Time Point`,
                                      sep = " ")
        rownames(casting.6.24) = paste(casting.6.24$Treatment,
                                       casting.6.24$`Cell Line`,
                                       casting.6.24$`Time Point`,
                                       sep = " ")
        rownames(casting.26.1) = paste(casting.26.1$Treatment,
                                       casting.26.1$`Cell Line`,
                                       casting.26.1$`Time Point`,
                                       sep = " ")
        rownames(casting.26.24) = paste(casting.26.24$Treatment,
                                        casting.26.24$`Cell Line`,
                                        casting.6.24$`Time Point`,
                                        sep = " ")

        casting.1[,c(1,2,3)] <- NULL
        casting.6[,c(1,2,3)] <- NULL
        casting.26[,c(1,2,3)] <- NULL
        casting.1h[,c(1,2,3)] <- NULL
        casting.24h[,c(1,2,3)] <- NULL
        casting.6.1[,c(1,2,3)] <- NULL
        casting.6.24[,c(1,2,3)] <- NULL
        casting.26.1[,c(1,2,3)] <- NULL
        casting.26.24[,c(1,2,3)] <- NULL
        
        casting.1.m <<- data.matrix(casting.1)
        casting.6.m <<- data.matrix(casting.6)
        casting.26.m <<- data.matrix(casting.26)
        casting.1h.m <<- data.matrix(casting.1h)
        casting.24h.m <<- data.matrix(casting.24h)
        casting.6.1.m <<- data.matrix(casting.6.1)
        casting.6.24.m <<- data.matrix(casting.6.24)
        casting.26.1.m <<- data.matrix(casting.26.1)
        casting.26.24.m <<- data.matrix(casting.26.24)
        
        both.both.tar <<- rcorr(casting.1.m)
        both.both.txt <<- rcorr(t(casting.1.m))
        GBM6.both.tar <<- rcorr(casting.6.m)
        GBM6.both.txt <<- rcorr(t(casting.6.m))
        GBM26.both.tar <<- rcorr(casting.26.m)
        GBM26.both.txt <<- rcorr(t(casting.26.m))
        both.1.tar <<- rcorr(casting.1h.m)
        both.1.txt <<- rcorr(t(casting.1h.m))
        both.24.tar <<- rcorr(casting.24h.m)
        both.24.txt <<- rcorr(t(casting.24h.m))
        GBM6.1.tar <<- rcorr(casting.6.1.m)
        GBM6.1.txt <<- rcorr(t(casting.6.1.m))
        GBM6.24.tar <<- rcorr(casting.6.24.m)
        GBM6.24.txt <<- rcorr(t(casting.6.24.m))
        GBM26.1.tar <<- rcorr(casting.26.1.m)
        GBM26.1.txt <<- rcorr(t(casting.26.1.m))
        GBM26.24.tar <<- rcorr(casting.26.24.m)
        GBM26.24.txt <<- rcorr(t(casting.26.24.m))
}

# Heatmaps
HeatPlots <- function(){
        # hmcol <- colorRampPalette(c("darkblue", "white", "darkred"))
        hmcol <- colorRampPalette(c("#3a5387", "#e8e8e8", "#b25752"))(256)
        
        png('Heatmap_All.png', width = 16, height = 12,
            units = "in", res = 100)
        heatmap.2(casting.1.m, scale = "col", col = hmcol, key = T, 
                  margins = c(20, 20), trace = "none",
                  main = "GBM 6 & 26 with 1h & 24 h")
        dev.off()
        
        png('Heatmap_GBM6.png', width = 16, height = 12,
            units = "in", res = 100)
        heatmap.2(casting.6.m, scale = "col", col = hmcol, key = T, 
                  margins = c(20, 20), trace = "none",
                  main = "GBM 6 with 1h & 24h")
        dev.off()
        
        png('Heatmap_GBM26.png', width = 16, height = 12,
            units = "in", res = 100)
        heatmap.2(casting.26.m, scale = "col", col = hmcol, key = T, 
                  margins = c(20, 20), trace = "none",
                  main = "GBM 26 with 1h & 24h")
        dev.off()
        
        png('Heatmap_1h_both.png', width = 16, height = 12,
            units = "in", res = 100)
        heatmap.2(casting.1h.m, scale = "col", col = hmcol, key = T, 
                  margins = c(20, 20), trace = "none",
                  main = "GBM 6 & 26 with 1h")
        dev.off()
        
        png('Heatmap_24h_both.png', width = 16, height = 12,
            units = "in", res = 100)
        heatmap.2(casting.24h.m, scale = "col", col = hmcol, key = T, 
                  margins = c(20, 20), trace = "none",
                  main = "GBM 6 & 26 with 24h")
        dev.off()
        
        png('Heatmap_GBM6_1h.png', width = 16, height = 12,
            units = "in", res = 100)
        heatmap.2(casting.6.1.m, scale = "col", col = hmcol, key = T, 
                  margins = c(20, 20), trace = "none",
                  main = "GBM 6 with 1h")
        dev.off()
        
        png('Heatmap_GBM6_24h.png', width = 16, height = 12,
            units = "in", res = 100)
        heatmap.2(casting.6.24.m, scale = "col", col = hmcol, key = T, 
                  margins = c(20, 20), trace = "none",
                  main = "GBM 6 with 24h")
        dev.off()
        
        png('Heatmap_GBM26_1h.png', width = 16, height = 12,
            units = "in", res = 100)
        heatmap.2(casting.26.1.m, scale = "col", col = hmcol, key = T, 
                  margins = c(20, 20), trace = "none",
                  main = "GBM 26 with 1h")
        dev.off()
        
        png('Heatmap_GBM26_24h.png', width = 16, height = 12,
            units = "in", res = 100)
        heatmap.2(casting.26.24.m, scale = "col", col = hmcol, key = T, 
                  margins = c(20, 20), trace = "none",
                  main = "GBM 26 with 24h")
        dev.off()
}

# Corrplots
CorrPlots <- function(){
        png('Corrplot_All_Targets.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(both.both.tar$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 1h, 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_All_Treatments.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(both.both.txt$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 1h, 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        # Time Points
        
        png('Corrplot_1h_Targets.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(both.1.tar$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 1h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_1h_Treatments.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(both.1.txt$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 1h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_24h_Targets.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(both.24.tar$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_24h_Treatments.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(both.24.txt$r,
                 title = "Cell Lines: GBM 6, 26 Time Points: 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        # GBM 6
        
        png('Corrplot_GBM6_Targets.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM6.both.tar$r,
                 title = "Cell Lines: GBM 6 Time Points: 1h, 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        
        png('Corrplot_GBM6_Treatments.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM6.both.txt$r,
                 title = "Cell Lines: GBM 6 Time Points: 1h, 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_GBM6_1h_Treatments.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM6.1.txt$r,
                 title = "Cell Lines: GBM 6 Time Points: 1h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_GBM6_1h_Targets.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM6.1.tar$r,
                 title = "Cell Lines: GBM 6 Time Points: 1h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_GBM6_24h_Targets.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM6.24.tar$r,
                 title = "Cell Lines: GBM 6 Time Points: 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_GBM6_24h_Treatments.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM6.24.txt$r,
                 title = "Cell Lines: GBM 6 Time Points: 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        # GBM 26
        
        png('Corrplot_GBM26_Targets.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM26.both.tar$r,
                 title = "Cell Lines: GBM 26 Time Points: 1h, 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        
        png('Corrplot_GBM26_Treatments.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM26.both.txt$r,
                 title = "Cell Lines: GBM 26 Time Points: 1h, 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_GBM26_1h_Treatments.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM26.1.txt$r,
                 title = "Cell Lines: GBM 26 Time Points: 1h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_GBM26_1h_Targets.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM26.1.tar$r,
                 title = "Cell Lines: GBM 26 Time Points: 1h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_GBM26_24h_Treatments.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM26.24.txt$r,
                 title = "Cell Lines: GBM 26 Time Points: 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
        
        png('Corrplot_GBM26_24h_Targets.png', width = 16, height = 12,
            units = "in", res = 100)
        corrplot(GBM26.24.tar$r,
                 title = "Cell Lines: GBM 26 Time Points: 24h",
                 method = "color",
                 diag = FALSE,
                 type = "lower",
                 tl.col = "black",
                 mar = c(1, 1, 1, 1))
        dev.off()
}