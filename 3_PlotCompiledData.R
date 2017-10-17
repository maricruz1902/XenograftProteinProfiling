## Plot all data combined
PlotAllData <- function(){
        g <- ggplot(dat, aes(x = Target, y = NetShift))
        
        all.point <- g + 
                geom_point(aes(color = Treatment), 
                           position = "jitter", alpha = 0.7) +
                facet_grid(CellLine ~ TimePoint) + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle("Full Dataset")
        
        
        ggsave(all.point, filename = "everything_point.png", 
               width = 12, height = 8)
        
        all.boxplot <- g + 
                geom_boxplot(aes(fill = Treatment)) +
                facet_grid(CellLine ~ TimePoint) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle("Full Dataset")
        
        ggsave(all.boxplot, filename = "everything_boxplot.png",
               width = 12, height = 8)
}

## Plot each target as Net Shift vs Treatment
PlotEachTarget <- function(){
        targetList <- unique(dat$Target)
        
        for(i in targetList) {
                
                dat.tar <- filter(dat, Target == i & Treatment != "(+)-Serum" & 
                                          Treatment != "(-)-Serum")
                
                g <- ggplot(dat.tar, aes(x = Treatment, y = NetShift))
                
                plotName <- unlist(strsplit(i, "/"))[1]
                
                target.point <- g + 
                        geom_point(aes(color = Treatment), 
                                   position = "jitter") + 
                        facet_grid(CellLine~TimePoint) +
                        ggtitle(paste0("Target: ", i)) + 
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                target.box <- g + 
                        geom_boxplot(aes(fill = Treatment)) + 
                        facet_grid(CellLine~TimePoint) +
                        ggtitle(paste0("Target: ", i)) + 
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                dir.create("Target Plots", showWarnings = FALSE)
                
                ggsave(target.point, 
                       filename = paste0("Target Plots/", 
                                         plotName, "_point.png"), 
                       width = 12, height = 8)
                
                ggsave(target.box, 
                       filename = paste0("Target Plots/", 
                                         plotName, "_boxplot.png"),
                       width = 12, height = 8)
        }
}

## Plot each treatment as Net Shift vs Target
PlotEachTreatment <- function(){
        treatmentList <- unique(dat$Treatment)
        
        for (i in treatmentList) {
                
                dat.rx <- filter(dat, Treatment == i)
                
                g <- ggplot(dat.rx, aes(Target, NetShift))
                
                fig4 <- g + 
                        geom_point(aes(color = Target), position = "jitter") + 
                        facet_grid(CellLine~TimePoint) +
                        ggtitle(paste0("Treatment: ", i)) +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                fig5 <- g + 
                        geom_boxplot(aes(fill = Target)) + 
                        facet_grid(CellLine~TimePoint) +
                        ggtitle(paste0("Treatment: ", i)) +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                dir.create("Treatment Plots", showWarnings = FALSE)
                
                ggsave(fig4, filename = paste0("Treatment Plots/",
                                               i, "_point.png"),
                       width = 12, height = 8)
                
                ggsave(fig5, filename = paste0("Treatment Plots/",
                                               i, "_boxplot.png"),
                       width = 12, height = 8)
        }
}

## Plot Treatments
PlotTreatment <- function(control, treatment, targets, cellLine){
        dat.cntl <- filter(dat, Treatment == control & TimePoint == "1 h")
        dat.all <- rbind(filter(dat, Treatment == treatment), dat.cntl)
        dat.all$Treatment <- factor(dat.all$Treatment, 
                                    levels = c("DMSO", "(-)-Serum", 
                                               "Apitolisib", "Erlotinib",
                                               "Palbociclib", "GNE-317",
                                               "(+)-Serum"))
        
        # Plot all treatments for treatment
        g.all <- ggplot(dat.all, 
                        aes(x = interaction(TimePoint, Treatment, Target),
                            y = NetShift, 
                            fill = Target)) +
                geom_boxplot() + 
                facet_grid(CellLine~.) +
                labs(fill = "") + 
                scale_x_discrete(labels = rep(c("0 h", "1 h", "24 h"), 
                                              length(unique(dat.all$Target)))) +
                xlab("Treatment Time") + 
                ggtitle(paste0("Treatment: ", treatment)) +
                theme(axis.text.x = 
                              element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(g.all, 
               filename = paste0("Treatment Plots/", treatment, ".png"), 
               width = 12, height = 8)
        
        # Plot select treatments for treatment
        dat.rx <- filter(dat.all, Target %in% targets & 
                                 CellLine == cellLine)
        
        g <- ggplot(dat.rx, 
                    aes(x = interaction(TimePoint, Treatment, 
                                        Target, CellLine), 
                        y = NetShift, 
                        fill = Target))
        
        txt <- g +
                geom_boxplot() + 
                labs(fill = "") + 
                scale_x_discrete(labels = rep(c("0 h", "1 h", "24 h"), 
                                              length(targets))) +
                xlab("Treatment Time") + 
                ggtitle(paste0("Treatment: ", treatment, ", GBM 6")) +
                theme(axis.text.x = 
                              element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(txt, 
               filename = paste0("Treatment Plots/", treatment, "_", 
                                 cellLine, ".png"),
               width = 8, height = 6)
}

## Pairwise treatment comparisons
TreatmentComp <- function(treatments, targets, cellLine){
        dat.rx <- filter(dat, Treatment == treatments & Target %in% targets &
                                 CellLine == cellLine)
        print(head(dat.rx))
        g <- ggplot(dat.rx, aes(x = Target, 
                                group = interaction(Treatment, 
                                                    TimePoint,
                                                    Target),
                                y = NetShift))
        
        fig_comp <- g + 
                geom_boxplot(aes(fill = interaction(Treatment, TimePoint))) + 
                labs(fill = "") + 
                xlab("Target") + 
                ggtitle(paste0(treatments[1], " vs ", treatments[2])) +
                theme(axis.text.x = 
                              element_text(angle = 45, hjust = 1))
        
        fig_comp.2 <- g + 
                geom_boxplot(aes(fill = Treatment)) + 
                facet_wrap(~TimePoint) + 
                labs(fill = "") + 
                xlab("Target") + 
                ggtitle(paste0(treatments[1], " vs ", treatments[2])) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        dir.create(path = "Treatment Comparisons", showWarnings = FALSE)
        
        ggsave(fig_comp,
               filename = paste0("Treatment Comparisons/", treatments[1], "_", 
                                 treatments[2], "_", cellLine, ".png"), 
               width = 8, height = 6)
        ggsave(fig_comp.2, 
               filename = paste0("Treatment Comparisons/", treatments[1], "_", 
                                 treatments[2], "_", cellLine, "_2.png"), 
               width = 8, height = 6)
}

## Run all of the above functions to generate plots
PlotData <- function(){
        setwd("D:/Box Sync/Data/")
        
        # Load libraries and set theme for all plots
        library(tidyverse)
        library(ggthemes)
        theme_set(theme_few(base_size = 16))
        
        # Load in data to make plots
        dat <<- read_csv("compiledLabeled.csv")
        dat$Treatment <- factor(dat$Treatment, 
                                levels = c("DMSO", "(-)-Serum", 
                                           "Apitolisib", "Erlotinib",
                                           "Palbociclib", "GNE-317",
                                           "(+)-Serum"))
        
        # Save current wd to return to later and setwd to plots folder
        directory <- getwd()
        setwd("../XPP_Plots/")
        
        PlotAllData()
        
        PlotEachTreatment()
        
        PlotEachTarget()
        
        control <- "DMSO"
        txtList <- unique(dat$Treatment)
        compTargets <- c("pAktSer473", "pS6Ser235/6", "pS6Ser240/4",
                         "pp70S6KThr389", "pRbSer780", "pRbSer807/11")
        
        lapply(txtList, function(i){
                PlotTreatment(control = control, treatment = i,
                              targets = compTargets, cellLine = "GBM6")
                PlotTreatment(control = control, treatment = i,
                              targets = compTargets, cellLine = "GBM26")
        })
        
        
        
        # Pairwise treatment list
        txtPairs <- combn(unique(dat$Treatment), 2, simplify = FALSE)
        
        # Run through pair-wise list to plot treatment comparisons
        lapply(txtPairs, function(i){
                TreatmentComp(treatments = as.vector(i), targets = compTargets,
                              cellLine = "GBM6")
                TreatmentComp(treatments = as.vector(i), targets = compTargets,
                              cellLine = "GBM26")
        })
        setwd(directory)
}
