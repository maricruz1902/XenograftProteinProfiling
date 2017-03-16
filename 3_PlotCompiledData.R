library(tidyverse)

# Check if data file exists in working directory
# Download file if not found

filename <- "compiledLabeled.csv"

if (!file.exists("compiledLabeled.csv")){
        url <- "https://raw.githubusercontent.com/JamesHWade/XenograftProteinProfiling/master/compiledLabeled.csv"
        filename <- basename(url)
        download.file(url, filename)
}


# read in data
dat <- read_csv("compiledLabeled.csv")


# save current wd to return to later and setwd to plots folder
directory <- getwd()
setwd("../XPP_Plots/")

# plotting all data
g <- ggplot(dat, aes(Target, `Net Shift`))

fig <- g + geom_point(aes(color = Treatment, shape = factor(Replicate)), 
               position = "jitter", alpha = 1/2, size = 3) +
        facet_grid(`Cell Line` ~ `Time Point`) + theme_bw() +
        theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(fig, filename = "everything_point.png", width = 8, height = 6)
ggsave(fig, filename = "everything_point.pdf", width = 8, height = 6)


fig1 <- g + geom_boxplot(aes(color = Treatment)) +
        facet_grid(`Cell Line` ~ `Time Point`) + theme_bw() +
        theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(fig1, filename = "everything_boxplot.png", width = 8, height = 6)
ggsave(fig1, filename = "everything_boxplot.pdf", width = 8, height = 6)

# write loop to go through each target and plot the result

targetList <- unique(dat$Target)

for(i in targetList) {
        
        dat.tar <- filter(dat, Target == i & Treatment != "(+)-Serum" & 
                        Treatment != "DMSO")
        
        dat.tar$Treatment <- factor(dat.tar$Treatment, level = c("(+)-Serum", 
                "(-)-Serum", "DMSO", "Erlotinib", "GNE-317", "Apitolisib",
                "Palbociclib"))
        
        g <- ggplot(dat.tar, aes(Treatment, `Net Shift`, color = Treatment))
        
        plotName <- unlist(strsplit(i, "/"))[1]
        
        fig2 <- g + geom_point(position = "jitter", 
                               aes(shape = factor(Replicate))) + 
                facet_grid(`Cell Line`~`Time Point`) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle(i) + theme_bw() +
                theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggsave(fig2, filename = paste0(plotName, "_point.png"), width = 8, height = 6)
        ggsave(fig2, filename = paste0(plotName, "_point.pdf"), width = 8, height = 6)
        
        fig3 <- g + geom_boxplot() + facet_grid(`Cell Line`~`Time Point`) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle(i) + theme_bw() +
                theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggsave(fig3, filename = paste0(plotName, "_boxplot.png"), width = 8, heigh = 6)
        ggsave(fig3, filename = paste0(plotName, "_boxplot.pdf"), width = 8, heigh = 6)
        
}


# write loop to go through each treatment and plot the result

treatmentList <- unique(dat$Treatment)

for (i in treatmentList) {
        
        dat.rx <- filter(dat, Treatment == i & !grepl("Abl|p53|HIF", Target))
        
        g <- ggplot(dat.rx, aes(Target, `Net Shift`, color = Target))
        
        fig4 <- g + geom_point(position = "jitter", aes(shape = factor(Replicate))) + 
                facet_grid(`Cell Line`~`Time Point`) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle(i) + theme_bw() +
                theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        fig4
        ggsave(fig4, filename = paste0(i, "_point.png"), width = 8, height = 6)
        ggsave(fig4, filename = paste0(i, "_point.pdf"), width = 8, height = 6)
        
        fig5 <- g + geom_boxplot() + facet_grid(`Cell Line`~`Time Point`) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle(i) + theme_bw() +
                theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggsave(fig5, filename = paste0(i, "_boxplot.png"), width = 8, heigh = 6)
        ggsave(fig5, filename = paste0(i, "_boxplot.pdf"), width = 8, heigh = 6)
}

# theme for treatments

theme <- theme_bw() + theme(panel.grid = element_blank(),
                            text = element_text(size = 14),
                            axis.title.x=element_blank(),
                            axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
                            axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank())

# Palbociclib

control <- "(-)-Serum"
treatment <- "Palbociclib"
targets <- c("pRb Ser780", "pRb Ser807/11", "pGSK3ß Ser9", "pPDK1 Ser341")

dat.cntl <- filter(dat, Treatment == control & Target %in% targets &
                           `Time Point` == 1)
dat.rx <- filter(dat, Treatment == treatment & Target %in% targets)

dat.rx2 <- rbind(dat.cntl, dat.rx)

dat.rx2.6 <- filter(dat.rx2, `Cell Line` == 6)

g <- ggplot(dat.rx2.6, aes(interaction(factor(`Time Point`),
                Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_palbo_6 <- g + geom_boxplot() + theme + labs(fill = "")

fig_palbo_6

ggsave(fig_palbo_6, filename = paste0('Palbociclib_GBM6.png'), width = 8, height = 6)
ggsave(fig_palbo_6, filename = paste0('Palbociclib_GBM6.pdf'), width = 8, height = 6)

dat.rx2.26 <- filter(dat.rx2, `Cell Line` == 26)

g <- ggplot(dat.rx2.26, aes(interaction(factor(`Time Point`),
                                       Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_palbo_26 <- g + geom_boxplot() + theme + labs(fill = "")

fig_palbo_26

ggsave(fig_palbo_26, filename = paste0('Palbociclib_GBM26.png'), width = 8, height = 6)
ggsave(fig_palbo_26, filename = paste0('Palbociclib_GBM26.pdf'), width = 8, height = 6)

# Erlotinib
control <- "(-)-Serum"
treatment <- "Erlotinib"
targets <- c("p-MAPK","p-Akt_308", "p-Akt_473", "p-S6_235", "p-S6_240", "p-mTOR")

dat.cntl <- filter(dat, Treatment == control & Target %in% targets &
                           `Time Point` == 1)
dat.rx <- filter(dat, Treatment == treatment & Target %in% targets)

dat.rx2 <- rbind(dat.cntl, dat.rx)

dat.rx2.6 <- filter(dat.rx2, `Cell Line` == 6)

g <- ggplot(dat.rx2.6, aes(interaction(factor(`Time Point`),
                Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_erlot_6 <- g + geom_boxplot() + theme + labs(fill = "")

fig_erlot_6

ggsave(fig_erlot_6, filename = paste0('Erlotinib_GBM6.png'), width = 8, height = 6)
ggsave(fig_erlot_6, filename = paste0('Erlotinib_GBM6.pdf'), width = 8, height = 6)

dat.rx2.26 <- filter(dat.rx2, `Cell Line` == 26)

g <- ggplot(dat.rx2.26, aes(interaction(factor(`Time Point`),
                Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_erlot_26 <- g + geom_boxplot() + theme + labs(fill = "")

fig_erlot_26

ggsave(fig_erlot_26, filename = paste0('Erlotinib_GBM26.png'), width = 8, height = 6)
ggsave(fig_erlot_26, filename = paste0('Erlotinib_GBM26.pdf'), width = 8, height = 6)

# GNE-317

control <- "(-)-Serum"
treatment <- "GNE-317"
targets <- c("p-Akt_308", "p-Akt_473", "p-S6_235", "p-S6_240", "p-mTOR", "p-GSK3b", "p-p70S6K")

dat.cntl <- filter(dat, Treatment == control & Target %in% targets &
                           `Time Point` == 1)
dat.rx <- filter(dat, Treatment == treatment & Target %in% targets)

dat.rx2 <- rbind(dat.cntl, dat.rx)

dat.rx2.6 <- filter(dat.rx2, `Cell Line` == 6)

g <- ggplot(dat.rx2.6, aes(interaction(factor(`Time Point`),
                                       Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_gne_6 <- g + geom_boxplot() + theme + labs(fill = "")

fig_gne_6

ggsave(fig_gne_6, filename = paste0('GNE-317_GBM6.png'), width = 8, height = 6)
ggsave(fig_gne_6, filename = paste0('GNE-317_GBM6.pdf'), width = 8, height = 6)

dat.rx2.26 <- filter(dat.rx2, `Cell Line` == 26)

g <- ggplot(dat.rx2.26, aes(interaction(factor(`Time Point`),
                                       Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_gne_26 <- g + geom_boxplot() + theme + labs(fill = "")

fig_gne_26

ggsave(fig_gne_26, filename = paste0('GNE-317_GBM26.png'), width = 8, height = 6)
ggsave(fig_gne_26, filename = paste0('GNE-317_GBM26.pdf'), width = 8, height = 6)


#GDC-0980

control <- "(-)-Serum"
treatment <- "GDC0980"
targets <- c("p-Akt_308", "p-Akt_473", "p-S6_235", "p-S6_240", 
             "p-mTOR", "p-GSK3b", "p-p70S6K")

dat.cntl <- filter(dat, Treatment == control & Target %in% targets &
                           `Time Point` == 1)
dat.rx <- filter(dat, Treatment == treatment & Target %in% targets)

dat.rx2 <- rbind(dat.cntl, dat.rx)

dat.rx2.6 <- filter(dat.rx2, `Cell Line` == 6)

g <- ggplot(dat.rx2.6, aes(interaction(factor(`Time Point`),
        Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_gdc_6 <- g + geom_boxplot() + theme + labs(fill = "")

fig_gdc_6

ggsave(fig_gdc_6, filename = paste0('GDC-0980_GBM6.png'), width = 8, height = 6)
ggsave(fig_gdc_6, filename = paste0('GDC-0980_GBM6.pdf'), width = 8, height = 6)

dat.rx2.26 <- filter(dat.rx2, `Cell Line` == 26)

g <- ggplot(dat.rx2.26, aes(interaction(factor(`Time Point`),
        Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_gdc_26 <- g + geom_boxplot() + theme + labs(fill = "")

fig_gdc_26

ggsave(fig_gdc_26, filename = paste0('GDC-0980_GBM26.png'), width = 8, height = 6)
ggsave(fig_gdc_26, filename = paste0('GDC-0980_GBM26.pdf'), width = 8, height = 6)


# GDC-0980 vs GNE-317

treatment <- c("GDC0980", "GNE-317")

target_treated <- c("p-Akt_308", "p-Akt_473", "p-S6_240", "p-S6_235", "p-GSK3b")

dat.rx2 <- filter(dat, Treatment == treatment & Target %in% target_treated)

g <- ggplot(dat.rx2, aes(interaction(factor(`Time Point`), Treatment, Target),
                         `Net Shift`, fill = factor(`Time Point`)))

fig_gvsg <- g + geom_boxplot() + facet_grid(`Cell Line`~.) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("GDC0980 vs GNE-317") + theme_bw() +
        theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(fill = "Time Point") +
        geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5), linetype = "longdash") +
        annotate("text", x = 2.5, y = 3500, label = "p-Akt Thr 308") +
        annotate("text", x = 6.5, y = 3500, label = "p-Akt Ser 473") +
        annotate("text", x = 10.5, y = 3500, label = "p-GSK-3b") +
        annotate("text", x = 14.5, y = 3500, label = "p-S6 Ser235") +
        annotate("text", x = 18.5, y = 3500, label = "p-S6 Ser240")

fig_gvsg

 ggsave(fig_gvsg, filename = paste0('GDCvsGNE.png'), width = 8, height = 6)
ggsave(fig_gvsg, filename = paste0('GDCvsGNE.pdf'), width = 8, height = 6)

# make some heatmaps
