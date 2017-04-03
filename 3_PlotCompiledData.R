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

dat$Treatment <- factor(dat$Treatment, level = c("(+)-Serum", "(-)-Serum", 
                "DMSO", "Erlotinib", "GNE-317", "Apitolisib", "Palbociclib"))


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
                            text = element_text(size = 22),
                            axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank())

# Palbociclib

control <- "DMSO"
treatment <- "Palbociclib"
targets <- c("pRb Ser780", "pRb Ser807/11", "pGSK3ß Ser9", "pPDK1 Ser341")

dat.cntl <- filter(dat, Treatment == control & Target %in% targets &
                           `Time Point` == 1)

dat.all <- rbind(filter(dat, Treatment == treatment), 
        filter(dat, Treatment == control & `Time Point` == 1))
dat.all$`Cell Line` <- factor(dat.all$`Cell Line`, labels = c("GBM 6", "GBM 26"))

g.all <- ggplot(dat.all,  aes(interaction(factor(`Time Point`),
                Treatment, Target), `Net Shift`, fill = Target)) +
        geom_boxplot() + theme + labs(fill = "") + 
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 15)) +
        xlab("Treatment Time") + facet_grid(`Cell Line`~.) +
        ggtitle("Treatment: Palbociclib") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

g.all
ggsave(g.all, filename = "Palbociclib.png", width = 12, height = 6)
ggsave(g.all, filename = "Palbociclib.pdf", width = 12, height = 6)

dat.rx <- filter(dat, Treatment == treatment & Target %in% targets)

dat.rx2 <- rbind(dat.cntl, dat.rx)

dat.rx2.6 <- filter(dat.rx2, `Cell Line` == 6)

g <- ggplot(dat.rx2.6, aes(interaction(factor(`Time Point`),
        Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_palbo_6 <- g + geom_boxplot() + theme + labs(fill = "") + 
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 4)) +
        xlab("Treatment Time") + theme(legend.position = c(0.2, 0.8)) +
        ggtitle("Treatment: Palbociclib, GBM 6")

fig_palbo_6

ggsave(fig_palbo_6, filename = "Palbociclib_GBM6.png", width = 8, height = 6)
ggsave(fig_palbo_6, filename = "Palbociclib_GBM6.pdf", width = 8, height = 6)

dat.rx2.26 <- filter(dat.rx2, `Cell Line` == 26)

g <- ggplot(dat.rx2.26, aes(interaction(factor(`Time Point`),
        Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_palbo_26 <- g + geom_boxplot() + theme + labs(fill = "") + 
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 3)) +
        xlab("Treatment Time") + theme(legend.position = c(0.8, 0.8)) +
        ggtitle("Treatment: Palbociclib, GBM 26")

fig_palbo_26

ggsave(fig_palbo_26, filename = paste0('Palbociclib_GBM26.png'), width = 8, height = 6)
ggsave(fig_palbo_26, filename = paste0('Palbociclib_GBM26.pdf'), width = 8, height = 6)

# Erlotinib

control <- "DMSO"
treatment <- "Erlotinib"
targets <- c("pMAPK Thr202/Tyr204","pAkt Thr308", "pAkt Ser473", 
             "pS6 Ser235/6", "pS6 Ser240/4", "pmTOR Ser2448")

dat.cntl <- filter(dat, Treatment == control & Target %in% targets &
                           `Time Point` == 1)

dat.all <- rbind(filter(dat, Treatment == treatment), 
                 filter(dat, Treatment == control & `Time Point` == 1))
dat.all$`Cell Line` <- factor(dat.all$`Cell Line`, labels = c("GBM 6", "GBM 26"))

g.all <- ggplot(dat.all,  aes(interaction(factor(`Time Point`),
                Treatment, Target), `Net Shift`, fill = Target)) +
        geom_boxplot() + theme + labs(fill = "") + 
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 15)) +
        xlab("Treatment Time") + facet_grid(`Cell Line`~.) +
        ggtitle("Treatment: Erlotinib") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

g.all

ggsave(g.all, filename = "Erlotinib.png", width = 12, height = 6)
ggsave(g.all, filename = "Erlotinib.pdf", width = 12, height = 6)

dat.rx <- filter(dat, Treatment == treatment & Target %in% targets)

dat.rx2 <- rbind(dat.cntl, dat.rx)

dat.rx2.6 <- filter(dat.rx2, `Cell Line` == 6)

g <- ggplot(dat.rx2.6, aes(interaction(factor(`Time Point`),
                Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_erlot_6 <- g + geom_boxplot() + theme + labs(fill = "") +
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 6)) +
        xlab("Treatment Time") + theme(legend.position = c(0.3, 0.7)) +
        ggtitle("Treatment: Erlotinib, GBM 6") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig_erlot_6

ggsave(fig_erlot_6, filename = paste0('Erlotinib_GBM6.png'), width = 8, height = 6)
ggsave(fig_erlot_6, filename = paste0('Erlotinib_GBM6.pdf'), width = 8, height = 6)

dat.rx2.26 <- filter(dat.rx2, `Cell Line` == 26)

g <- ggplot(dat.rx2.26, aes(interaction(factor(`Time Point`),
                Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_erlot_26 <- g + geom_boxplot() + theme + labs(fill = "") +
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 6)) +
        xlab("Treatment Time") +  theme(legend.position = c(0.3, 0.7)) +
        ggtitle("Treatment: Erlotinib, GBM 26") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig_erlot_26

ggsave(fig_erlot_26, filename = paste0('Erlotinib_GBM26.png'), width = 8, height = 6)
ggsave(fig_erlot_26, filename = paste0('Erlotinib_GBM26.pdf'), width = 8, height = 6)

# GNE-317

control <- "DMSO"
treatment <- "GNE-317"
targets <- c("pAkt Ser473", "pS6 Ser235/6", "pS6 Ser240/4", "pmTOR Ser2448", 
             "pGSK3β Ser9", "pp70S6K Thr389")

dat.cntl <- filter(dat, Treatment == control & Target %in% targets &
                           `Time Point` == 1)
dat.rx <- filter(dat, Treatment == treatment & Target %in% targets)

dat.all <- rbind(filter(dat, Treatment == treatment), 
                 filter(dat, Treatment == control & `Time Point` == 1))
dat.all$`Cell Line` <- factor(dat.all$`Cell Line`, labels = c("GBM 6", "GBM 26"))

g.all <- ggplot(dat.all,  aes(interaction(factor(`Time Point`),
                                          Treatment, Target), `Net Shift`, fill = Target)) +
        geom_boxplot() + theme + labs(fill = "") + 
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 15)) +
        xlab("Treatment Time") + 
        facet_grid(`Cell Line`~.) +
        ggtitle("Treatment: GNE-317") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

g.all

ggsave(g.all, filename = "GNE-317.png", width = 12, height = 6)
ggsave(g.all, filename = "GNE-317.pdf", width = 12, height = 6)

dat.rx2 <- rbind(dat.cntl, dat.rx)

dat.rx2.6 <- filter(dat.rx2, `Cell Line` == 6)

g <- ggplot(dat.rx2.6, aes(interaction(factor(`Time Point`),
        Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_gne_6 <- g + geom_boxplot() + theme + labs(fill = "") +
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 6)) +
        xlab("Treatment Time") + theme(legend.position = c(0.2, 0.8)) +
        ggtitle("Treatment: GNE-317, GBM 6") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig_gne_6

ggsave(fig_gne_6, filename = paste0('GNE-317_GBM6.png'), width = 8, height = 6)
ggsave(fig_gne_6, filename = paste0('GNE-317_GBM6.pdf'), width = 8, height = 6)

dat.rx2.26 <- filter(dat.rx2, `Cell Line` == 26)

g <- ggplot(dat.rx2.26, aes(interaction(factor(`Time Point`),
        Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_gne_26 <- g + geom_boxplot() + theme + labs(fill = "") +
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 6)) +
        xlab("Treatment Time") + theme(legend.position = c(0.2, 0.8)) +
        ggtitle("Treatment: GNE-317, GBM 26") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig_gne_26

ggsave(fig_gne_26, filename = paste0('GNE-317_GBM26.png'), width = 8, height = 6)
ggsave(fig_gne_26, filename = paste0('GNE-317_GBM26.pdf'), width = 8, height = 6)


#Apitolisib

control <- "DMSO"
treatment <- "Apitolisib"
targets <- c("pAkt Ser473", "pS6 Ser235/6", "pS6 Ser240/4", "pmTOR Ser2448", 
             "pGSK3β Ser9", "pp70S6K Thr389")

dat.cntl <- filter(dat, Treatment == control & Target %in% targets &
                           `Time Point` == 1)
dat.rx <- filter(dat, Treatment == treatment & Target %in% targets)

dat.all <- rbind(filter(dat, Treatment == treatment), 
                 filter(dat, Treatment == control & `Time Point` == 1))
dat.all$`Cell Line` <- factor(dat.all$`Cell Line`, labels = c("GBM 6", "GBM 26"))

g.all <- ggplot(dat.all,  aes(interaction(factor(`Time Point`),
                                          Treatment, Target), `Net Shift`, fill = Target)) +
        geom_boxplot() + theme + labs(fill = "") + 
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 15)) +
        xlab("Treatment Time") + 
        facet_grid(`Cell Line`~.) +
        ggtitle("Treatment: Apitolisib") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

g.all

ggsave(g.all, filename = "Apitolisib.png", width = 12, height = 6)
ggsave(g.all, filename = "Apitolisib.pdf", width = 12, height = 6)

dat.rx2 <- rbind(dat.cntl, dat.rx)

dat.rx2.6 <- filter(dat.rx2, `Cell Line` == 6)

g <- ggplot(dat.rx2.6, aes(interaction(factor(`Time Point`),
        Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_gdc_6 <- g + geom_boxplot() + theme + labs(fill = "") +
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 6)) +
        xlab("Treatment Time") + theme(legend.position = c(0.2, 0.8)) +
        ggtitle("Treatment: Apitolisib, GBM 6") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig_gdc_6

ggsave(fig_gdc_6, filename = paste0('Apitolisib_GBM6.png'), width = 8, height = 6)
ggsave(fig_gdc_6, filename = paste0('Apitolisib_GBM6.pdf'), width = 8, height = 6)

dat.rx2.26 <- filter(dat.rx2, `Cell Line` == 26)

g <- ggplot(dat.rx2.26, aes(interaction(factor(`Time Point`),
        Treatment, Target, factor(`Cell Line`)), `Net Shift`, fill = Target))

fig_gdc_26 <- g + geom_boxplot() + theme + labs(fill = "") +
        scale_x_discrete(labels = rep(c("Control", "1 h", "24 h"), 6)) +
        xlab("Treatment Time") + theme(legend.position = c(0.2, 0.8)) +
        ggtitle("Treatment: Apitolisib, GBM 26") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig_gdc_26

ggsave(fig_gdc_26, filename = paste0('Apitolisib_GBM26.png'), width = 8, height = 6)
ggsave(fig_gdc_26, filename = paste0('Apitolisib_GBM26.pdf'), width = 8, height = 6)


# Apitolisib vs GNE-317

treatment <- c("Apitolisib", "GNE-317")
target_treated <- c("pAkt Ser473", "pS6 Ser235/6", "pS6 Ser240/4", "pmTOR Ser2448", 
             "pGSK3β Ser9", "pp70S6K Thr389")

dat.rx2 <- filter(dat, Treatment == treatment & Target %in% target_treated)

g <- ggplot(dat.rx2, aes(interaction(factor(`Time Point`), Treatment, Target),
                         `Net Shift`, fill = Target))

fig_gvsg <- g + geom_boxplot() + facet_grid(.~`Cell Line`) + theme + labs(fill = "") +
        scale_x_discrete(labels = rep(c("GNE-317 1h", "GNE-317 24h", 
                "Apitolisib 1h", "Apitolisib 24h"), 6)) +
        xlab("Treatment Time") + ggtitle("Treatment: Apitolisib vs GBM 26") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

fig_gvsg

fig_gvsg.2 <- g + geom_boxplot() + facet_grid(`Time Point`~`Cell Line`) + theme + labs(fill = "") +
        scale_x_discrete(labels = rep(c("GNE-317","Apitolisib"), 6)) +
        xlab("Treatment Time") + ggtitle("Treatment: Apitolisib vs GBM 26") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fig_gvsg.2

ggsave(fig_gvsg, filename = paste0('ApitolisibvsGNE317.png'), width = 8, height = 6)
ggsave(fig_gvsg, filename = paste0('ApitolisibvsGNE317.pdf'), width = 8, height = 6)

# make some heatmaps
