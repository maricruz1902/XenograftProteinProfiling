library(readr)
library(ggplot2)
library(dplyr)

dat <- read_csv("compiledLabeled.csv")
g <- ggplot(dat, aes(Target, `Net Shift`))

fig <- g + geom_point(aes(color = Treatment, shape = factor(Replicate)), 
               position = "jitter", alpha = 1/2, size = 3) +
        facet_grid(`Cell Line` ~ `Time Point`) + theme_bw() +
        theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig
ggsave(fig, filename = "everything_point.png", width = 8, height = 6)
ggsave(fig, filename = "everything_point.pdf", width = 8, height = 6)


fig1 <- g + geom_boxplot(aes(color = Treatment)) +
        facet_grid(`Cell Line` ~ `Time Point`) + theme_bw() +
        theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig1

ggsave(fig1, filename = "everything_boxplot.png", width = 8, height = 6)
ggsave(fig1, filename = "everything_boxplot.pdf", width = 8, height = 6)

# write loop to go through each target and plot the result

targetList <- unique(dat$Target)

for(i in targetList) {
        dat.tar <- filter(dat, Target == i)
        
        g <- ggplot(dat.tar, aes(Treatment, `Net Shift`, color = Treatment))
        
        fig2 <- g + geom_point(position = "jitter", 
                               aes(shape = factor(Replicate))) + 
                facet_grid(`Cell Line`~`Time Point`) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle(i) + theme_bw() +
                theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggsave(fig2, filename = paste0(i, "_point.png"), width = 8, height = 6)
        ggsave(fig2, filename = paste0(i, "_point.pdf"), width = 8, height = 6)
        
        fig3 <- g + geom_boxplot() + facet_grid(`Cell Line`~`Time Point`) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                ggtitle(i) + theme_bw() +
                theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggsave(fig3, filename = paste0(i, "_boxplot.png"), width = 8, heigh = 6)
        ggsave(fig3, filename = paste0(i, "_boxplot.pdf"), width = 8, heigh = 6)
        
}


# write loop to go through each treatment and plot the result

treatmentList <- unique(dat$Treatment)

for (i in treatmentList) {
        
        dat.rx <- filter(dat, Treatment == i)
        
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

# Palbociclib

treatment <- c("Palbociclib", "(-)-Serum")

target_treated <- c("p-Rb_780", "p-Rb_807", "p-GSK3b")

dat.rx2 <- filter(dat, Treatment == treatment & Target %in% target_treated)
        
g <- ggplot(dat.rx2, aes(interaction(factor(`Time Point`), Treatment, Target),
                         `Net Shift`, fill = Target))

fig_palbo <- g + geom_boxplot() + facet_grid(`Cell Line`~.) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(treatment) + theme_bw() +
        theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(fill = "Target") +
        geom_vline(xintercept = c(4.5, 8.5), linetype = "longdash") +
        annotate("text", x = 2.5, y = 7500, label = "p-GSK-3b") +
        annotate("text", x = 6.5, y = 7500, label = "p-Rb Ser780") +
        annotate("text", x = 10.5, y = 7500, label = "p-Rb Ser807")

fig_palbo

ggsave(fig_palbo, filename = paste0('Palbociclib.png'), width = 8, height = 6)
ggsave(fig_palbo, filename = paste0('Palbociclib.pdf'), width = 8, height = 6)

# Erlotinib

treatment <- c("Erlotinib", "(-)-Serum")

target_treated <- c("p-Akt_473", "p-PDK1", "p-S6_235")

dat.rx2 <- filter(dat, Treatment == treatment & Target %in% target_treated)

g <- ggplot(dat.rx2, aes(interaction(factor(`Time Point`), Treatment, Target),
                         `Net Shift`, fill = Target))

fig_erlot <- g + geom_boxplot() + facet_grid(`Cell Line`~.) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(treatment) + theme_bw() +
        theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(fill = "Target") +
        geom_vline(xintercept = c(4.5, 8.5), linetype = "longdash") +
        annotate("text", x = 2.5, y = 3500, label = "p-Akt Ser473") +
        annotate("text", x = 6.5, y = 3500, label = "p-PDK1") +
        annotate("text", x = 10.5, y = 3500, label = "p-GSK-3b")

fig_erlot

ggsave(fig_erlot, filename = paste0('Erlotinib.png'), width = 8, height = 6)
ggsave(fig_erlot, filename = paste0('Erlotinib.pdf'), width = 8, height = 6)


# GNE-317

treatment <- c("GNE-317", "(-)-Serum")

target_treated <- c("p-Akt_308", "p-Akt_473", "p-S6_240", "p-S6_235", "p-p70S6K")

dat.rx2 <- filter(dat, Treatment == treatment & Target %in% target_treated)

g <- ggplot(dat.rx2, aes(interaction(factor(`Time Point`), Treatment, Target),
                         `Net Shift`, fill = factor(`Time Point`)))

fig_gne <- g + geom_boxplot() + facet_grid(`Cell Line`~.) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(treatment) + theme_bw() +
        theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(fill = "Time Point") +
        geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5), linetype = "longdash") +
        annotate("text", x = 2.5, y = 3250, label = "p-Akt Thr 308") +
        annotate("text", x = 6.5, y = 3250, label = "p-Akt Ser 473") +
        annotate("text", x = 10.5, y = 3250, label = "p-p70S6K") +
        annotate("text", x = 14.5, y = 3250, label = "p-S6 Ser235") +
        annotate("text", x = 18.5, y = 3250, label = "p-S6 Ser240")
        
fig_gne

ggsave(fig_gne, filename = paste0('GNE-317.png'), width = 8, height = 6)
ggsave(fig_gne, filename = paste0('GNE-317.pdf'), width = 8, height = 6)


#GDC-0980

treatment <- c("GDC0980", "(-)-Serum")

target_treated <- c("p-Akt_308", "p-Akt_473", "p-S6_240", "p-S6_235", "p-p70S6K")

dat.rx2 <- filter(dat, Treatment == treatment & Target %in% target_treated)

g <- ggplot(dat.rx2, aes(interaction(factor(`Time Point`), Treatment, Target),
                         `Net Shift`, fill = factor(`Time Point`)))

fig_gdc <- g + geom_boxplot() + facet_grid(`Cell Line`~.) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(treatment) + theme_bw() +
        theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(fill = "Time Point") +
        geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5), linetype = "longdash") +
        annotate("text", x = 2.5, y = 3500, label = "p-Akt Thr 308") +
        annotate("text", x = 6.5, y = 3500, label = "p-Akt Ser 473") +
        annotate("text", x = 10.5, y = 3500, label = "p-p70S6K") +
        annotate("text", x = 14.5, y = 3500, label = "p-S6 Ser235") +
        annotate("text", x = 18.5, y = 3500, label = "p-S6 Ser240")

fig_gdc

ggsave(fig_gdc, filename = paste0('GDC-0980.png'), width = 8, height = 6)
ggsave(fig_gdc, filename = paste0('GDC-0980.pdf'), width = 8, height = 6)


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
