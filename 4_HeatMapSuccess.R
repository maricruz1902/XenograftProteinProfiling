setwd("D:/Box Sync/Data/")
library(tidyverse)
library(RColorBrewer)
library(reshape2)

dat <- read_csv("compiledLabeled.csv")
dat <- select(dat, -Experiment, -Shift.1, -Shift.2, -Experiment, -Replicate)

dat.1 <- filter(dat, Treatment != "(+)-Serum" & Treatment != "DMSO")
dat.2 <- filter(dat.1, !grepl("Abl|p53|HIF", Target))

dat.6 <- filter(dat.2, `Cell Line` == 6)
dat.26 <- filter(dat.2, `Cell Line` == 26)

casting.1 <- dcast(dat.1, Treatment + `Cell Line` + `Time Point` ~ 
                        Target, mean, value.var = "Net Shift")

casting.2 <- dcast(dat.2, Treatment + `Cell Line` + `Time Point` ~ 
                        Target, mean, value.var = "Net Shift")

casting.6 <- dcast(dat.6, Treatment + `Cell Line` + `Time Point` ~ 
                           Target, mean, value.var = "Net Shift")

casting.26 <- dcast(dat.26, Treatment + `Cell Line` + `Time Point` ~ 
                           Target, mean, value.var = "Net Shift")

rownames(casting.1) = paste(casting.1$Treatment, casting.1$`Cell Line`, 
                        casting.1$`Time Point`, sep = " ")

rownames(casting.2) = paste(casting.2$Treatment, casting.2$`Cell Line`, 
                            casting.2$`Time Point`, sep = " ")

rownames(casting.6) = paste(casting.6$Treatment, casting.6$`Cell Line`, 
                            casting.6$`Time Point`, sep = " ")

rownames(casting.26) = paste(casting.26$Treatment, casting.26$`Cell Line`, 
                            casting.26$`Time Point`, sep = " ")

fig.1 <- ggplot(dat, aes(interaction(Treatment, factor(`Cell Line`), 
                        factor(`Time Point`)), Target, fill = `Net Shift`)) + 
                geom_tile() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))

casting.1[,c(1,2,3)] <- NULL
casting.2[,c(1,2,3)] <- NULL
casting.6[,c(1,2,3)] <- NULL
casting.26[,c(1,2,3)] <- NULL

casting.1.m <- data.matrix(casting.1)
casting.2.m <- data.matrix(casting.2)
casting.6.m <- data.matrix(casting.6)
casting.26.m <- data.matrix(casting.26)

hmcol <-colorRampPalette(c("darkblue", "white", "darkred"))(256)

heatmap(casting.1.m, scale = "col", col = hmcol)
dev.print(png,"heatmap_all.png",width=8,height=6,units="in",res=300)
dev.off()

heatmap(casting.1.m, scale = "col", col = hmcol)
dev.print(pdf,"heatmap_all.pdf",width=8,height=6)
dev.off()

heatmap(casting.2.m, scale = "col", col = hmcol)
dev.print(png,"heatmap_all_2.png",width=8,height=6,units="in",res=300)
dev.off()

heatmap(casting.2.m, scale = "col", col = hmcol)
dev.print(pdf,"heatmap_all_2.pdf",width=8,height=6)
dev.off()

heatmap(casting.6.m, scale = "col", col = hmcol)
dev.print(png,"heatmap_GBM6.png",width=8,height=6,units="in",res=300)
dev.off()

heatmap(casting.6.m, scale = "col", col = hmcol)
dev.print(pdf,"heatmap_GBM6.pdf",width=8,height=6)
dev.off()

heatmap(casting.26.m, scale = "col", col = hmcol)
dev.print(png,"heatmap_GBM26.png",width=8,height=6,units="in",res=300)
dev.off()

heatmap(casting.26.m, scale = "col", col = hmcol)
dev.print(pdf,"heatmap_GBM26.pdf",width=8,height=6)
dev.off()


