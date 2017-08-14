setwd("D:/Box Sync/Data/")
library(tidyverse)
library(RColorBrewer)
library(reshape2)

dat <- read_csv("compiledLabeled.csv")

dat.1 <- filter(dat, !grepl("p53|Abl", Target))

dat.6 <- filter(dat.1, `Cell Line` == "GBM 6")
dat.26 <- filter(dat.1, `Cell Line` == "GBM 26")

ggplot(dat.1, aes(x = Target, 
                y = interaction(Treatment, `Time Point`, `Cell Line`), 
                fill = `Net Shift`)) + geom_tile()

casting.1 <- dcast(dat.1, Treatment + `Cell Line` + `Time Point` ~ 
                        Target, mean, value.var = "Net Shift")

casting.6 <- dcast(dat.6, Treatment + `Cell Line` + `Time Point` ~ 
                           Target, mean, value.var = "Net Shift")

casting.26 <- dcast(dat.26, Treatment + `Cell Line` + `Time Point` ~ 
                           Target, mean, value.var = "Net Shift")

rownames(casting.1) = paste(casting.1$Treatment, casting.1$`Cell Line`, 
                        casting.1$`Time Point`, sep = " ")

rownames(casting.6) = paste(casting.6$Treatment, casting.6$`Cell Line`, 
                            casting.6$`Time Point`, sep = " ")

rownames(casting.26) = paste(casting.26$Treatment, casting.26$`Cell Line`, 
                            casting.26$`Time Point`, sep = " ")

fig.1 <- ggplot(dat, 
                aes(x = interaction(Treatment, `Cell Line`, `Time Point`), 
                    y = Target, 
                    fill = `Net Shift`)) + 
                geom_tile() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1))
fig.1


casting.1[,c(1,2,3)] <- NULL
casting.6[,c(1,2,3)] <- NULL
casting.26[,c(1,2,3)] <- NULL

casting.1.m <- data.matrix(casting.1)
casting.6.m <- data.matrix(casting.6)
casting.26.m <- data.matrix(casting.26)

hmcol <- colorRampPalette(c("#3a5387", "#e8e8e8", "#b25752"))(256)
# png('test.png', width = 1200, height = 800, res = 100)
heatmap(casting.1.m, scale = "col", col = hmcol, margins = c(10, 10))
# dev.off()

heatmap(cor(casting.1.m), scale = "none", col = hmcol, margins = c(12,8))

heatmap(casting.6.m, scale = "col", col = hmcol, margins = c(10, 10))
heatmap(cor(casting.6.m), scale = "none", col = hmcol, margins = c(12,8))

heatmap(casting.26.m, scale = "col", col = hmcol, margins = c(10, 10))
heatmap(cor(casting.26.m), scale = "none", col = hmcol, margins = c(12,8))