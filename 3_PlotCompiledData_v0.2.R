library(readr)
library(ggplot2)
library(dplyr)

dat <- read_csv("compiledLabeled.csv")
g <- ggplot(dat, aes(Target, `Net Shift`))
g + geom_point(aes(color = Treatment, shape = factor(Replicate)), 
               position = "jitter", alpha = 1/2, size = 3) +
        facet_grid(`Cell Line` ~ `Time Point`) # + theme_bw() +
        #theme(panel.grid = element_blank(), axis.title.x=element_blank()) +
        #theme(axis.text.x = element_text(angle = 45, hjust = 1))

dat.tar <- filter(dat, Target == "p-Akt_473")
g <- ggplot(dat.tar, aes(Treatment, `Net Shift`, color = Treatment))
g + geom_point(position = "jitter", aes(shape = factor(Replicate))) + facet_grid(`Cell Line`~`Time Point`)
g + geom_boxplot() + facet_grid(`Cell Line`~`Time Point`)
