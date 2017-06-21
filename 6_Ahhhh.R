setwd("D:/Box Sync/Data/")
library(tidyverse)
library(ggthemes)

dat <- read_csv("compiledLabeled.csv")
dat <- select(dat, -Shift.1, -Shift.2, - Group)
dat <- filter(dat, !grepl("p53|Abl", Target))
dat$Id <- dat$Ring %% 4

head(dat)

dat.dcast <<- dcast(dat,
                    Treatment + `Time Point` + `Cell Line` + Id ~ Target,
                    value.var = "Net Shift",
                    fun.aggregate = mean)
head(dat.dcast)
dat.dcast$`Cell Line` <- factor(dat.dcast$`Cell Line`)

plot <- ggplot(dat.dcast, 
               aes(x = `pRb Ser807/11`, 
                   y = `pS6 Ser235/6`, 
                   color = Treatment,
                   shape = interaction(`Cell Line`, `Time Point`)))

plot + geom_point(size = 6) + theme_few()



plot(dat.dcast$`hydroxy-HIF Pro564` ~ dat.dcast$`pAkt Ser473`)
