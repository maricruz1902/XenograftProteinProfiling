rm(list=ls())
setwd("D:/Box Sync/Data/")
library(tidyverse)
library(RColorBrewer)
library(reshape2)

dat <- read_csv("compiledLabeled.csv")
# dat.1 <- select(dat, -Experiment, -Shift.1, -Shift.2, -Experiment, -Replicate)

dat.6 <- filter(dat, `Cell.Line` == 6)
dat.10 <- filter(dat, `Cell.Line` == 10)

casting.1.alt <- dcast(dat.1, Target ~ Treatment + `Cell.Line` + `Time Point`, 
                       mean, value.var = "Net Shift")

casting.1 <- dcast(dat.1, Treatment + `Cell.Line` + `Time Point` ~ 
                           Target, mean, value.var = "Net Shift")
casting.6 <- dcast(dat.6, Treatment + `Cell.Line` + `Time Point` ~ 
                           Target, mean, value.var = "Net Shift")

casting.10 <- dcast(dat.10, Treatment + `Cell.Line` + `Time Point` ~ 
                            Target, mean, value.var = "Net Shift")

casting.1 <- na.omit(casting.1)
casting.6 <- na.omit(casting.6)
casting.10 <- na.omit(casting.10)
library("reshape2")
names(casting.1)[3]<-paste("Time.Point")

casting.1$Cell.Line <- as.character(casting.1$Cell.Line)
casting.1$Time.Point <- as.character(casting.1$Time.Point)
casting.1$Treatment <- as.character(casting.1$Treatment)

df <- casting.1
df$Treatment_ID <- paste(df$Treatment,df$Cell.Line, df$Time.Point)
df <- df[ -c(1, 3, 2) ]
df2 <- df[,c(16,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)]
casting.1 <- df2

GBM_all <- casting.1
GBM6 <- casting.6
GBM10 <- casting.10

names(GBM_all)

# MANOVA test
names(GBM_all) <- sub(" ", ".", names(GBM_all))
names(GBM_all) <- sub("/", ".", names(GBM_all))
names(GBM_all)[6]<-paste("pGSK3b.Ser9")
names(GBM_all)[2]<-paste("hHIF.Pro564")
names(GBM_all)[5]<-paste("pcAbl.Tyr204")
names(GBM_all)



Treatment_num <- 1:nrow(GBM_all)
GBM_all <- cbind(GBM_all, Treatment_num)

manova <- manova(cbind(hHIF.Pro564, pAkt.Ser473, pAkt.Thr308, pGSK3b.Ser9, pMAPK.Thr202.Tyr204, 
             pPDK1.Ser341, pRb.Ser780, pRb.Ser807.11,pS6.Ser235.6, pS6.Ser240.4, 
             pSrc.Tyr416, pcAbl.Tyr204, pmTOR.Ser2448, pp53.Ser15, pp70S6K.Thr389) 
       ~ Treatment_num, 
       data = GBM_all)

summary(manova)

summary.aov(manova)


heatmap(manova$residuals, scale = "none")

manova2 <- manova(cbind(hHIF.Pro564, pAkt.Ser473, pAkt.Thr308, pGSK3b.Ser9, pMAPK.Thr202.Tyr204, 
                       pPDK1.Ser341, pRb.Ser780, pRb.Ser807.11,pS6.Ser235.6, pS6.Ser240.4, 
                       pSrc.Tyr416, pcAbl.Tyr204, pmTOR.Ser2448, pp53.Ser15, pp70S6K.Thr389) 
                 ~ Treatment_num, 
                 data = GBM_all,
                 subset = Treatment_ID, c("pRb.Ser807.11" %in% "pAkt.Ser473"))
summary(manova2)

