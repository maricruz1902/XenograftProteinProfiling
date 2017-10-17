# Import and scale data ---------------------------------------------------

rm( list = ls() )
library(tidyverse)
library(reshape2)

setwd("D:/Box Sync/Data/")
dat <- read_csv("compiledLabeled.csv")
dat <- dat %>% mutate(n = Ring %% 4,
                      LogTransformed = log(`Net Shift`))
# dat <- filter(dat, !(Target %in% c("pp53 Ser15", "pc-Abl Tyr204",
#                                    "h-HIF1a Pro564")) & Replicate != 1)

datScaled <- dat %>%
        group_by(Target) %>%
        mutate(NormLog = scale(LogTransformed),
               Normalized = scale(`Net Shift`))

cv <- function(x){sd(x)/mean(x) * 100}
datSum <- datScaled %>%
        group_by(Target, `Cell Line`, Treatment, `Time Point`) %>%
        summarize_at(vars("Net Shift", "Normalized", 
                          "NormLog", "LogTransformed"),
                     funs(mean, sd, cv))


# Assess data normality ---------------------------------------------------

shapiro.test(log(datScaled$`Net Shift`))
shapiro.test(datScaled$Normalized)
shapiro.test(datScaled$NormLog)

par(mfrow=c(1,3))
hist(datScaled$NormLog, breaks = 100)
hist(datScaled$Normalized, breaks = 100)
hist(datScaled$LogTransformed, breaks = 100)

par(mfrow=c(2,2))
x = datScaled$Normalized
qqnorm(x, main = "Normalized"); qqline(x)
y = datScaled$LogTransformed
qqnorm(y, main = "Log Transformed"); qqline(y)
z = datScaled$NormLog
qqnorm(z, main = "Log Transformed & Normalized"); qqline(z)
qqnorm(rnorm(1000), main = "Normal Distribution"); qqline(z)


# Cast data for multivariate analysis -------------------------------------

targetList <- unique(dat$Target)

cast <- dcast(data = datScaled, Treatment + `Cell Line` + `Time Point` + 
                      Replicate + n ~ Target, mean, value.var = "NormLog")

man <- manova(cbind(`pc-Abl Tyr204`, `pPDK1 Ser341`,
                    `pGSK3b Ser9`, `pp70S6K Thr389`,
                    `pS6 Ser235/6`, `pRb Ser780`,
                    `pAkt Thr308`, `pS6 Ser240/4`,
                    `pRb Ser807/11`, `pmTOR Ser2448`,
                    `pMAPK Thr202/Tyr204`, `pAkt Ser473`,
                    `pp53 Ser15`, `pSrc Tyr416`,
                    `h-HIF1a Pro564`) ~ `Time Point` + `Cell Line` + Treatment,
              data = cast)
library(RVAideMemoire)
data(iris)

summary(manova(as.matrix(iris[,1:4])~iris$Species))

adonis(iris[,1:4]~Species, method = "euclidian")

scaleCols <- scale(select(cast, -Treatment, -`Cell Line`,
                           -`Time Point`, -Replicate, -n), center = T)
qqnorm(scaleCols); qqline(scaleCols)
labelCols <- select(cast, Treatment, `Cell Line`, `Time Point`, Replicate, n)
castScaled <- cbind(scaleCols, labelCols)

pairsList <- t(combn(targetList, 2))


cellLine <- unique(dat$`Cell Line`)
treatmentList <- unique(dat$Treatment)
timePoint <- unique(dat$`Time Point`)
combo <- expand.grid(Treatment = treatmentList, 
                     `Cell Line` = cellLine,
                     `Time Point` = timePoint)

combo2 <- expand.grid(Treatment = treatmentList, 
                      `Cell Line` = cellLine,
                      `Time Point` = timePoint,
                      Target = targetList)

# test <- list()
# for (i in seq_len(nrow(pairsList))){
#         pair <- pairsList[i,]
#         prot1 <- pair[1]
#         prot2 <- pair[2]
#         pairDat <- filter(dat, Target %in% pair)
#         temp <- list()
#         for(j in seq_len(nrow(labelCols))){
#                 txt <- labelCols[j, 1]
#                 cl <- labelCols[j, 2]
#                 tp <- labelCols[j, 3]
#                 p1 <- filter(pairDat, Treatment == txt & `Cell Line` == cl &
#                                     `Time Point` == tp & Target == prot1) %>%
#                         select(`Net Shift`) %>% unlist()
#                 p2 <- filter(pairDat, Treatment == txt & `Cell Line` == cl &
#                                      `Time Point` == tp & Target == prot2) %>%
#                         select(`Net Shift`) %>% unlist()
#                 temp[[j]] <- tidy(t.test(p1, p2))
#                 print(tidy(t.test(p1, p2)))
#         }
#         test[[pair]] <- temp
# }


pairMat <- matrix(nrow = nrow(combo), ncol = nrow(pairsList)*4)



# for(i in seq_len(nrow(combo))){
#         txt <- combo[i, 1]
#         cl <- combo[i, 2]
#         tp <- combo[i, 3]
#         pairDat <- filter(castScaled, Treatment == txt & `Cell Line` == cl &
#                                   `Time Point` == tp)
#         print(paste0("i = ", i, "/", nrow(combo)))
#         for(j in seq_len(nrow(pairsList))){
#                 print(paste0(j, "/", nrow(pairsList)))
#                 prot1 <- pairsList[j, 1]
#                 prot2 <- pairsList[j, 2]
#                 # print(paste(prot1, prot2))
#                 p1 <- pairDat %>% select(prot1) %>% unlist()
#                 p2 <- pairDat %>% select(prot2) %>% unlist()
#                 cNum <- (j-1)*4
#                 pairMat[i,cNum+1] <- mean(p1, na.rm = TRUE) - 
#                         mean(p2, na.rm = TRUE)
#                 pairMat[i,cNum+2] <- sqrt(sd(p1, na.rm = TRUE)^2 + 
#                         sd(p2, na.rm = TRUE)^2)
#                 pairMat[i,cNum+3] <- length(p1[!is.na(p1)])
#                 pairMat[i,cNum+4] <- length(p2[!is.na(p2)])
#         }
# }

a <- filter(castScaled, Treatment == "(-)-Serum" & `Cell Line` == "GBM 26" &
                    `Time Point` == "1 h") %>% select("pPDK1 Ser341")

b <- filter(datScaled, Treatment == "(-)-Serum" & `Cell Line` == "GBM 26" &
                    `Time Point` == "1 h" & Target == "pPDK1 Ser341")
b$Normalized

library(dplyr)
set.seed(12)
n = 1000
df <- data.frame(ID = sample(c("A","B","C","D","E"), size=n, replace=TRUE),
                 score = runif(n, 0, 10))

scaledByID <- 
        df %>%
        group_by(ID) %>%
        mutate(scaledScore = scale(score))

notScaledByID <- 
        df %>%
        mutate(scaledScore = scale(score))

scaledByID$scaledScore == notScaledByID$scale

file1 <- read_csv("file1.csv")
file1$Experiment <- "Experiment Name"

dat.all <- rbind(file1, file2, file2)