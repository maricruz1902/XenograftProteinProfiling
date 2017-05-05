setwd("D:/Box Sync/Data/")
library(tidyverse)
library(reshape2)
library(Hmisc)
library(corrplot)
library(igraph)


dat <- read_csv("compiledLabeled.csv")
dat <- select(dat, -Shift.1, -Shift.2, - Group)
dat <- filter(dat, !grepl("p53|Abl", Target))
dat$Id <- dat$Ring %% 4

PlotNetwork <- function(cellLine, treatment, timePoint){
        
        
        dat.1 <- filter(dat, Treatment %in% treatment & `Cell Line` %in% 
                                cellLine & `Time Point` %in% timePoint)
        

        dat.alt <- dcast(dat.1, `Cell Line` + Treatment + `Time Point` +
                Replicate + Id ~ Target, value.var = "Net Shift")

        rownames(dat.alt) = paste(dat.alt$`Cell Line`, dat.alt$Treatment,
                dat.alt$`Time Point`, dat.alt$Replicate, dat.alt$Id, sep = " ")

        dat.alt[,c(1,2,3,4,5)] <- NULL
        m.alt <<- as.matrix(dat.alt)
        holder.r <- rcorr(scale(m.alt))
        holder <- cor(m.alt)
        
        g <- graph.adjacency(holder.r$r, mode = "undirected", 
                             weighted = TRUE, diag = FALSE)
        
        custom.layout <- function(g, ...) {
                # layout.drl(g, weights = E(g)$weight, ...) # For bigger graphs
                layout.fruchterman.reingold(g, weights = E(g)$weight, ...)
        }
        
        # Format edges
        E(g)$cor <- E(g)$weight
        E(g)$weight <- abs(E(g)$cor)
        E(g)$color <- ifelse(E(g)$cor < 0, "blue", "darkred")
        E(g)$width <- atanh(E(g)$weight)
        
        # Format vertices
        #V(g)$size <- 3*abs(rowSums(holder.r$r))
        V(g)$label.color <- "black"
        V(g)$label.font <- 11
        
        # plot(g, layout = custom.layout)
        # plot(g, layout = layout.fruchterman.reingold)
        g <- delete.edges(g, E(g)[ abs(weight) < 0.9 ])
        plot(g, layout = layout_with_kk)
        
        #corrplot(holder.r$r, type = "upper")
        #corrplot(holder.r$P, type = "upper")
}

PlotNetwork(cellLine = c(6), treatment = c("Erlotinib"), timePoint = c(24))

m.control <- m.alt

m.alt <- rcorr(m.alt, m.control)
corrplot(m.alt$r, type = "upper")

hmcol <- colorRampPalette(c("darkblue", "white", "darkred"))(256)

heatmap(m.alt$r, col = hmcol)

# I need to find a way to get unique identifiers before dcast so that the 
# values don't get combined
