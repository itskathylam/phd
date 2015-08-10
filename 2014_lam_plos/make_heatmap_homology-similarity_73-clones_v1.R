#install.packages("gplots")
library(gplots)

setwd("D:\\Dropbox\\PhD\\Manuscripts\\Pooled Sequencing\\Resubmission\\Homology analysis")
x <- read.csv("out.csv", header=FALSE, sep=",", as.is=TRUE)
x

#get names
x.names <- sort(unique(x$V1))
x.names

# create a matrix of the right size and put names on it
x.sim <- matrix(0, length(x.names), length(x.names))
dimnames(x.sim) <- list(x.names, x.names)

# create indices by converting names to numbers and create the normal and reversed to fill in all the matrix
x.ind <- rbind(cbind(match(x[[1]], x.names), match(x[[2]], x.names)), cbind(match(x[[2]], x.names), match(x[[1]], x.names)))
x.sim[x.ind] <- rep(x[[3]], 2)
x.sim

#heatmap
par(mar=c(5,5))
x.sim.matrix = as.matrix(x.sim)
lwid=c(1.5,4.5)	
lhei = c(0.75,4.5)
heatmap.2(x.sim.matrix, trace="none", lwid=lwid, lhei=lhei, density.info=c("none"), cexRow=0.75, cexCol=0.75, margins=c(6.5,6.5))



