#Class-09_10_18

library(SpikeInSubset)

#Load data
data(spikein95)

#Normalize
rma.eset <- rma(spikein95)

#Expression levels put into rma.e
rma.e <- exprs(rma.eset)

d <- rowMeans(rma.e[,1:3]) - rowMeans(rma.e[,4:6])
a <- rowMeans(rma.e)

source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")
library("genefilter")
pData(rma.eset) <- pData(mas5.eset)
tt <- rowMeans(rma.e)
lod <- tt
plot(d, lod, cex = 0.25, main = "Volcano plot for MA", xlim = c(-2, 2), xlab = "M", ylab = "A", yaxt = "n")
axis(2, at = seq(0, 3, by = 1), labels = 10^(-seq(0, 3, by = 1)))
points(d[spikedIndex], lod[spikedIndex], pch = 19, col = "red")
abline(h = 2, v = c(-1, 1))

#cex =size of dots
#xlim= limits of x axes
#yaxt= "n" = do not use default axis next line you put your own axis in
#at = seq(0, 3, by = 1)  = devide sequence 0-3 by 1= 0,1,2,3
#at= -seq(0, 3, by = 1) = 0,-1,-2,-3

plot(d, lod, cex = 0.25, main = "Volcano plot for MA", xlim = c(-2, 2), xlab = "M", ylab = "A", yaxt = "n")
axis(2, at = seq(0, 3, by = 1), labels = 10^(-seq(0, 3, by = 1)))
axis(2, at = seq(0, 6, by = 1), labels = 10^(-seq(0, 6, by = 1)))

library(affy)

#set working directory
setwd("/Users/Amanda/Desktop/estrogen")

targetsFile <-"estrogen.txt"

#Read the file
pd <- read.AnnotatedDataFrame(targetsFile,header=TRUE,sep="",row.names=1)
#can just use estrogen.txt instead of target file
#header is true considers first row in dataset header if put as FALSE, first line is condsidered a sample and default headers are given, no seperator
#give row names to object pd

pData(pd)

#samples as rows and variables as columns =8 samples estrogen present and time

ER <- pData(pd)$estrogen
#in the data, pulling estrogen column 

Time <- factor(pData(pd)$time.h)
#in the data, pulling time column. Factor assigns levels (unique values) (10,48)
 
design <- model.matrix(~ER+Time)
design

#linear fit matrix that give binary values (0,1) to columns that have values (1) or do not (0)

design2 <- model.matrix(~ER*Time)
design2

#linear fit that multiplies binary values in columns 

raw <-ReadAffy(celfile.path = "/Users/Amanda/Desktop/estrogen", filenames=rownames(pData(pd)),phenoData = pd)
raw

# read chiptakes, path where files are, file names are row names,whatever is in pd is phenoData

boxplot(raw,col="red")

#Assignment 3 Page on what is linear regression/linear fit 