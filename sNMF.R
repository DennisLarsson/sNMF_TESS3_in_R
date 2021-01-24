# version 1.11

ck=require("LEA")
if (ck==FALSE) {
  install.packages("devtools")
  devtools::install_github("bcm-uga/LEA")
}
library(parallel)

setwd("/path/to/workDirectory/")  #give the path to the work directory where you files are and/or where you would like to output files to go
filename="spicatum.vcf"           #give the name of the file if it is in the work directory or the full path if it is somewhere else
outputname="spicatum"             #give the prefix for all of the output files
popmap="popmap_spicatum.txt"      #give the name of the popmap if it is in the work directory or the full path if it is somewhere else

NrCores=NULL # NULL can be replaced with a number of cores if you don't want to use all.
if (is.null(NrCores)) {NrCores=detectCores()} #detects how many cores are available

pop <- read.delim(popmap, header = FALSE)
pop_sorted<-pop[order(pop[,2]),]

vcf2geno(filename, output.file = paste(outputname,".geno",sep=""), force = TRUE) #converts the vcf to geno format 

obj.snmf <- snmf(paste(outputname,".geno",sep=""), K = 1:10, repetitions = 100, project = "new", CPU = NrCores, entropy = T, iterations = 2000) #runs the actual analysis

#if a projects has already been run, you can load it with the below command. Just uncomment it (remove the starting '#') and change the path to the file
#obj.snmf <- load.snmfProject("/path/to/workDirectory/spicatum.snmfProject") 

#this set up a list of cross entropy for each run and K.
CE_AllK <- list()
for (K in 1:10) {
  CE <- vector()
  x <- 1
  for (i in seq(K,length(obj.snmf@runs),10)){
    CE[x] <- obj.snmf@runs[[i]]@crossEntropy
    x=x+1
  }
  CE_AllK[[K]] <- CE
}

#this calculates the mean cross entropy
mean.ce = vector()
for (i in 1:10) {
  mean.ce[i] <- mean(unlist(CE_AllK[i]))
}
diff.mean.ce <- vector()
i=1
while (i <= length(mean.ce)-1){
  diff.mean.ce[i] <- mean.ce[i]-mean.ce[i+1]
  i = i +1
}

# Here the pdf is opened. everything that is plotted after here is plotted in the pdf file. the pdf is closed after dev.off(). skip this is dev.off() if you want to 
# display the plots in Rstudio
pdf(file=paste(outputname,"_snmf.pdf",sep=""), height = 5, width = 8, title = paste(outputname,"_snmf",sep=""))
plot(obj.snmf, pch = 19, col = "blue",main = "Cross-entropy")
plot(diff.mean.ce, xlab = "K value", ylab = "Difference in mean cross-entropy", main = "Difference in mean cross-entropy\nK(i)-K(i+1)")

boxplot(CE_AllK, main = "Cross-entropy box-plot", ylab = "Cross-entropy")

#sets up a loop to plot the admixture ratios for each K
for (K in 2:10) {
  bestRun <- match(min(unlist(CE_AllK[K])),unlist(CE_AllK[K])) #find the run number of the run with the lowest cross entropy

  qmatrix = Q(obj.snmf, K = K, run=bestRun) #extracts the best (lowest cross entropy) qmatrix (matrix of admixture ratios) from the full dataset.
  colorsPlot = c("red", "blue", "orange", "green","purple", "brown","darkgrey", "yellow", "darkgreen", "cyan")
  barplot(t(qmatrix), border = NA, space = 0, ylab = "Ancestry coefficients", col = colorsPlot, main = paste("Ancestry coefficients for K=",K,sep = ""))
  
  #plot lines and names of populations into the plot
  axis(1, tapply(1:nrow(pop), pop[,2], mean), unique(pop_sorted[,2]), las=2, cex.axis=1, tick = F, line = -0.8)
  abline(v=tapply(1:nrow(pop), pop[,2],max), lty=2, lwd=0.5)
  
  #write the qmatrix to a csv file that can be used in other applications.
  write.csv(qmatrix,paste(outputname,"_K",K,"_snmf.csv",sep=""), row.names = pop$V1)
}
dev.off()
