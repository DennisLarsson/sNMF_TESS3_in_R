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

vcf2geno(filename, output.file = paste(outputname,".geno",sep=""), force = TRUE) #converts the vcf to geno format 

obj.snmf <- snmf(paste(outputname,".geno",sep=""), K = 1:10, repetitions = 100, project = "new", CPU = NrCores, entropy = T, iterations = 2000) #runs the actual analysis

#if a projects has already been run, you can load it with the below command. Just uncomment it (remove the starting '#') and change the path to the file
#obj.snmf <- load.snmfProject("/path/to/workDirectory/spicatum.snmfProject") 

#this set up a database for storing the cross entropy values. if you want to run more than 10 Ks you need to add more, for example if you want  for K=11
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

mean.cverror = vector()
for (i in 1:10) {
  mean.cverror[i] <- mean(unlist(CE_AllK[i]))
}
diff.mean.cverror <- vector()
i=1
while (i <= length(mean.cverror)-1){
  diff.mean.cverror[i] <- mean.cverror[i]-mean.cverror[i+1]
  i = i +1
}

pdf(file=paste(outputname,"_snmf.pdf",sep=""), height = 5, width = 8, title = outputname)
plot(obj.snmf, pch = 19, col = "blue",main = "Cross-entropy")
plot(diff.mean.cverror, xlab = "K value", ylab = "Difference in mean cross-entropy", main = "Difference in mean cross-entropy\nK(i)-K(i+1)")
boxplot(CE_AllK, main = "Cross-entropy box-plot", ylab = "Cross-entropy")

for (K in 2:10) {
  bestRun <- match(min(unlist(CE_AllK[K])),unlist(CE_AllK[K]))

  qmatrix = Q(obj.snmf, K = K, run=bestRun)
  colorsPlot = c("red", "blue", "orange", "green","purple", "brown","darkgrey", "yellow", "darkgreen", "cyan")
  barplot(t(qmatrix), border = NA, space = 0, ylab = "Ancestry coefficients", col = colorsPlot, main = paste("Ancestry coefficients for K=",K,sep = ""))
  
  pop <- read.delim(popmap, header = FALSE)
  pop_sorted<-pop[order(pop[,2]),]
  axis(1, tapply(1:nrow(pop), pop[,2],mean),unique(pop_sorted[,2]),las=2, cex.axis=0.5,tick = F,line = -0.8)
  abline(v=tapply(1:nrow(pop), pop[,2],max), lty=2, lwd=0.5)
  K_list <- seq(1,K)
  i=1
  columnNames <- list()
  while (i <= length(K_list)) {
    columnNames[i] <- paste("K",K,"_",i,sep = "")
    i=i+1
  }
  
  write.csv(qmatrix,paste(outputname,"_K",K,"_snmf.csv",sep=""), row.names = pop$V1)
}


dev.off()
