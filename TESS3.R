#Version 1.21

ck=require("tess3r")
if (ck==FALSE) {
  install.packages("devtools")
  devtools::install_github("bcm-uga/TESS3_encho_sen")
}
ck=require("LEA")
if (ck==FALSE) {
  install.packages("devtools")
  devtools::install_github("bcm-uga/LEA")
}
library(parallel)

setwd("/path/to/workDirectory/")             #give the path to the work directory where you files are and/or where you would like to output files to go
filename="spicatum.vcf"                      #give the name of the file if it is in the work directory or the full path if it is somewhere else
outputname="spicatum"                        #give the prefix for all of the output files
popmap="popmap_spicatum.txt"                 #give the name of the popmap if it is in the work directory or the full path if it is somewhere else
coordinates_file="coordinates_spicatum.txt"  #give the name of the coordinates file if it is in the work directory or the full path if it is somewhere else

NrCores=NULL # NULL can be replaced with a number of cores if you don't want to use all.
if (is.null(NrCores)) {NrCores=detectCores()} #detects how many cores are available

coordinates <- read.delim(coordinates_file, sep="\t", header = FALSE) #loads in the coordinates file into a table
vcf2geno(filename, output.file = paste(outputname,".geno",sep=""), force = TRUE) #converts the vcf to geno format 
geno2lfmm(paste(outputname,".geno",sep=""), output.file = paste(outputname,".lfmm",sep=""), force = TRUE) #converts geno format into lfmm format
system(paste("sed -i 's/9/NA/g' ",outputname,".lfmm",sep = "")) # re formats missingness in the lfmm file from 9s to NAs
lfmm_file <- read.delim(paste(outputname,".lfmm",sep=""), sep=" ", header = FALSE) # reads in the lfmm file

#this runs the tess 3 analysis 
tess3.obj <- tess3(X = lfmm_file, coord = as.matrix(coordinates), K = 1:10, method = "projected.ls", ploidy = 2, openMP.core.num = NrCores, rep = 100)

#this calculates the mean cross validation score (rmse)
mean.cverror = vector()
for (i in 1:10) {
  mean.cverror[i] <- mean(tess3.obj[[i]][["rmse"]])
}
diff.mean.cverror <- vector()
i=1
while (i <= length(mean.cverror)-1){
  diff.mean.cverror[i] <- mean.cverror[i]-mean.cverror[i+1]
  i = i +1
}
K_range = seq(2,10)
# Here the pdf is opened. everything that is plotted after here is plotted in the pdf file. the pdf is closed after dev.off(). skip this is dev.off() if you want to 
# display the plots in Rstudio
pdf(file=paste(outputname,"_tess3.pdf",sep=""), height = 5, width = 12, title = outputname)
par(mfrow = c(1, 2))
plot(tess3.obj, pch = 19, col = "blue", xlab = "Number of ancestral populations", ylab = "Cross-validation score",main = "Cross-validation score for each K")
plot(diff.mean.cverror, xaxt="n", xlab = "K value", ylab = "Difference in mean CV error", main = "Difference in mean CV error\nK(i)-K(i+1)")
axis(side=1, at=1:9, labels = K_range)

#sets up a loop to plot the admixture ratios for each K
for (K in 2:10) {
  par(mfrow = c(1, 2)) #set up the window to have two parts, one for the barplot and one for the map
  q.matrix <- qmatrix(tess3.obj, K = K) #extracts the best (lowest cross validation score) qmatrix (matrix of admixture ratios) from the full dataset.
  
  my.colors <- c("red", "blue", "orange", "green","purple", "brown","darkgrey", "yellow", "darkgreen", "cyan")
  my.palette <- CreatePalette(my.colors, 9)
  
  #This plots the barplot with the admixture ratios
  barplot(q.matrix, border = NA, space = 0, sort.by.Q = F, ylab = "Ancestry coefficients", col = my.palette, main = paste("Ancestry coefficients for K=",K,sep="")) -> bp
  
  #plot lines and names of populations into the plot
  pop <- read.delim(popmap, header = FALSE)
  pop_sorted<-pop[order(pop[,2]),]
  axis(1, tapply(1:nrow(pop), pop[,2], mean), unique(pop_sorted[,2]), las=2, cex.axis=1, tick = F, line = -0.8)
  abline(v=tapply(1:nrow(pop), pop[,2],max), lty=2, lwd=0.5)
  
  #plot the diffusion map of ancestries
  plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10), main = paste("Diffusion map of ancestry for K=",K), 
       xlab = "Longitude", ylab = "Latitude", resolution = c(600,600), cex = .4, col.palette = my.palette)
  
  #writes the qmatrix to a csv file for future use
  write.csv(q.matrix,paste(outputname,"_K",K,"_tess3.csv",sep=""), row.names = pop$V1)
}
dev.off()

