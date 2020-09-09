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

NrCores=NULL
if (is.null(NrCores)) {NrCores=detectCores()}

coordinates <- read.delim(coordinates_file, sep="\t", header = FALSE)
vcf2geno(filename, output.file = paste(outputname,".geno",sep=""), force = TRUE)
geno2lfmm(paste(outputname,".geno",sep=""), output.file = paste(outputname,".lfmm",sep=""), force = TRUE)
system(paste("sed -i 's/9/NA/g' ",outputname,".lfmm",sep = ""))
lfmm_file <- read.delim(paste(outputname,".lfmm",sep=""), sep=" ", header = FALSE)

tess3.obj <- tess3(X = lfmm_file, coord = as.matrix(coordinates), K = 1:10, method = "projected.ls", ploidy = 2, openMP.core.num = NrCores, rep = 100)

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

pdf(file=paste(outputname,"_tess3.pdf",sep=""), height = 5, width = 12, title = outputname)
par(mfrow = c(1, 2))
plot(tess3.obj, pch = 19, col = "blue", xlab = "Number of ancestral populations", ylab = "Cross-validation score",main = "Cross-validation score for each K")
plot(diff.mean.cverror, xlab = "K value", ylab = "Difference in mean CV error", main = "Difference in mean CV error\nK(i)-K(i+1)")

r=seq(2,10)
for (K in r) { 
  par(mfrow = c(1, 2))
  q.matrix <- qmatrix(tess3.obj, K = K)
  
  my.colors <- c("red", "blue", "orange", "green","purple", "brown","darkgrey", "yellow", "darkgreen", "cyan")
  my.palette <- CreatePalette(my.colors, 9)
  
  barplot(q.matrix, border = NA, space = 0, sort.by.Q = F,
          ylab = "Ancestry coefficients", col = my.palette,
          main = paste("Ancestry coefficients for K=",K)) -> bp
  
  pop <- read.delim(popmap, header = FALSE)
  pop_sorted<-pop[order(pop[,2]),]
  axis(1, tapply(1:nrow(pop), pop[,2],mean),unique(pop_sorted[,2]),las=2, cex.axis=0.4,tick = F,line = -0.8)
  abline(v=tapply(1:nrow(pop), pop[,2],max), lty=2, lwd=0.5)
  
  plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),  
       main = paste("Diffusion map of ancestry for K=",K),
       xlab = "Longitude", ylab = "Latitude", 
       resolution = c(600,600), cex = .4, 
       col.palette = my.palette)
  write.csv(q.matrix,paste(outputname,"_K",K,"_tess3.csv",sep=""), row.names = pop$V1)
}
dev.off()

