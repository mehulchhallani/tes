### Image Pdf

############ MainFinalALLWIINDOWS
cat("\014")
rm(list = ls())

fileList <- list.files("E:/Classification/Dataset", pattern="\\.RData$", full.names=TRUE)
DATATmp <- c()
OJA <- c()

setwd("E:/SmoothingProj/")
source(paste(getwd(), '/functions_V0.R', sep=''))
source(paste(getwd(), '/Despikee.R', sep=''))
source(paste(getwd(), '/Despike.R', sep=''))
source(paste(getwd(), '/scanSmootherOriginal1.R', sep=''))
source(paste(getwd(), '/RDATALOAD.R', sep=''))
source(paste(getwd(), '/calibration.R', sep=''))
source(paste(getwd(), '/SpectralMeanFilter.R', sep=''))
source(paste(getwd(), '/SpectralMedianFilter.R', sep=''))
source(paste(getwd(), '/ApplySmoothing1.R', sep=''))
source(paste(getwd(), '/functions_V1.R', sep=''))

setwd("E:/SmoothingProj/")
require(modelTransfer)
require(fields)
require(oce)
require(raster)
require(rgdal)
require(imager)
require(despike)

DIM <- c()
dataTmp <- c()
MICE <- c()
dimm <- c()
MEDIAN <- list()
MEAN <- list()
MAINSCANS <- vector('list',
                    length(fileList))
M15 <- c()
labssMean <- c("0", "Mean3", "Mean5", "Mean7")
labssMedian <- c("0", "Median3", "Median5", "Median7")
SCANN = c()
# labs = c('UnSmoothed', 'Median', 'Mean', "SG.MEAN", "SG.MEDIAN")

pdf(file='C:/Users/da32sac/Desktop/AllWINDOWS.Condition_1.pdf',
    width=11, height=6, onefile = TRUE, bg = "transparent")

for (i in 1:length(fileList))
{
  SCANN <- c()
  cat("\014")
  load(fileList[[i]])
  
  if (unique(DATA$Properties[,20]) == c("biopsy")) {
    next
  }
  
  print(paste("Scan number= ",i))
  print(paste("The mice info is",unique(DATA$Properties[,11])))
  labs <- DATA$Bewertung[,1]
  M15 <- DATA$Properties[,15]
  dim <- c(DATA$Dim[[1, 1]],DATA$Dim[[1, 2]])
  MICE <- rbind(MICE,   unique(DATA$Properties[,11]))
  DIM <- rbind(DIM, dim)
  print(paste("The dim is", dim))
  dataTmp<- RAW$y
  wn <- RAW$x
  calispec <- colMeans(RAW$Standard)
  
  if (i==1){
    WNCONST <- wn
  }
  
  ResultLIST <- list()
  
  ######################DESPIKE
  # dataTmp <- oce::despike(dataTmp,reference = "reference",n=4,k =5)
  # dataTmp <- oce::despike(dataTmp,reference = "reference",n=1,k =3)
  
  dataTmp <- despike(wn = wn, Data=dataTmp,
                     Thresh=7,Background=TRUE,
                     N.Median=7,pSd=30)
  
  ################### Wavenumber cali #######################################
  dataTmp <- calibration(MAT = dataTmp, WN =  wn, calispec = calispec,
                         wnNew=seq(from=300, to=3150, by=3), degree=2, WHICH=TRUE,
                         tPeak=c(390.9,465.1,504.0,651.6,710.8,797.2,834.5,857.9,968.7,1168.6,1236.8,
                                 1278.5,1323.9,1371.5,1648.4,2931.1,3102.4),
                         sdPeak=5, Plot=FALSE,i = i)
  
  wnNew <- as.numeric(colnames(dataTmp))
  
  
  ###### Smooth
  Spec_BA <- dataTmp
  ResultLIST <- ApplySmoothing1(Spec_BA =dataTmp , DATA = DATA)
  
  ########### Baseline
  dataTmp <- removeBaseline(wnNew,dataTmp,
                            method = "SNIP", iterations = 20)
  ResultLIST$SpectralMean3 <- removeBaseline(wnNew, ResultLIST$SpectralMean3, 
                                             method = "SNIP", iterations = 20)
  ResultLIST$SpectralMedian3 <- removeBaseline(wnNew, ResultLIST$SpectralMedian3,
                                               method = "SNIP", iterations = 20)
  ResultLIST$SpectralMean5 <- removeBaseline(wnNew, ResultLIST$SpectralMean5,
                                             method = "SNIP", iterations = 20)
  ResultLIST$SpectralMedian5 <- removeBaseline(wnNew, ResultLIST$SpectralMedian5,
                                               method = "SNIP", iterations = 20)
  ResultLIST$SpectralMean7 <- removeBaseline(wnNew, ResultLIST$SpectralMean7,
                                             method = "SNIP", iterations = 20)
  ResultLIST$SpectralMedian7 <- removeBaseline(wnNew, ResultLIST$SpectralMedian7,
                                               method = "SNIP", iterations = 20)
  
  ResultLIST$SpatialMean3 <- removeBaseline(wnNew, ResultLIST$SpatialMean3, 
                                            method = "SNIP", iterations = 20)
  ResultLIST$SpatialMedian3 <- removeBaseline(wnNew, ResultLIST$SpatialMedian3,
                                              method = "SNIP", iterations = 20)
  ResultLIST$SpatialMean5 <- removeBaseline(wnNew, ResultLIST$SpatialMean5,
                                            method = "SNIP", iterations = 20)
  ResultLIST$SpatialMedian5 <- removeBaseline(wnNew, ResultLIST$SpatialMedian5,
                                              method = "SNIP", iterations = 20)
  ResultLIST$SpatialMean7 <- removeBaseline(wnNew, ResultLIST$SpatialMean7,
                                            method = "SNIP", iterations = 20)
  ResultLIST$SpatialMedian7 <- removeBaseline(wnNew, ResultLIST$SpatialMedian7,
                                              method = "SNIP", iterations = 20)
  
  # matplot(wnNew, t(Spec_BA), type='l')
  
  ### Normalization
  ix <- which(wnNew%in%c(300:1850, 2815:3150))
  Spec_NM <- dataTmp/sqrt(rowSums(dataTmp[,ix]^2))
  
  SpatialMedian3 <- ResultLIST$SpatialMedian3/sqrt(rowSums(ResultLIST$SpatialMedian3[,ix]^2))
  SpatialMean3 <- ResultLIST$SpatialMean3/sqrt(rowSums(ResultLIST$SpatialMean3[,ix]^2))
  
  SpatialMedian5 <- ResultLIST$SpatialMedian5/sqrt(rowSums(ResultLIST$SpatialMedian5[,ix]^2))
  SpatialMean5 <- ResultLIST$SpatialMean5/sqrt(rowSums(ResultLIST$SpatialMean5[,ix]^2))
  
  SpatialMedian7 <- ResultLIST$SpatialMedian7/sqrt(rowSums(ResultLIST$SpatialMedian7[,ix]^2))
  SpatialMean7 <- ResultLIST$SpectralMean7/sqrt(rowSums(ResultLIST$SpatialMean7[,ix]^2))
  
  SpectralMedian3 <- ResultLIST$SpectralMedian3/sqrt(rowSums(ResultLIST$SpectralMedian3[,ix]^2))
  SpectralMean3 <- ResultLIST$SpectralMean3/sqrt(rowSums(ResultLIST$SpectralMean3[,ix]^2))
  
  SpectralMedian5 <- ResultLIST$SpectralMedian5/sqrt(rowSums(ResultLIST$SpectralMedian5[,ix]^2))
  SpectralMean5 <- ResultLIST$SpectralMean5/sqrt(rowSums(ResultLIST$SpectralMean5[,ix]^2))
  
  SpectralMedian7 <- ResultLIST$SpectralMedian7/sqrt(rowSums(ResultLIST$SpectralMedian7[,ix]^2))
  SpectralMean7 <- ResultLIST$SpectralMean7/sqrt(rowSums(ResultLIST$SpectralMean7[,ix]^2))
  
  ######### Silent region remove
  ix <- which(wnNew%in%c(400:1900, 2800:3150))
  wnNew <- wnNew[ix]
  Spec_NM <- Spec_NM[,ix]
  
  SpatialMedian3 <- SpatialMedian3[,ix]
  SpatialMean3 <- SpatialMean3[,ix]
  
  SpatialMedian5 <- SpatialMedian5[,ix]
  SpatialMean5 <- SpatialMean5[,ix]
  
  SpatialMedian7 <- SpatialMedian7[,ix]
  SpatialMean7 <- SpatialMean7[,ix]
  
  SpectralMedian3 <- SpectralMedian3[,ix]
  SpectralMean3 <- SpectralMean3[,ix]
  
  SpectralMedian5 <- SpectralMedian5[,ix]
  SpectralMean5 <- SpectralMean5[,ix]
  
  SpectralMedian7 <- SpectralMedian7[,ix]
  SpectralMean7 <- SpectralMean7[,ix]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  ##### Image
  
  ########### SpatialMean
  par(mfrow=c(2,4))
  # layout(matrix(c(1,2,3,0,4,5), 2, 3, byrow = TRUE))
  
  ###### Before Smoothing
  OJA <-array(Spec_NM,dim=c(DATA$Dim[[1, 1]],
                            DATA$Dim[[1, 2]],ncol(Spec_NM)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Scan", i, "Smoothing 0"))
  
  #window3
  OJA <-array(SpatialMean3,dim=c(DATA$Dim[[1, 1]],
                            DATA$Dim[[1, 2]],ncol(SpatialMean3)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spatial Mean 3"))
  
  #window5
  OJA <-array(SpatialMean5,dim=c(DATA$Dim[[1, 1]],
                                 DATA$Dim[[1, 2]],ncol(SpatialMean5)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spatial Mean 5"))
  
  #window7
  OJA <-array(SpatialMean7,dim=c(DATA$Dim[[1, 1]],
                                 DATA$Dim[[1, 2]],ncol(SpatialMean7)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spatial Mean 7"))
  
  mtext("Condition 1 (Before smoothing)", side = 3, line = -21, outer = TRUE)
  ########### SpatialMedian
  # par(mfrow=c(2,2))
  # layout(matrix(c(1,2,3,0,4,5), 2, 3, byrow = TRUE))
  
  ###### Before Smoothing
  OJA <-array(Spec_NM,dim=c(DATA$Dim[[1, 1]],
                            DATA$Dim[[1, 2]],ncol(Spec_NM)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Scan", i, "Smoothing 0"))
  
  #window3
  OJA <-array(SpatialMedian3,dim=c(DATA$Dim[[1, 1]],
                                 DATA$Dim[[1, 2]],ncol(SpatialMedian3)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spatial Median 3"))
  
  #window5
  OJA <-array(SpatialMedian5,dim=c(DATA$Dim[[1, 1]],
                                 DATA$Dim[[1, 2]],ncol(SpatialMedian5)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spatial Median 5"))
  
  #window7
  OJA <-array(SpatialMedian7,dim=c(DATA$Dim[[1, 1]],
                                 DATA$Dim[[1, 2]],ncol(SpatialMedian7)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spatial Median 7"))
  
  
  
  
  ########## Spectral
  ########### SpatialMean
  par(mfrow=c(2,4))
  # layout(matrix(c(1,2,3,0,4,5), 2, 3, byrow = TRUE))
  
  ###### Before Smoothing
  OJA <-array(Spec_NM,dim=c(DATA$Dim[[1, 1]],
                            DATA$Dim[[1, 2]],ncol(Spec_NM)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Scan", i, "Smoothing 0"))
  
  #window3
  OJA <-array(SpectralMean3,dim=c(DATA$Dim[[1, 1]],
                                 DATA$Dim[[1, 2]],ncol(SpectralMean3)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spectral Mean 3"))
  
  #window5
  OJA <-array(SpectralMean5,dim=c(DATA$Dim[[1, 1]],
                                 DATA$Dim[[1, 2]],ncol(SpectralMean5)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spectral Mean 5"))
  
  #window7
  OJA <-array(SpectralMean7,dim=c(DATA$Dim[[1, 1]],
                                 DATA$Dim[[1, 2]],ncol(SpectralMean7)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spectral Mean 7"))
  
  
  
  mtext("Condition 1 (Before smoothing)", side = 3, line = -21, outer = TRUE)
  ########### SpatialMedian
  # par(mfrow=c(2,2))
  # layout(matrix(c(1,2,3,0,4,5), 2, 3, byrow = TRUE))
  
  ###### Before Smoothing
  OJA <-array(Spec_NM,dim=c(DATA$Dim[[1, 1]],
                            DATA$Dim[[1, 2]],ncol(Spec_NM)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Scan", i, "Smoothing 0"))
  
  #window3
  OJA <-array(SpectralMedian3,dim=c(DATA$Dim[[1, 1]],
                                   DATA$Dim[[1, 2]],ncol(SpectralMedian3)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spectral Median 3"))
  
  #window5
  OJA <-array(SpectralMedian5,dim=c(DATA$Dim[[1, 1]],
                                   DATA$Dim[[1, 2]],ncol(SpectralMedian5)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spectral Median 5"))
  
  #window7
  OJA <-array(SpectralMedian7,dim=c(DATA$Dim[[1, 1]],
                                   DATA$Dim[[1, 2]],ncol(SpectralMedian7)))
  X<-as.cimg(OJA[,,250])
  plot(X,axes=FALSE)
  title(paste("Spectral Median 7"))
  
  
  
  #Spatial
  par(mfrow=c(1,1))
  ####### Window3,5,7
  plot(wnNew, colMeans(Spec_NM, na.rm = TRUE), type ='l', lwd=1.5, col = 1,
       xlab=expression(Wavenumber / cm^-1), 
       ylab="Normalized Raman Intensity",
       ylim = range(colMeans((Spec_NM)-0.001, na.rm=TRUE),
                    colMeans(Spec_NM+0.08, na.rm=TRUE)))
  lines(wnNew, colMeans(SpatialMean3, na.rm = TRUE), type ='l', lwd=1.5, col = 2)
  lines(wnNew, colMeans(SpatialMean5, na.rm = TRUE), type ='l', lwd=1.5, col = 3)
  lines(wnNew, colMeans(SpatialMean7, na.rm = TRUE), type ='l', lwd=1.5, col = 4)
  legend('bottomright', legend=paste(as.character(labssMean), ''),
         fill=1:length(labssMean), bty="n", inset=0.09,
         cex = 0.8)
  
  
  lines(wnNew, colMeans(Spec_NM+0.07, na.rm = TRUE), type ='l', lwd=1.5, col = 1)
  lines(wnNew, colMeans(SpatialMedian3+0.07, na.rm = TRUE), type ='l', lwd=1.5, col = 2)
  lines(wnNew, colMeans(SpatialMedian5+0.07, na.rm = TRUE), type ='l', lwd=1.5, col = 3)
  lines(wnNew, colMeans(SpatialMedian7+0.07, na.rm = TRUE), type ='l', lwd=1.5, col = 4)
  legend('topright', legend=paste(as.character(labssMedian), ''),
         fill=1:length(labssMedian), bty="n",
         cex = 0.8)
  title(paste("Scan", i ,"Condition 1","for Spatial"))
  
  
  
  
  # Spectral
  ####### Window3,5,7
  par(mfrow=c(1,1))
  plot(wnNew, colMeans(Spec_NM, na.rm = TRUE), type ='l', col = 1,
       xlab=expression(Wavenumber / cm^-1),
       ylab="Normalized Raman Intensity", lwd=1.5,
       ylim = range(colMeans((Spec_NM)-0.001, na.rm=TRUE),
                    colMeans(Spec_NM+0.08, na.rm=TRUE)))
  lines(wnNew, colMeans(SpectralMean3, na.rm = TRUE), type ='l', col = 2, lwd=1.5)
  lines(wnNew, colMeans(SpectralMean5, na.rm = TRUE), type ='l', col = 3, lwd=1.5)
  lines(wnNew, colMeans(SpectralMean7, na.rm = TRUE), type ='l', col = 4, lwd=1.5)
  legend('bottomright', legend=paste(as.character(labssMean), ''),
         fill=1:length(labssMean), bty="n", inset=0.09,
         cex = 0.8)
  
  
  lines(wnNew, colMeans(Spec_NM+0.07, na.rm = TRUE), type ='l', col = 1, lwd=1.5)
  lines(wnNew, colMeans(SpectralMedian3+0.07, na.rm = TRUE), type ='l', col = 2, lwd=1.5)
  lines(wnNew, colMeans(SpectralMedian5+0.07, na.rm = TRUE), type ='l', col = 3, lwd=1.5)
  lines(wnNew, colMeans(SpectralMedian7+0.07, na.rm = TRUE), type ='l', col = 4, lwd=1.5)
  legend('topright', legend=paste(as.character(labssMedian), ''),
         fill=1:length(labssMedian), bty="n",
         cex = 0.8)
  title(paste("Scan", i ,"Condition 1","for Spectral"))
  
}

dev.off()
dev.off()
dev.off()

