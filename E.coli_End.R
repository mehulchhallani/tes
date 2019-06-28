!!!!!!!!!!!!!!!

rm(list = ls())
cat("\014")

rm(list=ls())

while(dev.cur()!=1)
{
dev.off()
}

forPath <- function(x){ x<-1 }
setwd(paste(getSrcDirectory(forPath), '/', sep=''))

setwd("E:/Regression analysis of phases CO2/E.coli/3.End")
source(paste(getwd(), '/All_Functions_V4_All.Phases.R', sep=''))
setwd("E:/Regression analysis of phases CO2/E.coli/3.End")

require(caret)
require(zoo)
require(pls)
require(e1071)
require(baseline)
require(oce)
require(MASS)
require(gsubfn)
require(stringr)
require(EMSC)
require(gtools)
require(randomForest)

caliSub <- '4-AAP'
fPattern <- '.txt'

paths <- dir(getwd(), full.names=TRUE, recursive=TRUE, pattern=fPattern)
paths <- unlist(lapply(paths, dirname))

ix <- lapply(caliSub, grep, x=paths, fixed=TRUE)[[1]]
dataPath <- mixedsort(sort(unique(paths[-ix], ".txt", sep = ""))) ## path of data
aapPath <- mixedsort(sort(unique(paths[ix], ".txt", sep = "")))   ## path of 4-AAP

### extract the time when 4AAP was measured
aapDate <- c()
for (j in 1:length(aapPath))
{
fileTmp <- dir(aapPath[j], pattern=fPattern)[1]
timeIndex <- gregexpr("[0-9]{6}_[0-9]{6}", fileTmp)[[1]]
timeCode <- substr(fileTmp, timeIndex, timeIndex + 12)
aapDate <- c(aapDate, as.POSIXct(strptime(timeCode, "%y%m%d_%H%M%S")))
}
aapDate <- as.POSIXlt(aapDate, tz = "", origin="1970-01-01")
aapDate <- aapDate[order(aapDate)]

batches <- c()
DATA <- c()
WN <- c()
Cons <- c()
BATCHES <- c()
PHASE <- c()

for(i in 1:length(dataPath))
{
dataTmp <- readData(path=dataPath[i], pattern=fPattern, 
                    excludeDirs='4-AAP',
                    batches=c('Batch1', 'Batch2', 'Batch3'))
wnRAW <- dataTmp$wn
batches <- c(batches, dataTmp$batches)
RAW <- dataTmp$data

dataTmp$data <- despike(wnRAW, Data=dataTmp$data,
                        Thresh=7,Background=TRUE,
                        N.Median=7,pSd=30)

dataTmp$data <- despike(wnRAW, Data=dataTmp$data,
                        Thresh=3,Background=TRUE,
                        N.Median=3,pSd=20)

### calibration
### decide time when data is measured
fileTmp <- dir(dataPath[i], pattern=fPattern)[1]
timeIndex <- gregexpr("[0-9]{6}_[0-9]{6}", fileTmp)[[1]]
timeCode <- substr(fileTmp, timeIndex, timeIndex + 12)
dataDate <- as.POSIXct(strptime(timeCode, "%y%m%d_%H%M%S"))

### decide which 4-AAP file to use for calibration according to time

ix <- which(aapDate <= dataDate)
if(length(ix)<1)
  ix <- which(appDate>dataDate)[1]
else
  ix <- ix[length(ix)]

#Wavenumber calibration
##### This is essentially the same function that you designed long back
dataTmp <- calibration(wn=dataTmp$wn, Data=dataTmp$data, 
                       caliPath=aapPath[ix], caliPattern=fPattern, caliSuffix='4-AAP',
                       wnNew=seq(from=300,to=3150,by=3), degree=2, WHICH=TRUE,
                       tPeak=c(390.9,465.1,504.0,651.6,710.8,797.2,834.5,857.9,968.7,1168.6,1236.8,
                               1278.5,1323.9,1371.5,1648.4,2931.1,3102.4),sdPeak=5, Plot=TRUE)

DATA <- rbind(DATA, dataTmp$data)
WN <- dataTmp$wn

tz <- strsplit(dataPath[i], "/")
BATCH <- tz[[1]][5]
tmpcon <- tz[[1]][6]
tmpcon <- sub("\\%.*", "", tmpcon)
Phase <- tz[[1]][8]

BATCHES <- c(BATCHES ,rep(BATCH, nrow(dataTmp$data)))
Cons <- c(Cons ,rep(tmpcon, nrow(dataTmp$data)))
PHASE <- c(PHASE ,rep(Phase, nrow(dataTmp$data)) )
}

dev.off()
Cons <- as.numeric(Cons)
unique(Cons)
unique(BATCHES)
str(DATA)

###### Baseline correction
Spectra_BA <- baseline.als(spectra = DATA,lambda = 4.5, p = 0.001)
pdf(file =paste(getwd(),"/Baseline.pdf",sep = ""),
    width=11, height=6, onefile = TRUE, bg = "transparent")
for (i in 1:nrow(DATA))
{
  plot(WN,DATA[i, ], type = 'l', col = 1)
  lines(WN, Spectra_BA$baseline[i,], type = 'l', col = 2)
  title(paste("Spec number", i, "belonging to concentration", Cons[i], "from batch", BATCHES[i]))
}
dev.off()
Spectra_BA <- Spectra_BA$corrected

##### Norma
ix <- which(WN%in%c(500:1800, 2800:3150))
Spectra_NM <- Spectra_BA/sqrt(rowSums(Spectra_BA[,ix]^2))

##### Mean Spec & SD & difference & Batch difference
pdf(file =paste(getwd(),"/1.MeanSpec.pdf",sep = ""),
    width=11, height=6, onefile = TRUE, bg = "transparent")
plot.mean.sd.diff(wn = WN, data = Spectra_NM,
                  labels = Cons, type = "mean",
                  strain = "E.coli end phase")
dev.off()

####### SD spec
pdf(file =paste(getwd(),"/2.SDSpec.pdf",sep = ""),
    width=11, height=6, onefile = TRUE, bg = "transparent")
plot.mean.sd.diff(wn = WN, data = Spectra_NM,
                  labels = Cons, type = "mean.sd",
                  strain = "E.coli end phase")
dev.off()

###### Diff spec between groups 0 and 10
pdf(file =paste(getwd(),"/3.DiffSpec.pdf",sep = ""),
    width=11, height=6, onefile = TRUE, bg = "transparent")
plot.mean.sd.diff(wn = WN, data = Spectra_NM,
                  labels = Cons, type = "difference",
                  diff.betwn = c(1,5),strain = "E.coli end phase")
dev.off()

###### Diff between batches
pdf(file =paste(getwd(),"/4.BatchDifference.pdf",sep = ""),
    width=11, height=6, onefile = TRUE, bg = "transparent")
plot.mean.sd.diff(wn = WN, data = Spectra_NM,
                  labels = Cons, batches = BATCHES,type = "BatchDifference",
                  strain = "E.coli end phase")
dev.off()

######### Silent region remove
ix <- which(WN%in%c(500:1800, 2800:3150))
Spectra_NM <- Spectra_NM[,ix]
WN <- WN[ix]

##### This gives us the concentrations of the dataset
##### We have to select the ones that we need to build the training model
unique(Cons)

save(Spectra_NM, file = paste(getwd(),"/Spectra_NM.RData",sep = ""))
save(Cons, file = paste(getwd(),"/Cons.RData",sep = ""))
save(BATCHES, file = paste(getwd(),"/BATCHES.RData",sep = ""))
save(WN, file = paste(getwd(),"/WN.RData",sep = ""))

load(paste(getwd(),"/Spectra_NM.RData",sep = ""))
load(paste(getwd(),"/Cons.RData",sep = ""))
load(paste(getwd(),"/BATCHES.RData",sep = ""))
load(paste(getwd(),"/WN.RData",sep = ""))

# ######## select which concentrations you want to train ideally to obtain
# ######## the sigmoidal function, it should be, 0 and 10 which corresponds to c(1,2)

############## LDA
pdf(file =paste(getwd(),"/LDA.pdf",sep = ""),
    width=12, height=7, onefile = TRUE, bg = "transparent")
tz <- LDA.Disc.Values.Prediction(Spec = Spectra_NM, Batches = BATCHES,
                                 labels = Cons, TRAIN = c(1,5), nPCs = 5:30,
                                 strain = "E.coli end phase")
dev.off()

############## PLS
pdf(file =paste(getwd(),"/PLS.pdf",sep = ""),
    width=12, height=7, onefile = TRUE, bg = "transparent")
tz <- PLS.Disc.Values.Prediction(Spec = Spectra_NM, Batches = BATCHES,
                                 labels = Cons, TRAIN = c(1,5), nPCs = 5:30,
                                 strain = "E.coli end phase")
dev.off()

############## Random Forest
pdf(file =paste(getwd(),"/RandomForest.pdf",sep = ""),
    width=12, height=7, onefile = TRUE, bg = "transparent")
tz <- RANDOMFOREST.Disc.Values.Prediction(Spec = Spectra_NM, Batches = BATCHES,
                                          labels = Cons, TRAIN = c(1,5), nPCs = 5:30,
                                          strain = "E.coli end phase")
dev.off()

############## SVM
pdf(file =paste(getwd(),"/SVM.pdf",sep = ""),
    width=12, height=7, onefile = TRUE, bg = "transparent")
tz <- SVM.Disc.Values.Prediction(Spec = Spectra_NM, Batches = BATCHES,
                                 labels = as.numeric(Cons), 
                                 TRAIN = c(1,5), nPCs = 5:30,
                                 strain = "E.coli end phase")
dev.off()

############## GLM
pdf(file =paste(getwd(),"/GLM.pdf",sep = ""),
    width=12, height=7, onefile = TRUE, bg = "transparent")
tz <- GLM.Disc.Values.Prediction(Spec = Spectra_NM, Batches = BATCHES,
                                 labels = as.numeric(Cons), 
                                 TRAIN = c(1,5), nPCs = 5:30,
                                 strain = "E.coli end phase")
dev.off()
