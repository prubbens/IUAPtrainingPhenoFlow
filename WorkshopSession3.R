#### only execute first time on your own machine ####
# install.packages(devtools)
# library(devtools)
# install_github("rprops/Phenoflow_package")


#### set data directory ####
dataloc <- "/media/projects2/IUAPworkshop" #folder containing all your fcs files

#### load required packages ####
library(Phenoflow)
library(ggplot2)
#### set RNG seed for reproducibility ####
set.seed(777)


#### load the data ####
flowData <- read.flowSet(path=dataloc,transformation = FALSE, pattern=".fcs")

#### check attributes ####
attributes(flowData)

#### transform ####
flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]
remove(flowData)

#### set gates and verify ####
### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(7.8,7.8,16,16,3,8,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check (check for 1 and 5)
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[10], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

#### normalization to 0,1 range ==> for bandwidth 0,001 in fingerprint calculation
summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = max(summary[,1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))

### Randomly resample to the lowest sample size
flowData_transformed <- FCS_resample(flowData_transformed, replace=FALSE)
### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

#### Calculate ecological parameters from normalized fingerprint ####
### Densities will be normalized to the interval [0,1]
### n = number of replicates
### d = rounding factor
Diversity.fbasis <- Diversity(fbasis,d=3,plot=TRUE, R=20) #boot kept low for workshop, should be 999 or more
Evenness.fbasis <- Evenness(fbasis,d=3,n=1,plot=TRUE) #Pareto Evenness (Wittebolle, 2009)
Structural.organization.fbasis <- So(fbasis,d=3,n=1,plot=TRUE) #Structural Organization parameter (Koch et al. (2014), Frontiers in Microbiology)
Coef.var.fbasis <- CV(fbasis,d=3,n=1,plot=TRUE) # Coefficient of Variation (CV) 


#### Export ecological data to .csv file in the chosen directory ####
write.csv2(file="results.metrics.csv",
           cbind(Diversity.fbasis, Evenness.fbasis,
                 Structural.organization.fbasis,
                 Coef.var.fbasis))
#### extract metadata ####
evennessclass <- sub("(^.+?)\\_.*$","\\1",as.character(Diversity.fbasis$Sample_name))
timepoint <- sub("(^.+?)\\_.*(t[15]).*$","\\2",as.character(Diversity.fbasis$Sample_name))

Diverstity.fbasis.meta <- data.frame(Diversity.fbasis,
                                     initial_even=evennessclass,
                                     time=timepoint)


#### Compare alpha diversity ####
alphacompgg.D0 <- ggplot(data=Diverstity.fbasis.meta,
                      aes(x=time,y=D0,fill=initial_even))

alphacompgg.D0 + geom_boxplot()

alphacompgg.D2 <- ggplot(data=Diverstity.fbasis.meta,
                         aes(x=time,y=D2,fill=initial_even))

alphacompgg.D2 + geom_boxplot()

#### Beta-diversity assessment of fingerprint ####
beta.div <- beta_div_fcm(fbasis, ord.type="NMDS", dist="jaccard", iter=20)

plot_beta_fcm(beta.div,color = evennessclass,shape = timepoint,
              labels=c(levels(factor(evennessclass)),levels(factor(timepoint))))


#### Creating a rectangle gate for counting HNA and LNA cells ####
# this also exports the exact cell densities
#only valid for data gathered on BD C6 accuri FC!
# HN/LNA according to PRest et al.
rGate_HNA <- rectangleGate("FL1-H"=c(asinh(20000), 20)/max,"FL3-H"=c(0,20)/max, 
                           filterId = "HNA bacteria")
### Normalize total cell gate
sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3)/max,ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

### Check if rectangle gate is correct, if not, adjust rGate_HNA
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=rGate_HNA,
       scales=list(y=list(limits=c(0,1)),
                   x=list(limits=c(0.4,1))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, 
                                                          cex=2), smooth=FALSE)
### Extract the cell counts
a <- filter(flowData_transformed, rGate_HNA) 
HNACount <- summary(a);HNACount <- toTable(HNACount)
s <- filter(flowData_transformed, polyGate1)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)

### Exporting cell counts to .csv file to working directory
write.csv2(file="results.counts.csv",
           data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                      Total.cells = TotalCount$true,HNA.cells = HNACount$true))
