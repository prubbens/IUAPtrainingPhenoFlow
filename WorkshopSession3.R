#### only execute first time on your own machine ####
# install.packages(devtools)
# library(devtools)
# install_github(rprops/Phenoflow_package)


#### set data directory ####
dataloc <- "/media/projects2/IUAPworkshop" #folder containing all your fcs files

#### load required packages ####
library(Phenoflow)

#### set RNG seed for reproducibility ####
set.seed(777)


#### load the data ####
flowData <- read.flowSet(path=dataloc,transformation = FALSE, pattern=".fcs")
