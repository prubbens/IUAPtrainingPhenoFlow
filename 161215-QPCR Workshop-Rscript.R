rm (list = ls())
library(ade4TkGUI)
library(grofit)
library(readxl)
library(xlsx) #requires Java on your PC!
#import data as dataframe
# ade4TkGUI()
#   ## "Read a data file"_no header
#   ## Copy "QPCR-time"-Import from clipboard in dataframe "time1"
#   ## Copy "QPCR-data"-Import from clipboard in dataframe "data1"
# time1
# data1

data1 <- readxl::read_excel("161216-Data sets 1 species.xls",sheet="QPCR-data",col_names = FALSE)
time1 <- readxl::read_excel("161216-Data sets 1 species.xls",sheet="QPCR-time",col_names = FALSE)


MyOptL<-grofit.control(fit.opt="m",model.type=c("gompertz"),log.y.gc=TRUE, interactive = FALSE)
tableresults<-grofit(time1,data1,FALSE,MyOptL)$gcFit$gcTable
print(tableresults)
write.csv2(x = tableresults,file = "tableresults.csv") #write.csv2 is for semicolon separated files (, as decimal separator)

#### if you have Java on your machine you can write directly to the excel file ####
write.xlsx(tableresults,"161216-Data sets 1 species.xls",sheetName = "Eff2",append = TRUE)