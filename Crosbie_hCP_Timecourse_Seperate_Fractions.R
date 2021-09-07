install.packages(c("gplots", "lme4", "reshape", "reshape2",
                   "ggplot2", "ggrepel", "data.table", "dplyr", "tidyr","survival", "doSNOW", "snow", "foreach", 'stringr',"readr","randomForest", "minpack.lm"),repos='http://cran.us.r-project.org')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos='http://cran.us.r-project.org')
BiocManager::install("MSstats")
a

rm(list = ls())
library(stringr)
library(readr)
library(MSstats)
library(plyr); library(dplyr)
library(tidyverse)

#A couple of notes before we get started. I hate how reliant on setwd() this is, I know it's inefficient setting absolute 
#file paths and that this will break if anyone else tries to run it on their computer without respecifying the path. 
#I hope to replace these with relative path files in the future, but this works for now. 
#I also have a lot of work to do with comparison figure generation. The QC files function extremely well, but I'm having 
#some trouble getting the comparisons to do what I'd like. This will likely come together once I figure out an efficient way
#to generate a consolidated summed fractions file to work from. 

#TO DO:
  #Change outputs to Results file. DONE
  #Create summed fractions consolidated protein level output WRITTEN BUT NOT TESTED
  #Create summed fractions consolidated peptide level output (MIGHT JUST IMPORT AND RENAME EACH PEPTIDE.TSV)
  #Generate comparison figures from protein level output WRITTEN BUT NOT TESTED
  #Integrate solubility analysis and generate solubility profile heatmap from that data
  #Integrate GSEA and generate Bubble Plots

#Cell Fraction processing
rootDir <- "C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Data/Cell/" # Specify the path of the directory containing MSstats.csv.
print(str_c("Using IonQuant's result from ", rootDir))

# Read MSstats.csv file.
raw_Cell <- read_csv(str_c(rootDir, "MSstats.csv"), na = c("", "NA", "0"))
raw_Cell$ProteinName <- factor(raw_Cell$ProteinName)
raw_Cell$PeptideSequence <- factor(raw_Cell$PeptideSequence)

# Change root directory for MSstats
print(str_c("Root DIR: ", rootDir))
setwd(rootDir)

# Processing the data using MSstats
processedData_Cell <- dataProcess(raw_Cell)

#Generating QC Plots of Cell Fraction
#Changing Working Directory to Results QC Plots folder
setwd("C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Results/QC Plots/Cell")
dataProcessPlots(data = processedData_Cell, type = "QCplot", ylimUp = 35, width = 5, height = 5, address = "Cell_Individual_")
dataProcessPlots(data = processedData_Cell, type = "QCplot", which.Protein = "allonly", ylimUp = 35, width = 5, height = 5, address = "Cell_Total_")
dataProcessPlots(data = processedData_Cell, type = "Profileplot", ylimUp = 35,
                 featureName="NA", width = 5, height = 5, address = "Cell_Profile_")
dataProcessPlots(data = processedData_Cell, type = "Conditionplot", ylimUp = 35,
                 width = 5, height = 5, address = "Cell_Condition_")

#Generating Cell Fraction Protein-Intensity Table 

print("Preparing a protein-intensity table for the Cellular Fraction.")
x <- processedData_Cell$RunlevelData
y <- data.frame(matrix(nrow=length(unique(x$Protein)), ncol=length(unique(x$SUBJECT_NESTED)) + 1))
colnames(y) <- c("Protein", unique(as.character(x$SUBJECT_NESTED)))
y$Protein <- unique(x$Protein)
for (protein in y$Protein) {
  for (subject in unique(x$SUBJECT_NESTED)) {
    if (length(x[x$Protein == protein & x$SUBJECT_NESTED == subject,]$LogIntensities) > 0) {
      y[y$Protein == protein, match(subject, colnames(y))] <- x[x$Protein == protein & x$SUBJECT_NESTED == subject,]$LogIntensities
    }
  }
}

#Write the .csv output to the Spreadsheets Output folder in Results
write_csv(y, str_c("C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Results/Spreadsheets Outputs/", "MSstats_protein_Cell.csv"), na = "")


#sECM Fraction processing
rootDir <- "C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Data/sECM/" # Specify the path of the directory containing MSstats.csv.
print(str_c("Using IonQuant's result from ", rootDir))

# Read MSstats.csv file.
raw_sECM <- read_csv(str_c(rootDir, "MSstats.csv"), na = c("", "NA", "0"))
raw_sECM$ProteinName <- factor(raw_sECM$ProteinName)
raw_sECM$PeptideSequence <- factor(raw_sECM$PeptideSequence)

# Change root directory for MSstats
print(str_c("Root DIR: ", rootDir))
setwd(rootDir)

# Processing the data using MSstats
processedData_sECM <- dataProcess(raw_sECM)

#Generating QC Plots of sECM Fraction
#Changing Working Directory to Results QC Plots folder
setwd("C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Results/QC Plots/sECM")
dataProcessPlots(data = processedData_sECM, type = "QCplot", ylimUp = 35, width = 5, height = 5, address = "sECM_Individual_")
dataProcessPlots(data = processedData_sECM, type = "QCplot", which.Protein = "allonly", ylimUp = 35, width = 5, height = 5, address = "sECM_Total_")
dataProcessPlots(data = processedData_sECM, type = "Profileplot", ylimUp = 35,
                 featureName="NA", width = 5, height = 5, address = "sECM_Profile_")
dataProcessPlots(data = processedData_sECM, type = "Conditionplot", ylimUp = 35,
                 width = 5, height = 5, address = "sECM_Condition_")

#Generating sECM Fraction Protein-Intensity Table 
print("Preparing a protein-intensity table for the sECM Fraction.")
x <- processedData_sECM$RunlevelData
y <- data.frame(matrix(nrow=length(unique(x$Protein)), ncol=length(unique(x$SUBJECT_NESTED)) + 1))
colnames(y) <- c("Protein", unique(as.character(x$SUBJECT_NESTED)))
y$Protein <- unique(x$Protein)
for (protein in y$Protein) {
  for (subject in unique(x$SUBJECT_NESTED)) {
    if (length(x[x$Protein == protein & x$SUBJECT_NESTED == subject,]$LogIntensities) > 0) {
      y[y$Protein == protein, match(subject, colnames(y))] <- x[x$Protein == protein & x$SUBJECT_NESTED == subject,]$LogIntensities
    }
  }
}

write_csv(y, str_c("C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Results/Spreadsheets Outputs/", "MSstats_protein_sECM.csv"), na = "")


#iECM Fraction processing
rootDir <- "C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Data/iECM/" # Specify the path of the directory containing MSstats.csv.
print(str_c("Using IonQuant's result from ", rootDir))

# Read MSstats.csv file.
raw_iECM <- read_csv(str_c(rootDir, "MSstats.csv"), na = c("", "NA", "0"))
raw_iECM$ProteinName <- factor(raw_iECM$ProteinName)
raw_iECM$PeptideSequence <- factor(raw_iECM$PeptideSequence)

# Change root directory for MSstats
print(str_c("Root DIR: ", rootDir))
setwd(rootDir)

# Processing the data using MSstats
processedData_iECM <- dataProcess(raw_iECM)

#Generating QC Plots of iECM Fraction
#Changing Working Directory to Results QC Plots folder
setwd("C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Results/QC Plots/iECM")
dataProcessPlots(data = processedData_iECM, type = "QCplot", ylimUp = 35, width = 5, height = 5, address = "iECM_Individual")
dataProcessPlots(data = processedData_iECM, type = "QCplot", which.Protein = "allonly", ylimUp = 35, width = 5, height = 5, address = "iECM_Total")
dataProcessPlots(data = processedData_iECM, type = "Profileplot", ylimUp = 35,
                 featureName="NA", width = 5, height = 5, address = "iECM_Profile")
dataProcessPlots(data = processedData_iECM, type = "Conditionplot", ylimUp = 35,
                 width = 5, height = 5, address = "iECM_Total")

#Generating iECM Fraction Protein-Intensity Table 
print("Preparing a protein-intensity table for the iECM Fraction.")
x <- processedData_iECM$RunlevelData
y <- data.frame(matrix(nrow=length(unique(x$Protein)), ncol=length(unique(x$SUBJECT_NESTED)) + 1))
colnames(y) <- c("Protein", unique(as.character(x$SUBJECT_NESTED)))
y$Protein <- unique(x$Protein)
for (protein in y$Protein) {
  for (subject in unique(x$SUBJECT_NESTED)) {
    if (length(x[x$Protein == protein & x$SUBJECT_NESTED == subject,]$LogIntensities) > 0) {
      y[y$Protein == protein, match(subject, colnames(y))] <- x[x$Protein == protein & x$SUBJECT_NESTED == subject,]$LogIntensities
    }
  }
}

write_csv(y, str_c("C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Results/Spreadsheets Outputs/", "MSstats_protein_iECM.csv"), na = "")

#Next need to consolidate all Fractional Intensity Sheets
#Read in all fractional MSStats_Protein.csv files

#Changing Working Directory to the Spreadsheets Output folder in Results
setwd("C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Results/Spreadsheets Outputs/")

#This next section merges all Fractional CSVs into a single file
MSStats_Protein_Combined_Fractions <- list.files(path = "C:/Users/maslankm/Documents/R/Crosbie_hCP_Timecourse/Results/Spreadsheets Outputs/",
                                                     pattern = "*.csv",full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows()
#I'm still not positive if the above section works, still need to test but don't have the data yet
#Next need to sum all identical proteins for each subject to consolidate fractions
MSStats_Protein_Consolidated_Fractions <- MSStats_Protein_Combined_Fractions %>%
  group_by(Protein,colnames())  #will likely need to figure out how to specify subject
  summarise(
    Protein_Sum = sum(colnames())
  )
#I don't think this will work. I think I need a way to specify each subject and iteratively perform this for each. 
  #a loop will be needed to perform that but we'll play around with this to see if it can work

#need to make consolidated processedData File
#Generating 1v1 comparisons
comparison1 <- matrix(c(-1,1,0),nrow = 1) #Control vs IPF
#row.names(comparison1) <- "Control-IPF"

#Control_vs_IPF_Comparison <- groupComparison(contrast.matrix = comparison1, data = processedData)

comparison2 <- matrix(c(-1,0,1),nrow = 1) #Control vs RA-ILD
#row.names(comparison2) <- "Control-RA-ILD"

#Control_vs_RA_ILD_Comparison <- groupComparison(contrast.matrix = comparison2, data = processedData)

comparison3 <- matrix(c(0,-1,1),nrow = 1) #IPF vs RA-ILD
#row.names(comparison3) <- "IPF-RA-ILD"

#IPF_vs_RA_ILD_Comparison <- groupComparison(contrast.matrix = comparison3, data = processedData)


#Lets try this to bring them all together
comparison <- rbind(comparison1,comparison2,comparison3)
rownames(comparison) <- c('Control-IPF','Control-RA_ILD','IPF-RA_ILD') #be sure to check that labels are ordered correctly

Results <- groupComparison(processedData, contrast.matrix = comparison)

StatisticsReport <- Results$ComparisonResults
write.csv(StatisticsReport, file = 'StatisticsReport.csv')

SignificantResults <- StatisticsReport[StatisticsReport$pvalue < 0.1, ]
write.csv(SignificantResults, file = 'SignificantResults')


#Still need to figure out how to compare all 3 for heatmaps + Comparison Plots
#Doubt this'll work but worth a try
#comparison4 <- matrix(c(-1,1,1),nrow = 1) #Control vs IPF & RA-ILD
#row.names(comparison4) <- "Control-IPF-RA_ILD"

#Control_vs_IPF_and_RA_ILD_Comparison <- groupComparison(contrast.matrix = comparison4, data = processedData)
#Shit didn't work

#Generate Plots
groupComparisonPlots(data = StatisticsReport, type = "VolcanoPlot")
groupComparisonPlots(data = StatisticsReport, type = "Heatmap")
groupComparisonPlots(data = StatisticsReport, type = "Heatmap", numProtein = 50)
groupComparisonPlots(data = StatisticsReport, type = "ComparisonPlot")

#Well shit, these didn't work
#groupComparisonPlots(data = Control_vs_IPF_and_RA_ILD_Comparison$ComparisonResult, type = "Heatmap", numProtein = 50)
#groupComparisonPlots(data = Control_vs_IPF_and_RA_ILD_Comparison$ComparisonResult, type = "ComparisonPlot")

#Quantification of Protein Abundances
ConsolidatedProcessedData <- rbind(processedData_Cell,processedData_sECM,processedData_iECM)
SampleQuant <- quantification(ConsolidatedProcessedData)
GroupQuant <- quantification(ConsolidatedProcessedData , type = 'group')
