#!/usr/bin/env Rscript


library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)

# print warnings as they happen instead of collecting them for after a loop ends
options(warn=1)

# define valid parameters
parameters <- list(
  ctResultsFilePath=list("ctResultsFilePath_argument", "file", "ctResultsFile_96_well.txt", T),
  setUpFilePath=list("setUpFilePath_argument", "file", "setUpFile_96_well.txt", T),
  outputFilePath=list("output_file_path_argument", "string", "/outputFolder", T),
  controlSample=list("control_sample_argument", "string", "control_sample", T),
  controlPrimer=list("control_primer_argument", "string", "control_primer", T),
  sampleOfInterest=list("sample_of_interest_argument", "string", "sample_of_interest", T),
  primerOfInterest=list("primer_of_interest_argument", "string", "primer_of_interest", T),
  CTvalueColumn=list("CT_value_Column_argument", "string", "Content")
)



# print help if necessary
args <- commandArgs(trailingOnly=T)
if (any(grepl("^--help", args)) || length(args) == 0) {
  usage <- "Usage: Quick_qPCR_Analysis.R"
  for (parameter in names(parameters)) {
    usage <- paste0(usage, " ")
    if (length(parameters[[parameter]]) <= 3 || !parameters[[parameter]][[4]])
      usage <- paste0(usage, "[")
    usage <- paste0(usage, "--", parameter, "=", parameters[[parameter]][[3]])
    if (length(parameters[[parameter]]) <= 3 || !parameters[[parameter]][[4]])
      usage <- paste0(usage, "]")
  }
  message(usage)
  quit("no", ifelse(length(args) == 0, 1, 0))
}


# make sure mandatory arguments are present
for (parameter in names(parameters))
  if (length(parameters[[parameter]]) > 3 && parameters[[parameter]][[4]])
    if (!any(grepl(paste0("^--", parameter, "="), args), perl=T))
      stop(paste0("Missing mandatory argument: --", parameter))

# set default values
for (parameter in names(parameters))
  assign(parameters[[parameter]][[1]], ifelse(parameters[[parameter]][[2]] == "file", "", parameters[[parameter]][[3]]))

# parse command-line parameters
for (arg in args) {
  argName <- sub("=.*", "", sub("^--", "", arg, perl=T), perl=T)
  argValue <- sub("^[^=]*=", "", arg, perl=T)
  if (!(argName %in% names(parameters)) || !grepl("^--", arg, perl=T))
    stop(paste("Unknown parameter:", arg))
  if (parameters[[argName]][[2]] == "bool") {
    if (argValue %in% c("TRUE", "T", "FALSE", "F")) {
      assign(parameters[[argName]][[1]], as.logical(argValue))
    } else {
      stop(paste0("Invalid argument to --", argName))
    }
  } else if (parameters[[argName]][[2]] == "string") {
    assign(parameters[[argName]][[1]], argValue)
  } else if (parameters[[argName]][[2]] == "numeric") {
    if (is.na(suppressWarnings(as.numeric(argValue))))
      stop(paste0("Invalid argument to --", argName))
    assign(parameters[[argName]][[1]], as.numeric(argValue))
  } else if (parameters[[argName]][[2]] == "file") {
    if (file.access(argValue) == -1)
      stop(paste("Cannot read file:", argValue))
    assign(parameters[[argName]][[1]], argValue)
  }
}

# check if required packages are installed
if (!suppressPackageStartupMessages(require(tidyverse)))
  warning("Package 'tidyverse' is not installed.")
if (!suppressPackageStartupMessages(require(ggplot2)))
  warning("Package 'ggplot2' is not installed.")
if (!suppressPackageStartupMessages(require(dplyr)))
  warning("Package 'dplyr' is not installed.")

#functions 

#to label CT values by merging the setUpFile and ctResultsFile together 
merge_qpcr_files <- function(setUpFilePath, ctResultsFilePath, output_file_path) {
  # Read in the two files
  read.table(setUpFilePath, header = TRUE, stringsAsFactors = FALSE) -> data1
  #reading in BioRad "Quantification Cq Results_0.txt" file which puts CQ values under Content Column when read in which is then renamed to CT. Adjust as need to put CT values under CT column 
  read.table(ctResultsFilePath, header = TRUE,stringsAsFactors = FALSE, fill = TRUE) %>% dplyr::rename(CT = CT_value_Column_argument) -> data2
  
  # Merge the files by the "well" column
  merged_data <- merge(data1[,c("Well", "Sample", "Primer", grep("CT", colnames(data1), value = TRUE))],
                       data2[,c("Well", grep("CT", colnames(data2), value = TRUE))],
                       by = "Well", all.x = TRUE)
  
  
  # Return the merged data with only the "ct", "gene", and "sample" columns
  merged_data[,c("CT", "Sample", "Primer")]   -> merged_data
  
  #eliminate rows that have blank as a value 
  
  merged_data <- subset(merged_data, Sample != "blank")
  
  file_name <- basename(ctResultsFilePath)
  
  file_name_CT <- paste(file_name, "RAW_CT_results_labeled.txt", sep = "_")
  write.table(merged_data, file = file.path(output_file_path, file_name_CT),  sep = "\t", col.names = TRUE, row.names = FALSE)
  
  return(merged_data)
}

#to generate standard devations and averages
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#to analyze the qPCR results 
analyze_qPCR <- function(merged_data, control_sample, control_primer, sample_of_interest, primer_of_interest, output_file_path) {
  # Subset the table to only include the control and sample of interest
  subset1 <- subset(merged_data, Sample == control_sample & Primer == control_primer)
  subset2 <- subset(merged_data, Sample == sample_of_interest & Primer == control_primer)
  subset3 <- subset(merged_data, Sample == control_sample & Primer == primer_of_interest)
  subset4 <- subset(merged_data, Sample == sample_of_interest & Primer == primer_of_interest)
  
  average_subset_1 <- mean(subset1$CT)
  average_subset_2 <- mean(subset2$CT)
  
  subset3 %>% mutate(dCT = subset3$CT - average_subset_1) -> subset3_dCT
  subset4 %>% mutate(dCT = subset4$CT - average_subset_2) -> subset4_dCT
  
  average_subset_3 <- mean(subset3_dCT$dCT)
  
  subset3_dCT %>% mutate(ddCT = subset3_dCT$dCT - average_subset_3) -> subset3_ddCT
  subset4_dCT %>% mutate(ddCT = subset4_dCT$dCT - average_subset_3) -> subset4_ddCT
  
  subset3_ddCT %>% mutate(ddCT_2 = 2^(-subset3_ddCT$ddCT)) -> subset3_ddCT
  subset4_ddCT %>% mutate(ddCT_2 = 2^(-subset4_ddCT$ddCT)) -> subset4_ddCT
  
  rbind(subset3_ddCT, subset4_ddCT) -> subset5
  
  df2 <- data_summary(subset5, varname="ddCT_2", 
                      groupnames=c("Sample", "Primer"))
  
  file_name <- paste(sample_of_interest, primer_of_interest, control_sample, control_primer, "ddCT_2_table.txt", sep = "_")
  write.table(subset5, file = file.path(output_file_path, file_name), sep = "\t", row.names = FALSE)
  
  file_name <- paste(sample_of_interest, primer_of_interest, control_sample, control_primer, "summary_stats.txt", sep = "_")
  write.table(df2, file = file.path(output_file_path, file_name), sep = "\t", row.names = FALSE)
  
  
  p<- ggplot(df2, aes(x=Sample, y=ddCT_2, fill=Primer)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=ddCT_2-sd, ymax=ddCT_2+sd), width=.2,
                  position=position_dodge(.9)) +
    geom_hline(yintercept = 1, linetype = "dashed")
  
  file_name <- paste(sample_of_interest, primer_of_interest, control_sample, control_primer, "plot.pdf", sep = "_")
  ggsave(filename = file.path(output_file_path, file_name), plot = p , width = 10, height = 5)
  
  return(list(ddCT_2_table = subset5, plot = p, summary_stats = df2))
}


#running the functions 

merged_data <- merge_qpcr_files(setUpFilePath_argument, ctResultsFilePath_argument, output_file_path_argument)

analyze_qPCR(merged_data, control_sample_argument, control_primer_argument, sample_of_interest_argument, primer_of_interest_argument, output_file_path_argument)





