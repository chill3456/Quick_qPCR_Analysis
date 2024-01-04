library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)





#define paths to results and set up documents 
ctResultsPlate1FilePath_argument <- "example_files/ctResultsFile_96_well.txt"
ctResultsPlate2FilePath_argument <- "example_files/ctResultsFile_384_well.txt"

setUpPlate1FilePath_argument<- "example_files/setUpFile_96_well.txt"
setUpPlate2FilePath_argument<- "example_files/setUpFile_384_well.txt"

output_file_path_argument <- "outputExamplesNoArgs/"

#define samples and primers to test
control_sample_argument <- "control_sample"

control_primer_argument <- "control_primer"

control_sample1_argument <- "control_sample1"

sample_of_interest1_argument <- "sample_of_interest1"

primer_of_interest1_argument <- "primer_of_interest"

sample_of_interest2_argument <- "sample_of_interest2"

primer_of_interest2_argument <- "primer_of_interest2"

CT_value_Column_argument <- "Content"

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
  ggsave(filename = file.path(output_file_path, file_name), plot = p, , width = 10, height = 5)
  
  return(list(ddCT_2_table = subset5, plot = p, summary_stats = df2))
}


#running the functions 
merged_data_plate_1 <- merge_qpcr_files(setUpPlate1FilePath_argument, ctResultsPlate1FilePath_argument, output_file_path_argument)

merged_data_plate_2 <- merge_qpcr_files(setUpPlate2FilePath_argument, ctResultsPlate2FilePath_argument, output_file_path_argument)

merged_data <- rbind(merged_data_plate_1,merged_data_plate_2)

sample_of_interest1_with_primer_of_interest1_analysis <- analyze_qPCR(merged_data, control_sample_argument, control_primer_argument, sample_of_interest1_argument, primer_of_interest1_argument, output_file_path_argument)

sample_of_interest2_with_primer_of_interest2_analysis <- analyze_qPCR(merged_data, control_sample1_argument, control_primer_argument, sample_of_interest2_argument, primer_of_interest2_argument, output_file_path_argument)


#can plot together 

objects <- mget(ls(pattern = "_analysis"))

all_summary_stats <- lapply(objects, function(x) x$summary_stats)

dplyr::bind_rows(all_summary_stats) -> combined_summary_stats

combined_summary_stats$Sample = forcats::fct_rev(factor(combined_summary_stats$Sample))

ggplot(combined_summary_stats, aes(x=Primer, y=ddCT_2, fill=Sample)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=ddCT_2-sd, ymax=ddCT_2+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 1, linetype = "dashed")

file_name_CT <-  "combined_summary_stats.pdf"

ggsave(file = file.path(output_file_path_argument, file_name_CT), width = 10, height = 5)




