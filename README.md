<p align="center">
  <a>
    <img width="30%" src="https://github.com/chill3456/Quick_qPCR_Analysis/blob/master/assets/graphic.png" alt="graphic describing Quick_qPCR_analysis.R png">
  </a>
  <h1 align="center">Quick qPCR Analysis</h1>
</p>




# Description 

Quick_qPCR_Analysis.R is an R script that allows for analysis of qPCR using a single R script that can be ran on the command line. 

# Required packages 

`install.packages("tidyverse")`
`install.packages("ggplot2")`
`install.packages("dplyr")`

# Required Files 
A set up txt file with plate position corresponding to sample and primer names. Plate positions are in "Well" column, samples under "Sample" column, and primers under "Primer" column. There are 96 and 384 examples. 

`setUpFile_96_well.txt`
`setUpFile_384_well.txt`

Set up file can be generated using a .xlsx file that allows you to input sample and primer names according to their positon on the plate. Plate positions are put in "Well" column, samples under "Sample" column, and primers under "Primer" column. There are 96 and 384 examples. 

`platesetup_to_setUpFile_96_well.xlsx`
`platesetup_to_setUpFile_384_well.xlsx`

It also requires a results txt file with plate position corresponding to CT or CQ values. Plate positions should be in the "Well" column. The example file is from the BioRad platform and when read in the CQ values are read in under the Content Column if CQ or CT values are under different column a different columnn can be specified as an argument. 

`ctResultsFile_96_well.txt`
`ctResultsFile_384_well.txt`


# Running the script 

The script can be run on the command line. 

`Rscript Quick_qPCR_Analysis.R --ctResultsFilePath=/example_files/ctResultsFile_96_well.txt --setUpFilePath=/example_files/setUpFile_96_well.txt --outputFilePath=/outputExamples --controlSample=control_sample --controlPrimer=control_primer --sampleOfInterest=sample_of_interest1 --primerOfInterest=primer_of_interest"`

# Outputs

The script outputs a table that labels the CT values with a primer and sample name in outputExamples/RAW_CT_results_labeled.txt 

| CT          | Sample              | Primer         |
| ----------- | ------------------- | -------------- |
| 32.11672579	| control_sample	    | control_primer |
| 31.69335192	| control_sample	    | control_primer |
| 30.86998086	| control_sample	    | control_primer |
| 31.65992152	| control_sample	    | control_primer |
| 18.09801643	| sample_of_interest1 |	control_primer |
| 18.09438969	| sample_of_interest1 |	control_primer |
| 18.08867095	| sample_of_interest1 |	control_primer |
| 18.18379642	| sample_of_interest1 |	control_primer |
| 17.79293803	| sample_of_interest2 |	control_primer |

Calculated ddCT values for each sample are output in outputExamples/sample_of_interest1_primer_of_interest_control_sample_control_primer_ddCT_2_table.txt

| CT	        | Sample              | Primer             | dCT         | ddCT         |	ddCT_2      |
| ----------- | ------------------- | ------------------ | ----------- | ------------ | ----------- |
| 33.45958998 |	control_sample      |	primer_of_interest | 1.874594956 | -0.025127303 |	1.017569478 |
| 33.37342965 |	control_sample      |	primer_of_interest | 1.788434629 | -0.111287629 |	1.080191896 |
| 33.6905051  |	control_sample      |	primer_of_interest | 2.105510076 | 0.205787818  |	0.867065076 |
| 33.41534439 |	control_sample      |	primer_of_interest | 1.830349372 | -0.069372886 |	1.04926049  |
| 22.39176874 |	sample_of_interest1 |	primer_of_interest | 4.275550373 | 2.375828114  |	0.19266573  |
| 22.30511258 |	sample_of_interest1 |	primer_of_interest | 4.188894204 | 2.289171946  |	0.20459291  |
| 22.30515157 |	sample_of_interest1 |	primer_of_interest | 4.1889332   | 2.289210941  |	0.20458738  |
| 22.32831807 |	sample_of_interest1 |	primer_of_interest | 4.212099702 | 2.312377443  | 0.201328394 |

Also outputs a summary file that shows average ddCT values plotted outputExamples/sample_of_interest1_primer_of_interest_control_sample_control_primer_summary_stats.txt

| Sample            	| Primer             | ddCT_2	     | sd          |
| ------------------- | ------------------ | ----------- | ----------- |
| control_sample	    | primer_of_interest | 1.003521735 | 0.094495337 |
| sample_of_interest1	| primer_of_interest | 0.200793603 | 0.005632518 |

A graph is output showing ddCT values 

<img src='https://github.com/chill3456/Quick_qPCR_Analysis/blob/master/assets/sample_of_interest1_primer_of_interest_control_sample_control_primer_plot.png' width = '50%' alt="graph showing results of Quick_qPCR_analysis.R png">


# Running with multiple plates

The script can also be run outside of the command line by assigning values to the arugment variables. 

`Quick_qPCR_analysis_no_arguments.R`

You can also assign diffent variables such as those needed to use multiple plates. 

`ctResultsPlate1FilePath_argument <- "example_files/ctResultsFile_96_well.txt"`

`ctResultsPlate2FilePath_argument <- "example_files/ctResultsFile_384_well.txt"`

`setUpPlate1FilePath_argument<- "example_files/setUpFile_96_well.txt"`

`setUpPlate2FilePath_argument<- "example_files/setUpFile_384_well.txt"`


`output_file_path_argument <- "outputExamplesNoArgs/"`

`CT_value_Column_argument <- "Content" `

Also Define samples and primers to test
`control_sample_argument <- "control_sample"`

`control_primer_argument <- "control_primer"`

`control_sample1_argument <- "control_sample1"`

`sample_of_interest1_argument <- "sample_of_interest1"`

`primer_of_interest1_argument <- "primer_of_interest"`

`sample_of_interest2_argument <- "sample_of_interest2"`

`primer_of_interest2_argument <- "primer_of_interest2"`


When running merge_qpcr_files assign the function to a variable which will contain the labeled CT values. 

`merged_data_plate_1 <- merge_qpcr_files(setUpPlate1FilePath_argument, ctResultsPlate1FilePath_argument, output_file_path_argument)`

`merged_data_plate_2 <- merge_qpcr_files(setUpPlate2FilePath_argument, ctResultsPlate2FilePath_argument, output_file_path_argument)`

Those variables can then be combined. 

`merged_data <- rbind(merged_data_plate_1,merged_data_plate_2)`

This allows you to run the analysis using data from multiple plates while also saving the results to a variable

`sample_of_interest1_with_primer_of_interest1_analysis <- analyze_qPCR(merged_data, control_sample_argument, control_primer_argument, sample_of_interest1_argument, primer_of_interest1_argument, output_file_path_argument)`

`sample_of_interest2_with_primer_of_interest2_analysis <- analyze_qPCR(merged_data, control_sample1_argument, control_primer_argument, sample_of_interest2_argument, primer_of_interest2_argument, output_file_path_argument)`

The summary stats can then be combined and used to plot multiple samples or primers together. 

`objects <- mget(ls(pattern = "_analysis"))`

`all_summary_stats <- lapply(objects, function(x) x$summary_stats)`

`dplyr::bind_rows(all_summary_stats) -> combined_summary_stats`

`combined_summary_stats$Sample = forcats::fct_rev(factor(combined_summary_stats$Sample))`

`ggplot(combined_summary_stats, aes(x=Primer, y=ddCT_2, fill=Sample)) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=ddCT_2-sd, ymax=ddCT_2+sd), width=.2, position=position_dodge(.9)) + geom_hline(yintercept = 1, linetype = "dashed")`

`file_name_CT <-  "combined_summary_stats.pdf"`

`ggsave(file = file.path(output_file_path_argument, file_name_CT))`

<img src='https://github.com/chill3456/Quick_qPCR_Analysis/blob/master/assets/combined_summary_stats.png' width = '50%' alt="graph showing results of Quick_qPCR_analysis_no_arguments.R png">
