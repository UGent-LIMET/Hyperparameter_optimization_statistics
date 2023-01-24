# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Part II: multivariate analysis


##########R Pipeline - Part II: multivariate analysis##########
print(Sys.time())
start_time <- Sys.time()
print("R pipeline - Part II: multivariate analysis - start!")
# Part II: multivariate analysis

## data_loading
setwd(path_data_in)
COLLUMN_NR_COMPARISON1 <- COLLUMN_NR_LAST_BEFORE_COMPARISONS + 1        #column number of the column 'Comparison1' for pairwise comparison (OPLSDA)
COLLUMN_NR_MULTIPLE_COMPARISON1 <- COLLUMN_NR_COMPARISON1 + AMOUNT_OF_COMPARISONS         #column number of the column 'MultipleComparison1' for multiple comparison (PLSDA, LIMMA)
COLLUMN_NR_PROJECTION1 <- COLLUMN_NR_MULTIPLE_COMPARISON1 + AMOUNT_OF_MULTIPLE_COMPARISONS
COLLUMN_NR_START_VARIABLES <- COLLUMN_NR_PROJECTION1 + AMOUNT_OF_PROJECTIONS

sampleMetadata <- load_sampleMetadata(INPUT_SAMPLES)
#Add "X" to sampleNames
sampleMetadata[,1] <- paste(replicate(nrow(sampleMetadata),"X"), sampleMetadata[,1], sep="")
check_nrow_SM <- nrow(sampleMetadata)

if (INPUT_VARIABLES == VARIABLEMETADATA_FROM_PIPELINE){
  INPUT_VARIABLES <- paste(name_project, '_variableMetadata.txt', sep="")
  if(exists("COLLUMN_NR_START_SAMPLES") == FALSE){ 		#if after merge, will be given value 21, so do not change
    COLLUMN_NR_START_SAMPLES <- 20  #always 20 (auto and manual must be same format)
  }
}
variableMetadata <- load_variableMetadata(INPUT_VARIABLES)
#note: automatically adds X to samples at this point


## set directory to output
setwd(path_data_out)

## merge variables from variableMetadata into the sampleMetadata
### make sampleMatrix
suppressMessages(library(data.table))
variableMetadata_from_start_samples <- subset(variableMetadata, select = -c(2:(COLLUMN_NR_START_SAMPLES-1))) #remove info exept CompID
sampleMatrix <- transpose(variableMetadata_from_start_samples)
colnames(sampleMatrix) <- variableMetadata_from_start_samples[ ,1]
sampleMatrix <- as.data.frame(sapply(sampleMatrix, as.numeric))
sampleMatrix$SampleName <- colnames(variableMetadata_from_start_samples)
sampleMatrix <- sampleMatrix[-1,] #remove first row (colnames compids + samplename; so 1 extra col than no of variables), is captured in colnames
#write_dataframe_as_txt_file(sampleMatrix, 'sampleMatrix.txt')

#merge variable intensities from compIDs to correct sample (if order not same), output order is sorted by samplename
sampleMetadata <- merge_accoding_to_SampleName(sampleMetadata, sampleMatrix)
write_dataframe_as_txt_file(sampleMetadata, 'sampleMetadata_variableMatrix_merged.txt')


##Check sampleMatrx 
#see if same number of rows (amount of sampleNames) as when loaded prev, so merge was succesfull
if(nrow(sampleMetadata) != check_nrow_SM){
  stop("ERROR: Part II: multivariate analysis stopped because SampleMetadata and VariableMetadata are incompatible to merge into correct sampleMatrix.")
} 

#########Variable filtering###########
## optional: Delete variables with coeficent of variation of QCs > 30%
if(COEFICIENT_OF_VARIANCE == COEFICIENT_OF_VARIANCE_30){
  stopifnot("QC" %in% sampleMetadata$Type |"IQC" %in% sampleMetadata$Type)
  QC_metadata <- sampleMetadata[sampleMetadata$Type == 'QC' | sampleMetadata$Type == 'IQC',] 
  
  CV_values <- apply(QC_metadata[,COLLUMN_NR_START_VARIABLES:length(QC_metadata)], 2, function(x) sd(x)/mean(x))
  matrix_with_CV_values <- data.frame(x=colnames(QC_metadata)[COLLUMN_NR_START_VARIABLES:length(QC_metadata)],
                                      y=CV_values)
  names(matrix_with_CV_values) <- c("CompID", "CV_values")
  retained_variables_CV <- matrix_with_CV_values[(matrix_with_CV_values[,2]<=0.3),]
  retained_variables_CV <- retained_variables_CV[,1]
  
  all_matrix <- from_df_to_matrix(sampleMetadata)
  all_matrix <- as.data.frame(all_matrix)
  
  matrix_filter_coeficient_variance <- all_matrix[, names(all_matrix) %in% (retained_variables_CV)]
  
  sampleMetadata_until_start_variables <- subset(sampleMetadata, select = c(1:(COLLUMN_NR_START_VARIABLES-1)))
  sampleMetadata_filter_coeficient_variance <- cbind(sampleMetadata_until_start_variables, matrix_filter_coeficient_variance)
  write_dataframe_as_txt_file(sampleMetadata_filter_coeficient_variance, 'sampleMetadata_filter_coeficient_variance.txt')
  
  sampleMetadata <- sampleMetadata_filter_coeficient_variance
}
if(COEFICIENT_OF_VARIANCE == COEFICIENT_OF_VARIANCE_0){
  sampleMetadata <- sampleMetadata
}


## optional: retain variable in present in 80% of the QCs samples
if(QC_PREFILTER == FILTER_QC_POOL_80){
  stopifnot("QC" %in% sampleMetadata$Type |"IQC" %in% sampleMetadata$Type)
  QC_metadata <- sampleMetadata[sampleMetadata$Type == 'QC' | sampleMetadata$Type == 'IQC',] 
  number_QCs <- nrow(QC_metadata)
  perc_present <- colSums(QC_metadata[,COLLUMN_NR_START_VARIABLES:length(QC_metadata)] != 0)/number_QCs
  matrix_with_perc_present <- data.frame(x=colnames(QC_metadata)[COLLUMN_NR_START_VARIABLES:length(QC_metadata)],
                                         y=perc_present)
  names(matrix_with_perc_present) <- c("CompID", "perc_present")
  retained_variables_QCpool <- matrix_with_perc_present[(matrix_with_perc_present[,2]>=0.8),]
  retained_variables_QCpool <- retained_variables_QCpool[,1] #compIDs of good ones
  
  all_matrix <- from_df_to_matrix(sampleMetadata)
  all_matrix <- as.data.frame(all_matrix)
  
  matrix_filter_QCpool <- all_matrix[, names(all_matrix) %in% (retained_variables_QCpool)]
  
  sampleMetadata_until_start_variables <- subset(sampleMetadata, select = c(1:(COLLUMN_NR_START_VARIABLES-1)))
  sampleMetadata_filter_QCpool <- cbind(sampleMetadata_until_start_variables, matrix_filter_QCpool)
  write_dataframe_as_txt_file(sampleMetadata_filter_QCpool, 'sampleMetadata_filter_QCpool.txt')
  
  sampleMetadata <- sampleMetadata_filter_QCpool
}
if(QC_PREFILTER == NO_FILTER_QC_POOL){
  sampleMetadata <- sampleMetadata
}

####################



## check if data has a normal distribution per Type(s)
### Check normality of all variables across all the samples (blanco, QC, STD, sample)
all_matrix <- from_df_to_matrix(sampleMetadata)

amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(all_matrix)
total_amount_of_variables <- ncol(all_matrix)
percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)

report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across all the samples (", percentage, "%).", sep="")
name_report_normality <- paste(name_project,'_skewkurtplot_normalities.txt', sep="")
append_result_to_report(report_check_normality, name_report_normality)

skewkurtplot <- plot_skewkurt(all_matrix)
#pdf(paste(name_project,'_skewkurtplot_all.pdf', sep=""),height=5,width=7)
#plot(skewkurtplot)
#dev.off()

png(paste(name_project,'_skewkurtplot_all.png', sep=""), width=7, height=5, units="in", res=150)
plot(skewkurtplot)
dev.off()

### Check normality of all variables across all biological samples + QC or IQC
samples_QC_metadata <- sampleMetadata[sampleMetadata$Type == 'Sample' | sampleMetadata$Type == 'QC' | sampleMetadata$Type == 'IQC',] 
samples_QC_matrix <- from_df_to_matrix(samples_QC_metadata)

amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(samples_QC_matrix)
total_amount_of_variables <- ncol(samples_QC_matrix)
percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)

report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the samples and QCs (", percentage, "%).", sep="")
append_result_to_report(report_check_normality, name_report_normality)

skewkurtplot <- plot_skewkurt(samples_QC_matrix)
#pdf(paste(name_project,'_skewkurtplot_samples_QC.pdf', sep=""),height=5,width=7)
#plot(skewkurtplot)
#dev.off()

png(paste(name_project,'_skewkurtplot_samples_QC.png', sep=""), width=7, height=5, units="in", res=150)
plot(skewkurtplot)
dev.off()

### Check normality of the biological samples
samples_metadata <- sampleMetadata[sampleMetadata$Type == 'Sample',] 
samples_matrix <- from_df_to_matrix(samples_metadata)

amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(samples_matrix)
total_amount_of_variables <- ncol(samples_matrix)
percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)

report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the biological samples (", percentage, "%).", sep="")
append_result_to_report(report_check_normality, name_report_normality)

skewkurtplot <- plot_skewkurt(samples_matrix)
#pdf(paste(name_project,'_skewkurtplot_samples.pdf', sep=""),height=5,width=7)
#plot(skewkurtplot)
#dev.off()

png(paste(name_project,'_skewkurtplot_samples.png', sep=""), width=7, height=5, units="in", res=150)
plot(skewkurtplot)
dev.off()

### Check ISTD = internal standards 
ISTD_metadata <- sampleMetadata[sampleMetadata$Type == 'ISTD',]
if(nrow(ISTD_metadata) == 0){
  report_check_normality <- paste("There are no ISTDs present to check the normality for.", sep="")
  append_result_to_report(report_check_normality, name_report_normality)
}
if(nrow(ISTD_metadata) > 0){
  ISTD_matrix <- from_df_to_matrix(ISTD_metadata)
  
  amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(ISTD_matrix)
  total_amount_of_variables <- ncol(ISTD_matrix)
  percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)
  
  report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the ISTDs (", percentage, "%).", sep="")
  append_result_to_report(report_check_normality, name_report_normality)
  
  skewkurtplot <- plot_skewkurt(ISTD_matrix)
  #pdf(paste(name_project,'_skewkurtplot_ISTD.pdf', sep=""),height=5,width=7)
  #plot(skewkurtplot)
  #dev.off()
  
  png(paste(name_project,'_skewkurtplot_ISTD.png', sep=""), width=7, height=5, units="in", res=150)
  plot(skewkurtplot)
  dev.off()
}

### Check STD = internal standard mixture 
STD_metadata <- sampleMetadata[sampleMetadata$Type == 'STD',]
if(nrow(STD_metadata) == 0){
  report_check_normality <- paste("There are no STDs present to check the normality for.", sep="")
  append_result_to_report(report_check_normality, name_report_normality)
}
if(nrow(STD_metadata) > 0){
  STD_matrix <- from_df_to_matrix(STD_metadata)
  
  amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(STD_matrix)
  total_amount_of_variables <- ncol(STD_matrix)
  percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)
  
  report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the STDs (", percentage, "%).", sep="")
  append_result_to_report(report_check_normality, name_report_normality)
  
  skewkurtplot <- plot_skewkurt(STD_matrix)
  #pdf(paste(name_project,'_skewkurtplot_STD.pdf', sep=""),height=5,width=7)
  #plot(skewkurtplot)
  #dev.off()
  
  png(paste(name_project,'_skewkurtplot_STD.png', sep=""), width=7, height=5, units="in", res=150)
  plot(skewkurtplot)
  dev.off()
}

### Check normality of QCs or IQCs
QC_metadata <- sampleMetadata[sampleMetadata$Type == 'QC' | sampleMetadata$Type == 'IQC',] 
QC_matrix <- from_df_to_matrix(QC_metadata)

if(nrow(QC_matrix) == 0){
  report_check_normality <- paste("There are no (I)QCs present to check the normality for.", sep="")
  append_result_to_report(report_check_normality, name_report_normality)
}
if(nrow(QC_matrix) > 0){
  amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(QC_matrix)
  total_amount_of_variables <- ncol(QC_matrix)
  percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)
  
  report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the QCs (", percentage, "%).", sep="")
  append_result_to_report(report_check_normality, name_report_normality)
  
  skewkurtplot <- plot_skewkurt(QC_matrix)
  #pdf(paste(name_project,'_skewkurtplot_QC.pdf', sep=""),height=5,width=7)
  #plot(skewkurtplot)
  #dev.off()
  
  png(paste(name_project,'_skewkurtplot_QC.png', sep=""), width=7, height=5, units="in", res=150)
  plot(skewkurtplot)
  dev.off()
}

##########Normalisation##########
## Perform 0th normalisation with IS (optional)
if(NORMALIZE_METHOD0 == NORMALIZE_WITH_IS){ 
  #IS or internal standard (either spicked in sample or added in solvent flow instrument
  #= intensity compID 1 / intenstity Comp IS for 1 sample/spectrum
  
  #check needed input (mz value IS) ok
  stopifnot(class(IS_MZ) == "numeric") #only if present and number
  stopifnot(class(PPM) == "numeric") #only if present and number
  
  #find CompID nr in VM to use below in SM
  ppm <- PPM
  standard <- IS_MZ
  standard_down <- standard - ((ppm*standard)/1000000) 
  standard_up <- standard + ((ppm*standard)/1000000)
  variableMetadata_standard <- variableMetadata[variableMetadata$MZ >= standard_down & variableMetadata$MZ <= standard_up, ]
  
  #if >=1 found, 1select closest to mz
  if(nrow(variableMetadata_standard) > 1){  
    ppm <- ppm/2 #half ppm and retry
    standard_down <- standard - ((ppm*standard)/1000000) 
    standard_up <- standard + ((ppm*standard)/1000000)
    variableMetadata_standard <- variableMetadata[variableMetadata$MZ >= standard_down & variableMetadata$MZ <= standard_up, ]
  }
  #if no found, stop script
  if(nrow(variableMetadata_standard) == 0){  
    stop("Error: IS was not found in variableMetadata. Check presence and adjust ppm")
  }
  compid_standard <- as.numeric(variableMetadata_standard$CompID)
  
  #at end all Intensties from IS will be =1 (devide by itself per sample) = 1#col in sm with IS
  averageIS_before <- mean(as.matrix(sampleMetadata[,(COLLUMN_NR_START_VARIABLES+compid_standard-1)])) #nr compid = 1-2-3, so same + info from Sm before
  
  #devide by vector with IS intensities
  IS_intensities <- sampleMetadata[,(COLLUMN_NR_START_VARIABLES+compid_standard-1)]
  
  #slow but sure correct way devide I/I from IS per samplerow
  matrix <- NULL
  for(rownr in 1:nrow(sampleMetadata)){
    rowvector <- sampleMetadata[rownr,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] / IS_intensities[rownr]
    matrix <- rbind(matrix, rowvector)
  }
  sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- matrix
  
  #must be re-scaled to same order of scale (afrondingsfouten bij kleine getallen tijdens berek PCA, daarom all getallen maal zelfde factor)
  averageIS_after <- mean(as.matrix(sampleMetadata[,(COLLUMN_NR_START_VARIABLES+compid_standard-1)]))
  factor <- averageIS_before/averageIS_after #average after is always 1, so just *average before (but good checkpoint)
  sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- factor * sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)]
}
if(NORMALIZE_METHOD0 == NORMALIZE_NOT){
  sampleMetadata <- sampleMetadata #nothing changes
}
#write_dataframe_as_txt_file(sampleMetadata, 'IS_normalized_sampleMetadata.txt')


## Perform 1th normalisation with TIC (optional)
if(NORMALIZE_METHOD1 == NORMALIZE_WITH_TIC){ 
  #TIC or total ion count normalisation accord to literature
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3124646/
  #this normalization will be based on the sum of all intensity values in the mass spectrum (i.e., TIC)
  #= intensity compID 1 / sum(intenstity alls compIDs) for 1 sample/spectrum
  averageTIC_before <- mean(as.matrix(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)])) #sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)])
  
  #devide by sum per sample
  normalized <-  apply(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)], 1, function(x){t(x/sum(x)) }) #iterate over rows =1 (samples), but saves results as column!!! 
  sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- t(normalized) #needs to be transposed for correct tic-norm
  
  #must be re-scaled to same order of scale (afrondingsfouten bij kleine getallen tijdens berek PCA, daarom all getallen maal zelfde factor)
  averageTIC_after <- mean(as.matrix(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)]))
  factor <- averageTIC_before/averageTIC_after
  sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- factor * sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)]
}
if(NORMALIZE_METHOD1 == NORMALIZE_WITH_ITIC){
  #generaly known as 'standard scaler'
  #= intensity sample 1 / sum(intenstity all samples) for 1 metabolite
  sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- apply(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)], 2, function(x) x/sum(x)) #iterate over cols =2 (compids)
}
if(NORMALIZE_METHOD1 == NORMALIZE_WITH_MAX){
  #'sumarizednormalized' named at lca
  #= intensity sample 1 / max(intenstity all samples) for 1 metabolite
  sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- apply(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)], 2, function(x) x/max(x)) #iterate over cols =2 (compids)
}
if(NORMALIZE_METHOD1 == NORMALIZE_NOT){
  sampleMetadata <- sampleMetadata #nothing changes
}
#write_dataframe_as_txt_file(sampleMetadata, 'tic_normalized_sampleMetadata.txt')


### Normalize samples
if(NORMALIZE_METHOD == NORMALIZE_WITH_QCs){
  stopifnot("QC" %in% sampleMetadata$Type | "IQC" %in% sampleMetadata$Type)
  samples_metadata <- sampleMetadata[sampleMetadata$Type == 'Sample',]
  QC_metadata <- sampleMetadata[sampleMetadata$Type == 'QC'|sampleMetadata$Type == 'IQC',] ##also allow IQC
  normalized_samples_metadata <- normalize_with_average_QCs(samples_metadata, QC_metadata)
}
if(NORMALIZE_METHOD == NORMALIZE_WITH_IQCs){
  stopifnot("QC" %in% sampleMetadata$Type | "IQC" %in% sampleMetadata$Type)
  samples_IQC_metadata <- sampleMetadata[sampleMetadata$Type == 'IQC' | sampleMetadata$Type == 'Sample' |sampleMetadata$Type == 'QC' ,] ##also allow QC
  normalized_samples_metadata <- normalize_with_IQCs(samples_IQC_metadata)  #no info, except sampleName as last collumn
}
if(NORMALIZE_METHOD == NORMALIZE_QC_RLSC){
  stopifnot("QC" %in% sampleMetadata$Type | "IQC" %in% sampleMetadata$Type)
  #order qc == 1 in function normalize_QC_RLSC
  normalized_samples_metadata <- normalize_QC_RLSC(sampleMetadata) #need to be in same format = matrix w/o info + no inf/na
}
if(NORMALIZE_METHOD == NORMALIZE_NOT){
  samples_metadata <- sampleMetadata[sampleMetadata$Type == 'Sample',]
  normalized_samples_metadata <- normalize_not(samples_metadata) #change nothing but need to be in same format = matrix w/o info + no inf/na
}
#write_dataframe_as_txt_file(normalized_samples_metadata, 'normalized_samples_matrix.txt')

### Normalize (I)QCs (for plotting in pca-x instead of OPLSDA) 
if ("IQC" %in% sampleMetadata$Type | "QC" %in% sampleMetadata$Type){
  QC_metadata <- sampleMetadata[sampleMetadata$Type == 'QC' | sampleMetadata$Type == 'IQC',]
  if(NORMALIZE_METHOD == NORMALIZE_WITH_QCs){
    normalized_QC_metadata <- normalize_with_average_QCs(QC_metadata, QC_metadata)
  }
  if(NORMALIZE_METHOD == NORMALIZE_WITH_IQCs){
    normalized_QC_metadata <- normalize_with_average_QCs(QC_metadata, QC_metadata)
  }
  if(NORMALIZE_METHOD == NORMALIZE_QC_RLSC){
    normalized_QC_metadata <- NULL  #empty because all samples+qc + std+blanks included in loess output samples
  }
  if(NORMALIZE_METHOD == NORMALIZE_NOT){
    normalized_QC_metadata <- normalize_not(QC_metadata)
  }
}
if ("IQC" %in% sampleMetadata$Type == FALSE & "QC" %in% sampleMetadata$Type == FALSE){
  normalized_QC_metadata <- normalize_not(QC_metadata)
}

### replace intensities with normalized values (both samples + (I)QCs)
normalized_samples_QC_metadata <- rbind(normalized_samples_metadata, normalized_QC_metadata)
#write_dataframe_as_txt_file(normalized_samples_QC_metadata, 'normalized_samples_QC_matrix.txt')

sampleMetadata_until_start_variables <- subset(sampleMetadata, select = c(1:(COLLUMN_NR_START_VARIABLES-1))) 
sampleMetadata_samples_QC_until_start_variables <- sampleMetadata_until_start_variables[sampleMetadata$Type == 'Sample' | sampleMetadata$Type == 'QC' | sampleMetadata$Type == 'IQC',]

normalizedMetadata <- merge_accoding_to_SampleName(sampleMetadata_samples_QC_until_start_variables, normalized_samples_QC_metadata)
write_dataframe_as_txt_file(normalizedMetadata, 'normalizedMetadata.txt')

####################



## Check normality of biological samples after normalization
nsamples_metadata <- normalizedMetadata[normalizedMetadata$Type == 'Sample',] 
nsamples_matrix <- from_df_to_matrix(nsamples_metadata)

amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(nsamples_matrix)
total_amount_of_variables <- ncol(nsamples_matrix)
percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)

report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the samples after normalisation (", percentage, "%).", sep="")
append_result_to_report(report_check_normality, name_report_normality)

skewkurtplot <- plot_skewkurt(nsamples_matrix)
#pdf(paste(name_project,'_skewkurtplot_normalized_samples.pdf', sep=""),height=5,width=7)
#plot(skewkurtplot)
#dev.off()

png(paste(name_project,'_skewkurtplot_normalized_samples.png', sep=""), width=7, height=5, units="in", res=150)
plot(skewkurtplot)
dev.off()

## Log transformation and Pareto scaling for samples + QCs toghether
lognormalizedMetadata_samples_QC <- log_transformation(normalizedMetadata)
scalednormalizedMetadata_samples_QC <- Pareto_scaling(lognormalizedMetadata_samples_QC)
name_df <- 'normalized, log transformed and pareto scaled samples and QC Metadata.txt'
write_dataframe_as_txt_file(scalednormalizedMetadata_samples_QC, name_df)

#write transponed txt too
TscalednormalizedMetadata_samples_QC <- t(scalednormalizedMetadata_samples_QC[COLLUMN_NR_START_VARIABLES:ncol(scalednormalizedMetadata_samples_QC)])
TscalednormalizedMetadata_samples_QC <- as.data.frame(TscalednormalizedMetadata_samples_QC)
colnames(TscalednormalizedMetadata_samples_QC) <- scalednormalizedMetadata_samples_QC$SampleName
TscalednormalizedMetadata_samples_QC$CompID <- rownames(TscalednormalizedMetadata_samples_QC)
order <- c(ncol(TscalednormalizedMetadata_samples_QC), 1:(ncol(TscalednormalizedMetadata_samples_QC)-1))
TscalednormalizedMetadata_samples_QC <- TscalednormalizedMetadata_samples_QC[, order] #reorder, compid as 1st col
name_df <- 'normalized, log transformed and pareto scaled samples and QC Metadata_T.txt'
write_dataframe_as_txt_file(TscalednormalizedMetadata_samples_QC, name_df)

# cols with value O everywhere before log/pareto are now NA. NA -> 0 for ok sharipo, delete not needed at this time (below @pca are deleted)
scaledmatrix_samples_QC <- from_df_to_matrix(scalednormalizedMetadata_samples_QC)
scaledmatrix_samples_QC <- data.frame(sapply(scaledmatrix_samples_QC, function(x) ifelse(is.na(x), 0, x)))

#check normal
amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(scaledmatrix_samples_QC)
total_amount_of_variables <- ncol(scaledmatrix_samples_QC)
percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)
report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the samples after normalisation, log transformation and pareto scaling of both samples and QCs (", percentage, "%).", sep="")
append_result_to_report(report_check_normality, name_report_normality)

skewkurtplot <- plot_skewkurt(scaledmatrix_samples_QC)
#pdf(paste(name_project,'_skewkurtplot_samples_QCs_scaled_and_log-transformed.pdf', sep=""),height=5,width=7)
#plot(skewkurtplot)
#dev.off()

png(paste(name_project,'_skewkurtplot_samples_QCs_scaled_and_log-transformed.png', sep=""), width=7, height=5, units="in", res=150)
plot(skewkurtplot)
dev.off()

## Log transformation and Pareto scaling for only samples
nsamples_metadata <- normalizedMetadata[normalizedMetadata$Type == 'Sample',]
lognormalizedMetadata_samples <- log_transformation(nsamples_metadata)
scalednormalizedMetadata_samples <- Pareto_scaling(lognormalizedMetadata_samples)
name_df <- 'normalized, log transformed and pareto scaled biological samples Metadata.txt'
write_dataframe_as_txt_file(scalednormalizedMetadata_samples, name_df)

#write transponed txt too
TscalednormalizedMetadata_samples <- t(scalednormalizedMetadata_samples[COLLUMN_NR_START_VARIABLES:ncol(scalednormalizedMetadata_samples)])
TscalednormalizedMetadata_samples <- as.data.frame(TscalednormalizedMetadata_samples)
colnames(TscalednormalizedMetadata_samples) <- scalednormalizedMetadata_samples$SampleName
TscalednormalizedMetadata_samples$CompID <- rownames(TscalednormalizedMetadata_samples)
order <- c(ncol(TscalednormalizedMetadata_samples), 1:(ncol(TscalednormalizedMetadata_samples)-1))
TscalednormalizedMetadata_samples <- TscalednormalizedMetadata_samples[, order] #reorder, compid as 1st col
name_df <- 'normalized, log transformed and pareto scaled samples Metadata_T.txt'
write_dataframe_as_txt_file(TscalednormalizedMetadata_samples, name_df)

# cols with value O everywhere before log/pareto are now NA. NA -> 0 for ok sharipo, delete not needed at this time (below @pca are deleted)
scaledmatrix_samples <- from_df_to_matrix(scalednormalizedMetadata_samples)
scaledmatrix_samples <- data.frame(sapply(scaledmatrix_samples, function(x) ifelse(is.na(x), 0, x)))

#check normal
amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(scaledmatrix_samples)
total_amount_of_variables <- ncol(scaledmatrix_samples)
percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)
report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the samples after normalisation, log transformation and pareto scaling of samples (", percentage, "%).", sep="")
append_result_to_report(report_check_normality, name_report_normality)

skewkurtplot <- plot_skewkurt(scaledmatrix_samples)
#pdf(paste(name_project,'_skewkurtplot_samples_scaled_and_log-transformed.pdf', sep=""),height=5,width=7)
#plot(skewkurtplot)
#dev.off()

png(paste(name_project,'_skewkurtplot_samples_scaled_and_log-transformed.png', sep=""), width=7, height=5, units="in", res=150)
plot(skewkurtplot)
dev.off()

##########PCA1##########
## PCA plots (unsupervised)
### PCA across all biological samples+ QCs after normalization
#snelste voor pca:
#http://www.gastonsanchez.com/visually-enforced/how-to/2012/06/17/PCA-in-R/
#https://cran.r-project.org/web/packages/svd/index.html
#https://stackoverflow.com/questions/8299460/what-is-the-fastest-way-to-calculate-first-two-principal-components-in-r
#nsamples_QC_matrix <- from_df_to_matrix(normalizedMetadata) #== zelfde name as voor scaling
nsamples_QC_matrix <- scaledmatrix_samples_QC

####remove vars that are 0 and NA(N) for all samples 
#which(apply(nsamples_QC_matrix, 2, var)==0) #print bad_compIDs
if(length(which(apply(nsamples_QC_matrix, 2, var)==0)) >1){
  nsamples_QC_matrix <- nsamples_QC_matrix[ , apply(nsamples_QC_matrix, 2, var) != 0]
}
nsamples_QC_matrix <- data.frame(sapply(nsamples_QC_matrix, function(x) ifelse(is.nan(x), NA, x)))
nsamples_QC_matrix <- nsamples_QC_matrix[,which(unlist(lapply(nsamples_QC_matrix, function(x)!all(is.na(x)))))]

row.names(nsamples_QC_matrix) <- scalednormalizedMetadata_samples_QC$SampleName 
pca_samples <- prcomp(nsamples_QC_matrix)   #if not scaled yet  , scale.=TRUE
#for interpretation see:
#https://www.analyticsvidhya.com/blog/2016/03/pca-practical-guide-principal-component-analysis-python/
scores = as.data.frame(pca_samples$x) # x = 'scores' with row: sample, PC: column

PCAplot <- plot_pca(scores, nsamples_QC_matrix)
#pdf(paste(name_project,'_PCA across samples and QCs.pdf', sep=""),height=7,width=7)
#plot(PCAplot)
#dev.off()

png(paste(name_project,'_PCA across samples and QCs.png', sep=""), width=7, height=5, units="in", res=150)
plot(PCAplot)
dev.off()

try({
  png(paste(name_project,'_PCA biplot across samples and QCs.png', sep=""), width=7, height=5, units="in", res=150)
  biplot(pca_samples, scale = 0, xlabs=rep(".", ncol(pca_samples$rotation)), ylabs=rep(".", nrow(pca_samples$rotation)))
  dev.off()
})


### PCA across all biological samples after normalization
####remove vars that are 0 for all samples
nsamples_metadata <- scalednormalizedMetadata_samples[scalednormalizedMetadata_samples$Type == 'Sample',] 
nsamples_matrix <- from_df_to_matrix(nsamples_metadata)

which(apply(nsamples_matrix, 2, var)==0) #print bad_compIDs
if(length(which(apply(nsamples_matrix, 2, var)==0)) >1){
  nsamples_matrix <- nsamples_matrix[ , apply(nsamples_matrix, 2, var) != 0]
}

nsamples_matrix <- data.frame(sapply(nsamples_matrix, function(x) ifelse(is.nan(x), NA, x)))
nsamples_matrix <- nsamples_matrix[,which(unlist(lapply(nsamples_matrix, function(x)!all(is.na(x)))))]

row.names(nsamples_matrix) <- nsamples_metadata$SampleName 
pca_samples <- prcomp(nsamples_matrix)  #, scale.=TRUE
scores = as.data.frame(pca_samples$x) # x = 'scores' with row: sample, PC: column

PCAplot <- plot_pca(scores, nsamples_matrix)
#pdf(paste(name_project,'_PCA across samples.pdf', sep=""),height=7,width=7)
#plot(PCAplot)
#dev.off()

png(paste(name_project,'_PCA across samples.png', sep=""), width=7, height=5, units="in", res=150)
plot(PCAplot)
dev.off()

try({
  png(paste(name_project,'_PCA biplot across samples.png', sep=""), width=7, height=5, units="in", res=150)
  biplot(pca_samples, scale = 0, xlabs=rep(".", ncol(pca_samples$rotation)), ylabs=rep(".", nrow(pca_samples$rotation)))
  dev.off()
})


#PCAplot
#####################


##########OPLS-DA##########
nsamples_metadata <- scalednormalizedMetadata_samples[scalednormalizedMetadata_samples$Type == 'Sample',]

#report summary Q2, R2X, R2Y, etc. for all pairwise comps
oplsda_cum_scores <- NULL
oplsda_cum_score <- NULL #1st time so empty and no bug

## OPLS-DA for each pairwise comparison
if(AMOUNT_OF_COMPARISONS >= 1){
  for(pairwise_comparison in 1:AMOUNT_OF_COMPARISONS){
    
    #pairwise_comparison <- 1
    print(paste("start calculation comparison ", pairwise_comparison, " using OPLS", sep=""))
    
    ## OPLS-DA for each condition
    suppressMessages(library(ropls))
    #https://www.bioconductor.org/packages/devel/bioc/vignettes/ropls/inst/doc/ropls-vignette.html
    #https://rdrr.io/bioc/ropls/man/opls.html
    # pca: pred=NA, ortho=NA, no grouping; 
    # plsda: pred=NA, ortho=NA, groups(>2) converted to{0,1} "PLS2 implementation"; 
    # oplsda: pred=1, ortho=NA group{0,1} => but if set ortho==1 makes moden and otherwise ERROR so put 1+1 input for first round of analysis
    
    ### remove non-essential data for comparison
    #remove samples that are not included in comp (Na instead of 0/1)
    comp_ <- nsamples_metadata[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1]
    stopifnot(class(comp_)=='integer') #need only 2 possible {-1,1} arguments, rest is NA and not ""!
    samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    
    
    #remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
    samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
    #samples_matrix_comp_no0 <- samples_matrix_comp[, !apply(samples_matrix_comp == 0, 2, all)]  #remove 0
    samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x)))
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))]
    
    comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1])
    
    samplenames <- samples_metadata_comp$SampleName
    samplenames <- substr(samplenames, 2, nchar(samplenames)) 
    
    #optional change {0,1} to {-1,1} and run with comp as number, less strikt selection so more metabolites then with factor
    #minimal restriction with {0,1} but less correct (model not centered), even more metabolitse obtained
    #comp <- samples_metadata_comp[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1] 
    
    ### OPLS-DA plot
    try({
      max_crossval <- length(comp)
      if(max_crossval <= 7){
        max_crossval <- max_crossval
      }
      if(max_crossval > 7){
        max_crossval <- 7
      }  
      
      
      name_plot <- paste(name_project, "_OPLSDA_comparison_" , pairwise_comparison, '.pdf', sep="") #only works in development mode, in autorun is saved as Rplot.pdf but often more than 1pdf is made
      #pdf(name_plot, height=5,width=7)
      opls_comp <- opls(samples_matrix_comp_no0, comp, 
                        predI=1, orthoI=ORT, algoC="nipals", crossvalI=max_crossval,
                        permI=100,log10L=FALSE,scaleC="none", fig.pdfC=name_plot ,plotL=TRUE,printL=TRUE) #plotL=T prints plot to rplot.pdf #printL=T shows scores in console
      plot(opls_comp)
      #dev.off()
      # remark: in terminal saves under rplots.pdf, in rbox rstudio: no plots save but see in plots window. see plot_bis for results
      
      #report summary Q2, R2X, R2Y, etc.
      oplsda_cum_score <- as.data.frame(opls_comp@summaryDF)
      rownames(oplsda_cum_score) <- paste0("pairwise comparison ", pairwise_comparison)
      #name_df <- paste(name_project,'_OPLSDA_scores_comparison_', pairwise_comparison, '.txt', sep="")
      #write_dataframe_as_txt_file(oplsda_cum_score, name_df)
      
      #retain only vars that were not less than 2.2e-16 @opls-da 
      retained_variables_oplsda <- attributes(opls_comp@vipVn) 
      retained_variables_oplsda <- unlist(retained_variables_oplsda, use.names=FALSE)
      #removed_variables_oplsda <- attributes(opls_comp@xZeroVarVi)
      samples_matrix_comp_no0 <- samples_matrix_comp_no0[ , retained_variables_oplsda] #colnames retain if in list
      
      
      ### permutation test
      oplsda_scores <- data.frame(getScoreMN(opls_comp))
      permutation_comp <- data.frame(opls_comp@suppLs$permMN) #Info on permutations
      
      #t-test based on permutations: compare if model RMSEE differs significantly from 100 permuted RMSEEs
      normality_RMSEE <- shapiro.test(permutation_comp[,4]) #Test normality of the permuted RMSEEs: p < 0.05 --> yes
      #name_report_permutation <- paste(name_project,'_Report_OPLS_perm_tests_pairw_comp.txt', sep="")
      if (normality_RMSEE[2] < 0.05){
        #name_df <- paste('Permutation test for comparison ', pairwise_comparison, ".", sep="")
        #append_result_to_report(name_df, name_report_permutation)
        ttest_permutation_comp <- t.test(permutation_comp[,4],alternative="greater", 
                                         mu=opls_comp@summaryDF$RMSEE,conf.level=0.95) #p<<<0.05 --> RMSEE used model differs significantly from the permuted ones
        # model w slechter met permutatie want random metab wisselen aan stalen, als zelfde zou zijn als echt model, zou model via toeval evengoed zijn
        print_ttest <- capture.output(print(ttest_permutation_comp))
        #append_result_to_report(print_ttest, name_report_permutation)
        
        oplsda_cum_score$permutation_p_value <- ttest_permutation_comp$p.value
      }
      
      if (normality_RMSEE[2] >= 0.05){
        #name_df <- paste('Permutation test for comparison ', pairwise_comparison, ".", sep="")
        #append_result_to_report(name_df, name_report_permutation)
        print_ttest <- "The permuted RMSEEs do not have a normal distribution, hence no t-test may be performed."
        #append_result_to_report(print_ttest, name_report_permutation)
        
        oplsda_cum_score$permutation_p_value <- print_ttest
      }  
      
      
      #oplsda score plot weight values
      PC1 <- opls_comp@scoreMN[ ,1] #PC1 is kept (pre=1)
      PC2 <- opls_comp@orthoScoreMN[ ,1] #PC2 is kept (ort =1)
      score_weights_oplsda_comp <- cbind(PC1, PC2)
      score_weights_oplsda_comp <- as.data.frame(score_weights_oplsda_comp)
      rownames(score_weights_oplsda_comp) <- samples_metadata_comp$SampleName
      name_df <- paste(name_project, '_OPLSDA_weights_pairwise_comparison_', pairwise_comparison, '.txt', sep="")
      write_matrix_as_txt_file(score_weights_oplsda_comp, name_df) 
      #plot(score_weights_oplsda_comp)
      
      #nice ggplot with 2*PCs (pre1 and ort1)
      N <- nrow(score_weights_oplsda_comp)
      HotellingellipseOPLSDA<-Hotellingellipse(N,score_weights_oplsda_comp)
      dev.off()
      OPLSDAplot <- plot_oplsda(score_weights_oplsda_comp, comp, opls_comp, HotellingellipseOPLSDA)
      #pdf(paste(name_project, "_OPLSDA pairwise comparison " , pairwise_comparison, "_bis", ".pdf", sep=""), height=5, width=7)
      #plot(OPLSDAplot)
      #dev.off()
      
      png(paste(name_project, "_OPLSDA pairwise comparison " , pairwise_comparison, "_bis", ".png", sep=""), width=7, height=5, units="in", res=150)
      plot(OPLSDAplot)
      dev.off()
      
      
      ### selection metabolites
      #make outputfile with mz + time info 
      #retain only vars that were not less than 2.2e-16 @opls-da => done above      
      #retain only detected variables in comparison              => done above
      variableMetadata_info <- variableMetadata #subset(variableMetadata, select = c(1:(COLLUMN_NR_START_SAMPLES-1)))
      row.names(variableMetadata_info) <- paste(variableMetadata_info[ ,1])
      
      #retain only compoundIDs if in comparison (eg after filtering some removed) 
      variablenames <- colnames(samples_matrix_comp_no0)
      variablenames <- substr(variablenames, 2,nchar(variablenames))
      output_comp <- variableMetadata_info[rownames(variableMetadata_info) %in% variablenames, ]
      output_comp <- output_comp[ ,1:(COLLUMN_NR_START_SAMPLES-1)]

      #keep only vm intensities in in samplenames of comparison + only only compoundIDs if in comparison (eg after filtering some removed) 
      samplenames_w_X <- paste0("X", samplenames)
      output_intensities <- variableMetadata_info[ ,samplenames_w_X] 
      output_intensities <- output_intensities[rownames(output_intensities) %in% variablenames, ]
      
      
      ### VIP plot
      # first remove variables with all 0 values in collumn
      vip_comp <- getVipVn(opls_comp, orthoL = FALSE) #(default is FALSE and the predictive VIP is returned)
      comp <- as.numeric(comp)
      pvip <- apply(samples_matrix_comp_no0, 2, function(samples_matrix_comp_no0) cor.test(samples_matrix_comp_no0, comp)[["p.value"]])
      vipdf <- as.data.frame(cbind(vip_comp,pvip))
      
      VIPplot_comp <- plot_vip(vipdf)      
      #name_plot <- paste(name_project, "_OPLSDA VIP-plot Comparison " , pairwise_comparison, ".pdf", sep="")
      #pdf(name_plot, height=5,width=7)
      #plot(VIPplot_comp)
      #dev.off()
      
      png(paste(name_project, "_OPLSDA VIP-plot Comparison " , pairwise_comparison, ".png", sep=""), width=7, height=5, units="in", res=150)
      plot(VIPplot_comp)
      dev.off()
      
      # add vip scores to output file
      Variableselection_comp <-as.data.frame(matrix(ncol=9,nrow=(length(samples_matrix_comp_no0))))
      names(Variableselection_comp)<-c("Covariance","Correlation", "Corr_cutoff",
                                       "VIP_values","VIP_1","CI_center","CI_low","CI_high","CI_No_0")
      Variableselection_comp$VIP_values<-opls_comp@vipVn
      Variableselection_comp$VIP_1<-as.factor((opls_comp@vipVn>1)*1) #1 if vip > 1
      
      Variableselection_write <- cbind(output_comp, Variableselection_comp, output_intensities)
      name_df <- paste(name_project, '_OPLSDA_Variables_before_selection_comparison_', pairwise_comparison, '.txt', sep="")
      write_dataframe_as_txt_file(Variableselection_write, name_df) 
      
      
      
      ### S-plot
      #Calculate S-plot according to Wiklund et al. (2008), same technique as used in SIMCA. 
      #Description is in Supporting Information,as the NIPALS algorithm is used.
      
      #Calculate covariance ('p' in SIMCA) and correlation ('p(corr)' in SIMCA)
      Splotframe_comp<-as.data.frame(matrix(ncol=2,nrow=(length(samples_matrix_comp_no0)-0)))
      names(Splotframe_comp)<-c("Covariance","Correlation")
      for(i in 1:(length(samples_matrix_comp_no0)-0)){
        Splotframe_comp[i,1]<-(t(opls_comp@scoreMN)%*%samples_matrix_comp_no0[,i])/(t(opls_comp@scoreMN)%*%opls_comp@scoreMN)   #Covariance calculation for each variable
        Splotframe_comp[i,2]<-(t(opls_comp@scoreMN)%*%samples_matrix_comp_no0[,i])/(sqrt(sum(as.numeric(opls_comp@scoreMN)^2))*sqrt(sum(as.numeric(samples_matrix_comp_no0[,i])^2))) #Correlation calculation. Norm used is Frobenius norm (the Euclidean norm of x treated as if it were a vector)
      }
      
      #choose/calculate correlation coefficent (Y-axis)
      if(CUTOFF_CORRELATION == FIXED_CORRELATION_COEFFICIENT){
        Cutoffvalue_corr <- FIXED_CORRELATION_COEFFICIENT
      }
      if(CUTOFF_CORRELATION == CRITICAL_CORRELATION_COEFFICIENT){
        #Find critical cut-off value for correlation: two-tailed test
        # Critical correlation coefficient at sample size N = number of samples
        N <- nrow(samples_matrix_comp_no0)
        Cutoffvalue_corr<-critical.r(N,alpha=0.05) #(function from Functions_R_Script_multivariate.R file)
      }
      
      Splotframe_comp$Corr_cutoff <-(Splotframe_comp[,2]>Cutoffvalue_corr)*1 | (Splotframe_comp[,2]<(-Cutoffvalue_corr)*1)
      
      Splot_comp <- plot_Splot(Splotframe_comp)
      #name_plot <- paste(name_project, "_OPLSDA S-plot Comparison " , pairwise_comparison, ".pdf", sep="")
      #pdf(name_plot, height=5,width=7)
      #plot(Splot_comp)
      #dev.off()
      
      png(paste(name_project, "_OPLSDA S-plot Comparison " , pairwise_comparison, ".png", sep=""), width=7, height=5, units="in", res=150)
      plot(Splot_comp)
      dev.off()
      
      # add S-plot scores to output file
      Variableselection_comp$Covariance<-Splotframe_comp$Covariance
      Variableselection_comp$Correlation<-Splotframe_comp$Correlation
      Variableselection_comp$Corr_cutoff<-as.factor(Splotframe_comp$Corr_cutoff) #TRUE if corr(x) < -/+cutoff < corr(x).
      
      Variableselection_write <- cbind(output_comp, Variableselection_comp, output_intensities)
      name_df <- paste(name_project, '_OPLSDA_Variables_before_selection_comparison_', pairwise_comparison, '.txt', sep="")
      write_dataframe_as_txt_file(Variableselection_write, name_df) 
      
      
      ###Loading plot 
      #for calculation CIs using bootstrap
      suppressMessages(library(boot))
      
      #Put same amount of orthogonal components found to be optimal in previous models
      optimal_predI <- as.numeric(as.character(opls_comp@summaryDF[5]))
      optimal_orthoI <- as.numeric(as.character(opls_comp@summaryDF[6]))
      comp <- as.numeric(comp)
      df_loading <- cbind(comp, samples_matrix_comp_no0)
      
      loadingboot <- function(data, indices, select_predI, select_orthoI, select_crossval){
        
        loadingduringboot1 <- opls(data[indices,2:length(data)], data$comp[indices],predI=select_predI, orthoI=select_orthoI, 
                                   algoC="nipals", crossvalI=select_crossval,permI=5,log10L=FALSE,scaleC="none",plotL=FALSE,printL=FALSE)
        
        variablenames <- colnames(data)[-1]             #all names before OPLSDA
        #variablenames <- substr(variablenames, 2,nchar(variablenames)) #only if "X1" => when x instead of "1", check @samplesnoO!
        colnames(data)[-1] <- variablenames
        df <- matrix(nrow=length(data[-1]), ncol=1) #make matrix of all names with NA
        #df[,] <- 0                                  #replace all NA to 0
        df <- as.data.frame(df)
        df$names <- colnames(data[-1])
        
        retained_variables_oplsda <- attributes(loadingduringboot1@vipVn)       #names kept after oplsda
        retained_variables_oplsda <- unlist(retained_variables_oplsda, use.names=FALSE)
        variablenames2 <- retained_variables_oplsda
        #variablenames2 <- substr(variablenames2, 2,nchar(variablenames2)) #only if "X1" =>
        retained_loadings <- loadingduringboot1@loadingMN                       #loading values of names kept after oplsda
        retained_loadings <- data.frame(retained_loadings)
        retained_loadings$names <- variablenames2
        
        all_loadings <- merge(df, retained_loadings, by.x="names", all.x = TRUE, sort=FALSE) #merge full (empty) df before oplsda + retained with values
        all_loadings <- all_loadings[ ,-2]                                      #remove empty col from all_variables needed for creation df 
        all_loadings <- all_loadings[match(df$names, all_loadings$names),]      # correct order names
        all_loadings <- all_loadings[ ,-1]                                      #delete col names
        all_loadings <- as.vector(all_loadings)
        
        return(all_loadings)
      }
      
      #resampling
      bootresults <- boot(data = df_loading, statistic = loadingboot, R=200,    
                          select_predI=optimal_predI, select_orthoI=optimal_orthoI, select_crossval=max_crossval)
      #min 200x resampling, else error too little bootstrap for CI calculation
      #write_dataframe_as_txt_file(bootresults$t, 'bootresultsT.txt')
      #write_dataframe_as_txt_file(bootresults$t0, 'bootresultsT0.txt')
      
      #Put CIs in dataframe
      myCIsframe_comp <-as.data.frame(matrix(ncol=3,nrow=(length(samples_matrix_comp_no0)-0))) #-0 is info colums vb 'comp'
      names(myCIsframe_comp)<-c("Center_CI","Lower_limit","Upper_limit")
      
      for(i in 1:length(samples_matrix_comp_no0)){
        try({
          #i <- 1
          CI<-boot.ci(bootresults,conf=0.95,type="bca",index=i)  # 'adjusted bootstrap percentile' (BCA) method
          myCIsframe_comp[i,1]<-CI$t0 #Center of CI
          myCIsframe_comp[i,2]<-CI$bca[4] #Lower limit of CI
          myCIsframe_comp[i,3]<-CI$bca[5] #Upper limit of CI
        })
      }
      #write_dataframe_as_txt_file(myCIsframe_comp, "CIresults.txt") 
      
      #Loading plot OPLS-DA with 95% CI
      OPLSDAloadingplot_comp <- plot_loading(myCIsframe_comp) 
      #name_plot <- paste(name_project, "_OPLSDA loadingplot Comparison " , pairwise_comparison, ".pdf", sep="")
      #pdf(name_plot, height=5,width=7)
      #plot(OPLSDAloadingplot_comp)
      #dev.off()
      
      png(paste(name_project, "_OPLSDA loadingplot Comparison " , pairwise_comparison, ".png", sep=""), width=7, height=5, units="in", res=150)  
      plot(OPLSDAloadingplot_comp)
      dev.off()
      
      # add loading scores to output file
      Variableselection_comp$CI_center<-myCIsframe_comp$Center_CI
      Variableselection_comp$CI_low<-myCIsframe_comp$Lower_limit
      Variableselection_comp$CI_high<-myCIsframe_comp$Upper_limit
      #Variableselection_comp$CI_0 <-(myCIsframe_comp[,2] <= 0 & myCIsframe_comp[,3] >= 0) #returns true if CI in 0
      Variableselection_comp$CI_No_0 <-((myCIsframe_comp[,2] > 0 & myCIsframe_comp[,3] > 0) | (myCIsframe_comp[,2] < 0 & myCIsframe_comp[,3] < 0)) #returns true if CI NOT in 0, wanted!
      Variableselection_comp$CI_No_0[is.na(Variableselection_comp$CI_No_0)] <- "NO_CI" #all values no CI was calculated NA -> 0, become boolean elsewise
      
      Variableselection_write <- cbind(output_comp, Variableselection_comp, output_intensities)
      name_df <- paste(name_project, '_OPLSDA_Variables_before_selection_comparison_', pairwise_comparison, '.txt', sep="")
      write_dataframe_as_txt_file(Variableselection_write, name_df) 
      
      # selection based on vip, CI and corr(x)
      Variableselection_comp <- cbind(output_comp, Variableselection_comp, output_intensities) #one time merge for selection. above 3times only for write away to file and replace if extra info
      Variables_VIP_1_comp<-Variableselection_comp[Variableselection_comp$VIP_1==1,]
      Variables_corr_cutoff_comp<-Variableselection_comp[Variableselection_comp$Corr_cutoff==TRUE,] 
      Variables_CIno0_comp<-Variableselection_comp[Variableselection_comp$CI_No_0==TRUE,] #Select only variables with CI that do not contain 0
      Variables_VIP_1_CIno0_corr_cutoff_comp<-Variableselection_comp[Variableselection_comp$CI_No_0==TRUE & Variableselection_comp$Corr_cutoff==TRUE & Variableselection_comp$VIP_1==1,]
      
      name_df <- paste(name_project, '_OPLSDA_Variables_after_selection_VIP1_CIno0_Corr_cutoff_comparison_', pairwise_comparison, '.txt', sep="")
      try(write_dataframe_as_txt_file(Variables_VIP_1_CIno0_corr_cutoff_comp, name_df)) #error if empty df: nrow(dataframe) >= 1 is not TRUE
      
      total_amount_of_variables <- nrow(Variableselection_comp)
      amount_variables_after_selection <- nrow(Variables_VIP_1_CIno0_corr_cutoff_comp)
      report_variables_after_selection <- paste("There are " , amount_variables_after_selection, " significant variables retrieved from the total of " , total_amount_of_variables , " variables obtained in the OPLSDA model of comparison ", pairwise_comparison, ".", sep="")
      name_report_variables <- paste(name_project,'_OPLSDA_nr_variables_after_selection_pairwise_comparisons.txt', sep="")
      append_result_to_report(report_variables_after_selection, name_report_variables) 
      
      
      
      ### Heatmap
      suppressMessages(library(mixOmics)) #not in R3.6.1, R3.5.3
      suppressMessages(library(RColorBrewer))
      colors<-brewer.pal(11,name="RdYlBu")
      pal<-colorRampPalette(colors)
      
      #Heat map of the data with VIP > 1, no 0 in the CI and corr > cutoff
      Variables_VIP_1_CIno0_corr_cutoff_comp$CompID <- as.character(Variables_VIP_1_CIno0_corr_cutoff_comp$CompID)
      Tsamples_matrix_comp_no0 <- t(samples_matrix_comp_no0)
      rownames(Tsamples_matrix_comp_no0) <- substring(rownames(Tsamples_matrix_comp_no0) , 2)  #remove "X" before compid to ifnd in selection df
      Tsamples_matrix_comp_no0_selection <- Tsamples_matrix_comp_no0[rownames(Tsamples_matrix_comp_no0) %in% Variables_VIP_1_CIno0_corr_cutoff_comp$CompID, ]
      
      colnames(Tsamples_matrix_comp_no0_selection) <- samples_metadata_comp$SampleName
      name_plot <- paste(name_project, "_OPLSDA Heatmap comparison " , pairwise_comparison, sep="")
      
      try({  
        colours <-pal(100) #default scheme
        #colours <-brewer.pal(n=9, name="PuBuGn") #colorblind scheme
        
        #error if empty df after selection: must have n >= 2 objects to cluster
        heatmap_comp <- cim(data.matrix(Tsamples_matrix_comp_no0_selection), color=colours, save='pdf', name.save = name_plot, 
                            row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
        
        cim(data.matrix(Tsamples_matrix_comp_no0_selection), color=colours, save='png', name.save = name_plot,
            row.names=TRUE,col.names=TRUE, cluster = "both",center=F,scale=F, margins = c(10,10))
      })
      detach(package:mixOmics)
      
    })
    
    #append oplsda score to 1df
    if(length(oplsda_cum_score) > 0){ #not from prev opls paste, only add found
      oplsda_cum_scores <- rbind(oplsda_cum_scores, oplsda_cum_score)
    }
    oplsda_cum_score <- NULL
    
    print(paste("finished calculation comparison ", pairwise_comparison, " using OPLS", sep=""))
    
  }
  
  #report summary Q2, R2X, R2Y, etc. for all pairwise comps
  oplsda_cum_scores$Comparison <- rownames(oplsda_cum_scores)
  oplsda_cum_scores <- oplsda_cum_scores[ , c(10, 1, 2, 3, 4, 5, 6, 7, 8, 9)] #reorder cols     
  name_df <- paste(name_project,'_OPLSDA_scores_pairwise_comparisons.txt', sep="")
  write_dataframe_as_txt_file(oplsda_cum_scores, name_df) 
}
###################

##########RF##########
nsamples_metadata <- scalednormalizedMetadata_samples[scalednormalizedMetadata_samples$Type == 'Sample',]

## OPLS-DA for each pairwise comparison
if(AMOUNT_OF_COMPARISONS >= 1){
  for(pairwise_comparison in 1:AMOUNT_OF_COMPARISONS){
    
    #pairwise_comparison <- 1
    print(paste("start calculation comparison ", pairwise_comparison, " using RF", sep=""))
    
    ## RF for each comparison
    suppressMessages(library(randomForest))
    #https://www.r-bloggers.com/2018/01/how-to-implement-random-forests-in-r/
    
    ### remove non-essential data for comparison
    #remove samples that are not included in comp (Na instead of 0/1)
    comp_ <- nsamples_metadata[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1]
    stopifnot(class(comp_)=='integer') #need only 2 possible {-1,1} arguments, rest is NA and not ""!
    samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    
    
    #remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
    samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
    #samples_matrix_comp_no0 <- samples_matrix_comp[, !apply(samples_matrix_comp == 0, 2, all)]  #remove 0
    samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x)))
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))]
    
    comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1])
    
    
    #column with labels cbind to matrix for RF
    samples_matrix_comp_no0 <- cbind(samples_matrix_comp_no0, comp)
    
    
    # Split into Train and Validation sets
    # Training Set : Validation Set = 70 : 30 (random)
    set.seed(100)
    train <- sample(nrow(samples_matrix_comp_no0), 0.7*nrow(samples_matrix_comp_no0), replace = FALSE)
    TrainSet <- samples_matrix_comp_no0[train,]
    ValidSet <- samples_matrix_comp_no0[-train,]
    #summary(TrainSet)
    #summary(ValidSet)
    
    
    #if fails:
    try({
      # Create a Random Forest model with parameters == 500trees + nr variabl tried at each split 10
      model1 <- randomForest(comp ~ ., data = TrainSet, importance = TRUE, ntree = 500, mtry = 10)
      #model1
      
      # Using For loop to identify the right mtry for model => manual chack
      #a=c()
      #i=5
      #for (i in 3:8) {
      #  model3 <- randomForest(comp ~ ., data = TrainSet, ntree = 500, mtry = i, importance = TRUE)
      #  predValid <- predict(model3, ValidSet, type = "class")
      #  a[i-2] = mean(predValid == ValidSet$comp)
      #}
      #a
      #plot(3:8,a)
      
      #write summary model with param and score arror rate
      name_report_rf <- paste(name_project,'_RF_scores_pairw_comp.txt', sep="")
      title <- paste('RF_comparison ', pairwise_comparison, ".", sep="")
      append_result_to_report(title, name_report_rf)
      print_rf <- capture.output(print(model1))
      append_result_to_report(print_rf, name_report_rf)
      print_space <- ""
      append_result_to_report(print_space, name_report_rf)
      
      
      ## Predicting on Validation set
      predValid <- predict(model1, ValidSet, type = "class")
      
      # Checking classification accuracy
      append_result_to_report("Predicting on Validation set", name_report_rf)
      accuracy <- mean(predValid == ValidSet$comp)     
      append_result_to_report(paste0("accuracy model: ", accuracy), name_report_rf)
      
      #write confusion matrix
      conf_matrix <- table(predValid,ValidSet$comp)
      append_result_to_report(capture.output(print(conf_matrix)), name_report_rf)
      print_space <- ""
      append_result_to_report(print_space, name_report_rf)
      
      #plot confus matrix nice (long format)
      library(reshape2)
      confusion_matrix_long <- melt(conf_matrix)
      colnames(confusion_matrix_long)[2] <- "ActualValid"
      confusion_matrix_long$predValid <- as.factor(confusion_matrix_long$predValid) #as factor for plot
      confusion_matrix_long$ActualValid <- as.factor(confusion_matrix_long$ActualValid)
      
      cm_plot <- plot_ConfusionMatrix_comp(confusion_matrix_long)
      
      png(paste(name_project, "_RF_confusion_matrix_comparison_" , pairwise_comparison, ".png", sep=""), width=7, height=5, units="in", res=150)
      plot(cm_plot)
      dev.off()
      
      
      ## importance_features
      importance_features <- as.data.frame(importance(model1))
      importance_features$CompID <- as.character(rownames(importance_features))
      name_df <- paste(name_project, '_RF_importance_features_comparison_', pairwise_comparison, '.txt', sep="")
      write_dataframe_as_txt_file(importance_features, name_df) 
      
      #plot show the drop in mean accuracy for each of the variables.
      png(paste(name_project, "_RF_importance_features_comparison_" , pairwise_comparison, ".png", sep=""), width=7, height=5, units="in", res=150)
      varImpPlot(model1)
      dev.off()
      
      print(paste("finished calculation comparison ", pairwise_comparison, " using RF", sep=""))
    })
  }
}
###################


##########LIMMA##########
nsamples_metadata <- scalednormalizedMetadata_samples[scalednormalizedMetadata_samples$Type == 'Sample',]

## OPLS-DA for each pairwise comparison
if(AMOUNT_OF_COMPARISONS >= 1){
  for(pairwise_comparison in 1:AMOUNT_OF_COMPARISONS){
    
    #pairwise_comparison <- 1
    print(paste("start calculation comparison ", pairwise_comparison, " using LIMMA", sep=""))
    
    ### remove non-essential data for comparison
    #remove samples that are not included in comp (Na instead of 0/1)
    comp_ <- nsamples_metadata[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1]
    stopifnot(class(comp_)=='integer') #need only 2 possible {-1,1} arguments, rest is NA and not ""!
    samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    
    
    #remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
    samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
    #samples_matrix_comp_no0 <- samples_matrix_comp[, !apply(samples_matrix_comp == 0, 2, all)]  #remove 0
    samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x)))
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))]
    
    comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1])

    suppressMessages(library(limma))
    
    ### input file rename
    Tsamples_matrix_comp_no0 <- t(samples_matrix_comp_no0)
    expressionset <- Tsamples_matrix_comp_no0
    
    #optional: add block eg batch effect contructed in blocks
    #block <- factor(samples_metadata_comp$Kooi)
    
    ### design matrix
    f <- factor(comp) 
    design_comp_matrix <- model.matrix(~0+f) #make dummy matrix with only 0/1
    #in case paired etc, add block eg batch effect contructed in blocks: design_comp_matrix <- model.matrix(~block+f)  
    colnames(design_comp_matrix) <- c("V1", "V2")    #name of each group in comparison -1, 1
    #write_dataframe_as_txt_file(design_comp_matrix, 'design_comp_matrix.txt') 
    
    ### fit linear model for each group for a series of variables
    fit <- lmFit(expressionset, design=design_comp_matrix)
    #write_dataframe_as_txt_file(fit, 'fit.txt') 
    
    ###add contract (1 combi since pairwise) to find differentially present features
    contrasts <- "V1-V2"
    contrastMatrix_pairwise <- makeContrasts(contrasts=contrasts, levels=colnames(design_comp_matrix))
    contrast.matrix <- contrastMatrix_pairwise
    #write_dataframe_as_txt_file(contrastMatrix_pairwise, 'contrastMatrix_pairwise.txt') 
    
    
    ### fit model with contrasts
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2) #empirical Bayes moderated t-statistics and associated Benjamini Hochberg adjusted p-values 
    #write_dataframe_as_txt_file(fit2, 'fit2.txt') 
    
    
    ### selection metabolites
    #CompID names not as "X1"
    variablenames <- colnames(samples_matrix_comp_no0)
    variablenames <- substr(variablenames, 2,nchar(variablenames))
    colnames(samples_matrix_comp_no0) <- variablenames
    Tsamples_matrix_comp_no0 <- t(samples_matrix_comp_no0) #again for beneath correct name compID
    
    #make outputfile with mz + time info
    variableMetadata_info <- variableMetadata #subset(variableMetadata, select = c(1:(COLLUMN_NR_START_SAMPLES-1)))
    row.names(variableMetadata_info) <- paste(variableMetadata_info[ ,1])
    
    output_comp <- merge(variableMetadata_info, Tsamples_matrix_comp_no0, by="row.names", all.y = TRUE, sort=FALSE) #retain only if in comp
    output_comp <- output_comp[ , 1:ncol(variableMetadata_info)]
    #write_dataframe_as_txt_file(output_comp, 'empty_output_comp.txt') 
    
    #Export top variables for each combination whitin comparison
    #topTable(fit2, coef=1, adjust="BH") #coef1 = 1e vgl.
    output_coef <- topTable(fit2, coef=1, adjust="BH", number=nrow(fit2)) 
    output_coef$CompID <- rownames(output_comp)
    output_coef_with_info <- merge(output_coef, output_comp, by="CompID", sort=FALSE)
    name_df <- paste(name_project, '_LIMMA_Toptable_output_comp_', pairwise_comparison, '_coef_', colnames(fit2), '.txt', sep="")
    write_dataframe_as_txt_file(output_coef_with_info, name_df)

    # Export 1 result {-1,0,1} matrix for pairwise comparison
    results <- decideTests(fit2)
    results <- as.data.frame(results)
    results$CompID <- rownames(output_comp)
    results_with_info <- merge(results, output_comp, by="CompID", sort=FALSE)
    name_df <- paste(name_project, '_LIMMA_ResultsTable_output_comparison', pairwise_comparison, '.txt', sep="")
    write_dataframe_as_txt_file(results_with_info, name_df)
    
    
    print(paste("finished calculation comparison ", pairwise_comparison, " using LIMMA", sep=""))
  }
}
###################    
    

##########PLS-DA##########
#report summary Q2, R2X, R2Y, etc. for all multi comps
plsda_cum_scores <- NULL
plsda_cum_score <- NULL

## PLS-DA for each multiple comparison
if(AMOUNT_OF_MULTIPLE_COMPARISONS >= 1){
  for(multiple_comparison in 1:AMOUNT_OF_MULTIPLE_COMPARISONS){
    
    print(paste("start calculation multiple comparison ", multiple_comparison, " using PLS", sep=""))
    
    ## PLS-DA for each condition
    suppressMessages(library(ropls))
    #https://www.bioconductor.org/packages/devel/bioc/vignettes/ropls/inst/doc/ropls-vignette.html
    #https://rdrr.io/bioc/ropls/man/opls.html
    # pca: pred=NA, ortho=NA, no grouping; 
    # plsda: pred=NA, ortho=NA, groups(>2) depicted with numbers! converted to{0,1} "PLS2 implementation"; 
    # oplsda: pred=1, ortho=NA group{0,1} => but if set ortho==1 makes moden and otherwise ERROR so put 1+1 input for first round of analysis
    
    ### remove non-essential data for comparison
    #remove samples that are not included in comp (Na instead of 0/1)
    nsamples_metadata <- scalednormalizedMetadata_samples[scalednormalizedMetadata_samples$Type == 'Sample',] 
    comp_ <- nsamples_metadata[,COLLUMN_NR_MULTIPLE_COMPARISON1+multiple_comparison-1]
    stopifnot(class(comp_)=='integer') #need {0,1,2, ...} arguments, rest is NA and not ""!
    samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    
    
    #remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
    samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
    #samples_matrix_comp_no0 <- samples_matrix_comp[, !apply(samples_matrix_comp == 0, 2, all)]  #remove 0
    samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x)))
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))]
    
    comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_MULTIPLE_COMPARISON1+multiple_comparison-1])
    
    
    ### PLS-DA plot
    try({
      max_crossval <- length(comp)
      if(max_crossval <= 7){
        max_crossval <- max_crossval
      }
      if(max_crossval > 7){
        max_crossval <- 7
      }  
      
      #pls-da: ortho =0 (want geen ortholoog model) en pred= NA (niet op 1 want min. 2 nodig voor plot PC1+PC2)!
      name_plot <- paste(name_project, "_PLSDA_multiple_comparison_" , multiple_comparison, '.pdf', sep="") #only works in development mode, in autorun is saved as Rplot.pdf and more>1 is made
      #pdf(name_plot, height=5,width=7) 
      pls_comp <- opls(samples_matrix_comp_no0, comp,
                       predI=PRE, orthoI=0, algoC="nipals", crossvalI=max_crossval,permI=100, 
                       log10L=FALSE,scaleC="none", fig.pdfC=name_plot, plotL=TRUE,printL=TRUE)
      plot(pls_comp)
      #dev.off()
      # remark: in terminal saves under rplots.pdf, in rbox rstudio: no plots save but see in plots window. see plot_bis for results
      
      #report summary Q2, R2X, R2Y, etc.
      plsda_cum_score <- as.data.frame(pls_comp@summaryDF)
      rownames(plsda_cum_score) <- paste0("multiple comparison ", multiple_comparison)
      #name_df <- paste(name_project, '_PLSDA_scores_multiple_comparison_', multiple_comparison, '.txt', sep="")
      #write_dataframe_as_txt_file(plsda_cum_score, name_df) 
      
      
      #retain only vars that were not less than 2.2e-16 @pls-da 
      retained_variables_plsda <- attributes(pls_comp@vipVn) 
      retained_variables_plsda <- unlist(retained_variables_plsda, use.names=FALSE)
      #removed_variables_plsda <- attributes(pls_comp@xZeroVarVi)
      samples_matrix_comp_no0 <- samples_matrix_comp_no0[ , retained_variables_plsda] #colnames retain if in list
      
      
      ### permutation test
      plsda_scores <- data.frame(getScoreMN(pls_comp))
      permutation_comp <- data.frame(pls_comp@suppLs$permMN) #Info on permutations
      
      #t-test based on permutations: compare if model RMSEE differs significantly from 100 permuted RMSEEs
      normality_RMSEE <- shapiro.test(permutation_comp[,4]) #Test normality of the permuted RMSEEs: p < 0.05 --> yes
      #name_df <- paste('Permutation test for multiple comparison ', multiple_comparison, ".", sep="")
      #name_report_permutation2 <- paste(name_project,'_Report_PLS_permutation_tests_multiple_comparisons.txt', sep="")
      if (normality_RMSEE[2] < 0.05){
        #append_result_to_report(name_df, name_report_permutation2)
        ttest_permutation_comp <- t.test(permutation_comp[,4],alternative="greater", 
                                         mu=pls_comp@summaryDF$RMSEE,conf.level=0.95) #p<<<0.05 --> RMSEE used model differs significantly from the permuted ones
        print_ttest <- capture.output(print(ttest_permutation_comp))
        #append_result_to_report(print_ttest, name_report_permutation2)
        
        plsda_cum_score$permutation_p_value <- ttest_permutation_comp$p.value
      }
      
      if (normality_RMSEE[2] >= 0.05){
        #name_df <- paste('Permutation test for multiple comparison ', multiple_comparison, ".", sep="")
        #append_result_to_report(name_df, name_report_permutation2)
        print_ttest <- "The permuted RMSEEs do not have a normal distribution, hence no t-test may be performed."
        #append_result_to_report(print_ttest, name_report_permutation2)
        
        plsda_cum_score$permutation_p_value <-print_ttest
      }  
      
      #nice ggplot: pls plot_bis 
      if (ncol(plsda_scores) > 1){
        N <- nrow(samples_matrix_comp_no0)
        HotellingellipsePLSDA<-Hotellingellipse(N,plsda_scores)
        dev.off()
        PLSDAplot <- plot_plsda(plsda_scores, comp, pls_comp, HotellingellipsePLSDA)
        
        #pdf(paste(name_project, "_PLSDA multiple comparison " , multiple_comparison, "_bis", ".pdf", sep=""), height=5, width=7)
        #plot(PLSDAplot)
        #dev.off()
        
        png(paste(name_project, "_PLSDA multiple comparison " , multiple_comparison, "_bis", ".png", sep=""), width=7, height=5, units="in", res=150)
        plot(PLSDAplot)
        dev.off()
        
        #plsda score plot weight values
        plsda_score_comp <- pls_comp@scoreMN
        rownames(plsda_score_comp) <- samples_metadata_comp$SampleName
        score_weights_plsda_comp <- plsda_score_comp[ ,1:2] #PC1 and PC2 are kept
        name_df <- paste(name_project, '_PLSDA_weights_multiple_comparison_', multiple_comparison, '.txt', sep="")
        write_matrix_as_txt_file(score_weights_plsda_comp, name_df) 
        #plot(score_weights_plsda_comp)
      }
      if (ncol(plsda_scores) <= 1){
        print("No plot plsda_bis is possibe, only 1 PC was found.")
        remark_plsda_plot <- "A pls model is made but no plsda_plot_bis could be plotted because only 1 PC was retained after dimension reduction."
        name_df <- paste(name_project, '_PLSDA_remark_plot_multiple_comparison_', multiple_comparison, '.txt', sep="")
        write_dataframe_as_txt_file(remark_plsda_plot, name_df) 
        
        #plsda score plot weight values
        plsda_score_comp <- plsda_scores #== pls_comp@scoreMN is >1 dim this object exists
        rownames(plsda_score_comp) <- samples_metadata_comp$SampleName
        score_weights_plsda_comp <- plsda_score_comp[ ,1] #PC1 is kept
        name_df <- paste(name_project, '_PLSDA_weights_multiple_comparison_', multiple_comparison, '.txt', sep="")
        write_matrix_as_txt_file(score_weights_plsda_comp, name_df) 
        #plot(score_weights_plsda_comp)
        
        #pdf(paste(name_project, "_PLSDA multiple comparison " , multiple_comparison, "_bis", ".pdf", sep=""), height=5, width=7)
        #plot(score_weights_plsda_comp) #todo new function with nic ggplot in future if needed
        #dev.off()
        
        png(paste(name_project, "_PLSDA multiple comparison " , multiple_comparison, "_bis", ".png", sep=""), width=7, height=5, units="in", res=150)
        plot(score_weights_plsda_comp)
        dev.off()
      }
      
    })
    
    #append oplsda score to 1df
    if(length(plsda_cum_score) > 0){ #not from prev opls paste, only add found
      plsda_cum_scores <- rbind(plsda_cum_scores, plsda_cum_score)
    }
    plsda_cum_score <- NULL
    print(paste("finished calculation multiple comparison ", multiple_comparison, " using PLS", sep=""))
  }
  
  #report summary Q2, R2X, R2Y, etc. for all pairwise comps
  plsda_cum_scores$Comparison <- rownames(plsda_cum_scores)
  plsda_cum_scores <- plsda_cum_scores[ , c(10, 1, 2, 3, 4, 5, 6, 7, 8, 9)] #reorder cols     
  name_df <- paste(name_project,'_PLSDA_scores_multiple_comparisons.txt', sep="")
  write_dataframe_as_txt_file(plsda_cum_scores, name_df) 
}    
###################


##########Limma##########
## multiple comp with Limma = Linear Models for microarray data
#limmaUsersGuide(view=TRUE)
if(AMOUNT_OF_MULTIPLE_COMPARISONS >= 1){
  for(multiple_comparison in 1:AMOUNT_OF_MULTIPLE_COMPARISONS){
    
    print(paste("start calculation multiple comparison ", multiple_comparison, " using LIMMA", sep=""))
    
    ### remove non-essential data for comparison
    #remove samples that are not included in comp (Na instead of 0/1)
    nsamples_metadata <- scalednormalizedMetadata_samples[scalednormalizedMetadata_samples$Type == 'Sample',] 
    comp_ <- nsamples_metadata[,COLLUMN_NR_MULTIPLE_COMPARISON1+multiple_comparison-1]
    stopifnot(class(comp_)=='integer') #need {0,1,2, ...} arguments, rest is NA and not ""!
    samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    
    
    #remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
    samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
    #samples_matrix_comp_no0 <- samples_matrix_comp[, !apply(samples_matrix_comp == 0, 2, all)]  #remove 0
    samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x)))
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))]
    
    comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_MULTIPLE_COMPARISON1+multiple_comparison-1])
    
    suppressMessages(library(limma))
    ### input file rename
    Tsamples_matrix_comp_no0 <- t(samples_matrix_comp_no0)
    expressionset <- Tsamples_matrix_comp_no0
    
    #optional: add block eg batch effect contructed in blocks
    #block <- factor(samples_metadata_comp$Kooi)
    
    ### design matrix
    f <- factor(comp) 
    design_comp_matrix <- model.matrix(~0+f) #make dummy matrix with only 0/1
    #in case paired etc, add block eg batch effect contructed in blocks: design_comp_matrix <- model.matrix(~block+f)  
    colnames(design_comp_matrix) <- paste("V", levels(f), sep="")    #name of each group in multiple comparison
    #write_dataframe_as_txt_file(design_comp_matrix, 'design_comp_matrix.txt') 
    
    ### fit linear model for each group for a series of variables
    fit <- lmFit(expressionset, design=design_comp_matrix)
    #write_dataframe_as_txt_file(fit, 'fit.txt') 
    
    ### add contrast
    nr_of_groups <- ncol(design_comp_matrix)
    
    
    #make input for makecontrast() with pairwise comparisons
    contrasts1 <-c()
    
    for(i in 1:(length(as.character(levels(f))))){
      #i <- 9
      group <- levels(f)[i]
      group <- paste("V", group, sep="")
      
      for(j in 1:(length(as.character(levels(f))))-1){
        group2 <- levels(f)[j+1]
        group2 <- paste("V", group2, sep="")
        
        if(group2 > group){
          contrast1 <- paste(group, "-", group2, sep="")
          contrasts1 <- paste(contrasts1, contrast1, sep=", ") 
        }
      }
    }
    
    contrasts <- substr(contrasts1, 3, (nchar(contrasts1)))
    contrasts <- unlist(strsplit(contrasts, ", ", fixed=FALSE))
    #contrasts <- c(contrasts, "V1+V2-V3-V4", "V1-V2+V3-V4", "V1+V2+V3+V4","V1-V2-V3+V4", "(V1-V2)-(V3-V4)")
    #with "(V1-V2)-(V3-V4)" to check difference, extract the comparisons of interest as contrasts vb Dis treatm/untreatm and WT treatm/untreatm 
    #contrasts <- c(contrasts, "(V1-V2)-(V3-V4-V5-V6-V7-V8)") #kooi effect
    contrasts <- dput(as.character(contrasts))
    
    contrastMatrix_pairwise <- makeContrasts(contrasts=contrasts, levels=colnames(design_comp_matrix))
    contrast.matrix <- contrastMatrix_pairwise
    #write_dataframe_as_txt_file(contrastMatrix_pairwise, 'contrastMatrix_pairwise.txt') 
    
    
    ### fit model with contrasts
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2) #empirical Bayes moderated t-statistics and associated Benjamini Hochberg adjusted p-values 
    #write_dataframe_as_txt_file(fit2, 'fit2.txt') 
    
    
    ### selection metabolites
    #CompID names not as "X1"
    variablenames <- colnames(samples_matrix_comp_no0)
    variablenames <- substr(variablenames, 2,nchar(variablenames))
    colnames(samples_matrix_comp_no0) <- variablenames
    Tsamples_matrix_comp_no0 <- t(samples_matrix_comp_no0) #again for beneath correct name compID
    
    #make outputfile with mz + time info
    variableMetadata_info <- variableMetadata #subset(variableMetadata, select = c(1:(COLLUMN_NR_START_SAMPLES-1)))
    row.names(variableMetadata_info) <- paste(variableMetadata_info[ ,1])
    
    output_comp <- merge(variableMetadata_info, Tsamples_matrix_comp_no0, by="row.names", all.y = TRUE, sort=FALSE) #retain only if in comp
    output_comp <- output_comp[ , 1:ncol(variableMetadata_info)]
    #write_dataframe_as_txt_file(output_comp, 'empty_output_comp.txt') 
    
    #Export top variables for each combination whitin comparison
    #topTable(fit2, coef=1, adjust="BH") #coef1 = 1e vgl.
    results <- decideTests(fit2)
    for(coeficent in 1:ncol(results)){
      #coeficent <- 1
      output_coef <- topTable(fit2, coef=coeficent, adjust="BH", number=nrow(fit2)) 
      output_coef$CompID <- rownames(output_comp)
      output_coef_with_info <- merge(output_coef, output_comp, by="CompID", sort=FALSE)
      name_df <- paste(name_project, '_LIMMA_Toptable_output_multiple_comp_', multiple_comparison, '_coef_', colnames(fit2)[coeficent], '.txt', sep="")
      write_dataframe_as_txt_file(output_coef_with_info, name_df)
    }
    
    # Export 1 result {-1,0,1} matrix for multiple comparison
    results <- decideTests(fit2)
    results <- as.data.frame(results)
    results$CompID <- rownames(output_comp)
    results_with_info <- merge(results, output_comp, by="CompID", sort=FALSE)
    name_df <- paste(name_project, '_LIMMA_ResultsTable_output_multiple_comparison', multiple_comparison, '.txt', sep="")
    write_dataframe_as_txt_file(results_with_info, name_df)
    
    # Venn diagram
    results <- decideTests(fit2)
    if (ncol(results) <=5){   
      png(paste(name_project, '_LIMMA_Venn_diagram_multiple_comparison_', multiple_comparison, '.png', sep=""), width=7, height=5, units="in", res=150)
      vennDiagram(results)
      dev.off()
    }
    
    # Export 1 log2 FC matrix for multiple comparison
    #The statistic fit2$F and the corresponding fit2$F.p.value combine the comparisons into one F-test. 
    topTableF_comp <- topTableF(fit2, number=nrow(fit2)) #print for each variable score, no selection adj p <0.05
    variablenames <- rownames(topTableF_comp)
    variablenames <- substr(variablenames, 2,nchar(variablenames))
    topTableF_comp$CompID <- variablenames
    topTableF_comp_with_info <- merge(topTableF_comp, output_comp, by="CompID", sort=FALSE)
    name_df <- paste(name_project, '_LIMMA_ToptableF_output_multiple_comparison', multiple_comparison, '.txt', sep="")
    write_dataframe_as_txt_file(topTableF_comp_with_info, name_df)
    
    
    print(paste("finished calculation multiple comparison ", multiple_comparison, " using LIMMA", sep=""))
  }
} 
###################


##########RF##########
## RF for each multiple comparison
if(AMOUNT_OF_MULTIPLE_COMPARISONS >= 1){
  for(multiple_comparison in 1:AMOUNT_OF_MULTIPLE_COMPARISONS){
    
    #multiple_comparison <- 1
    print(paste("start calculation multiple comparison ", multiple_comparison, " using RF", sep=""))
    
    ## PLS-DA for each condition
    suppressMessages(library(randomForest))
    #https://www.r-bloggers.com/2018/01/how-to-implement-random-forests-in-r/
    
    ### remove non-essential data for comparison
    #remove samples that are not included in comp (Na instead of 0/1)
    nsamples_metadata <- scalednormalizedMetadata_samples[scalednormalizedMetadata_samples$Type == 'Sample',] 
    comp_ <- nsamples_metadata[,COLLUMN_NR_MULTIPLE_COMPARISON1+multiple_comparison-1]
    stopifnot(class(comp_)=='integer') #need {0,1,2, ...} arguments, rest is NA and not ""!
    samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    
    
    #remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
    samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
    #samples_matrix_comp_no0 <- samples_matrix_comp[, !apply(samples_matrix_comp == 0, 2, all)]  #remove 0
    samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x)))
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))]
    
    comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_MULTIPLE_COMPARISON1+multiple_comparison-1])
    
    
    #column with labels cbind to matrix for RF
    samples_matrix_comp_no0 <- cbind(samples_matrix_comp_no0, comp)
    
    
    # Split into Train and Validation sets
    # Training Set : Validation Set = 70 : 30 (random)
    set.seed(100)
    train <- sample(nrow(samples_matrix_comp_no0), 0.7*nrow(samples_matrix_comp_no0), replace = FALSE)
    TrainSet <- samples_matrix_comp_no0[train,]
    ValidSet <- samples_matrix_comp_no0[-train,]
    #summary(TrainSet)
    #summary(ValidSet)
    
    
    #if fails:
    try({
      # Create a Random Forest model with parameters == 500trees + nr variabl tried at each split 10
      model1 <- randomForest(comp ~ ., data = TrainSet, importance = TRUE, ntree = 500, mtry = 10)
      #model1
      
      # Using For loop to identify the right mtry for model => manual chack
      #a=c()
      #i=5
      #for (i in 3:8) {
      #  model3 <- randomForest(comp ~ ., data = TrainSet, ntree = 500, mtry = i, importance = TRUE)
      #  predValid <- predict(model3, ValidSet, type = "class")
      #  a[i-2] = mean(predValid == ValidSet$comp)
      #}
      #a
      #plot(3:8,a)
      
      #write summary model with param and score arror rate
      name_report_rf <- paste(name_project,'_RF_scores_multiple_comp.txt', sep="")
      title <- paste('RF_multiple_comparison ', multiple_comparison, ".", sep="")
      append_result_to_report(title, name_report_rf)
      print_rf <- capture.output(print(model1))
      append_result_to_report(print_rf, name_report_rf)
      print_space <- ""
      append_result_to_report(print_space, name_report_rf)
      
      
      ## Predicting on Validation set
      predValid <- predict(model1, ValidSet, type = "class")
      
      # Checking classification accuracy
      append_result_to_report("Predicting on Validation set", name_report_rf)
      accuracy <- mean(predValid == ValidSet$comp)     
      append_result_to_report(paste0("accuracy model: ", accuracy), name_report_rf)
      
      #write confusion matrix
      conf_matrix <- table(predValid,ValidSet$comp)
      append_result_to_report(capture.output(print(conf_matrix)), name_report_rf)
      print_space <- ""
      append_result_to_report(print_space, name_report_rf)
      
      #plot confus matrix nice (long format)
      library(reshape2)
      confusion_matrix_long <- melt(conf_matrix)
      colnames(confusion_matrix_long)[2] <- "ActualValid"
      confusion_matrix_long$predValid <- as.factor(confusion_matrix_long$predValid) #as factor for plot
      confusion_matrix_long$ActualValid <- as.factor(confusion_matrix_long$ActualValid)
      
      cm_plot <- plot_ConfusionMatrix_mcomp(confusion_matrix_long)
      
      png(paste(name_project, "_RF_confusion_matrix_multiple_comparison_" , multiple_comparison, ".png", sep=""), width=7, height=5, units="in", res=150)
      plot(cm_plot)
      dev.off()
      
      
      ## importance_features
      importance_features <- as.data.frame(importance(model1))
      importance_features$CompID <- as.character(rownames(importance_features))
      name_df <- paste(name_project, '_RF_importance_features_multiple_comparison_', multiple_comparison, '.txt', sep="")
      write_dataframe_as_txt_file(importance_features, name_df) 
      
      #plot show the drop in mean accuracy for each of the variables.
      png(paste(name_project, "_RF_importance_features_multiple_comparison_" , multiple_comparison, ".png", sep=""), width=7, height=5, units="in", res=150)
      varImpPlot(model1)
      dev.off()
      
      print(paste("finished calculation multiple comparison ", multiple_comparison, " using RF", sep=""))
    })
  }
}    
###################


##########PCA##########
## PCA for projection 
if(AMOUNT_OF_PROJECTIONS >= 1){
  for(projection in 1:AMOUNT_OF_PROJECTIONS){
    
    #projection <- 1 #for testing
    
    print(paste("start calculation PCA ", projection, " using PCA", sep=""))
    
    ### remove non-essential data for comparison
    #remove samples that are not included in comp (Na instead of 0/1)
    nsamples_metadata <- scalednormalizedMetadata_samples_QC #also QC can be in PCA
    comp_ <- nsamples_metadata[,COLLUMN_NR_PROJECTION1+projection-1]
    stopifnot(class(comp_)=='integer') #need {0,1,2, ...} arguments, rest is NA and not ""!
    samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    
    
    #remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
    samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
    #samples_matrix_comp_no0 <- samples_matrix_comp[, !apply(samples_matrix_comp == 0, 2, all)]  #remove 0
    samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x)))
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))]
    
    comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_PROJECTION1+projection-1])
    
    ### PCA for each projection
    suppressMessages(library(pcaMethods))
    if(PC_AMOUNT == DEFAULT_PC_AMOUNT){
      number_samples <- nrow(samples_matrix_comp_no0)     #Amount of samples
      number_samples <- number_samples-1
    }
    if(PC_AMOUNT == FIXED_PC_AMOUNT){
      number_samples = FIXED_PC_AMOUNT
      if(number_samples >= nrow(samples_matrix_comp_no0)){
        stop("ERROR: Part II: multivariate analysis stopped because number of PCs exceeds number of samples-1.")
      }
    }
    PCA_comp <- pca(samples_matrix_comp_no0,nPcs=number_samples,scale="none",completeObs=TRUE,subset=NULL,cv="q2",center=FALSE) #Calculate full PCA model (max PCs is # of samples - 1)
    #plot(PCA_comp)   #Rule of thumb: amount of PCs to choose: Where Q2 (calcul by cross validation) is max. => use this (Other rule of thumb: total var. explained +-80%)
    #name_plot <- paste(name_project, "_PCA projection " , projection, ".pdf", sep="")
    #pdf(name_plot, height=5,width=7)
    #plot(PCA_comp)
    #dev.off()
    
    png( paste(name_project, "_PCA projection " , projection, ".png", sep=""), width=7, height=5, units="in", res=150)
    plot(PCA_comp)
    dev.off()
    
    #set.seed(42) #Set seed of random number generator for reproducible results
    #Q2CV <- Q2(PCA_comp,originalData = completeObs(PCA_comp),fold=7,nruncv = 5,type="krzanowski") #Choose optimal number of components from this
    #Q2max <- max(Q2CV[,5])
    #nPCsopt <- which(grepl(Q2max, Q2CV[,5])) #nr of row (row1 == pc1)
    #if(nPCsopt == 1){
    nPCsopt <- 3
    #}
    PCAoptnumber <- pca(samples_matrix_comp_no0,nPcs=nPCsopt,scale="none",completeObs=TRUE,subset=NULL,cv="q2",center=FALSE)
    PCAscores_comp <- as.data.frame(PCAoptnumber@scores)
    
    HotellingellipsePCA_comp <- Hotellingellipse(number_samples,PCAscores_comp)  #Calculate Hotelling's T2 ellipse for the PCA, only use PC1 and PC2 for this in the func
    dev.off()
    PCAplot_comp <- plot_pcascores(PCAscores_comp, comp)
    #PCAplot_comp
    #name_plot <- paste(name_project, "_PCA_scoreplot " , projection, ".pdf", sep="")
    #pdf(name_plot) #, height=5,width=7)
    #plot(PCAplot_comp)
    #dev.off()
    
    png(paste(name_project, "_PCA_scoreplot " , projection, ".png", sep=""), width=7, height=5, units="in", res=150)
    plot(PCAplot_comp)
    dev.off()
    
    #pca score plot weight values
    score_weights_pca_comp <- PCAscores_comp[ ,1:3] #PC1, PC2 and PC3 are kept
    rownames(score_weights_pca_comp) <- samples_metadata_comp$SampleName
    name_df <- paste(name_project, '_PCA_weights_projection_', projection, '.txt', sep="")
    write_matrix_as_txt_file(score_weights_pca_comp, name_df) 
    #plot(score_weights_pca_comp)
    
    print(paste("finished calculation PCA ", projection, " using PCA", sep=""))
  }
}    
###################


print("R pipeline - Part II: multivariate analysis - done!")
print(Sys.time())
end_time <- Sys.time()
print(end_time - start_time)
#
#####################

