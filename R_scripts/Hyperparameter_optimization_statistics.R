# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: hyperparameter optimization statistics


##########R Pipeline - Part II: Hyperparameter optimization statistics##########
print(Sys.time())
start_time <- Sys.time()
print("R pipeline - Part II: Hyperparameter optimization statistics - start!")
# Part II: Hyperparameter optimization statistics



#todo exapnd to pareto, log in final stat analysis...
TRANSFORMATION_NOT <- 'no transformation'
TRANSFORMATION_LOG <- 'Log transformation' #natural log with (x+1) so no inf

SCALING_NOT <- 'no scaling'
SCALING_PARETO <- 'Pareto scaling'

TRANSFORMATION <- c(TRANSFORMATION_LOG) #TRANSFORMATION_NOT, 
SCALING <- c(SCALING_PARETO) #SCALING_NOT, 



###repeat necessary steps of part 2: statisitcal analysis (incl. variable filtering)
########merge, skew, variable selection############
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


## set directory to output/param 
dir.create(file.path(path_data_out, 'param'))
path_data_out_param <- file.path(path_data_out, 'param')
setwd(path_data_out_param)


## merge variables from variableMetadata into the sampleMetadata
### make sampleMatrix
library(data.table)
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
  stop("ERROR: Part II: statistical analysis stopped because SampleMetadata and VariableMetadata are incompatible to merge into correct sampleMatrix.")
} 


### start hyperparameter optimization
#########hyperparameter###########
# method: gridsearch method (compare all combinations)
# todo: if needed in future, do with random optimization instead of all combi's

#need to start from original SM every loop:
save_sampleMetadata <- sampleMetadata

#skip combi's QC normalisation and QC-related filters instead of retrieving error if no QC in SM:
if(!("QC" %in% sampleMetadata$Type | "IQC" %in% sampleMetadata$Type)){
  COEFICIENT_OF_VARIANCE <- c(COEFICIENT_OF_VARIANCE_0) #instead of vector with options described above
  QC_PREFILTER <- c(NO_FILTER_QC_POOL) #instead of vector with options described above
  NORMALIZE_METHOD <- NORMALIZE_NOT #instead of vector with options described above
}

## list combinations of hyperparameters to optimize
report_hyperparameter_test <- NULL
combi_nr <-1
amount_of_filter1 <- length(COEFICIENT_OF_VARIANCE)
amount_of_filter2 <- length(QC_PREFILTER)
amount_of_is_param <- length(NORMALIZE_METHOD0)
amount_of_tic_param <- length(NORMALIZE_METHOD1)
amount_of_qc_param <- length(NORMALIZE_METHOD)
amount_of_transf_param <- length(TRANSFORMATION)
amount_of_scaling_param <- length(SCALING)


for(filter1 in 1:amount_of_filter1){
  for(filter2 in 1:amount_of_filter2){
  
    #need to start from original SM every loop: #start from merged SM
    sampleMetadata <- save_sampleMetadata
    
    #filter1 <- 1
    #filter2 <- 1
    #combi_nr <-1
    #print(COEFICIENT_OF_VARIANCE[filter1])
    #print(QC_PREFILTER[filter2])
    
    
    ### variable filtering
    ## optional: Delete variables with coeficent of variation of QCs > 30%
    if(COEFICIENT_OF_VARIANCE[filter1] == COEFICIENT_OF_VARIANCE_30){
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
      #write_dataframe_as_txt_file(sampleMetadata_filter_coeficient_variance, paste0('sampleMetadata_filter_coeficient_variance', combi_nr, '.txt'))
      
      sampleMetadata <- sampleMetadata_filter_coeficient_variance
    }
    if(COEFICIENT_OF_VARIANCE[filter1] == COEFICIENT_OF_VARIANCE_0){
      sampleMetadata <- sampleMetadata
    }
    
    
    ## optional: retain variable in present in 80% of the QCs samples
    if(QC_PREFILTER[filter2] == FILTER_QC_POOL_80){
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
      #write_dataframe_as_txt_file(sampleMetadata_filter_QCpool, paste0('sampleMetadata_filter_QCpool', combi_nr, '.txt'))
      
      sampleMetadata <- sampleMetadata_filter_QCpool
    }
    if(QC_PREFILTER[filter2] == NO_FILTER_QC_POOL){
      sampleMetadata <- sampleMetadata
    }
   
    
    #########skew###########
    ## check if data has a normal distribution per Type(s)
    ### Check normality of all variables across all the samples (blanco, QC, STD, sample)
    all_matrix <- from_df_to_matrix(sampleMetadata)
    
    amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(all_matrix)
    total_amount_of_variables <- ncol(all_matrix)
    percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)
    
    report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across all the samples (", percentage, "%).", sep="")
    name_report_normality <- paste(name_project,'_Report_normalities.txt', sep="")
    append_result_to_report(report_check_normality, name_report_normality)
    
    skewkurtplot <- plot_skewkurt(all_matrix)
    #pdf(paste(name_project,'_skewkurtplot_all', combi_nr, '.pdf', sep=""),height=5,width=7)
    #plot(skewkurtplot)
    #dev.off()
    
    png(paste(name_project,'_skewkurtplot_all', combi_nr, '.png', sep=""), width=7, height=5, units="in", res=150)
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
    #pdf(paste(name_project,'_skewkurtplot_samples_QC', combi_nr, '.pdf', sep=""),height=5,width=7)
    #plot(skewkurtplot)
    #dev.off()
    
    png(paste(name_project,'_skewkurtplot_samples_QC', combi_nr, '.png', sep=""), width=7, height=5, units="in", res=150)
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
    #pdf(paste(name_project,'_skewkurtplot_samples', combi_nr, '.pdf', sep=""),height=5,width=7)
    #plot(skewkurtplot)
    #dev.off()
    
    png(paste(name_project,'_skewkurtplot_samples', combi_nr, '.png', sep=""), width=7, height=5, units="in", res=150)
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
      #pdf(paste(name_project,'_skewkurtplot_ISTD', combi_nr, '.pdf', sep=""),height=5,width=7)
      #plot(skewkurtplot)
      #dev.off()
      
      png(paste(name_project,'_skewkurtplot_ISTD', combi_nr, '.png', sep=""), width=7, height=5, units="in", res=150)
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
      #pdf(paste(name_project,'_skewkurtplot_STD', combi_nr, '.pdf', sep=""),height=5,width=7)
      #plot(skewkurtplot)
      #dev.off()
      
      png(paste(name_project,'_skewkurtplot_STD', combi_nr, '.png', sep=""), width=7, height=5, units="in", res=150)
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
      #pdf(paste(name_project,'_skewkurtplot_QC', combi_nr, '.pdf', sep=""),height=5,width=7)
      #plot(skewkurtplot)
      #dev.off()
      
      png(paste(name_project,'_skewkurtplot_QC', combi_nr, '.png', sep=""), width=7, height=5, units="in", res=150)
      plot(skewkurtplot)
      dev.off()
    }
    ####################
    
    
    
    #########loop tic, loop qc###########
    #need to start from original SM every loop: #changed since log/pareto so re-save for next loop
    save_sampleMetadataLogPar <- sampleMetadata

    for(isparam in 1:amount_of_is_param){
      for(ticparam in 1:amount_of_tic_param){
        for(param in 1:amount_of_qc_param){
          
          #need to start from original SM every loop: #changed since log/pareto
          sampleMetadata <- save_sampleMetadataLogPar
          
          #isparam <- 1
          #ticparam <- 1  #for testing
          #param <- 1   #for testing
          #print(combi_nr)
          #print(NORMALIZE_METHOD0[isparam])
          #print(NORMALIZE_METHOD1[ticparam])
          #print(NORMALIZE_METHOD[param])
            
          # individual loops per type of hyperparameter (is, tic, qc)
          
          ##########Normalisation##########
          
          ## Perform 0th normalisation with IS (optional)
          if(NORMALIZE_METHOD0[isparam] == NORMALIZE_WITH_IS){ 
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
          if(NORMALIZE_METHOD0[isparam] == NORMALIZE_NOT0){
            sampleMetadata <- sampleMetadata #nothing changes
          }
          #write_dataframe_as_txt_file(sampleMetadata, paste0('NORMALIZE_METHOD0_sampleMetadata', combi_nr, '.txt'))
          
          
          ## Perform 1th normalisation with TIC (optional)
          if(NORMALIZE_METHOD1[ticparam] == NORMALIZE_WITH_TIC){ 
            #TIC or total ion count normalisation accord to literature
            #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3124646/
            #this normalization will be based on the sum of all intensity values in the mass spectrum (i.e., TIC)
            #= intensity compID 1 / sum(intenstity alls compIDs) for 1 sample/spectrum
            averageTIC_before <- mean(as.matrix(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)])) #sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)])
            
            #devide by sum per sample
            normalized <-  apply(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)], 1, function(x){t(x/sum(x)) }) #iterate over rows =1 (samples), but saves results as column!!! 
            if(is.nan(mean(normalized))){
              #in case (when working with sparse data: 1 sample with for all features area=0, x devide sum(x) = /0 = NaN
              normalized <- data.frame(sapply(as.data.frame(normalized), function(x) ifelse(is.nan(x), 0, x)))
            }
            sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- t(normalized) #needs to be transposed for correct tic-norm
            
            #must be re-scaled to same order of scale (afrondingsfouten bij kleine getallen tijdens berek PCA, daarom all getallen maal zelfde factor)
            averageTIC_after <- mean(as.matrix(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)]))
            factor <- averageTIC_before/averageTIC_after
            sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- factor * sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)]
          }
          if(NORMALIZE_METHOD1[ticparam] == NORMALIZE_WITH_ITIC){
            #done by progenisis as 'tic normalisation', iterate over samples
            #= intensity sample 1 / sum(intenstity all samples) for 1 metabolite
            sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- apply(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)], 2, function(x) x/sum(x)) #iterate over cols =2 (compids)
          }
          if(NORMALIZE_METHOD1[ticparam] == NORMALIZE_WITH_MAX){
            #'sumarizednormalized' named at lca
            #= intensity sample 1 / max(intenstity all samples) for 1 metabolite
            sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- apply(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)], 2, function(x) x/max(x)) #iterate over cols =2 (compids)
          }
                                                                                        if(NORMALIZE_METHOD1[ticparam] == NORMALIZE_WITH_MEDIAN){
            #this normalization will be based on the median of all intensity values in the mass spectrum, medI 
            #= intensity compID 1 / median(intenstity alls compIDs) for 1 sample/spectrum
            average_before <- mean(as.matrix(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)])) #sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)])
            
            #devide by sum per sample
            normalized <-  apply(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)], 1, function(x){t(x/median(x)) }) #iterate over rows =1 (samples), but saves results as column!!! 
            if(is.nan(mean(normalized))){
              #in case (when working with sparse data: 1 sample with for all features area=0, x devide sum(x) = /0 = NaN OR Inf*
              #*) "Depending on the programming environment and the type of number (e.g. floating point, integer) being divided by zero, it may generate positive or negative infinity by the IEEE 754 floating point standard [...]"
              normalized <- data.frame(sapply(as.data.frame(normalized), function(x) ifelse(is.nan(x), 0, x)))
              normalized <- data.frame(sapply(as.data.frame(normalized), function(x) ifelse(is.infinite(x), 0, x)))
            }
            if(min(normalized) == max(normalized)){
              #in case sparse data, can be ALL data becomes 0 values... so best errormessage
              stop("ERROR: Part II: multivariate analysis stopped because median normalisation is not applicable on sparse data.")
            }
            sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- t(normalized) #needs to be transposed for correct tic-norm

            #must be re-scaled to same order of scale (afrondingsfouten bij kleine getallen tijdens berek PCA, daarom all getallen maal zelfde factor)
            average_after <- mean(as.matrix(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)]))
            factor <- average_before/average_after
            sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- factor * sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)]
          }
          if(NORMALIZE_METHOD1[ticparam] == NORMALIZE_WITH_QUANTILE){
            #https://en.wikipedia.org/wiki/Quantile_normalization
            #https://www.statology.org/quantile-normalization-in-r/ 
            average_before <- mean(as.matrix(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)])) #sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)])
            
            #perform quantile normalization
            library(preprocessCore)
            sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- as.data.frame(normalize.quantiles(as.matrix(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)])))
            
            #check if correct see quantiles same, ok
            #sapply(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)], function(x) quantile(x, probs = seq(0, 1, 1/4)))
            
            #must be re-scaled to same order of scale (afrondingsfouten bij kleine getallen tijdens berek PCA, daarom all getallen maal zelfde factor)
            average_after <- mean(as.matrix(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)]))
            factor <- average_before/average_after
            sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)] <- factor * sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)]
          }
          if(NORMALIZE_METHOD1[ticparam] == NORMALIZE_NOT1){
            sampleMetadata <- sampleMetadata #nothing changes
          }
          #write_dataframe_as_txt_file(sampleMetadata, paste0('NORMALIZE_METHOD1_sampleMetadata', combi_nr, '.txt'))
          
          
          ### Normalize samples
          if(NORMALIZE_METHOD[param] == NORMALIZE_WITH_QCs){
            stopifnot("QC" %in% sampleMetadata$Type | "IQC" %in% sampleMetadata$Type)
            samples_metadata <- sampleMetadata[sampleMetadata$Type == 'Sample',]
            QC_metadata <- sampleMetadata[sampleMetadata$Type == 'QC' |sampleMetadata$Type == 'IQC',] ##also allow IQC
            normalized_samples_metadata <- normalize_with_average_QCs(samples_metadata, QC_metadata)
          }
          if(NORMALIZE_METHOD[param] == NORMALIZE_WITH_IQCs){
            stopifnot("QC" %in% sampleMetadata$Type | "IQC" %in% sampleMetadata$Type)
            samples_IQC_metadata <- sampleMetadata[sampleMetadata$Type == 'IQC' | sampleMetadata$Type == 'Sample' |sampleMetadata$Type == 'QC' ,] ##also allow QC
            normalized_samples_metadata <- normalize_with_IQCs(samples_IQC_metadata)  #no info, except sampleName as last collumn
          }
          if(NORMALIZE_METHOD[param] == NORMALIZE_QC_RLSC){
            stopifnot("QC" %in% sampleMetadata$Type | "IQC" %in% sampleMetadata$Type)
            #order qc == 1 in function normalize_QC_RLSC
            normalized_samples_metadata <- normalize_QC_RLSC(sampleMetadata) #need to be in same format = matrix w/o info + no inf/na
          }
          if(NORMALIZE_METHOD[param] == NORMALIZE_NOT){
            samples_metadata <- sampleMetadata[sampleMetadata$Type == 'Sample',]
            normalized_samples_metadata <- normalize_not(samples_metadata) #change nothing but need to be in same format = matrix w/o info + no inf/na
          }
          #write_dataframe_as_txt_file(normalized_samples_metadata, paste0('normalized_samples_matrix', combi_nr, '.txt'))
          
          ### Normalize (I)QCs (for plotting in pca-x instead of OPLSDA) 
          if ("IQC" %in% sampleMetadata$Type | "QC" %in% sampleMetadata$Type){
            QC_metadata <- sampleMetadata[sampleMetadata$Type == 'QC' | sampleMetadata$Type == 'IQC',]
            if(NORMALIZE_METHOD[param] == NORMALIZE_WITH_QCs){
              normalized_QC_metadata <- normalize_with_average_QCs(QC_metadata, QC_metadata)
            }
            if(NORMALIZE_METHOD[param] == NORMALIZE_WITH_IQCs){
              normalized_QC_metadata <- normalize_with_average_QCs(QC_metadata, QC_metadata) #this uses average from all QCs 
              #QC_metadata_use <- QC_metadata[order(QC_metadata$Order),]  ##if keeping QCs and wanting to divide each qc in a pair (based on order) by the average across the pair
              #QC_only <- data.table(QC_metadata_use[,c(1, 8:ncol(QC_metadata_use))]) 
              #g <- QC_only[, rep(1:.N, each=2, length=.N)]
              #nms <- names(QC_only)[-1]
              #divbyMean <- function(x) replace(x, 1:2, x[1:2]/mean(x[1:2]))
              #out <- as.data.frame(QC_only[, c(nms) := lapply(.SD[,-1], divbyMean), by=g])
              #normalized_QC_metadata <- dplyr::relocate(out, 1, .after=last_col())    ##reposition the sample names to last col        
            }
            if(NORMALIZE_METHOD[param] == NORMALIZE_QC_RLSC){
              normalized_QC_metadata <- NULL  #empty because all samples+qc + std+blanks included in loess output samples
            }
            if(NORMALIZE_METHOD[param] == NORMALIZE_NOT){
              normalized_QC_metadata <- normalize_not(QC_metadata)
            }
          }
          if ("IQC" %in% sampleMetadata$Type == FALSE & "QC" %in% sampleMetadata$Type == FALSE){
            normalized_QC_metadata <- normalize_not(QC_metadata)
          }
          
          ### replace intensities with normalized values (both samples + (I)QCs)
          normalized_samples_QC_metadata <- rbind(normalized_samples_metadata, normalized_QC_metadata)
          #write_dataframe_as_txt_file(normalized_samples_QC_metadata, paste0('normalized_samples_QC_matrix', combi_nr, '.txt'))
          
          sampleMetadata_until_start_variables <- subset(sampleMetadata, select = c(1:(COLLUMN_NR_START_VARIABLES-1))) 
          sampleMetadata_samples_QC_until_start_variables <- sampleMetadata_until_start_variables[sampleMetadata$Type == 'Sample' | sampleMetadata$Type == 'QC' | sampleMetadata$Type == 'IQC',]
          
          normalizedMetadata <- merge_accoding_to_SampleName(sampleMetadata_samples_QC_until_start_variables, normalized_samples_QC_metadata)
          #write_dataframe_as_txt_file(normalizedMetadata, paste0('NORMALIZE_METHOD_normalizedMetadata', combi_nr, '.txt'))
          ####################
          
          
          #########skew###########
          ## Check normality of biological samples after normalization
          nsamples_metadata <- normalizedMetadata[normalizedMetadata$Type == 'Sample',] 
          nsamples_matrix <- from_df_to_matrix(nsamples_metadata)
          
          amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(nsamples_matrix)
          total_amount_of_variables <- ncol(nsamples_matrix)
          percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)
          
          report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the samples after normalisation (", percentage, "%). ", combi_nr, sep="")
          append_result_to_report(report_check_normality, name_report_normality)
          
          skewkurtplot <- plot_skewkurt(nsamples_matrix)
          #pdf(paste(name_project,'_skewkurtplot_normalized_samples', combi_nr, '.pdf', sep=""),height=5,width=7)
          #plot(skewkurtplot)
          #dev.off()
          
          png(paste(name_project,'_skewkurtplot_normalized_samples', combi_nr, '.png', sep=""), width=7, height=5, units="in", res=150)
          plot(skewkurtplot)
          dev.off()
          ####################
         
           

          #############loops log, pareto###################
          #save for each combi start with chosen normalisation
          save_normalizedMetadata <- normalizedMetadata

          
          # individual loops per type of hyperparameter (transf, scaling)
          for(transf_param in 1:amount_of_transf_param){
            for(scaling_param in 1:amount_of_scaling_param){
              
              #transf_param <- 1
              #scaling_param <- 1
              #print(combi_nr)
              #print(TRANSFORMATION[transf_param])
              #print(SCALING[scaling_param])
              
              #load from saved 
              normalizedMetadata <- save_normalizedMetadata
              
              ## Log transformation and Pareto scaling for samples + QCs toghether
              if(TRANSFORMATION[transf_param] == TRANSFORMATION_NOT){
                lognormalizedMetadata_samples_QC <- normalizedMetadata #nothing changes
              }
              if(TRANSFORMATION[transf_param] == TRANSFORMATION_LOG){
                lognormalizedMetadata_samples_QC <- log_transformation(normalizedMetadata)
              }
              if(SCALING[scaling_param] == SCALING_NOT){
                scalednormalizedMetadata_samples_QC <- normalizedMetadata #nothing changes
              }
              if(SCALING[scaling_param] == SCALING_PARETO){
                scalednormalizedMetadata_samples_QC <- Pareto_scaling(lognormalizedMetadata_samples_QC)
              }
              #name_df <- paste(name_project,'_normalized, log transformed and pareto scaled samples and QC Metadata', combi_nr, '.txt', sep="")
              #write_dataframe_as_txt_file(scalednormalizedMetadata_samples_QC, name_df)
              
              # cols with value O everywhere before log/pareto are now NA. NA -> 0 for ok sharipo, delete not needed at this time (below @pca are deleted)
              scaledmatrix_samples_QC <- from_df_to_matrix(scalednormalizedMetadata_samples_QC)
              scaledmatrix_samples_QC <- data.frame(sapply(scaledmatrix_samples_QC, function(x) ifelse(is.na(x), 0, x)))
              
              #check normal
              amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(scaledmatrix_samples_QC)
              total_amount_of_variables <- ncol(scaledmatrix_samples_QC)
              percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)
              report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the samples after normalisation, log transformation and pareto scaling of both samples and QCs (", percentage, "%). ", combi_nr, sep="")
              append_result_to_report(report_check_normality, name_report_normality)
              
              skewkurtplot <- plot_skewkurt(scaledmatrix_samples_QC)
              #pdf(paste(name_project,'_skewkurtplot_samples_QCs_scaled_and_log-transformed', combi_nr, '.pdf', sep=""),height=5,width=7)
              #plot(skewkurtplot)
              #dev.off()
              
              png(paste(name_project,'_skewkurtplot_samples_QCs_scaled_and_log-transformed', combi_nr, '.png', sep=""), width=7, height=5, units="in", res=150)
              plot(skewkurtplot)
              dev.off()
              
              ## Log transformation and Pareto scaling for only samples
              nsamples_metadata <- normalizedMetadata[normalizedMetadata$Type == 'Sample',]
              if(TRANSFORMATION[transf_param] == TRANSFORMATION_NOT){
                lognormalizedMetadata_samples <- nsamples_metadata #nothing changes
              }
              if(TRANSFORMATION[transf_param] == TRANSFORMATION_LOG){
                lognormalizedMetadata_samples <- log_transformation(nsamples_metadata)
              }
              if(SCALING[scaling_param] == SCALING_NOT){
                scalednormalizedMetadata_samples <- nsamples_metadata #nothing changes
              }
              if(SCALING[scaling_param] == SCALING_PARETO){
                scalednormalizedMetadata_samples <- Pareto_scaling(lognormalizedMetadata_samples)
              }
              #name_df <- paste(name_project,'_normalized, log transformed and pareto scaled biological samples Metadata', combi_nr, '.txt', sep="")
              #write_dataframe_as_txt_file(scalednormalizedMetadata_samples, name_df)

              
              # cols with value O everywhere before log/pareto are now NA. NA -> 0 for ok sharipo, delete not needed at this time (below @pca are deleted)
              scaledmatrix_samples <- from_df_to_matrix(scalednormalizedMetadata_samples)
              scaledmatrix_samples <- data.frame(sapply(scaledmatrix_samples, function(x) ifelse(is.na(x), 0, x)))
              
              #check normal
              amount_of_variables_with_normal_distribution <- Shapiro_Wilk_test(scaledmatrix_samples)
              total_amount_of_variables <- ncol(scaledmatrix_samples)
              percentage <- round(((amount_of_variables_with_normal_distribution/total_amount_of_variables)*100), digits=2)
              report_check_normality <- paste("There are " , amount_of_variables_with_normal_distribution, " variables retrieved with a normal distribution from the total of " , total_amount_of_variables , " variables across the samples after normalisation, log transformation and pareto scaling of samples (", percentage, "%). ", combi_nr, sep="")
              append_result_to_report(report_check_normality, name_report_normality)
              
              skewkurtplot <- plot_skewkurt(scaledmatrix_samples)
              #pdf(paste(name_project,'_skewkurtplot_samples_scaled_and_log-transformed', combi_nr, '.pdf', sep=""),height=5,width=7)
              #plot(skewkurtplot)
              #dev.off()
              
              png(paste(name_project,'_skewkurtplot_samples_scaled_and_log-transformed', combi_nr, '.png', sep=""), width=7, height=5, units="in", res=150)
              plot(skewkurtplot)
              dev.off()
              ####################
              
              
              
              ##########PCA 1 projection##########
              ## PCA for projection 
              projection <- 1   #!! always calculate KNN with 1st projection
            
              print(paste("start calculation PCA ", projection, "-", combi_nr, sep=""))
              
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
              
              #max 10pcs, even if default:
              if(number_samples > 30){
                number_samples <- 30
              }
      
              try({ 
                PCA_comp <- pca(samples_matrix_comp_no0,nPcs=number_samples,scale="none",completeObs=TRUE,subset=NULL,cv="q2",center=FALSE) #Calculate full PCA model (max PCs is # of samples - 1)
                #plot(PCA_comp)   #Rule of thumb: amount of PCs to choose: Where Q2 (calcul by cross validation) is max. => use this (Other rule of thumb: total var. explained +-80%)
                #name_plot <- paste(name_project, "_PCA projection " , projection, combi_nr, ".pdf", sep="")
                #pdf(name_plot, height=5,width=7)
                #plot(PCA_comp)
                #dev.off()
                
                png( paste(name_project, "_PCA projection " , projection, "-", combi_nr, ".png", sep=""), width=7, height=5, units="in", res=150)
                plot(PCA_comp)
                dev.off()
                
                ##############faster, do with first 2 PCs!!!!##############
                #set.seed(42) #Set seed of random number generator for reproducible results
                #Q2CV <- Q2(PCA_comp,originalData = completeObs(PCA_comp),fold=7,nruncv = 5,type="krzanowski") #Choose optimal number of components from this
                #Q2max <- max(Q2CV[,5])
                #nPCsopt <- which(grepl(Q2max, Q2CV[,5])) #nr of row (row1 == pc1)    #faster, do with first 2 PCs!!!!
                #if(nPCsopt == 1){
                nPCsopt <- 2
                #}
                #################
                  
                PCAoptnumber <- pca(samples_matrix_comp_no0,nPcs=nPCsopt,scale="none",completeObs=TRUE,subset=NULL,cv="q2",center=FALSE)
                PCAscores_comp <- as.data.frame(PCAoptnumber@scores)
                
                HotellingellipsePCA_comp <- Hotellingellipse(number_samples,PCAscores_comp)  #Calculate Hotellis ng'T2 ellipse for the PCA
                dev.off()
                PCAplot_comp <- plot_pcascores(PCAscores_comp, comp) 
                #PCAplot_comp
                #name_plot <- paste(name_project, "_PCA_scoreplot " , projection, combi_nr, ".pdf", sep="")
                #pdf(name_plot) #, height=5,width=7)
                #plot(PCAplot_comp)
                #dev.off()
                
                png(paste(name_project, "_PCA_scoreplot " , projection, "-", combi_nr, ".png", sep="")) #, width=7, height=5, units="in", res=150)
                plot(PCAplot_comp)
                dev.off()
                
                #pca score plot weight values
                score_weights_pca_comp <- PCAscores_comp[ ,1:2] #PC1 and PC2 are kept
                rownames(score_weights_pca_comp) <- samples_metadata_comp$SampleName
                name_df <- paste(name_project, '_score_weights_projection_', projection, combi_nr, '.txt', sep="")
                write_matrix_as_txt_file(score_weights_pca_comp, name_df) 
                #plot(score_weights_pca_comp)
              })
                
              print(paste("finished calculation PCA ", projection, "-", combi_nr, sep=""))
            
              ###################
              
                
              #########Score with KNN##########
              ### score each combination using KNN (! << testset kleinste groep # neighbors):
              #https://towardsdatascience.com/k-nearest-neighbors-algorithm-with-examples-in-r-simply-explained-knn-1f2c88da405c
          
              #score_weights_pca_comp <- PCAscores_comp[ ,1:2] #for testing
              
              #pca score plot weight values (start with this from projections write to txt scores, only 2PCs) #todo: now scores + comp info from last PCA projection in loop above...
              score_weights_pca_comp$Group <- comp
              score_weights_pca_comp$Group <- as.numeric(score_weights_pca_comp$Group)
              score_weights_pca_comp$Group[is.na(score_weights_pca_comp$Group)] <- 0 #no NAs in knn allowed
              
              ##cross-validata 7times:
              scores <- NULL
              max_crossval <- 7
              
              for(crossval in 1:max_crossval){
                #Generate a random number that is 70% of the total number of rows in dataset.
                ran <- sample(1:nrow(score_weights_pca_comp), 0.7 * nrow(score_weights_pca_comp)) 
                
                train <- score_weights_pca_comp[ran,1:2] #col of 2 pcs, part 1-100
                test <- score_weights_pca_comp[-ran,1:2] #92 of part 2
                train_category <- score_weights_pca_comp$Group[ran] #part 1-100 of train
                test_category <- score_weights_pca_comp$Group[-ran]
                
                ##PCA-KNN
                library(class) ##run knn function
                KNN_NR <- min(summary(comp))-1 #amount of smallest group in comp -1 so odd number
                if(KNN_NR %% 2 == 0){
                  KNN_NR <- KNN_NR-1 #if even: -2 so odd number
                } 
                pred <- knn(as.matrix(train),as.matrix(test),cl=as.matrix(train_category),k=KNN_NR) #/search 4 nearest neighors
                pred <- as.character(pred)
                #length(pred)
                
                
                ##PCA-SVM 
                #library(e1071)
                #https://www.datacamp.com/community/tutorials/support-vector-machines-r
                #x <- as.matrix(train)
                #y <- as.matrix(train_category)
                #y <- as.factor(y)
                #model <- svm(x, y, kernel = "linear") 
                #print(model)
                #summary(model)
                #pred <- predict(model, as.matrix(test))
                #pred <- as.character(pred)

                
                ##PCA-LDA
                #library(MASS)
                #y <- as.matrix(train_category)
                #y <- as.factor(y)
                #model <- lda(y ~ ., data = train) #train is df
                #pred <- predict(model, test) #test is df
                #pred <- as.character(pred$class)
                #length(pred)
                
                
                ##make confusion matrix
                tab <- table(pred,test_category)
                
                
                
                ##score using WITH qc's, calc accuracy only for QC samples in Projection1
                if("QC" %in% samples_metadata_comp$Type |"IQC" %in% samples_metadata_comp$Type){
                  
                  #find number in projection1 of the QCs
                  subset_QCS <- samples_metadata_comp[samples_metadata_comp$Type == "QC" |samples_metadata_comp$Type == "IQC" ,]
                  nr_QC_projection1 <- subset_QCS$Projection1
                  
                  #stop if QCs in more than 1 group in Projection1
                  if(sd(nr_QC_projection1) > 0){
                    stop("Error: you can only have 1 group defined as QCs in Projection1")
                  }
                  
                  #QC is nr in projection1:
                  nr_QC_projection1 <- mean(nr_QC_projection1) #only 1 number needed
                  #in not all groups in test_category and/or pred: incomplete matrix 5*6 bv => find based on row and col name by using "3" ipv 3
                  nr_QC_projection1 <- as.character(nr_QC_projection1) #as string for column name in tab
                  
                  #if no pred/test_category contains the QC as group, next crossval + print (ok?):
                  if(!nr_QC_projection1 %in% rownames(tab)){
                    print(paste0("Remark: no prediction possible with group of QCs present in crossval ", crossval, " of ", max_crossval, " for combi_nr ", combi_nr))
                    next
                  }
                  
                  ##change confusion matrix with adjusted weigths
                  tab_weighted <- as.data.frame.matrix(tab) 
                  
                  #on row of QC group predicted: if good: leave, if bad each bad*F
                  value_QC_corr <- tab_weighted[nr_QC_projection1,nr_QC_projection1] #row first, quoted column name second, eg 3,3
                  
                  tab_weighted <- as.matrix(tab_weighted[nr_QC_projection1,])
                  score <- value_QC_corr/sum(tab_weighted)*100
                  
                }  
                
                ##score using crossval WO qc's, so same weight for each group in Projection1
                if(!("QC" %in% samples_metadata_comp$Type |"IQC" %in% samples_metadata_comp$Type)){
                  
                  #calc score
                  accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
                  #for df add score
                  score <- accuracy(tab) #% correct
                  #print(paste0(KNN_NR, " knn's -> accuracy: ", score))
                }
                
                #PCs_QCs <- score_weights_pca_comp[score_weights_pca_comp$Group == nr_QC_projection1,]
                #PC1_QCs <- mean(PCs_QCs$PC1)
                #PC2_QCs <- mean(PCs_QCs$PC2)
                #score <- PC2_QCs
                #since have to rewrite all, work with score op pc1 (number, so averae cv etc.)
                scores <- c(scores, score)
              }
              
              #score of knn (average)
              average_score <- mean(scores)
              
              #append df
              report_hyperparameter_test_combi <- NULL
              report_hyperparameter_test_combi$combination <- combi_nr
              report_hyperparameter_test_combi$QC_variance <- as.character(COEFICIENT_OF_VARIANCE[filter1])
              report_hyperparameter_test_combi$QC_presence <- as.character(QC_PREFILTER[filter2])
              report_hyperparameter_test_combi$IS_normalisation <- as.character(NORMALIZE_METHOD0[isparam])
              report_hyperparameter_test_combi$TIC_normalisation <- as.character(NORMALIZE_METHOD1[ticparam])
              report_hyperparameter_test_combi$QC_normalisation <- as.character(NORMALIZE_METHOD[param])
              report_hyperparameter_test_combi$Transformation <- as.character(TRANSFORMATION[transf_param])
              report_hyperparameter_test_combi$Scaling <- as.character(SCALING[scaling_param])
              #report_hyperparameter_test_combi$neighbours_nr <- KNN_NR
              report_hyperparameter_test_combi$score <- average_score
              report_hyperparameter_test_combi <- as.data.frame(report_hyperparameter_test_combi)
              
              report_hyperparameter_test <- rbind(report_hyperparameter_test, report_hyperparameter_test_combi)
              ###################
              
              
              combi_nr <- combi_nr + 1
            }
          }
        }
      }
    }
  }
}

#write resutls KNN scoring
write_dataframe_as_txt_file(report_hyperparameter_test, paste0(name_project,'_report_hyperparameter_test.txt'))

#best hyperparameter combination
optimal_combination <- report_hyperparameter_test[which.max(report_hyperparameter_test$score),] #row select with max score = accuraty correct/false
write_dataframe_as_txt_file(optimal_combination, paste0(name_project,'_report_optimal_combination.txt'))

#select best and save for part2: statistical analysis
COEFICIENT_OF_VARIANCE <- optimal_combination[2] 
QC_PREFILTER <- optimal_combination[3]
NORMALIZE_METHOD0 <- optimal_combination[4] 
NORMALIZE_METHOD1 <- optimal_combination[5] #4e col is tic 
NORMALIZE_METHOD <- optimal_combination[6] #5e col is qc
TRANSFORMATION <- optimal_combination[7]
SCALING <- optimal_combination[8]


### zip results folder
if (CODE_RUN_MODE == CODE_AUTORUN){
  files2zip <- dir(path_data_out_param, full.names = FALSE)
  zip(zipfile = 'paramZip', files = files2zip)
  file.remove(files2zip)
  write_dataframe_as_txt_file(report_hyperparameter_test, paste0(name_project,'_report_hyperparameter_test.txt')) #copy outside zip
  setwd(path_data_out)
  write_dataframe_as_txt_file(optimal_combination, paste0(name_project,'_report_optimal_combination.txt')) #set copy in output for reporting
}


print("R pipeline - hyperparameter optimization statistics - done!")
print(Sys.time())
end_time <- Sys.time()
print(end_time - start_time)
#
#####################
