## @knitr INFO
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Configuration


##########Global_settings##########

## @knitr settings
## options
RUN_CODE <- 'run this part of the pipeline'
DONT_RUN_CODE <- 'skip this part of the pipeline, keeps order: pre-processing - targeted analysis - statistical analysis - annotation'

## Adjustments
#Project name:
EXPERIMENT <- 'example' #structured and short, see READ_ME
POLARITY <- "positive" #{"positive", "negative"} #needed this format for annotation to work
# file_conversion and pre-processing need to be performed seperate
# only choose "both" when no pre-processing needed (eg. merge both ionisationmodes)
USER_COMMENT <- "Tutorial comment" #Add info about experiment, eg. explain (Multiple)Comparisons, to include in reports

RUN_PART_STATISTICAL <- RUN_CODE

#
#####################



##########Statistical_analysis##########
if(RUN_PART_STATISTICAL == RUN_CODE){

  ## @knitr statistics
  ## options
  #Source of file1 variableMetadata:
  VARIABLEMETADATA_FROM_PIPELINE <- 'automatically finds variableMetadata from R pipeline code (based on projectname), present in input folder'
  VARIABLEMETADATA_EXTERN <- 'you choose name of variableMetadata.txt below, present in input folder' #see READ_ME for format file

  #pretreatment options:
  HYPERPARAMETERS_AUTO <- 'automatically finds the best combination for normalisation (TIC, QC) parameters of the dataset, using info from "Projection1" (with all groups + QCs)'
  HYPERPARAMETERS_MANUAL <- 'mannually select parameters from options described below'

  #variable filtering options:
  COEFICIENT_OF_VARIANCE_30 <- 'QC variance < 30%'      #'retain only metabolietes with a variation of the coeficients (QCs) is below 30%'
  COEFICIENT_OF_VARIANCE_0 <- 'No QC variance filter'   #'retain all metabolites in analysis'

  FILTER_QC_POOL_80 <- 'QC presence > 80%'              #'retain metabolites if present in at least 80% of the QCs'
  NO_FILTER_QC_POOL <- 'No QC presence filter'          #'retain all metabolites, independent of the presente in the QCs'

  #normalisation options:
  #NORMALIZE_METHOD0 options
  NORMALIZE_WITH_IS <- 'IS normalisation'               #'devides the intensity of each variable (metabolite) by the intensity of the IS of this variable per sample/spectrum'
  NORMALIZE_NOT0 <- 'No IS normalisation'                #'if you do not want to perform normalisation, eg. when only coeficient of variance <30% instead'

  #NORMALIZE_METHOD1 options
  NORMALIZE_WITH_TIC <- 'TIC normalisation'             #'devides the intensity of each variable (metabolite) by the sum of all metabolotite intensities per sample/spectrum (TIC)'
  NORMALIZE_WITH_ITIC <- 'Autoscaling (ITIC)'           #'devides the intensity of each variable (metabolite) by the sum of all metabolotite intensities per CompID (inverse TIC)'
  NORMALIZE_WITH_MAX <- 'Max normalisation'             #'devides the intensity of each variable (metabolite) by the max metabolotite intensity per CompID'
  NORMALIZE_WITH_MEDIAN <- 'Median normalisation'       #'devides the intensity of each variable (metabolite) by the median of all metabolotite intensities per sample/spectrum (MedI)'
  NORMALIZE_WITH_QUANTILE <- 'Quantile normalisation'   #'makes the distributions identical in statistical properties across samples'
  NORMALIZE_NOT1 <- 'No TIC normalisation' 				#'if you do not want to perform normalisation, eg. when only coeficient of variance <30% instead'

  #NORMALIZE_METHOD options
  NORMALIZE_WITH_IQCs <- 'IQC normalisation'            #'for the intensities of each variable (metabolite): devides the samples with the mean intenstity of IQC present below the samples. So need to have IQC after samples'
  NORMALIZE_QC_RLSC <- 'LOESS normalisation'            #'based on the quality control sample based robust LOESS (locally estimated scatterplot smoothing) signal correction (QC-RLSC) method as described by Dunn et al. (2011) and impletemented statTarget (Luan 2017).'
  NORMALIZE_NOT <- 'No QC normalisation'                #'if you do not want to perform normalisation, eg. when only coeficient of variance <30% instead'

  #OPLS-DA: S-plot option:
  FIXED_CORRELATION_COEFFICIENT <- 'S-plot: selection correlation (Y-axis) set at fixed value eg. +/-0.5. Type below which fixed value'
  CRITICAL_CORRELATION_COEFFICIENT <- 'S-plot: selection correlation (Y-axis) automatically calculated based on 95% distribution'

  #PCA option:
  DEFAULT_PC_AMOUNT <- "calculates PCs depending on sample amount (n-1)"
  FIXED_PC_AMOUNT <- "set amount of PCs to calculate to certain integer, ! less than amount of samples present in projections"



  ## Adjustments
  #If you choose 'VARIABLEMETADATA_EXTERN', add additional info here:
  VARIABLEMETADATA_EXTERN <- 'EXAMPLE_variableMetadata.txt'   #'name.txt' of file. Ignore if file created from pipeline
  COLLUMN_NR_START_SAMPLES <- 20  #always 20 (auto and manual must be same format); unless extra col merged!

  #Source of variableMetadata:
  INPUT_VARIABLES <- VARIABLEMETADATA_EXTERN

  #input file2 sampleMetadata:
  INPUT_SAMPLES <- 'EXAMPLE_sampleMetadata.txt' #don't forget .txt
  COLLUMN_ORDER <- 6
  COLLUMN_NR_TYPE <- 7          #column number of the column 'Type'
  ORDER_NR_OF_FIRST_QC <- 1
  COLLUMN_NR_LAST_BEFORE_COMPARISONS <- 9
  AMOUNT_OF_COMPARISONS <- 1      #amount of pairwise comparisons, =0 if none
  AMOUNT_OF_MULTIPLE_COMPARISONS <- 1         #amount of multiple comparisons, =0 if none present
  AMOUNT_OF_PROJECTIONS <- 1     #amount of PCAs (or alternative projection methods in future), =0 if none present


  #choose method of pretreatment of dataset:
  HYPERPARAMETERS_SELECTION <- HYPERPARAMETERS_AUTO
  HYPERPARAMETERS_SELECTION <- HYPERPARAMETERS_MANUAL

  #If you choose 'HYPERPARAMETERS_MANUAL', add additional info here:
  if(HYPERPARAMETERS_SELECTION == HYPERPARAMETERS_MANUAL){
    #variable filtering method:
    COEFICIENT_OF_VARIANCE <- COEFICIENT_OF_VARIANCE_0
    QC_PREFILTER <- NO_FILTER_QC_POOL

    #normalisation method:
    NORMALIZE_METHOD0 <- NORMALIZE_NOT0
    NORMALIZE_METHOD1 <- NORMALIZE_NOT1
    NORMALIZE_METHOD <- NORMALIZE_WITH_IQCs   #check data contains correct annotation: QC or IQC !!

    #in case NORMALIZE_WITH_IS is chosen: add info about IS
    IS_MZ <- 554.2615 #give float of m/z of IS eg. 554.2615 for Leu-Enk in negative mode [M-H]- OR 556.2771 for Leu-Enk in positive mode [M+H]+
    PPM <- 100 #eg. 100ppm reims, check in VM
    #! please check POLARITY so you don't search mz in oppsite ionisation mode

    #in case NORMALIZE_QC_RLSC is chosen, set span width:
    SPAN <- 0.75 #default 0.75, if error during NORMALIZE_QC_RLSC: set to 1 (or test until best span)
  }

  #If you choose 'HYPERPARAMETERS_AUTO', which param do you want to include:
  if(HYPERPARAMETERS_SELECTION == HYPERPARAMETERS_AUTO){
    #syntax: c(xxx, xxxx, xxxx) or c(xxxx) with options from above
    #variable filtering method:
    COEFICIENT_OF_VARIANCE <- c(COEFICIENT_OF_VARIANCE_0)
    QC_PREFILTER <- c(NO_FILTER_QC_POOL)

    #normalisation method:
    NORMALIZE_METHOD0 <- c(NORMALIZE_NOT0)
    NORMALIZE_METHOD1 <- c(NORMALIZE_NOT1, NORMALIZE_WITH_TIC)
    NORMALIZE_METHOD <- c(NORMALIZE_NOT, NORMALIZE_WITH_IQCs, NORMALIZE_QC_RLSC)
    #remark: if no QC/IQC present, changed to NORMALIZE_NOT before start loops
    #reamrk: reason order start w not: when tie, choose first in order combinr's. So choose combi w the least amount of transform as possible
    #remark: calculation time increases expon with extra option

    #in case NORMALIZE_WITH_IS is chosen: add info about IS
    IS_MZ <- 554.2615 #give float of m/z of IS eg. 554.2615 for Leu-Enk in negative mode [M-H]- OR 556.2771 for Leu-Enk in positive mode [M+H]+
    PPM <- 100 #eg. 100ppm reims, check in VM
    #! please check POLARITY so you don't search mz in oppsite ionisation mode

    #in case NORMALIZE_QC_RLSC is chosen, set span width:
    SPAN <- 0.75 #default 0.75, if error during NORMALIZE_QC_RLSC: set to 1 (or test until best span)
  }

  #OPLS-DA setting:
  #calculation y-axis S-plot of OPLS-DA:
  FIXED_CORRELATION_COEFFICIENT <- 0.5 #type absolute value of cutoff
  CUTOFF_CORRELATION <- FIXED_CORRELATION_COEFFICIENT

  #force code to make OPLS-DA model:
  ORT <- NA    #default should be NA, integers eg. 1

  #force code to make PLS-DA model:
  PRE <- NA   #default should be NA, integers eg. 1

  #PCA settings:
  #PCA amount of PCs in projection:
  FIXED_PC_AMOUNT <- 20   #integer, less than amount of samples present in projections
  PC_AMOUNT <- FIXED_PC_AMOUNT

  #style
  SIZE_POINTS <- 1   #size samplepoints on scoreplots PCA, (O)PLS, depending on amount of variables. default: 5, big dataset: 3

}
#
#####################
