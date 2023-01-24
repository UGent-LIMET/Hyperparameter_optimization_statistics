Info:
*****
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: READ_ME

get started
***********

### requirements files ###
"variableMetadata" = output peak calling, is made automatically in part 1: pre-processing.
Alternatively, this can be imported directly here if the format is respected.
For this, the following rules apply:
File1 = "variableMetadata" = output peak calling (~sieve/CD)
	- need to be exact same colnames (order/name): 
		"CompID"	"MZ"	"Time"	"isotopes"	"adduct"	"pcgroup"	"name"	"fold"	"tstat"	"pvalue"	"mzmed"	"mzmin"	"mzmax"	"rtmed"	"rtmin"	"rtmax"	"npeaks"	"bio"	"blank"	
	- at the end of columns the intenstities from each unique samplename in a column eg.
		"19007s05"	"19007s15"	"19007s16"	"19007s70"	"19007s01"	"19007s06"
	- intensities are positive numbers, in case of targeted analysis replace < LOD, NF, ... by "0"
	- of which the solvent/blanks samples need to be the last (two) columns, needed for blank pre-filter in targeted analyse
	- All columns need to be present
	- CompIDs are integers eg. 1, 2, 3
	- do not use "#" in the file, because R recognizes this as 2 columns instead of 1	(eg. #Name for wrong function)
	- optional: add "MetaboliteName" column (columnnr set to 21 instead of 20) with names of standards (targeted analysis)

Needed for Part2 Statistical analysis:	
File2 = "sampleMetadata" = input sample information 
	- unique column names
	- SampleName = all sample names corresponding with the variableMetadata
		need to be the first column
		need to be unique names
		used as output on plots (pca, heatmap) so choose short name
		must begin with number (006 instead of S006)
		not spaces in name (R adds "." so does not recognize samplename with names from variableMetadata)
		same amount as samples in variableMetadata
		name can be same as sieve output,  eg. II_sample001 (including blank, QCs, ...)
	- all samples need to be present that were imported into R pipeline (same SampleNames as file.mzXML, same amount)
	- optional: Name = short easy samplename which you want displayed on heatmap
	- Batch = {1, 2, 3, ...} number which run samples were analysed (1 for only one analysis day)
	- Order = {1, 2, 3, ...} number according to runorder on machine (smallest number is 1, no skipped numbers); when using LOESS method: 1st QC = 1 (remove all before) todo in code!!
	- Type = {Sample = samples of interest; Blank = solvent control; QC = quality control OR IQC = internal quality control; STD = internal standard; ISTD = ?} ! no other syntax
	- other features may be added eg. group (before comparisons!)
	  for example:
		- Omics = {Metabolomics, Lipidomics, ...}
		- IonisationMode = {+, -, Both}
	- Comparison1, for the pairwise comparison between two setups with integer values {-1, 1 or empty value} as input cells
	- Comparison2, ... as many as needed, each in one column
	- MultipleComparison1, for the multiple comparison between more then two setups with integer values {1, 2, 3,... or empty cell => NA} as input cells
	- MultipleComparison2, ... as many as needed, each in one column
	- Projection1, for creation of labeled PCA (or other projection methods in future)
	  input are integer values {1, 2, 3,... or empty cell => NA}
	- Projection2, ... as many as needed, each PCA in one column
	- the last projection is the last column in table, no other column added here
	- do not use "#" in the file, because R recognizes this as 2 collumns instead of 1

- convert files to .txt file extension with tab as separator
- Paste files in Data/Input folder

!when exporting from excel: separator = dec (no 7815,5) and numbers are defined as 'number' (no 1,000,000)
!close excel if opened



### Adjust settings for your analysis ###
"configuration file" = all settings for the code to run.
- copy Configuration(reserve, see input).R into you input
- rename to Configuration.R
- Open Configuration.R 
- adjust settings under 'Configuration > Adjustments' according to the options given under 'Configuration > Options'
- projectname is structured and as short as possible: YYYYMMDD_initials_description_ionisationMode
  eg. 20191031_MDG_exp1_pos
  seperate with "_", don't use symbols (like +,-), spaces, ... 
  description part is optional, max 1 word
- polarity is needed for the pre-processing modules, options are "positive", "negative"
- user_comment is an optional, short comment that will be displayed on the first page of report; eg  explain (Multiple)Comparisons
- when only part of the pipeline is run: adjust the global settings (eg. projectname) and the selected parts
- Save Configuration.R
- close Configuration.R
