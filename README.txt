In this repository you will find all code necessary to reproduce the analyses and figures presented in 'Advancing AOPs for genotoxicity through quantifying KERs linking DNA alkylation to increase in mutations and chromosome aberrations' by Winkelman et al., 2026. 

**rawData** 					Contains all 'raw' data including extracted concentration response data and overviews of the manuscripts retrieved through systematic literature search

**processedData**				Contains processed data, generally meant to be used as input for subsequent scripts. This includes:
							-wrangled (postQC) data
							-data subsets based on experimental design
							-stored output from bootstrap analysis
							-overviews of data quantity per assay and compound 
							-overview of extracted manuscripts

**Scripts**					Contains:
							-.R files storing custom functions along with documentation on how to use those functions 
							-.Rmd files in which the actual analysis is performed (each Rmd starts with a short description of the purpose of that scripts)

**Figures**
	**BMD**
		**bias_summaries**		Shows the distribution of means and SEMs per concentration after bootstrap sampling and compares them to the values of the original data.
		**BMD_comparison**		Shows side by side BMDs along with CIs for different assays 
		**BMD_preparation**		Shows how data merging and log-transformation was handled
		**bootplots**			Shows BMD model fit on (log-transformed) original data along with BMD model fits on each of several bootstrap datasets
		**fitgrids**			Shows fits of 3 parameter model and 4 parameter model using 4 model families using log-transformed data	
	**ConcRes_comparison**			Shows comparisons of concentration-response curves for different studies while stratifying for experimental design
	**PieCharts_litSearch**			Shows the distribution of manuscripts retrieved from PubMed vs PubMed central and manuscripts for which data was / was not extracted
	**Response-Response**			Shows fitted response-response curves
	**VarianceDecomposition**	
		**Fig3				Shows the outcome of variance decomposition
		**FigS2				Shows the outcome of variance component analysis

-------------------------------------------------------------------------------------------
Overview of what script to run to reproduce figures and tables presented in the manuscript:
-------------------------------------------------------------------------------------------

		Derived from script		(sub)figure / Table stored under
Fig. 1		1_Examine_extracted_data	Figures/PieChartslitSearch/Fig1
Fig. 2		4_Compare_ConcRes		Figures/ConcRes_comparison/Fig2
Fig. 3		5_VarianceDecomposition		Figures/VarianceDecomposition/Fig3
Fig. 4		6_BMD_full_analysis		Figures/BMD_comparison/Fig4
Fig. 5		8_Response_Response		Figures/Response-Response/Fig5
Fig. S1		2_Data_wrangling_EMS/MMS/ENU	Figures/ConcRes_comparison/FigS1
Fig. S2		5_VarianceDecomposition		Figures/VarianceDecomposition/FigS2
Fig. S3		7_Visualize_BMD_pipeline	Figures/BMD/BMD_preparation/FigS3
Fig. S4		BMD_full_analysis.Rmd		Figures/BMD/BMD_comparison/FigS4
Fig. S5		Response_Response.Rmd		Figures/Response-Response/FigS5
Fig. S6		Response_Response.Rmd		Figures/Response-Response/FigS5
Table S1	-
Table S2	-
Table S3	1_Examine_extracted_data	processedData/S3.xlsx
Table S4	3_Count_studies			processedData/SummaryS4.rds
Table S5	3_Count_studies			processedData/SummaryS5.rds
Table S6	3_Count_studies			processedData/SummaryS6.rds
Table S7	3_Count_studies			processedData/SummaryS7.rds
Table S8	3_Count_studies			processedData/SummaryS8.rds
Table S9	3_Count_studies			processedData/SummaryS9.rds
Table S10	BMD_full_analysis.Rmd		processedData/BMD50ress.csv, processedData/BMDSDress.csv
Table S11	Response_Response.Rmd		
Table S12	-


Please note that some scripts use input that is generated from prior scripts, it is therefore advisable to run scripts in the numbered order in which they appear in the 'Scripts' folder.