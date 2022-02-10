# nasopharyngeal-carcinoma
M6A-mediated molecular patterns and tumor microenvironment infiltration characterization in nasopharyngeal carcinoma

Data downloading and processing 
mRNA expression data of NPC patients of two datasets, GSE102349, and GSE53815 datasets, were obtained from Gene Expression Omnibus (GEO, https://www.ncbi.nlm.nih.gov/geo/). 
The mRNA expression data of the GSE102349 dataset was generated using Illumina HiSeq 2000, this dataset includes mRNA expression data and clinical information from 113 NPC 
patients, and 88 cases with progression-free survival (PFS) involved in this study.The GSE102349 dataset was to re-run reads quality control, alignment, and quantitative gene expression analysis.

Univariate Cox analysis was performed with survival R package;
The consensus clustering was performed using the consensusClusterPlus package in R ;
GSVA using GSVA package in R to investigate the various biological function among distinct m6A subclasses;
ESTIMATE algorithm to estimate the composition of the tumor microenvironment;
the chemotherapy response rate using pRRophetic package in R. According to the Genomics of Drug Sensitivity in Cancer (GDSC) database (https://www.cancerrxgene.org/) 
The differentially expressed genes (DEGs) among NPC subclasses were identified using edgeR package with parameters of |log (fold change, FC) | > 1 and FDR < 0.05;
the m6A-related genes related to prognosis were identified using LASSO regression analysis using glmet package in R;
the time-dependent receiver operating characteristic (ROC) curve was constructed using the survivalROC R package;
A nomogram predicted model was constructed based on the risk factors using rms package in R;
