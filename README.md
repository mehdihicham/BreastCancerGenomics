# BreastCancerGenomics

This is the script that is used to reproduce the results obtained in the paper : "Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer", by Ciriello et al. published in Cell in 2015.

Dataset can be downloaded in cBioPortal : https://www.cbioportal.org/study/summary?id=brca_tcga_pub2015 

ILC_subtype.R is the main script. It requires gene_length.txt, produced by gene_length.R. I really recommend you not to run the latter, as it takes time to generated. gene_lenght.txt is store in this github. 
In order to run the script, you have to modify the first line of the code, in order to change the path of the brca_tcga_pub2015 directory.
