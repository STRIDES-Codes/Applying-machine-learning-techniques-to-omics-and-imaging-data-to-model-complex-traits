# Machine Learning Technique for Omics and Imaging Data: An Application in TCGA-GBM
## Abstract  
Clinical designs in biomedical research increasingly require the implementation of high-throughput multiomics approaches to characterize patient samples (e.g. RNAseq, ATACseq and imaging). Therefore, there is an urgent need to develop data analysis pipelines that are able to integrate such data in an efficient and generalizable manner. The canonical correlation analysis (CCA), which is rooted in the Pearson r correlation, is a multivariate technique that limits the probability of committing a Type I error (e.g., finding a difference, effect, or relationship when it really does not exist ) and may improve on many commonly-used parametric GLM variations (e.g., ANOVA, MANOVA, multiple regression, Pearson correlation, t test, discriminant analysis). Thus, the core objective of our project is to provide a CCA-based data analysis framework that underpins an integrative multimodal data analysis strategy and improves on alternative methods. Here, we aim to validate our approach using a glioma dataset from the Cancer Genome Atlas (TCGA), which contains histological images, transcriptomic, and methylation status data with ground truth survival labels.

![alt text](https://github.com/STRIDES-Codes/Applying-machine-learning-techniques-to-omics-and-imaging-data-to-model-complex-traits/blob/main/Images/CSHL2020_readme_images_workflow.png?raw=true)

## Team Members  
- **Md Ashad Alam**
- **Soner Koc**
- **Chengzhe Tian**
- **KuanJui "Ray" Su**
- **Rebecca "Becky" Waugh**
- **Udana Torian**
- **Amin Cheikhi**


## Dependencies:computer:
Google Cloud Platform(GCP) VM  
[National Cancer Institute GDC Data Portal](https://portal.gdc.cancer.gov/)  
Rstudio
RGCCA, * 2.1.2
CCA, * 1.2.0
rcc, * 1.0.0
dplyr, * 1.0.2
TCGAbiolinks, * 2.14.1