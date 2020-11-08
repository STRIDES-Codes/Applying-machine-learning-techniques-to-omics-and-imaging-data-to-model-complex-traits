rm(list=ls())
library("tidyverse")
library("devtools")
# install_github("isb-cgc/examples-R", build_vignettes=TRUE)
# library(ISBCGCExamples)
# 
# bigrquery::list_tables("isb-cgc", "tcga_201510_alpha")
# BiocManager::install("TCGAbiolinks")
library("TCGAbiolinks")
query.met <- GDCquery(project = "TCGA-GBM",
                      data.category = "DNA Methylation",
                      legacy = FALSE)
# GDCdownload(query.met)
# GBM.met <- GDCprepare(query = query.met,
#                     save = TRUE, 
#                     save.filename = "GBM.met.rda")

query.exp <- GDCquery(project = "TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification")
# GDCdownload(query.exp)
# GBM.exp <- GDCprepare(query = query.met,
#                       save = TRUE, 
#                       save.filename = "GBM_exp.rda")
query.CNV <- GDCquery(project = "TCGA-GBM",
                      data.category = "Copy Number Variation",
                      data.type = "Copy Number Segment")
# GDCdownload(query.CNV)
# GBM.CNV <- GDCprepare(query = query.met,
#                       save = TRUE, 
#                       save.filename = "GBM_CNV.rda")
query.SlideImage <- GDCquery(project = "TCGA-GBM",
                      data.category = "Biospecimen",
                      data.type = "Slide Image")
# GDCdownload(query.SlideImage)
# GBM.SlideImage <- GDCprepare(query = query.SlideImage,
#                       save = TRUE, 
#                       save.filename = "GBM_SlideImage.rda")

# Get all patients that have DNA methylation and gene expression.
common.patients <- intersect(substr(getResults(query.met, cols = "cases"), 1, 12),
                             substr(getResults(query.exp, cols = "cases"), 1, 12))
common.patients <- intersect(common.patients,substr(getResults(query.CNV, cols = "cases"), 1, 12)) %>% 
  intersect(.,substr(getResults(query.SlideImage, cols = "cases"), 1, 12))

mytab <-  getSampleFilesSummary(project = "TCGA-GBM")
# filter by common.patients
mytab <- mytab %>% 
       filter(.id%in%common.patients) %>% 
       select(colnames(mytab)[-grep(pattern = "Sequencing|Simple|Transcriptome Profiling_Isoform Expression Quantification_miRNA-Seq|Transcriptome Profiling_miRNA Expression Quantification_miRNA-Seq",x = colnames(mytab))])
View(mytab)
# exclude patients with zero Diagnositc slide
mytab.Dslide0 <- mytab[which(mytab$`Biospecimen_Slide Image_Diagnostic Slide`!=0),]


table(PID%in%mytab$.id)
# FALSE  TRUE 
# 16    73 
table(PID%in%mytab.Dslide0$.id)
# FALSE  TRUE 
# 42    47 
`%!in%` = Negate(`%in%`)
PID[PID%!in%mytab$.id]

# prepare barcode
query.CNV$barcode

# Only seelct the first 5 patients
query.met <- GDCquery(project = "TCGA-GBM",
                      data.category = "DNA Methylation",
                      legacy = FALSE,
                      barcode = common.patients[1:5])
query.exp <- GDCquery(project = "TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ",
                      barcode = common.patients[1:5])
DT::datatable(getResults(query.met, cols = c("data_type","cases")),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
DT::datatable(getResults(query.exp, cols = c("data_type","cases")), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)





library("tidyverse")

PID <- "TCGA-32-2616, TCGA-19-1389, TCGA-19-1787, TCGA-19-0957, TCGA-19-1390, TCGA-14-0813, TCGA-19-2629, TCGA-27-1838, TCGA-32-2632, TCGA-06-0145, TCGA-06-0137, TCGA-19-2623, TCGA-06-2558, TCGA-16-1045, TCGA-32-1982, TCGA-14-1034, TCGA-28-2513, TCGA-12-3649, TCGA-12-0821, TCGA-06-0876, TCGA-14-2554, TCGA-14-0789, TCGA-16-0846, TCGA-32-1970, TCGA-06-2563, TCGA-06-0877, TCGA-14-0817, TCGA-27-1833, TCGA-06-0875, TCGA-41-3392, TCGA-14-0790, TCGA-19-2620, TCGA-12-1597, TCGA-12-3652, TCGA-06-0211, TCGA-06-0125, TCGA-06-0152, TCGA-28-2509, TCGA-32-2638, TCGA-06-2562, TCGA-14-0736, TCGA-06-0879, TCGA-06-2567, TCGA-06-2559, TCGA-02-2485, TCGA-06-0210, TCGA-06-0190, TCGA-27-2521, TCGA-06-0124, TCGA-02-2486, TCGA-19-2625, TCGA-28-2514, TCGA-27-1834, TCGA-41-2571, TCGA-19-2619, TCGA-14-1823, TCGA-14-1825, TCGA-27-2523, TCGA-06-0171, TCGA-19-4065, TCGA-27-1831, TCGA-14-1829, TCGA-14-0871, TCGA-27-2524, TCGA-12-3650, TCGA-06-0878, TCGA-32-2634, TCGA-27-2528, TCGA-27-2526, TCGA-28-1753, TCGA-41-2573, TCGA-27-1830, TCGA-28-1747, TCGA-27-1836, TCGA-14-0787, TCGA-27-1832, TCGA-06-0129, TCGA-19-2624, TCGA-27-2519, TCGA-27-1837, TCGA-06-0882, TCGA-06-2561, TCGA-06-2569, TCGA-06-0221, TCGA-06-0130, TCGA-12-3653, TCGA-06-0141, TCGA-06-0881, TCGA-06-0139"
PID<- PID %>% strsplit(split = ",") %>% unlist
PID <- gsub(pattern = " ",replacement = "",x = PID)
mydata <- read.csv(file = "all_dataset.csv",stringsAsFactors = F)

(mydata$TCGA.ID%in%PID)
PID[PID%in%mydata$TCGA.ID]
table(PID%in%mydata$TCGA.ID)
# selected data
s.mydata <- mydata %>% filter(TCGA.ID%in%PID)