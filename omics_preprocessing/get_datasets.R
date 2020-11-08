setwd('/Users/chengzhetian/Desktop/hackathon/');
library(RTCGA);
releaseDate <- "2016-01-28";
cohorts <- 'GBM';

downloadTCGA( cancerTypes = cohorts, 
              dataSet = "GBM.Merge_methylation__humanmethylation27",
              destDir = "data2", 
              date = releaseDate );
downloadTCGA( cancerTypes = cohorts, 
              dataSet = "GBM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz",
              destDir = "data2", 
              date = releaseDate );
