setwd('/Users/chengzhetian/Desktop/hackathon/');

# set valid id
valid_patient_id <- c('TCGA-19-1389', 'TCGA-19-1787', 'TCGA-19-0957', 'TCGA-19-1390', 'TCGA-14-0813', 'TCGA-19-2629', 'TCGA-27-1838', 'TCGA-32-2632', 'TCGA-06-0145', 'TCGA-06-0137', 'TCGA-19-2623', 'TCGA-06-2558', 'TCGA-16-1045', 'TCGA-32-1982', 'TCGA-14-1034', 'TCGA-28-2513', 'TCGA-12-3649', 'TCGA-12-0821', 'TCGA-06-0876', 'TCGA-14-2554', 'TCGA-14-0789', 'TCGA-16-0846', 'TCGA-32-1970', 'TCGA-06-2563', 'TCGA-06-0877', 'TCGA-14-0817', 'TCGA-27-1833', 'TCGA-06-0875', 'TCGA-41-3392', 'TCGA-14-0790', 'TCGA-19-2620', 'TCGA-12-1597', 'TCGA-12-3652', 'TCGA-06-0211', 'TCGA-06-0125', 'TCGA-06-0152', 'TCGA-28-2509', 'TCGA-32-2638', 'TCGA-06-2562', 'TCGA-14-0736', 'TCGA-06-0879', 'TCGA-06-2567', 'TCGA-06-2559', 'TCGA-02-2485', 'TCGA-06-0210', 'TCGA-06-0190', 'TCGA-27-2521', 'TCGA-06-0124', 'TCGA-02-2486', 'TCGA-19-2625', 'TCGA-28-2514', 'TCGA-27-1834', 'TCGA-41-2571', 'TCGA-19-2619', 'TCGA-14-1823', 'TCGA-14-1825', 'TCGA-27-2523', 'TCGA-06-0171', 'TCGA-19-4065', 'TCGA-27-1831', 'TCGA-14-1829', 'TCGA-14-0871', 'TCGA-27-2524', 'TCGA-12-3650', 'TCGA-06-0878', 'TCGA-32-2634', 'TCGA-27-2528', 'TCGA-27-2526', 'TCGA-28-1753', 'TCGA-41-2573', 'TCGA-27-1830', 'TCGA-28-1747', 'TCGA-27-1836', 'TCGA-14-0787', 'TCGA-27-1832', 'TCGA-06-0129', 'TCGA-19-2624', 'TCGA-27-2519', 'TCGA-27-1837', 'TCGA-06-0882', 'TCGA-06-2561', 'TCGA-06-2569', 'TCGA-06-0221', 'TCGA-06-0130', 'TCGA-12-3653', 'TCGA-06-0141', 'TCGA-06-0881', 'TCGA-06-0139');
valid_gene_id <- as.character(read.csv('Multiomicsgenepanel.csv')$Symbol);

# get mRNA data
mRNA_dataset <- read.table('data2/RNAseq.txt', sep = '\t'); 
mRNA_dataset <- mRNA_dataset[-2, ];
mRNA_col_entry <- sapply(valid_patient_id, function(x) {
  temp <- sapply(1:ncol(mRNA_dataset), function(y) {
    startsWith(as.character(mRNA_dataset[1, y]), x);
  });
  return(which(temp)[1]);
});
mRNA_row_entry <- sapply(valid_gene_id, function(x) {
  return(which(startsWith(as.character(mRNA_dataset[, 1]), paste(x, '|', sep='')))[1]);
});

# get methylation data
methyl_dataset <- read.table('data2/Methyl.txt', sep = '\t');
methyl_col_entry <- sapply(valid_patient_id, function(x) {
  temp <- sapply(1:ncol(methyl_dataset), function(y) {
    startsWith(as.character(methyl_dataset[1, y]), x);
  });
  return(which(temp)[1]);
});
methyl_row_entry <- sapply(valid_gene_id, function(x) {
  temp1 <- as.character(methyl_dataset[, 3]);
  temp2 <- (temp1 == x) | 
    startsWith(temp1, paste(x, ';', sep='')) | 
    endsWith(temp1, paste(';', x, sep='')) |
    grepl(paste(';', x, ';', sep=''), temp1);
  return(which(temp2)[1]);
});

# finalize everything, part 1
row_to_keep <- !is.na(mRNA_row_entry) & !is.na(methyl_row_entry);
col_to_keep <- !is.na(mRNA_col_entry) & !is.na(methyl_col_entry);

mRNA_finalize <- mRNA_dataset[mRNA_row_entry[row_to_keep], mRNA_col_entry[col_to_keep]];
rownames(mRNA_finalize) <- valid_gene_id[row_to_keep];
colnames(mRNA_finalize) <- valid_patient_id[col_to_keep];

methyl_finalize <- methyl_dataset[methyl_row_entry[row_to_keep], methyl_col_entry[col_to_keep]];
rownames(methyl_finalize) <- valid_gene_id[row_to_keep];
colnames(methyl_finalize) <- valid_patient_id[col_to_keep];

# finalize everything, part 2 (found NA values)
row_to_keep2 <- apply(mRNA_finalize, 1, function(x) sum(is.na(x)) == 0) &
  apply(methyl_finalize, 1, function(x) sum(is.na(x)) == 0); 
mRNA_finalize <- mRNA_finalize[row_to_keep2, ];
methyl_finalize <- methyl_finalize[row_to_keep2, ];
write.table(mRNA_finalize, 'mRNA_selected.csv', sep = ',', quote = FALSE);
write.table(methyl_finalize, 'methyl_selected.csv', sep = ',', quote = FALSE);

# PCA
temp <- t(mRNA_finalize); class(temp) <- 'numeric'; temp <- t(temp);
mRNA_pca <- prcomp(temp, center=TRUE, scale.=TRUE);
write.table(mRNA_pca$rotation[, 1:10], 'mRNA_PCA.csv', sep = ',', quote = FALSE);

temp <- t(methyl_finalize); class(temp) <- 'numeric'; temp <- t(temp);
methyl_pca <- prcomp(temp, center=TRUE, scale.=TRUE);
write.table(methyl_pca$rotation[, 1:10], 'methyl_PCA.csv', sep = ',', quote = FALSE);

