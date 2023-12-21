#@author:chenbin
#@created time:2023/7/14
#@modified time: 2023/12/21
#@run time: 2023/12/21
#@comment: Script to calculate the PVE and draw the density plot of PVE.
#@计算每个eQTL对基因表达的解释方差R2
library(data.table)
library(gghalves)
library(ggplot2)

workdir <- "/home/chenb/Project/Plants/Wheat/nwafu_triticum_aestivum_reseq_20211113/05_Mutation/s6.chrs_part2whole.gatk/eGWAS/eQTL_summary"
setwd(workdir)
##################################################
##### set the input/output files' locations
outdir = "/home/chenb/Project/Plants/Wheat/nwafu_triticum_aestivum_reseq_20211113/05_Mutation/s6.chrs_part2whole.gatk/eGWAS/eQTL_summary"
outprefix = "WMG199_CSDV418_eQTL.uniq.QTLname.update_R2"
geno_infile = "/home/chenb/Project/Plants/Wheat/nwafu_triticum_aestivum_reseq_20211113/05_Mutation/s6.chrs_part2whole.gatk/WMG199_CSDV418_snp.biallelic.miss_maf.filter.id.CS.GT.VCFconverter.MatrixEQTL.vcf"
WW_gene_expr_infile = "/home/chenb/Project/Plants/Wheat/nwafu_triticum_aestivum_reseq_20211113/06_Expression/02_Expression_for_eQTL/PEER_estimation/gene.exp.wmg.WW.tpm0.1.Q5_95th.log2_qqnorm.filter.PEER_residuals.transpose.xls"
DS_gene_expr_infile = "/home/chenb/Project/Plants/Wheat/nwafu_triticum_aestivum_reseq_20211113/06_Expression/02_Expression_for_eQTL/PEER_estimation/gene.exp.wmg.DS.tpm0.1.Q5_95th.log2_qqnorm.filter.PEER_residuals.transpose.xls"
WW_eQTL.uniq.QTLname_infile = "/home/chenb/Project/Plants/Wheat/nwafu_triticum_aestivum_reseq_20211113/05_Mutation/s6.chrs_part2whole.gatk/eGWAS/normal_eQTL_WMG199_CSDV418/eQTL_analysis_by_MatrixEQTL_PEER/add_QTL_name/WMG199_CSDV418_WW_eQTL.uniq.QTLname.xls"
DS_eQTL.uniq.QTLname_infile = "/home/chenb/Project/Plants/Wheat/nwafu_triticum_aestivum_reseq_20211113/05_Mutation/s6.chrs_part2whole.gatk/eGWAS/drought_eQTL_WMG199_CSDV418/eQTL_analysis_by_MatrixEQTL_PEER/add_QTL_name/WMG199_CSDV418_DS_eQTL.uniq.QTLname.xls"
pca3_infile = "/home/chenb/Project/Plants/Wheat/nwafu_triticum_aestivum_reseq_20211113/05_Mutation/s6.chrs_part2whole.gatk/Population_genetics/WMG199_CSDV418_snp.biallelic.miss_maf.filter.id.CS.recodechrid.gcta.MatrixEQTL.pc3.xls"

##### Read in the data
## the genotype data
geno_dat <- fread(file = geno_infile, header = T, sep = "\t")
dim(geno_dat)
geno_dat[1:10, 1:10]

## the eQTLs data which have been assigned QTL names
# WW
WW_eQTL.uniq.QTLname_dat <- fread(file = WW_eQTL.uniq.QTLname_infile, header = T, sep = "\t")
dim(WW_eQTL.uniq.QTLname_dat)
WW_eQTL.uniq.QTLname_dat
# DS
DS_eQTL.uniq.QTLname_dat <- fread(file = DS_eQTL.uniq.QTLname_infile, header = T, sep = "\t")
dim(DS_eQTL.uniq.QTLname_dat)
DS_eQTL.uniq.QTLname_dat

## the gene expression data
# WW
WW_gene_expr_dat <- fread(file = WW_gene_expr_infile, header = T, sep = "\t")
names(WW_gene_expr_dat)[1] <- "ID"
dim(WW_gene_expr_dat)
WW_gene_expr_dat[1:10, 1:10]
# DS
DS_gene_expr_dat <- fread(file = DS_gene_expr_infile, header = T, sep = "\t")
names(DS_gene_expr_dat)[1] <- "ID"
dim(DS_gene_expr_dat)
DS_gene_expr_dat[1:10, 1:10]

## The top three genotyping principle components
covariate_pca3 = read.table(file = pca3_infile, header = T, row.names = 1, sep = '\t')
covariate_pca3[, 1:10]
covariate_pca3_transpose <- transpose(covariate_pca3)
head(covariate_pca3_transpose)
dim(covariate_pca3_transpose)
covariate_pca3_transpose_mat <- as.matrix(covariate_pca3_transpose)

##### Get the unique lead SNPs' genotypes across 200 accessions under WW and DS conditions
## the unique SNPs from WW and DS
WW_DS_unique_SNPs_ID <- union(WW_eQTL.uniq.QTLname_dat[, lead_snp], DS_eQTL.uniq.QTLname_dat[, lead_snp])
WW_DS_unique_SNPs_ID[1:20]
length(WW_DS_unique_SNPs_ID)
## the genotype data only containing  the unique SNPs
uniq_snpid_geno_dat <- geno_dat[which(ID %in% WW_DS_unique_SNPs_ID)]
dim(uniq_snpid_geno_dat)
uniq_snpid_geno_dat[1:10, 1:10]
## transpose
uniq_snpid_geno_dat_transpose <- data.table::transpose(uniq_snpid_geno_dat, make.names="ID", keep.names="ID")
dim(uniq_snpid_geno_dat_transpose)
uniq_snpid_geno_dat_transpose[1:10, 1:10]

##### Get the target genes expression PEER residuals under WW and DS conditions, respectively
### WW
## the unique Gene IDs from WW
WW_unique_geneid <- unique(WW_eQTL.uniq.QTLname_dat[, gene_id])
cat('The number of unique target Genes associated with markers:', length(WW_unique_geneid), '\n')
## the gene expression data from WW only containing the unique Gene IDs
WW_unique_gene_expr_dat <- WW_gene_expr_dat[which(ID %in% WW_unique_geneid)]
dim(WW_unique_gene_expr_dat)
WW_unique_gene_expr_dat[1:10, 1:10]
## transpose
WW_unique_gene_expr_dat_transpose <- data.table::transpose(WW_unique_gene_expr_dat, make.names="ID", keep.names="ID")
dim(WW_unique_gene_expr_dat_transpose)
WW_unique_gene_expr_dat_transpose[1:10, 1:8]

### DS
## the unique Gene IDs from DS
DS_unique_geneid <- unique(DS_eQTL.uniq.QTLname_dat[, gene_id])
cat('The number of unique target Genes associated with markers:', length(DS_unique_geneid), '\n')
## the gene expression data from DS only containing the unique Gene IDs
DS_unique_gene_expr_dat <- DS_gene_expr_dat[which(ID %in% DS_unique_geneid)]
dim(DS_unique_gene_expr_dat)
DS_unique_gene_expr_dat[1:10, 1:10]
## transpose
DS_unique_gene_expr_dat_transpose <- data.table::transpose(DS_unique_gene_expr_dat, make.names="ID", keep.names="ID")
dim(DS_unique_gene_expr_dat_transpose)
DS_unique_gene_expr_dat_transpose[1:10, 1:8]

##### Calculate the R2 using a two-step regression (full model & reduced model)
##### calculate the PVE for eQTLs from WW
## initialization
WW_eQTL.uniq.QTLname_dat[, SST := numeric()]
WW_eQTL.uniq.QTLname_dat[, R2_partial := numeric()]
WW_eQTL.uniq.QTLname_dat_df <- as.data.frame(WW_eQTL.uniq.QTLname_dat)
dim(WW_eQTL.uniq.QTLname_dat_df)
head(WW_eQTL.uniq.QTLname_dat_df)
###
SS2 = sum( (covariate_pca3_transpose[, 1] - mean(covariate_pca3_transpose[, 1]))^2 )
SS3 = sum( (covariate_pca3_transpose[, 2] - mean(covariate_pca3_transpose[, 2]))^2 )
SS4 = sum( (covariate_pca3_transpose[, 3] - mean(covariate_pca3_transpose[, 3]))^2 )
SP23 = sum( (covariate_pca3_transpose[, 1] - mean(covariate_pca3_transpose[, 1])) * (covariate_pca3_transpose[, 2] - mean(covariate_pca3_transpose[, 2])) )
SP24 = sum( (covariate_pca3_transpose[, 1] - mean(covariate_pca3_transpose[, 1])) * (covariate_pca3_transpose[, 3] - mean(covariate_pca3_transpose[, 3])) )
SP34 = sum( (covariate_pca3_transpose[, 2] - mean(covariate_pca3_transpose[, 2])) * (covariate_pca3_transpose[, 3] - mean(covariate_pca3_transpose[, 3])) )
for (i in 1:nrow(WW_eQTL.uniq.QTLname_dat_df)){
  geneid = WW_eQTL.uniq.QTLname_dat_df[i, "gene_id"]
  snpid = WW_eQTL.uniq.QTLname_dat_df[i, "lead_snp"]
  lead_snp_beta = WW_eQTL.uniq.QTLname_dat_df[i, "lead_snp_beta"]
  datx = data.frame(snpidx = uniq_snpid_geno_dat_transpose[, ..snpid], geneidx = WW_unique_gene_expr_dat_transpose[, ..geneid])
  # modelx <- lm(datx[,2] ~ datx[,1] + covariate_pca3_transpose_mat, na.action = na.omit)
  # summary(modelx)
  
  ### Calculate the inverse matrix of the coefficient matrix within the normal equation group
  SS1 = sum( (datx[,1] - mean(datx[,1]))^2 )
  SP12 = sum( (datx[,1] - mean(datx[,1])) * (covariate_pca3_transpose[, 1] - mean(covariate_pca3_transpose[, 1])) )
  SP13 = sum( (datx[,1] - mean(datx[,1])) * (covariate_pca3_transpose[, 2] - mean(covariate_pca3_transpose[, 2])) )
  SP14 = sum( (datx[,1] - mean(datx[,1])) * (covariate_pca3_transpose[, 3] - mean(covariate_pca3_transpose[, 3])) )
  # the coefficient matrix
  coef_matrix = as.matrix(data.frame(X1 = c(SS1, SP12, SP13, SP14),
                                     X2 = c(SP12, SS2, SP23, SP24),
                                     X3 = c(SP13, SP23, SS3, SP34), 
                                     X4 = c(SP14, SP24, SP34, SS4)))
  # the inverse matrix of the coefficient matrix
  rev_coef_matrix = solve(coef_matrix)
  ###
  SST = sum((datx[,2]-mean(datx[,2]))^2)
  U1 = (lead_snp_beta)^2 / rev_coef_matrix[1,1]
  R2_partial = U1 / SST
  WW_eQTL.uniq.QTLname_dat_df[i, "SST"] = SST
  WW_eQTL.uniq.QTLname_dat_df[i, "R2_partial"] = R2_partial
}

#check
WW_eQTL.uniq.QTLname_dat_df[c(1:5, 100:105, 1000:1005),]

##### calculate the PVE for eQTLs from DS
## initialization
DS_eQTL.uniq.QTLname_dat[, SST := numeric()]
DS_eQTL.uniq.QTLname_dat[, R2_partial := numeric()]
DS_eQTL.uniq.QTLname_dat_df <- as.data.frame(DS_eQTL.uniq.QTLname_dat)
dim(DS_eQTL.uniq.QTLname_dat_df)
head(DS_eQTL.uniq.QTLname_dat_df)
###
SS2 = sum( (covariate_pca3_transpose[, 1] - mean(covariate_pca3_transpose[, 1]))^2 )
SS3 = sum( (covariate_pca3_transpose[, 2] - mean(covariate_pca3_transpose[, 2]))^2 )
SS4 = sum( (covariate_pca3_transpose[, 3] - mean(covariate_pca3_transpose[, 3]))^2 )
SP23 = sum( (covariate_pca3_transpose[, 1] - mean(covariate_pca3_transpose[, 1])) * (covariate_pca3_transpose[, 2] - mean(covariate_pca3_transpose[, 2])) )
SP24 = sum( (covariate_pca3_transpose[, 1] - mean(covariate_pca3_transpose[, 1])) * (covariate_pca3_transpose[, 3] - mean(covariate_pca3_transpose[, 3])) )
SP34 = sum( (covariate_pca3_transpose[, 2] - mean(covariate_pca3_transpose[, 2])) * (covariate_pca3_transpose[, 3] - mean(covariate_pca3_transpose[, 3])) )
for (i in 1:nrow(DS_eQTL.uniq.QTLname_dat_df)){
  geneid = DS_eQTL.uniq.QTLname_dat_df[i, "gene_id"]
  snpid = DS_eQTL.uniq.QTLname_dat_df[i, "lead_snp"]
  lead_snp_beta = DS_eQTL.uniq.QTLname_dat_df[i, "lead_snp_beta"]
  datx = data.frame(snpidx = uniq_snpid_geno_dat_transpose[, ..snpid], geneidx = DS_unique_gene_expr_dat_transpose[, ..geneid])
  # modelx <- lm(datx[,2] ~ datx[,1] + covariate_pca3_transpose_mat, na.action = na.omit)
  # summary(modelx)
  
  ### Calculate the inverse matrix of the coefficient matrix within the normal equation group
  SS1 = sum( (datx[,1] - mean(datx[,1]))^2 )
  SP12 = sum( (datx[,1] - mean(datx[,1])) * (covariate_pca3_transpose[, 1] - mean(covariate_pca3_transpose[, 1])) )
  SP13 = sum( (datx[,1] - mean(datx[,1])) * (covariate_pca3_transpose[, 2] - mean(covariate_pca3_transpose[, 2])) )
  SP14 = sum( (datx[,1] - mean(datx[,1])) * (covariate_pca3_transpose[, 3] - mean(covariate_pca3_transpose[, 3])) )
  # the coefficient matrix
  coef_matrix = as.matrix(data.frame(X1 = c(SS1, SP12, SP13, SP14),
                                     X2 = c(SP12, SS2, SP23, SP24),
                                     X3 = c(SP13, SP23, SS3, SP34), 
                                     X4 = c(SP14, SP24, SP34, SS4)))
  # the inverse matrix of the coefficient matrix
  rev_coef_matrix = solve(coef_matrix)
  ###
  SST = sum((datx[,2]-mean(datx[,2]))^2)
  U1 = (lead_snp_beta)^2 / rev_coef_matrix[1,1]
  R2_partial = U1 / SST
  DS_eQTL.uniq.QTLname_dat_df[i, "SST"] = SST
  DS_eQTL.uniq.QTLname_dat_df[i, "R2_partial"] = R2_partial
}

#check
DS_eQTL.uniq.QTLname_dat_df[c(1:5, 100:105, 1000:1005),]

##### write out
write.table(WW_eQTL.uniq.QTLname_dat_df, file =  file.path(outdir, paste0(outprefix, ".WW.xls")),
            quote = F, sep = '\t', row.names = F, col.names = T)
write.table(DS_eQTL.uniq.QTLname_dat_df, file =  file.path(outdir, paste0(outprefix, ".DS.xls")),
            quote = F, sep = '\t', row.names = F, col.names = T)


########## Draw the violin plot for comparison of the explained variance (R2_partial) of local and distant eQTLs for expression at WW and DS condtions
# Local
WW.local.PVE.pdat <- WW_eQTL.uniq.QTLname_dat_df[which(WW_eQTL.uniq.QTLname_dat_df$eqtl_type == "Local"), c("eqtl_type", "R2_partial")]
DS.local.PVE.pdat <- DS_eQTL.uniq.QTLname_dat_df[which(DS_eQTL.uniq.QTLname_dat_df$eqtl_type == "Local"), c("eqtl_type", "R2_partial")]
# Distant
WW.distant.PVE.pdat <- WW_eQTL.uniq.QTLname_dat_df[which(WW_eQTL.uniq.QTLname_dat_df$eqtl_type == "Distant"), c("eqtl_type", "R2_partial")]
DS.distant.PVE.pdat <- DS_eQTL.uniq.QTLname_dat_df[which(DS_eQTL.uniq.QTLname_dat_df$eqtl_type == "Distant"), c("eqtl_type", "R2_partial")]

### use geom_half_violin() to draw the violin plot
p <- ggplot()+
  geom_half_violin(data = WW.local.PVE.pdat,
                   aes(x = eqtl_type, y = R2_partial, fill = "WW"),
                   side = "l", trim = T, color = "grey40", alpha = 0.75,
                   draw_quantiles = c(0.25, 0.5, 0.75),
                   position = position_nudge(x = 0, y = 0))+
  geom_half_violin(data = DS.local.PVE.pdat,
                   aes(x = eqtl_type, y = R2_partial, fill = "DS"),
                   side = "r", trim = T, color = "grey40", alpha = 0.75,
                   draw_quantiles = c(0.25, 0.5, 0.75),
                   position = position_nudge(x = 0, y = 0))+
  geom_half_violin(data = WW.distant.PVE.pdat,
                   aes(x = eqtl_type, y = R2_partial, fill = "WW"),
                   side = "l", trim = T, color = "grey40", alpha = 0.75,
                   draw_quantiles = c(0.25, 0.5, 0.75),
                   position = position_nudge(x = 0, y = 0))+
  geom_half_violin(data = DS.distant.PVE.pdat,
                   aes(x = eqtl_type, y = R2_partial, fill = "DS"),
                   side = "r", trim = T, color = "grey40", alpha = 0.75,
                   draw_quantiles = c(0.25, 0.5, 0.75),
                   position = position_nudge(x = 0, y = 0))+
  scale_x_discrete(limits = c("Local", "Distant"))+
  scale_fill_manual(limits = c("WW", "DS"),
                    name = "Treatment",
                    values = c("Dark Turquoise", "LightCoral"))+
  xlab("")+
  ylab(expression(R^2))+
  ylim(c(0,1))+
  theme_classic()+
  theme(legend.position = c(0.45, 0.95),
        legend.direction = "horizontal")
tiff(file = file.path(outdir,  paste0(outprefix, ".R2_partial_density.tiff")), width = 1200, height = 900)
print(p)
dev.off()
pdf(file = file.path(outdir,  paste0(outprefix, ".R2_partial_density.pdf")), width = 12, height = 9)
print(p)
dev.off()
