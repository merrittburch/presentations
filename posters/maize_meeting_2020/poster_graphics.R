# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-03-02
#
# Description 
#   - Script to make maize meeting figures
# ---------------------------------------------------------------


# Load in source scripts
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/ames2any_matrix.R')
source('~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/local_window_funs.R')
source("~/git_projects/haplotype_v_snp_gwas/src/R/pca_testing/count_gwas_overlaps.R")

# Load packages
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ape)
library(dplyr)
library(broom)
library(magrittr)


# ------------------
#   Load test data
# ------------------

# Load in Beagle imputed SNPs that were filtered with bcftools
ames.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/genos/hapmap_v321_snps/ames_for_pca_merged.vcf"
ames_vcf <- snpgdsVCF2GDS(ames.vcf, "ames.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("ames.gds")
genofile_ames <- snpgdsOpen(ames_vcf)

# Load in NAM SNPs that are shared (intersect) with the Ames SNPs
# nam.vcf <- "/workdir/mbb262/genotypes/nam/imputed/combined_ibm_nam_all_imputed_bcftools/filtered_for_pca"
nam.vcf <- "~/Box Sync/Cornell_PhD/labProjects/hap_gwa/data/genos/hapmap_v321_snps/nam_by_ames_sites_for_pca_allChroms.vcf"
nam_vcf <- snpgdsVCF2GDS(nam.vcf, "nam.gds", method="copy.num.of.ref", snpfirstdim=TRUE, ignore.chr.prefix = "S")
snpgdsSummary("nam.gds")
genofile_nam <- snpgdsOpen(nam_vcf)


# ------------------------------------------
#     Calculate ames2nam global PCs
#     uses ames2any_matrix functions
# ------------------------------------------

# Combine genofile objects
X_allNAM <- combine_genofile_matrix(genofile_ames = genofile_ames, genofile_any = genofile_nam)

# ------
# 3 gPCs
# ------ 

# Calculate Ames global PCs, get loadings
global_ames_BV_3gPCs <- get_ames_BV(X_ames = genofile_ames, Q = NULL, num_PCs = 3)

# Transfer loadings and coefficients over to all X (Ames, NAM) to get adjusted global PCs
# Get Q matrix (population structure)
ames2nam_gPCs_3gPCs_allNAM <- ames2any(X_all = X_allNAM$X_all,
                                       Q = NULL, B = global_ames_BV_3gPCs$B, V = global_ames_BV_3gPCs$V)
Q_all_3gPCs_allNAM <- cbind(1, ames2nam_gPCs_3gPCs_allNAM)


# PCA on NAM alone
global_nam_BV <- get_ames_BV(X_ames = genofile_nam, Q = NULL, num_PCs = 3)
# Make NAM X matrix and do formatting
X_all_nam_geno <- snpgdsGetGeno(genofile_nam, snpfirstdim = FALSE, with.id = TRUE, verbose = FALSE)
X_all_nam <- 2-X_all_nam_geno$genotype
colnames(X_all_nam) <- X_all_nam_geno$snp.id
rownames(X_all_nam) <- c(X_all_nam_geno$sample.id)
X_all_nam <- list(X_all = X_all_nam, snp_info = snpgdsSNPList(genofile_nam))

# Calculate PCs
nam2nam_gPCs <- ames2any(X_all = X_all_nam$X_all,
                         Q = NULL, B = global_nam_BV$B, V = global_nam_BV$V)

nam_only_pcs <- prcomp(X_all_nam$X_all, center = FALSE, rank. = 3)$rotation
temp_pcs <- prcomp(X_all_nam$X_all, center = FALSE, rank. = 3)

test_woIBMs <- data.frame(temp_pcs$x)
test_woIBMs$ids <- row.names(test_woIBMs)
test_woIBMs$ids <- gsub("M[0-9]{4}", "M0", test_woIBMs$ids)
test_woIBMs$ids <- gsub("E[0-9]{4}", "", test_woIBMs$ids)
test_woIBMs <- test_woIBMs[-which(test_woIBMs$ids == "M0"),]
plot(test_woIBMs$PC1, test_woIBMs$PC3)

# ------------------
#   Visulizations
# ------------------

# Load in heterotic group information
nam_het <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/nam_subpops.csv")
ames_het <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/ames_heterotic_groups.csv")

# Merge heterotic group information to NAM-only PCs
nam_pc_hets <- merge(test_woIBMs, nam_het, by.x = "ids", by.y = "taxa")

# Merge Ames heterotic group to Ames subset
ames_pcs <- data.frame(ames2nam_gPCs_3gPCs_allNAM)
ames_pcs <- ames_pcs[1:3545,]
ames_pcs$taxa <- rownames(ames_pcs)
ames_pcs$taxa <- gsub(":[0-9]{9}", "", ames_pcs$taxa)
ames_info <- merge(ames_pcs, ames_het, by.x = "taxa", by.y = "Accesion_N")

# Merge heteroitc group info to Ames 2 NAM PCs
ames2nam_pcs <- data.frame(ames2nam_gPCs_3gPCs_allNAM)
ames2nam_pcs <- ames2nam_pcs[3546:nrow(ames2nam_pcs),]
ames2nam_pcs$taxa <- rownames(ames2nam_pcs)
ames2nam_pcs$taxa <- gsub("E[0-9]{4}", "", ames2nam_pcs$taxa)
nam_info <- merge(ames2nam_pcs, nam_het, by.x = "taxa", by.y = "taxa")

# NAM  only
c <- ggplot(nam_pc_hets, aes(PC1,PC2)) +
  geom_point(aes(colour = factor(subpop2))) + 
  labs(x ="Principal Component 1", y = "Principal Component 2") +
  theme(axis.text = element_text(size = 13),
        legend.position="bottom")

ggsave(filename = "nam_only_pcs.png", plot = c, width = 8, units = "in")

# Ames PCs
b <- ggplot(ames_info, aes(PC1,PC2)) +
  geom_point(aes(colour = factor(Pop_structure))) + 
  labs(x ="Principal Component 1", y = "Principal Component 2") +
  theme(axis.text = element_text(size = 13),
        legend.position="bottom")

ggsave(filename = "ames_only_pcs.png", plot = b, width = 8, units = "in")

# Ames 2 NAM PCs
a <- ggplot(nam_info, aes(PC1,PC2)) +
  geom_point(aes(colour = factor(subpop2))) + 
  labs(x ="Principal Component 1", y = "Principal Component 2") +
  theme(axis.text = element_text(size = 13),
        legend.position="bottom")

ggsave(filename = "ames2nam_pcs.png", plot = a, width = 8, units = "in")



