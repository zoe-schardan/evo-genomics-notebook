#Q1. How would you characterize these populations at the state/regional level, in terms of Ho, He, and Fis?
## While all three Ho (observed heterozygosity) values seem low, WV has the highest Ho of the three states, indicating it has the greatest level of genetic variation.
## VA has the highest He (expected heterozygosity).
## OH and VA have a similar low number of heterozygotes than expected with HW Equilibrium (Fis/heterozygote deficit), with OH having the slightly larger deficit. WV has the most heterozygotes/least deficit.

#Q2. Can you come up with some R code, on your own (use any resources at your disposal), to calculate the ratio of Ho/He? 
#Can you then use an R command to add those values as an extra column in the table below? (table: print(diversity_summary) in week 3 tutorial)
##


ratio_diversity_summary <- data.frame(
  State = names(colMeans(pop_stats$Ho, na.rm = TRUE)),
  Ho    = colMeans(pop_stats$Ho, na.rm = TRUE),
  He    = colMeans(pop_stats$Hs, na.rm = TRUE),
  Fis   = colMeans(pop_stats$Fis, na.rm = TRUE),
`Ho/He` = colMeans(pop_stats$Ho, na.rm = TRUE)/colMeans(pop_stats$Hs, na.rm = TRUE)
)

#Q3. Examine the PCA output. How would you describe the population structure of S. saxatilis, if there is any? 
#Are there any patterns of distinctness among states or localities (subpopulations)?
##There appears to be some structure within these populations, with most of the WV subpopulations clustered together but still clustered with a decent number of 
##subpopulations from OH and VA. However, there is a secondary spread of primarily OH-LV and -SR-NP on the middle/middle-left side of the PCA plot, and a spread 
##of VA-BHG-DS, -ES, and -U2 below the OH secondary spread. All together, it would seem the subpopulations are largely the same across the three states, but several OH
##and VA subpopulations are more similar within the state than across states.

###code from week 3 assignment, had to edit some to get to run? worked the first time (9/2/26) so unsure what happened
suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(hierfstat)
  library(tidyverse)
  library(ggsci)
})

meta <- read.csv("master_metadata_scut_2026_07_27.csv", stringsAsFactors = FALSE)
vcf <- read.vcfR("fixed_Tassel_filtered_2023_12_06_noindels.vcf", verbose = FALSE)
genind_obj <- vcfR2genind(vcf)
sample_names <- indNames(genind_obj)
stopifnot(nrow(meta) == length(sample_names))
strata(genind_obj) <- data.frame(
  State = meta$state,
  Pop = meta$pop,
  Subpop = meta$subpop
)
pop(genind_obj) <- as.numeric(factor(meta$state))

hf_obj <- genind2hierfstat(genind_obj)
pop_stats <- basic.stats(hf_obj)
diversity_summary <- data.frame(
  State = names(colMeans(pop_stats$Ho, na.rm = TRUE)),
  Ho    = colMeans(pop_stats$Ho, na.rm = TRUE),
  He    = colMeans(pop_stats$Hs, na.rm = TRUE),
  Fis   = colMeans(pop_stats$Fis, na.rm = TRUE)
)
print(diversity_summary)

gen_imputed <- tab(genind_obj, NA.method = "mean")
pca_pass1 <- dudi.pca(gen_imputed, scannf = FALSE, nf = 2)
outliers <- rownames(pca_pass1$li)[pca_pass1$li$Axis2 > 25 | pca_pass1$li$Axis2 < -10]

genind_clean <- genind_obj[!indNames(genind_obj) %in% outliers, ]
gen_imputed_pass2 <- tab(genind_clean, NA.method = "mean")
pca_pass2 <- dudi.pca(gen_imputed_pass2, scannf = FALSE, nf = 2)
outliers_pass2 <- rownames(pca_pass2$li)[pca_pass2$li$Axis2 > 25]
all_outliers <- unique(c(outliers, outliers_pass2))

genind_clean <- genind_obj[!indNames(genind_obj) %in% all_outliers, ]
meta_clean   <- meta[!meta$name_in_vcf %in% all_outliers, ]
gen_imputed_final <- tab(genind_clean, NA.method = "mean")
pca_res <- dudi.pca(gen_imputed_final, scannf = FALSE, nf = 4)
percent_var <- (pca_res$eig / sum(pca_res$eig)) * 100

#meta_clean$PC1 <- pca_res$li$Axis1 ## 57 rows after excluding outliers; data has 59 rows. need to have same number of rows
#meta_clean$PC2 <- pca_res$li$Axis2

id="31764"
meta_clean$PC1 <- NA_real_
meta_clean$PC2 <- NA_real_
#makes vector of NAs so all 59 values are accounted for
pca_idx <- match(rownames(pca_res$li), indNames(genind_obj))
meta_clean$PC1[pca_idx] <- pca_res$li$Axis1
meta_clean$PC2[pca_idx] <- pca_res$li$Axis2
#reassigns calculated PC1/2 scores to their respective 57 remaining individuals
#marks outliers (34 and 47) as NA 

cat("Removed Outlier Samples:", paste(all_outliers, collapse = ", "), "\n")

p_pca <- ggplot(meta_clean, aes(x = PC1, y = PC2, fill = subpop, shape = state)) +
  geom_point(size = 3.5, stroke = 0.3, color = "black") +
  scale_shape_manual(values = c("OH" = 21, "VA" = 22, "WV" = 24)) +
  scale_fill_igv() +
  labs(
    title = "Principal Component Analysis (PCA - Outliers Removed)",
    x = paste0("PC1 (", round(percent_var[1], 1), "% Var)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "% Var)"),
    fill = "Subpopulation",
    shape = "State"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 3.5, color = "black"))
  )
print(p_pca)
ggsave("Figure1_PCA_Ordination_Cleaned.pdf", p_pca, width = 8, height = 6)