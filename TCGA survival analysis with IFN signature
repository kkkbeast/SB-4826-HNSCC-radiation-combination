BiocManager::install(c("TCGAbiolinks", "survival", "survminer", "SummarizedExperiment"))

# Load the libraries
library(ggplot2)
library(survminer)
library(TCGAbiolinks)
library(survival)
library(SummarizedExperiment)
library(dplyr)


# 1. Define your new 2-cancer cohort
my_projects_h <- c("TCGA-HNSC")
my_projects_R <- c("TCGA-READ")

my_projects_L <- c("TCGA-SARC")


# 2. Re-download Clinical Data
clin_query_pan <- GDCquery(project = my_projects_L,
                           data.category = "Clinical",
                           data.type = "Clinical Supplement",
                           data.format = "BCR Biotab")
GDCdownload(clin_query_pan)
clin_list <- GDCprepare(clin_query_pan)

# 3. Re-download Expression Data (This will take a few minutes!)
exp_query_pan <- GDCquery(project = my_projects_L,
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts")
GDCdownload(exp_query_pan)
exp_data_pan <- GDCprepare(exp_query_pan)



# Extract the raw matrix
expr_matrix_pan <- assay(exp_data_pan, "tpm_unstrand")
gene_info_pan <- rowData(exp_data_pan)

# FILTER FOR PRIMARY TUMORS ONLY (01)
sample_codes <- substr(colnames(expr_matrix_pan), 14, 15)
expr_matrix_tumor <- expr_matrix_pan[, sample_codes == "01"]
# Define target genes and find IDs
target_genes <- c("APOL6","ARID5B","ARL4A","AUTS2","B2M","BANK1","BATF2","BPGM","BST2","BTG1","C1R","C1S","CASP1","CASP3","CASP4","CASP7","CASP8","CCL2","CCL5","CCL7","CD274","CD38","CD40","CD69","CD74","CD86","CDKN1A","CFB","CFH","CIITA","CMKLR1","CMPK2","CSF2RB","CXCL10","CXCL11","CXCL9","RIGI","DDX60","DHX58","EIF2AK2","EIF4E3","EPSTI1","FAS","FCGR1A","FGL2","FPR1","CMTR1","GBP4","GBP6","GCH1","GPR18","GZMA","HERC6","HIF1A","HLA-A","HLA-B","HLA-DMA","HLA-DQA1","HLA-DRB1","HLA-G","ICAM1","IDO1","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","IFITM2","IFITM3","IFNAR2","IL10RA","IL15","IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7","IRF1","IRF2","IRF4","IRF5","IRF7","IRF8","IRF9","ISG15","ISG20","ISOC1","ITGB7","JAK2","KLRK1","LAP3","LATS2","LCP2","LGALS3BP","LY6E","LYSMD2","MARCHF1","TMT1B","MT2A","MTHFD2","MVP","MX1","MX2","MYD88","NAMPT","NCOA3","NFKB1","NFKBIA","NLRC5","NMI","NOD1","NUP93","OAS2","OAS3","OASL","OGFR","P2RY14","PARP12","PARP14","PDE4B","PELI1","PFKP","PIM1","PLA2G4A","PLSCR1","PML","PNP","PNPT1","HELZ2","PSMA2","PSMA3","PSMB10","PSMB2","PSMB8","PSMB9","PSME1","PSME2","PTGS2","PTPN1","PTPN2","PTPN6","RAPGEF6","RBCK1","RIPK1","RIPK2","RNF213","RNF31","RSAD2","RTP4","SAMD9L","SAMHD1","SECTM1","SELP","SERPING1","SLAMF7","SLC25A28","SOCS1","SOCS3","SOD2","SP110","SPPL2A","SRI","SSPN","ST3GAL5","ST8SIA4","STAT1","STAT2","STAT3","STAT4","TAP1","TAPBP","TDRD7","TNFAIP2","TNFAIP3","TNFAIP6","TNFSF10","TOR1B","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TXNIP","UBE2L6","UPP1","USP18","VAMP5","VAMP8","VCAM1","WARS1","XAF1","XCL1","ZBP1", "ZNFX1"
)
target_genes2 <- c("ADAR","B2M","BATF2","BST2","C1S","CASP1","CASP8","CCRL2","CD47","CD74","CMPK2","CNP","CSF1","CXCL10","CXCL11","DDX60","DHX58","EIF2AK2","ELF1","EPSTI1","MVB12A","TENT5A","CMTR1","GBP2","GBP4","GMPR","HERC6","HLA-C","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IL15","IL4R","IL7","IRF1","IRF2","IRF7","IRF9","ISG15","ISG20","LAMP3","LAP3","LGALS3BP","LPAR6","LY6E","MOV10","MX1","NCOA7","NMI","NUB1","OAS1","OASL","OGFR","PARP12","PARP14","PARP9","PLSCR1","PNPT1","HELZ2","PROCR","PSMA3","PSMB8","PSMB9","PSME1","PSME2","RIPK2","RNF31","RSAD2","RTP4","SAMD9","SAMD9L","SELL","SLC25A28","SP110","STAT2","TAP1","TDRD7","TMEM140","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TRIM5","TXNIP","UBA7","UBE2L6","USP18","WARS1")
target_gene_comb = union(target_genes, target_genes2)

gene_ids_pan <- rownames(gene_info_pan)[gene_info_pan$gene_name %in% target_gene_comb]

# Extract, transpose, and scale
sig_matrix_pan <- expr_matrix_tumor[gene_ids_pan, ]
sig_matrix_t_pan <- t(sig_matrix_pan)
scaled_matrix_pan <- scale(sig_matrix_t_pan)

# Calculate average signature score per patient
patient_scores_pan <- rowMeans(scaled_matrix_pan, na.rm = TRUE)

# Create the expression dataframe
sig_df_pan <- data.frame(
  bcr_patient_barcode = substr(rownames(scaled_matrix_pan), 1, 12),
  Signature_Score = as.numeric(patient_scores_pan)
)
sig_df_pan <- sig_df_pan[!duplicated(sig_df_pan$bcr_patient_barcode), ]


# Bind and clean using coalesce 
# 1. Find the exact names of the patient tables inside clin_list
patient_table_names <- grep("^clinical_patient", names(clin_list), value = TRUE)

# 2. Extract just those specific dataframes (HNSC and DLBCL)
patient_tables <- clin_list[patient_table_names]

# 3. Bind them together safely using dplyr
clean_clin_pan <- as.data.frame(bind_rows(patient_tables))

# Print the row count to verify it worked! (Should be over 500)
print(paste("Rows successfully extracted:", nrow(clean_clin_pan)))

# 4. Replace text NAs
clean_clin_pan[clean_clin_pan == "[Not Applicable]"] <- NA
clean_clin_pan[clean_clin_pan == "[Not Available]"] <- NA
clean_clin_pan[clean_clin_pan == "[Unknown]"] <- NA
clean_clin_pan[clean_clin_pan == ""] <- NA

# 5. Create dummy columns to prevent missing-column crashes
if(!"days_to_death" %in% colnames(clean_clin_pan)) clean_clin_pan$days_to_death <- NA
if(!"death_days_to" %in% colnames(clean_clin_pan)) clean_clin_pan$death_days_to <- NA
if(!"days_to_last_followup" %in% colnames(clean_clin_pan)) clean_clin_pan$days_to_last_followup <- NA
if(!"last_contact_days_to" %in% colnames(clean_clin_pan)) clean_clin_pan$last_contact_days_to <- NA

# 6. Coalesce into unified survival timelines
# (Note: "NAs introduced by coercion" warnings here are completely normal!)
clean_clin_pan$unified_death <- coalesce(as.numeric(clean_clin_pan$death_days_to),
                                         as.numeric(clean_clin_pan$days_to_death))
clean_clin_pan$unified_alive <- coalesce(as.numeric(clean_clin_pan$last_contact_days_to),
                                         as.numeric(clean_clin_pan$days_to_last_followup))

# 7. Calculate Survival Variables
clean_clin_pan$OS_status <- ifelse(toupper(clean_clin_pan$vital_status) == "DEAD", 1, 0)
clean_clin_pan$OS_time <- ifelse(clean_clin_pan$OS_status == 1,
                                 clean_clin_pan$unified_death,
                                 clean_clin_pan$unified_alive)

# 8. Filter out missing data and 0-day follow-ups
clean_clin_pan <- clean_clin_pan[!is.na(clean_clin_pan$OS_time) & clean_clin_pan$OS_time > 0, ]

print(paste("Final valid patients ready for plotting:", nrow(clean_clin_pan)))

# Merge clinical and expression
df_pan_final <- merge(clean_clin_pan, sig_df_pan, by = "bcr_patient_barcode")
df_pan_final <- df_pan_final[!is.na(df_pan_final$Signature_Score), ]

# ---------------------------------------------------------
# FILTER OUT 'NA' FROM TUMOR STATUS
# ---------------------------------------------------------

# Remove rows where tumor_status is NA
df_pan_final <- df_pan_final[!is.na(df_pan_final$tumor_status), ]

# Optional: If there are other weird TCGA placeholders like "[Unknown]" or "[Not Available]", 
# we can clean those out too just to be safe:
df_pan_final <- df_pan_final[df_pan_final$tumor_status != "[Unkown]" & 
                               df_pan_final$tumor_status != "[Not Available]", ]

# Print the results to verify it worked perfectly
print(paste("Total patients remaining with known tumor status:", nrow(df_pan_final)))

print("Patient counts by tumor status:")
print(table(df_pan_final$tumor_status))


print(paste("Total HNSC:", nrow(df_pan_final)))

# Create dummy columns to prevent errors if your specific project is missing one
if(!"new_tumor_event_dx_days_to" %in% colnames(df_pan_final)) df_pan_final$new_tumor_event_dx_days_to <- NA
if(!"days_to_new_tumor_event_after_initial_treatment" %in% colnames(df_pan_final)) df_pan_final$days_to_new_tumor_event_after_initial_treatment <- NA

# Unify the time-to-recurrence data
df_pan_final$time_to_new_tumor <- coalesce(
  as.numeric(df_pan_final$new_tumor_event_dx_days_to),
  as.numeric(df_pan_final$days_to_new_tumor_event_after_initial_treatment)
)

# Calculate PFS Status (1 = Progression or Death, 0 = Alive and Tumor-Free)
df_pan_final$PFS_status <- ifelse(!is.na(df_pan_final$time_to_new_tumor) | df_pan_final$OS_status == 1, 1, 0)

# Calculate PFS Time 
df_pan_final$PFS_time <- ifelse(
  df_pan_final$PFS_status == 1,
  pmin(df_pan_final$time_to_new_tumor, df_pan_final$unified_death, na.rm = TRUE),
  df_pan_final$unified_alive
)

# Clean out any missing or invalid data (negative days)
df_pfs_clean <- df_pan_final[!is.na(df_pan_final$PFS_time) & df_pan_final$PFS_time > 0, ]

grep("malignancy", colnames(df_pfs_clean), value = TRUE, ignore.case = TRUE)

# Step 2: Count how many patients fall into each category
# ('useNA = "always"' ensures you see missing data blank spaces)
table(df_pfs_clean$history_other_malignancy, useNA = "always")

# Step 3: Get just the absolute count of patients with a history 
num_with_history <- sum(df_pfs_clean$history_other_malignancy == "Yes", na.rm = TRUE)
cat("Total patients with a history of other malignancy:", num_with_history, "\n")
df_pfs_clean <- df_pfs_clean %>%
  filter(history_other_malignancy == "No")
table(df_pfs_clean$history_other_malignancy, useNA = "always")


# ---------------------------------------------------------
# 5-YEAR OS FOR THE ENTIRE COHORT 
# ---------------------------------------------------------

# 1. Define 5 years in days
five_years <- 1826



# 2. Start with a clean copy of the full dataset (removing invalid/negative days)
df_5yr_all <- df_pfs_clean[!is.na(df_pfs_clean$OS_time) & df_pfs_clean$OS_time > 0, ]

# 3. Apply the 5-Year Truncation (Censoring)
# If their time is greater than 1826 days, lock their time at 1826. Otherwise, keep real time.
df_5yr_all$OS_time_5yr <- ifelse(df_5yr_all$OS_time > five_years, five_years, df_5yr_all$OS_time)

# If their time was greater than 1826 days, they "survived" the 5 years, so status becomes 0 (Alive).
df_5yr_all$OS_status_5yr <- ifelse(df_5yr_all$OS_time > five_years, 0, df_5yr_all$OS_status)

# 4. Find the Mathematically Optimal Cutpoint for this specific 5-year window
res.cut.5yr <- surv_cutpoint(df_5yr_all, 
                             time = "OS_time_5yr", 
                             event = "OS_status_5yr", 
                             variables = "Signature_Score")

print("The mathematically optimal score threshold for 5-Year OS is:")
print(summary(res.cut.5yr))

# 5. Categorize the patients based on the algorithm's chosen cutpoint
df_5yr_opt <- surv_categorize(res.cut.5yr)

# Format the text so it looks perfect in the plot legend
df_5yr_opt$Group <- ifelse(df_5yr_opt$Signature_Score == "high", "High IFN Signature", "Low IFN Signature")
df_5yr_opt$Group <- factor(df_5yr_opt$Group, levels = c("Low IFN Signature", "High IFN Signature"))

# 6. Fit the 5-Year Survival Model
fit_5yr_all <- survfit(Surv(OS_time_5yr, OS_status_5yr) ~ Group, data = df_5yr_opt)

# 7. Generate the Beautiful Plot
plot_5yr_all <- ggsurvplot(fit_5yr_all, 
                           data = df_5yr_opt, 
                           pval = TRUE,              # P-value for the 5-year difference
                           conf.int = FALSE,
                           risk.table = TRUE,
                           palette = c("blue", "red"),
                           title = "5-Year LUSC Overall Survival",
                           xlab = "Time (Days)", 
                           ylab = "5-Year Survival Probability",
                           xlim = c(0, five_years),  # Lock the graph strictly to 1826 days
                           break.time.by = 365,      # Add clean tick marks every 1 year (365 days)
                           legend.title = "Algorithm-Defined Score",
                           legend.labs = c("Low IFN Signature", "High IFN Signature"))

print(plot_5yr_all)




