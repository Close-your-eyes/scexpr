SO_path <- "/Volumes/CMS_SSD_2TB/R_scRNAseq/2019_SLE_LN/data/SO_processed"
SO <- readRDS(file.path(SO_path, "SO_urine_blood_hg38_CR7_non_contam_rep1_wo_Introns_strict_QC_patient_integration_SCT_harmony_1_400_15_230515-113449.rds"))
SO_CD8_eff <- readRDS("/Volumes/CMS_SSD_2TB/R_scRNAseq/2019_SLE_LN/data/SO_processed/SO_urine_blood_hg38_CR7_non_contam_rep1_wo_Introns_strict_QC_patient_integration_CD8_eff_SCT_harmony_1_300_12_231016-163715.rds")
