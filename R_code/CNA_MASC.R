# load packages and functions
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(rcna))
suppressMessages(library(Seurat))
suppressMessages(library(stringr))

source("/home/jinamo@xsede.org/scripts/cytof/Seurat_functions.R")
source("/home/jinamo@xsede.org/scripts/cytof/MASC.R")

# set parameters
prop = commandArgs(trailingOnly=TRUE)[1] %>% as.numeric()
n_min = commandArgs(trailingOnly=TRUE)[2] %>% as.integer()
frac = commandArgs(trailingOnly=TRUE)[3] %>% as.numeric()

n_neighbors = commandArgs(trailingOnly=TRUE)[4] %>% as.integer()
min_dist = commandArgs(trailingOnly=TRUE)[5] %>% as.numeric()
resolution_list = commandArgs(trailingOnly=TRUE)[6] %>% as.numeric()

obj = readRDS(paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/SeuratObj_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))
obj@meta.data$subject_id <- stringr::str_split(obj@meta.data$OmiqFileIndex, pattern="_", simplify=TRUE)[,3] %>% 
  stringr::str_split(., pattern="\\-V", simplify=TRUE) %>%
  as.data.frame() %>%
  .[,1] %>%
  gsub("-","_",.)
umap_res <- readRDS(file=paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/umap_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))
# remove clusters less than cut off
min_cell_cluster = 30
clu = umap_res %>%
  dplyr::group_by(res_cell) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  dplyr::arrange(count) %>%
  dplyr::filter(count > min_cell_cluster) %>%
  .$res_cell %>%
  unique()
clu_logi = umap_res$res_cell %in% clu
umap_res = umap_res[clu_logi,]

batch_fl = readRDS(file=paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/meta_Prop",prop,"_Nmin",n_min,".rds"))
batch_fl$subject_id <- stringr::str_split(batch_fl$OmiqFileIndex, pattern="_", simplify=TRUE)[,3] %>% 
  stringr::str_split(., pattern="\\-V", simplify=TRUE) %>%
  as.data.frame() %>%
  .[,1] %>%
  gsub("-","_",.)
batch_fl = batch_fl[clu_logi,]
umap_res$subject_id = batch_fl$subject_id
umap_res$time_point = batch_fl$time_point

# meta data
meta<- read.table("/projects/jinamo@xsede.org/cytof/data/RA_AtRiskRA_Control_meta.txt",sep="\t",header=TRUE, na.strings=c("","NA")) %>%
  dplyr::mutate(age_at_biopsy = as.integer(age_at_biopsy),
                diagnosis = factor(diagnosis, levels=c("Control","AtRiskRA","RA")))
meta = meta[order(match(meta$subject_id, obj@meta.data$subject_id)),]
meta = merge(meta,obj@meta.data,by="subject_id")

meta_add = read.table("/projects/jinamo@xsede.org/cytof/data/meta_clinical_RA.txt", header = TRUE, sep = "\t") %>%
  dplyr::mutate(CDAI = as.numeric(CDAI),
                subject_id = gsub("-","_",subject_id),
                treatment = dplyr::case_when(
                  treatment == "1" ~ "naive",
                  treatment == "2" ~ "MTX failure",
                  treatment == "3" ~ "TNF failure",
                  treatment == "4" ~ "OA"
                )) %>%
  dplyr::select(-c(sex))
meta_add = dplyr::left_join(obj@meta.data[,!grepl("CDAI",colnames(obj@meta.data))],meta_add,by="subject_id")

meta_add2 = read.table("/projects/jinamo@xsede.org/cytof/data/clinical_data_for_FAN_from_Kat.txt", header = TRUE, sep = "\t") %>%
  dplyr::select(-c(sex))
meta_add2 = dplyr::left_join(obj@meta.data[,!grepl("CDAI",colnames(obj@meta.data))],meta_add2,by="subject_id")

meta_add3 = read.table("/projects/jinamo@xsede.org/cytof/data/RA_Cleaned_aggregate_vars.txt", header = TRUE, sep = "\t") %>%
  dplyr::select(c(subject_id,HAQ,das28_crp3,bmi,ra_group,pathotype_str,sdai,ccp_type,ccp_range,rf_range,mdtjc28_sum,mdsjc28_sum,diabetes,Race_white,Race_black,Race_multiple,Race_other,MTX,SSZ,HCQ,LEF,TNFi))
meta_add3 = dplyr::left_join(obj@meta.data,meta_add3,by="subject_id")

# merge meta data to Seurat object 
obj@meta.data$AtRisk = dplyr::case_when(
  meta$AMP_Group == "FDR" & meta$CCP == "Positive" ~ "FDR(+)/ACPA(+)",
  meta$AMP_Group == "GP" & meta$CCP == "Positive" ~ "FDR(-)/ACPA(+)",
  meta$AMP_Group == "FDR" & meta$CCP == "Negative" ~ "FDR(+)/ACPA(-)",
  meta$AMP_Group == "GP" & meta$CCP == "Negative" ~ "FDR(-)/ACPA(-)",
  obj@meta.data$disease == "Control" ~ "Control",
  obj@meta.data$disease == "RA" ~ "RA")
obj@meta.data$age_at_biopsy <- as.numeric(meta$age_at_biopsy)
obj@meta.data$sex <- as.numeric(factor(meta$sex, c('female', 'male')))
obj@meta.data$ethnicity <- as.numeric(factor(meta$ethnicity, c('not_hispanic', 'hispanic_latino')))
obj@meta.data$AMP_Group <- as.numeric(factor(meta$AMP_Group, c('GP', 'FDR')))
obj@meta.data$CCP <- as.numeric(factor(meta$CCP, c('Negative', 'Positive')))
obj@meta.data$CCP30_titer <- as.numeric(meta$CCP30_titer)
obj@meta.data$CCP31_titer <- as.numeric(meta$CCP31_titer)
obj@meta.data$RF <- as.numeric(factor(meta$RF, c('Negative', 'Positive')))
obj@meta.data$RF_IgM_titer <- as.numeric(meta$RF_IgM_titer)
obj@meta.data$RF_IgG_titer <- as.numeric(meta$RF_IgG_titer)
obj@meta.data$RF_IgA_titer <- as.numeric(meta$RF_IgA_titer)
obj@meta.data$batch <- as.numeric(factor(meta$batch))
obj@meta.data$AtRisk = factor(obj@meta.data$AtRisk,levels=c("Control","FDR(+)/ACPA(+)","FDR(-)/ACPA(+)","FDR(+)/ACPA(-)","FDR(-)/ACPA(-)","RA"))
obj@meta.data$alt_AtRisk = dplyr::case_when(
  obj@meta.data$AtRisk == "FDR(-)/ACPA(-)" ~ "FDR(-)/ACPA(-)",
  obj@meta.data$AtRisk == "FDR(+)/ACPA(-)" | obj@meta.data$AtRisk == "FDR(-)/ACPA(+)" | obj@meta.data$AtRisk == "FDR(+)/ACPA(+)" ~ "AtRiskRA",
  TRUE ~ "others")
obj@meta.data$ctap <- dplyr::case_when(
  meta$diagnosis == "AtRiskRA" ~ "AtRiskRA",
  meta$diagnosis == "Control" ~ "Control")
obj@meta.data$ctap <- ifelse(is.na(obj@meta.data$ctap),meta$ctap,obj@meta.data$ctap)
obj@meta.data$CDAI <- as.numeric(meta_add$CDAI)
obj@meta.data$krenn_lining <- as.numeric(meta_add$krenn_lining)
obj@meta.data$krenn_inflammation <- as.numeric(meta_add$krenn_inflammation)
obj@meta.data$treatment <- dplyr::case_when(
  meta$diagnosis == "AtRiskRA" ~ "AtRiskRA",
  meta$diagnosis == "Control" ~ "Control")
obj@meta.data$treatment <- ifelse(is.na(obj@meta.data$treatment),meta_add$treatment,obj@meta.data$treatment)
obj@meta.data$disease = dplyr::case_when(
  obj@meta.data$AtRisk == "FDR(+)/ACPA(+)" | obj@meta.data$AtRisk == "FDR(-)/ACPA(+)" | obj@meta.data$AtRisk == "FDR(+)/ACPA(-)" ~ "AtRiskRA",
  obj@meta.data$AtRisk == "Control" | obj@meta.data$AtRisk == "FDR(-)/ACPA(-)" ~ "Control",
  obj@meta.data$AtRisk == "RA" ~ "RA")
obj@meta.data$condition = obj@meta.data$disease
obj@meta.data$condition = factor(obj@meta.data$condition, levels = c('Control', 'AtRiskRA', "RA"))
obj@meta.data$CRP <- as.numeric(meta_add2$crp_result)
obj@meta.data$ESR <- as.numeric(meta_add2$esr_result)
obj@meta.data$RA_duration_years <- as.numeric(meta_add2$RA_duration_years)
obj@meta.data$CDAI <- ifelse(is.na(obj@meta.data$CDAI), as.numeric(meta_add2$cdai), obj@meta.data$CDAI)
obj@meta.data$SDAI <- as.numeric(meta_add3$sdai)
obj@meta.data$DAS28_ESR <- as.numeric(meta_add2$das28_esr3)
obj@meta.data$DAS28_CRP <- as.numeric(meta_add3$das28_crp3)
obj@meta.data$CCP31_titer <- ifelse(obj@meta.data$condition == "RA", as.numeric(meta_add2$ccp_result), obj@meta.data$CCP31_titer)
obj@meta.data$CCP31_cutoff <- as.numeric(meta_add3$ccp_range)
obj@meta.data$RF_IgM_titer <- ifelse(obj@meta.data$condition == "RA", as.numeric(meta_add2$rf_result), obj@meta.data$RF_IgM_titer)
obj@meta.data$RF_IgM_cutoff <- as.numeric(meta_add3$rf_range)
obj@meta.data$smoke_hx <- as.factor(meta_add2$smoke_hx)
obj@meta.data$smoke_current <- as.factor(meta_add2$smoke_current)
obj@meta.data$Methotrexate_only <- meta_add2$Methotrexate_only
obj@meta.data$TNFi_possible_nbDMARD <- meta_add2$TNFi_possible_nbDMARD
obj@meta.data$other_bDMARD <- meta_add2$other_bDMARD
obj@meta.data$other_nbDMARD_no_bDMARD <- meta_add2$other_nbDMARD_no_bDMARD
obj@meta.data$no_DMARD <- meta_add2$no_DMARD
obj@meta.data$bDMARD <- meta_add2$bDMARD
obj@meta.data$nbDMARD <- meta_add2$nbDMARD
obj@meta.data$Prednisone <- meta_add2$Prednisone
obj@meta.data$Prednisone_GE_7.5 <- ifelse(obj@meta.data$Prednisone, meta_add2$Prednisone_GE_7.5, FALSE)
obj@meta.data$time_point <- batch_fl$time_point
obj@meta.data$site <- ifelse(obj@meta.data$condition == "RA", meta_add$site, "University of Colorado") 
obj@meta.data$site <- factor(obj@meta.data$site, levels = c("University of Colorado","University of Rochester","Hospital for Special Surgery","Feinstein Institute","Cedars","UC San Diego","Northwestern","UK Birmingham","UK London","UAB","Columbia University")) 
obj@meta.data$site = as.integer(obj@meta.data$site)

# subset Control and AtRiskRA 
logi = obj$condition %in% c("Control","AtRiskRA")
obj_sub <- subset(x = obj, condition == "Control" | condition == "AtRiskRA")
obj_sub@graphs$RNA_nn = as.Graph(obj@graphs$RNA_nn[logi,logi])
obj_sub@graphs$RNA_snn = as.Graph(obj@graphs$RNA_snn[logi,logi])
obj_sub@meta.data$condition <- as.character(obj_sub@meta.data$condition)
obj_sub@meta.data$condition <- factor(obj_sub@meta.data$condition, c('Control', 'AtRiskRA'))
obj_sub@meta.data$condition_val <- as.numeric(factor(obj_sub@meta.data$condition, c('Control', 'AtRiskRA')))

obj_sub <- association.Seurat(
  seurat_object = obj_sub, 
  test_var = 'condition_val', 
  samplem_key = 'subject_id', 
  graph_use = 'RNA_nn', 
  verbose = TRUE,
  batches = NULL,
  covs = c("age_at_biopsy","sex")
)

saveRDS(obj_sub, paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/SeuratObj_Control_AtRiskRA_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))

# MASC: Mixed-effects association testing for single cells
## Sci Transl Med. 2018 Oct 17;10(463):eaaq0305.

df = cbind(umap_res[logi,],obj_sub@meta.data)

masc = MASC(dataset = df,
            cluster = df$res_cell,
            contrast = "condition",
            random_effects = c("subject_id"),
            fixed_effects = c("age_at_biopsy","sex"),
            verbose = TRUE, save_models = FALSE)

saveRDS(masc, paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/MASC_Control_AtRiskRA_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))

# subset AtRiskRA and RA
logi = obj$condition %in% c("RA","AtRiskRA")
obj_sub <- subset(x = obj, condition == "RA" | condition == "AtRiskRA")
obj_sub@graphs$RNA_nn = as.Graph(obj@graphs$RNA_nn[logi,logi])
obj_sub@graphs$RNA_snn = as.Graph(obj@graphs$RNA_snn[logi,logi])
obj_sub@meta.data$condition <- as.character(obj_sub@meta.data$condition)
obj_sub@meta.data$condition <- factor(obj_sub@meta.data$condition, c('RA', 'AtRiskRA'))
obj_sub@meta.data$condition_val <- as.numeric(factor(obj_sub@meta.data$condition, c('RA', 'AtRiskRA')))

obj_sub <- association.Seurat(
  seurat_object = obj_sub, 
  test_var = 'condition_val', 
  samplem_key = 'subject_id', 
  graph_use = 'RNA_nn', 
  verbose = TRUE,
  batches = NULL,
  covs = c("age_at_biopsy","sex")
)

saveRDS(obj_sub, paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/SeuratObj_AtRiskRA_RA_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))

df = cbind(umap_res[logi,],obj_sub@meta.data)

masc = MASC(dataset = df,
            cluster = df$res_cell,
            contrast = "condition",
            random_effects = c("subject_id"),
            fixed_effects = c("sex","age_at_biopsy"),
            verbose = TRUE, save_models = FALSE)

saveRDS(masc, paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/MASC_AtRiskRA_RA_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))

# subset Control and RA
logi = obj$condition %in% c("RA","Control")
obj_sub <- subset(x = obj, condition == "RA" | condition == "Control")
obj_sub@graphs$RNA_nn = as.Graph(obj@graphs$RNA_nn[logi,logi])
obj_sub@graphs$RNA_snn = as.Graph(obj@graphs$RNA_snn[logi,logi])
obj_sub@meta.data$condition <- as.character(obj_sub@meta.data$condition)
obj_sub@meta.data$condition <- factor(obj_sub@meta.data$condition, c('Control', 'RA'))
obj_sub@meta.data$condition_val <- as.numeric(factor(obj_sub@meta.data$condition, c('Control', 'RA')))

obj_sub <- association.Seurat(
  seurat_object = obj_sub, 
  test_var = 'condition_val', 
  samplem_key = 'subject_id', 
  graph_use = 'RNA_nn', 
  verbose = TRUE,
  batches = NULL,
  covs = c("age_at_biopsy","sex")
)

saveRDS(obj_sub, paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/SeuratObj_Control_RA_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))

df = cbind(umap_res[logi,],obj_sub@meta.data)

masc = MASC(dataset = df,
            cluster = df$res_cell,
            contrast = "condition",
            random_effects = c("subject_id"),
            fixed_effects = c("sex","age_at_biopsy"),
            verbose = TRUE, save_models = FALSE)

saveRDS(masc, paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/MASC_Control_RA_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))

# subset Control and RA(V0)
logi = (obj$condition %in% c("RA","Control")) & (obj$time_point == "V0")
obj_sub <- subset(x = obj, time_point == "V0" & (condition == "RA" | condition == "Control"))
obj_sub@graphs$RNA_nn = as.Graph(obj@graphs$RNA_nn[logi,logi])
obj_sub@graphs$RNA_snn = as.Graph(obj@graphs$RNA_snn[logi,logi])
obj_sub@meta.data$condition <- as.character(obj_sub@meta.data$condition)
obj_sub@meta.data$condition <- factor(obj_sub@meta.data$condition, c('Control', 'RA'))
obj_sub@meta.data$condition_val <- as.numeric(factor(obj_sub@meta.data$condition, c('Control', 'RA')))

obj_sub <- association.Seurat(
  seurat_object = obj_sub, 
  test_var = 'condition_val', 
  samplem_key = 'subject_id', 
  graph_use = 'RNA_nn', 
  verbose = TRUE,
  batches = NULL,
  covs = c("age_at_biopsy","sex")
)

saveRDS(obj_sub, paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/SeuratObj_Control_RAV0_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))

df = cbind(umap_res[logi,],obj_sub@meta.data)

masc = MASC(dataset = df,
            cluster = df$res_cell,
            contrast = "condition",
            random_effects = c("subject_id"),
            fixed_effects = c("sex","age_at_biopsy"),
            verbose = TRUE, save_models = FALSE)

saveRDS(masc, paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/MASC_Control_RAV0_Prop",prop,"_Nmin",n_min,"_topVar",frac,"_nneighbors",n_neighbors,"_mindist",min_dist,"_res",resolution_list,".rds"))
