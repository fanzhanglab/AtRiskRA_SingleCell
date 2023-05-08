#!/projects/jinamo@xsede.org/software/anaconda/envs/Renv/bin/Rscript

# load packages and functions
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(harmony))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))

prop = commandArgs(trailingOnly=TRUE)[1] %>% as.numeric() #1番目の引数を取得
n_min = commandArgs(trailingOnly=TRUE)[2] %>% as.integer() #2番目の引数を取得
frac = commandArgs(trailingOnly=TRUE)[3] %>% as.numeric() #3番目の引数を取得

pca_res <- readRDS(file=paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/PCA_embeddings_Prop",prop,"_Nmin",n_min,"_topVar",frac,".rds"))
batch_fl = readRDS(file=paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/meta_Prop",prop,"_Nmin",n_min,".rds"))

# batch correction by Harmony
print(Sys.time())
harmony_embeddings_all <- HarmonyMatrix(pca_res$x, 
                                        batch_fl, 
                                        c('batch','OmiqFileIndex'), # OmiqFileIndex: sample ID
                                        do_pca=FALSE)
print(Sys.time())

print(head(harmony_embeddings_all))
print(dim(harmony_embeddings_all))

saveRDS(harmony_embeddings_all, paste0("/projects/jinamo@xsede.org/cytof/data/B_cells/harmony_embeddings_Prop",prop,"_Nmin",n_min,"_topVar",frac,".rds"))


