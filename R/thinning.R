# Template thinning and predictive performance
# Code by Andrew Vigotsky
# Last edited: November 23, 2019
# 
# For "The hard limits of decoding mental states" 
# by Jabakhanji et al.
#
# This script builds templates using X random, nonzero
# voxels from each original template. It tests the
# performance of these thinned templates on studies
# 1 through 4.

library(RNifti)
library(pROC)
library(pbmcapply)
library(Matrix)
library(qlcMatrix)

setwd("/projects/b1090/Users/vigotsky/NPS_NSy")

## Load templates and convert them to sparse matrices
nps_filtered <- Matrix( apply(readNifti("Templates/NPS_filtered.nii"),  4, c), sparse = T)
pv_filtered  <- Matrix( apply(readNifti("Templates/PPV_filtered.nii"),  4, c), sparse = T)
nsy_filtered <- Matrix( apply(readNifti("Templates/NSyP_filtered.nii"), 4, c), sparse = T)
gpt_filtered <- Matrix( apply(readNifti("Templates/GPT_filtered.nii"),  4, c), sparse = T)

## Get the mask of each template
nps_vec <- nps_filtered[,22]
pv_vec  <-  pv_filtered[,22]
nsy_vec <- nsy_filtered[,22]
gpt_vec <- gpt_filtered[,22]

## Get the nonzero voxel indices. These are candidates from which to sample.
nps_candidates <- which(nps_vec != 0)
pv_candidates  <- which(pv_vec  != 0)
nsy_candidates <- which(nsy_vec != 0)
gpt_candidates <- which(gpt_vec != 0)

## Study combinations
studies <- data.frame(studies = c(1,2,3,3,3,4),
                      comparator = c(2,2,2,3,4,2))

## Organize templates
maps <- list(nps_filtered, pv_filtered, nsy_filtered, gpt_filtered)
map_names <- c("NPS","pPV","pNsy","GPT")
map_candidates <- list(nps_candidates, pv_candidates, nsy_candidates, gpt_candidates)

## We will loop over the studies and save their output one at a time
for(study in 1:nrow(studies)) {
  study_results <- c()
  std <- studies[study,"studies"]
  comp <- studies[study,"comparator"]
  stm1 <- readNifti(paste0("BetaMaps/Study_",
                           std,
                           "_Stm_1.nii"))
  stm2 <- readNifti(paste0("BetaMaps/Study_",
                           std,
                           "_Stm_",
                           comp,
                           ".nii"))
  
  stm1_vec <- Matrix(apply(stm1, 4, c), sparse=T)
  stm2_vec <- Matrix(apply(stm2, 4, c), sparse=T)
  
  for(map in 1:length(maps)) {
    print(map_names[map])
    set.seed(0)
    results <- c()
    nsim <- 500 # number of replicates
    levels <- c(10,20,25,50,100,250,500,1000,5000,7500,10000) # voxel counts
    results <- pbmclapply(levels, function(x) {
      
      ## Generate a data frame of voxels to sample for each replicate
      mat_vox <- lapply(1:nsim, function(i) {
        data.frame(replicate = i,
                   voxel = 
                        sample(map_candidates[[map]], 
                               min(x,length(map_candidates[[map]])))
        )
      })
      mat_vox <- do.call(rbind.data.frame, mat_vox)
      
      ## Create a mask of the sampled voxels for each replicate
      map_mask <- sparseMatrix(i = mat_vox$voxel, 
                               j = mat_vox$replicate, 
                               x = 1, 
                               dims = c(nrow(maps[[map]]),nsim) 
                               )
      
      ## Apply that mask to the current template
      mat_tmp <- maps[[map]][,1] * map_mask
      mat_tmp_binary <- map_mask
      mat_tmp_sign <- sign(mat_tmp)
      
      ## Calculate normalized dot products and AUCs for the unfiltered map
      nofilt_dps <- list(cosSparse(stm1_vec,mat_tmp),
                         cosSparse(stm2_vec,mat_tmp))
      nofilt.results <- sapply(1:nsim, function(i) { 
        df.nofilt <- rbind(data.frame(dp = nofilt_dps[[1]][,i],
                                      condition = 1),
                           data.frame(dp = nofilt_dps[[2]][,i],
                                      condition = 0))
        auc(roc(condition ~ dp, df.nofilt, direction = "<", quiet = T, algorithm = 3))[[1]]
      })
      
      ## Calculate normalized dot products and AUCs for the infinitely filtered map
      filt_dps <- list(cosSparse(stm1_vec,mat_tmp_binary),
                       cosSparse(stm2_vec,mat_tmp_binary))
      filt.results <- sapply(1:500, function(i) { 
        df.nofilt <- rbind(data.frame(dp = filt_dps[[1]][,i],
                                      condition = 1),
                           data.frame(dp = filt_dps[[2]][,i],
                                      condition = 0))
        auc(roc(condition ~ dp, df.nofilt, direction = "<", quiet = T, algorithm = 3))[[1]]
      })
      
      ## Calculate normalized dot products and AUCs for sign(map)
      sign_dps <- list(cosSparse(stm1_vec,mat_tmp_binary),
                       cosSparse(stm2_vec,mat_tmp_binary))
      sign.results <- sapply(1:500, function(i) { 
        df.nofilt <- rbind(data.frame(dp = sign_dps[[1]][,i],
                                      condition = 1),
                           data.frame(dp = sign_dps[[2]][,i],
                                      condition = 0))
        auc(roc(condition ~ dp, df.nofilt, direction = "<", quiet = T, algorithm = 3))[[1]]
      })
      
      data.frame(unfiltered = nofilt.results,
                 filtered = filt.results,
                 sign = sign.results)
    }, mc.cores = detectCores()-1)
    
    ## Organize results into a data frame
    results_array <- simplify2array(results)
    results.filtered <- c()
    for(x in 1:11) {
      results.filtered <- rbind(results.filtered,
                                data.frame(voxels = levels[x],
                                           map = map_names[map],
                                           study = std,
                                           comparator = comp,
                                           auc.unfiltered = as.vector(results_array[1,x])[[1]],
                                           auc.filtered = as.vector(results_array[2,x])[[1]],
                                           auc.sign = as.vector(results_array[3,x])[[1]]))
    }
    
    ## Add results from this map to the study results data frame
    study_results <- rbind(study_results, results.filtered)
  }
  
  ## Save results
  saveRDS(study_results, paste0("Thinning_Output/Thinning_Study_",std,"_Stm_",comp,".rds"))
}
