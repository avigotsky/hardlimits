# Evaluate predictive performance of filtered templates
# Code by Andrew Vigotsky
# Last edited: November 23, 2019
# 
# For "The hard limits of decoding mental states" 
# by Jabakhanji et al.
#
# This assesses the predictive performance of each template
# at different filtering levels. It tests the performance
# of these filtered templates on studies 1 through 4.

library(RNifti)
library(pROC)
library(pbmcapply)
library(Matrix)
library(qlcMatrix)
library(tidyr)
library(boot)

setwd("/projects/b1090/Users/vigotsky/NPS_NSy")

## Load templates and convert them to sparse matrices
nps_filtered <- Matrix( apply(readNifti("Templates/NPS_filtered.nii"),  4, c), sparse = T)
pv_filtered  <- Matrix( apply(readNifti("Templates/PPV_filtered.nii"),  4, c), sparse = T)
nsy_filtered <- Matrix( apply(readNifti("Templates/NSyP_filtered.nii"), 4, c), sparse = T)
gpt_filtered <- Matrix( apply(readNifti("Templates/GPT_filtered.nii"),  4, c), sparse = T)

## Organize templates
maps <- list(nps_filtered, pv_filtered, nsy_filtered, gpt_filtered)
map_names <- c("NPS","pPV","pNsy","GPT")
map_candidates <- list(nps_candidates, pv_candidates, nsy_candidates, gpt_candidates)

studies <- data.frame(studies = c(1,2,3,3,3,4),
                      comparator = c(2,2,2,3,4,2))
levels <- 1:22

## Function to calculate the AUCs for the bootstrap
auc_calc <- function(df.wide, index) {
  df.tmp <- df.wide[index,]
  df <- gather(df.tmp, condition, dp, 2:3)
  auc(roc(condition ~ dp, df, direction = "<", quiet = T, algorithm = 3))[[1]]
}

for(study in 1:nrow(studies)) {
  print(study)
  
  ## Import beta maps
  stm1 <- readNifti(paste0("BetaMaps/Study_",studies[study,"studies"],"_Stm_1.nii"))
  stm2 <- readNifti(paste0("BetaMaps/Study_",studies[study,"studies"],"_Stm_",studies[study,"comparator"],".nii"))
  
  stm1_vec <- Matrix(apply(stm1, 4, c), sparse=T)
  stm2_vec <- Matrix(apply(stm2, 4, c), sparse=T)
  
  results <- c()
  for(i in 1:length(maps)) {
    print(i)
    map <- maps[[i]]
    
    ## Calculate dot products
    dps <- list(cosSparse(stm1_vec,map),
                cosSparse(stm2_vec,map))
    sgn_dps <- list(cosSparse(stm1_vec,sign(map)),
                    cosSparse(stm2_vec,sign(map)))
    
    set.seed(0)
    results_tmp <- pbmclapply(levels, function(level) {
      ## Create data frame for standard filtering levels
      df <- rbind(data.frame(id = 1:ncol(stm1_vec),
                             dp = dps[[1]][,level],
                             condition = 1),
                  data.frame(id = 1:ncol(stm2_vec),
                             dp = dps[[2]][,level],
                             condition = 0))
      
      ## Prepare and run bootstrap for standard filtering levels
      df.wide <- spread(df, condition, dp)
      boots <- boot(df.wide, auc_calc, 10000)
      cis <- boot.ci(boots, type = "bca")$bca[4:5]
      ## To ensure we don't run into errors with bootstrapping when AUC = 1 ∀ bootstraps
      if(is.null(cis) && auc(roc(condition ~ dp, df, direction = "<", quiet = T, algorithm = 3))[[1]] == 1)
        cis <- c(1,1)
      
      ## Create data frame for sign(template)
      sgn.df <- rbind(data.frame(id = 1:ncol(stm1_vec),
                                 dp = sgn_dps[[1]][,level],
                                 condition = 1),
                      data.frame(id = 1:ncol(stm2_vec),
                                 dp = sgn_dps[[2]][,level],
                                 condition = 0))
      
      ## Prepare and run bootstrap for sign(template)
      df.wide <- spread(sgn.df, condition, dp)
      boots <- boot(df.wide, auc_calc, 10000)
      sgn.cis <- boot.ci(boots, type = "bca")$bca[4:5]
      ## To ensure we don't run into errors with bootstrapping when AUC = 1 ∀ bootstraps
      if(is.null(sgn.cis) && auc(roc(condition ~ dp, sgn.df, direction = "<", quiet = T, algorithm = 3))[[1]] == 1)
        sgn.cis <- c(1,1)
      
      
      ## Final data frame for this filtering level, template, and study
      data.frame(auc = auc(roc(condition ~ dp, df, direction = "<", quiet = T, algorithm = 3))[[1]],
                 sgn.auc = auc(roc(condition ~ dp, sgn.df, direction = "<", quiet = T, algorithm = 3))[[1]],
                 level = level,
                 map = map_names[i],
                 lb = cis[1],
                 ub = cis[2],
                 sgn.lb = sgn.cis[1],
                 sgn.ub = sgn.cis[2])
    }, mc.cores = detectCores()-1)
    
    ## Add template to study results
    results <- rbind(results,
                     do.call(rbind.data.frame,results_tmp))
  }
  
  ## Save study results
  saveRDS(results, paste0("Filtering_Output/Filtering_Study_",
                          studies[study,"studies"],"_Stm_",studies[study,"comparator"],".rds"))
}