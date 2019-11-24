library(ggplot2)

studies <- data.frame(studies = c(1,2,3,3,3,4),
                      comparator = c(2,2,2,3,4,2))

unfiltered_all <- c()
filtered_all <- c()
sign_all <- c()

for(study in 1:nrow(studies)) {
  maps <- readRDS(paste0("~/Dropbox/Rami- Brain Signature/Andrew's analysis & results/Thinning/Thinning_Study_",
                         studies[study,"studies"],
                         "_Stm_",
                         studies[study,"comparator"],
                         ".rds"))
  
  unfiltered <- aggregate(auc.unfiltered ~ voxels + map, maps, function(x) c(mean = mean(x), sd = sd(x)))
  unfiltered$auc <- unfiltered$auc.unfiltered[,1]
  unfiltered$std <- unfiltered$auc.unfiltered[,2]
  filtered <- aggregate(auc.filtered ~ voxels + map, maps, function(x) c(mean = mean(x), sd = sd(x)))
  filtered$auc <- filtered$auc.filtered[,1]
  filtered$std <- filtered$auc.filtered[,2]
  sgn <- aggregate(auc.sign ~ voxels + map, maps, function(x) c(mean = mean(x), sd = sd(x)))
  sgn$auc <- sgn$auc.sign[,1]
  sgn$std <- sgn$auc.sign[,2]
  
  unfiltered_all <- rbind(unfiltered_all,
                          cbind(unfiltered,
                                study = paste0("Study_",studies[study,"studies"],"_Stm_",studies[study,"comparator"])))
  filtered_all <- rbind(filtered_all,
                        cbind(filtered,
                              study = paste0("Study_",studies[study,"studies"],"_Stm_",studies[study,"comparator"])))
  sign_all <- rbind(sign_all,
                    cbind(sgn,
                          study = paste0("Study_",studies[study,"studies"],"_Stm_",studies[study,"comparator"])))
}

ggplot(data = unfiltered_all[unfiltered_all$std != 0,]) +
  geom_ribbon(aes(ymin = auc - std,
                  ymax = auc + std,
                  x = as.numeric(voxels),
                  fill = map),
              alpha = 0.2
  ) +
  geom_line(aes(x=as.numeric(voxels), 
                y=auc,
                color = map),
            size = 1) +
  coord_cartesian(ylim = c(0.4,1)) +
  ylab("AUC") +
  xlab("Number of voxels") +
  facet_grid(study~map)

ggsave(paste0("~/Dropbox/Rami- Brain Signature/Andrew's analysis & results/Thinning/Thinning_Unfiltered_ALL.pdf"))


ggplot(data = filtered_all[filtered_all$std != 0,]) +
  geom_ribbon(aes(ymin = auc - std,
                  ymax = auc + std,
                  x = as.numeric(voxels),
                  fill = map),
              alpha = 0.2
  ) +
  geom_line(aes(x=as.numeric(voxels), 
                y=auc,
                color = map),
            size = 1) +
  coord_cartesian(ylim = c(0.4,1)) +
  ylab("AUC") +
  xlab("Number of voxels") +
  facet_grid(study~map)

ggsave(paste0("~/Dropbox/Rami- Brain Signature/Andrew's analysis & results/Thinning/Thinning_Infinite_ALL.pdf"))

ggplot(data = sign_all[sign_all$std != 0,]) +
  geom_ribbon(aes(ymin = auc - std,
                  ymax = auc + std,
                  x = as.numeric(voxels),
                  fill = map),
              alpha = 0.2
  ) +
  geom_line(aes(x=as.numeric(voxels), 
                y=auc,
                color = map),
            size = 1) +
  coord_cartesian(ylim = c(0.4,1)) +
  ylab("AUC") +
  xlab("Number of voxels") +
  facet_grid(study~map)

ggsave(paste0("~/Dropbox/Rami- Brain Signature/Andrew's analysis & results/Thinning/Thinning_Sign_ALL.pdf"))
