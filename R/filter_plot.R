library(ggplot2)

studies <- data.frame(studies = c(1,2,3,3,3,4),
                      comparator = c(2,2,2,3,4,2))

results <- c()

for(study in 1:nrow(studies)) {
  maps <- readRDS(paste0("~/Dropbox/Rami- Brain Signature/Andrew's analysis & results/Filtering/Filtering_Study_",
                         studies[study,"studies"],
                         "_Stm_",
                         studies[study,"comparator"],".rds"))
  maps$study <- paste0("Study_",studies[study,"studies"],"_Stm_",studies[study,"comparator"])
  
  results <- rbind(results, maps)
}

  
lim = c(0.4,1)

ggplot() +
  geom_point(data = results[results$level == 1,],
             aes(x = as.numeric(level)-1,
                 y = auc,
                 color = map),
             position = position_dodge(width = 2),
             size = 2
  ) +
  geom_errorbar(data = results[results$level == 1,],
                aes(x = as.numeric(level)-1,
                    ymin = lb,
                    ymax = ub,
                    color = map),
                position = position_dodge(width = 2),
                size = .3
  ) +
  geom_ribbon(data = results[-which(results$level %in% c(1,22)),],
              aes(ymin = lb,
                  ymax = ub,
                  x = as.numeric(level)-1,
                  fill = map),
              alpha = 0.2
  ) +
  geom_line(data = results[-which(results$level %in% c(1,22)),],
            aes(x=as.numeric(level)-1, 
                y=auc,
                color = map),
            size = 1) +
  geom_point(data = results[results$level == 22,],
             aes(x = as.numeric(level),
                 y = auc,
                 color = map),
             position = position_dodge(width = 2),
             size = 2
  ) +
  geom_errorbar(data = results[results$level == 22,],
                aes(x = as.numeric(level),
                    ymin = lb,
                    ymax = ub,
                    color = map),
                position = position_dodge(width = 2),
                size = .3
  ) +
  ylab("AUC") +
  xlab("Filtering (mm)") +
  coord_cartesian(ylim = c(0.4,1)) +
  facet_grid(study~map)

ggsave("~/Dropbox/Rami- Brain Signature/Andrew's analysis & results/Filtering/Filtering_ALL.pdf")



ggplot() +
  geom_point(data = results[results$level == 1,],
             aes(x = as.numeric(level)-1,
                 y = sgn.auc,
                 color = map),
             position = position_dodge(width = 2),
             size = 2
  ) +
  geom_errorbar(data = results[results$level == 1,],
                aes(x = as.numeric(level)-1,
                    ymin = sgn.lb,
                    ymax = sgn.ub,
                    color = map),
                position = position_dodge(width = 2),
                size = .3
  ) +
  geom_ribbon(data = results[-which(results$level %in% c(1,22)),],
              aes(ymin = sgn.lb,
                  ymax = sgn.ub,
                  x = as.numeric(level)-1,
                  fill = map),
              alpha = 0.2
  ) +
  geom_line(data = results[-which(results$level %in% c(1,22)),],
            aes(x=as.numeric(level)-1, 
                y=sgn.auc,
                color = map),
            size = 1) +
  geom_point(data = results[results$level == 22,],
             aes(x = as.numeric(level),
                 y = sgn.auc,
                 color = map),
             position = position_dodge(width = 2),
             size = 2
  ) +
  geom_errorbar(data = results[results$level == 22,],
                aes(x = as.numeric(level),
                    ymin = sgn.lb,
                    ymax = sgn.ub,
                    color = map),
                position = position_dodge(width = 2),
                size = .3
  ) +
  ylab("AUC") +
  xlab("Filtering (mm)") +
  coord_cartesian(ylim = c(0.4,1)) +
  facet_grid(study~map)

ggsave("~/Dropbox/Rami- Brain Signature/Andrew's analysis & results/Filtering/Sign_Filtering_ALL.pdf")
