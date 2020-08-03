library(dplyr)
library(ggplot2)
library(tidyr)


corrs <- read.table(file = here("model_train_test/chronic","res_allmodels/corrs.txt"), header = TRUE)
corrs_transformed <- gather(corrs)
model_names <- c("full", "self_foreign", "self_binding", "foreign_binding", "self", "foreign", "binding")
label_names <- c("full", "self\n&\nforeign", "self\n&\nbinding", "foreign\n&\nbinding", "self", "foreign", "binding")

corrs_transformed$key <- factor(corrs_transformed$key, levels = model_names)
ggplot(corrs_transformed, aes(key, value))+
  stat_summary(geom = "bar", fun.y = mean, position = "dodge", width=0.7)+
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=0.15, size = 1.5)+
  scale_x_discrete(name = NULL, breaks=model_names, labels=label_names)+
  scale_y_continuous(name = "Spearman Correlation Coefficient", expand = c(0,0))+
  coord_cartesian(ylim = c(0, 0.38))+ 
  theme(text = element_text(size=15))
ggsave(here("model_train_test/chronic", "res_allmodels/Comp_spearman_chronic.pdf"), width = 5, height = 5)

aucs <- read.table(file = here("model_train_test/chronic", "res_allmodels/aucs0.25.txt"), header = TRUE)
aucs_transformed <- gather(aucs)
aucs_transformed$key <- factor(aucs_transformed$key, levels = model_names)
ggplot(aucs_transformed, aes(key, value))+
  stat_summary(geom = "bar", fun.y = mean, position = "dodge", width=0.7)+
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=0.15, size = 1.5)+
  scale_x_discrete(name = NULL, breaks=model_names, labels=label_names)+
  scale_y_continuous(name = "AUC", expand = c(0,0))+
  coord_cartesian(ylim = c(0.5, 0.73))+ 
#  coord_cartesian(ylim = c(0.5, 0.86))+
  theme(text = element_text(size=15))
ggsave(here("model_train_test/chronic", "res_allmodels/Comp_auc_chronic0.25.pdf"), width = 5, height = 5)
