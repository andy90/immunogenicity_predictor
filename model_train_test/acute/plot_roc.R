rm(list = ls())
library(matrixStats)
library(tidyverse)
dom_nondom_crits <- c(0.25)

for (dom_crit in dom_nondom_crits){
  df_tpr <- read.table(file = here("model_train_test/acute", paste("res_allmodels/tprs_full", dom_crit, ".txt", sep = "")), stringsAsFactors = FALSE)
  df_fpr <- read.table(file = here("model_train_test/acute", paste("res_allmodels/fpr_full", dom_crit, ".txt", sep = "")), stringsAsFactors = FALSE)
  
  tprs_full <- as.matrix(df_tpr)
  fpr <- df_fpr$V1
  
  tpr_mean <- rowMeans(tprs_full)
  tpr_mean <- c(0, tpr_mean,1)
  
  tpr_std <- rowSds(tprs_full)
  tpr_std <- c(0, tpr_std,1)
  

  fpr_ticks <- c(0, fpr,1)

  p <- ggplot(data = data.frame(x = fpr_ticks, y = tpr_mean))
  p <- p + geom_line(aes(x = x, y = y),size = 2,color= "red")
  p <- p + coord_equal()
  p <- p + geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed", size=0.5)
  p <- p + theme(text = element_text(size = 22))
  p <- p + scale_x_continuous(name = "False Positive Rate", limits = c(0,1), expand = c(0,0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"))
  p <- p + scale_y_continuous(name = "True Positive Rate", limits = c(0,1), expand = c(0,0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"))
  p
  ggsave(filename = here("model_train_test/acute", paste("res_allmodels/roc_full", dom_crit, ".pdf", sep = "")), width = 6.1, height = 6)
}