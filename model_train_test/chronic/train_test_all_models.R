# train and test all the models: full model and the models with 1 or 2 missing terms
get_auc <- function(experiment, prediction){
  crit0 <- c(0, sort(prediction))
  crit1 <- c(crit0[2:length(crit0)], crit0[length(crit0)]+1)
  crit_seq <- (crit0 + crit1)/2
  
  df_TP_FP <- numeric()
  for (i_crit in crit_seq){ # it is important to choose a good span of i_crit
    TP <-  sum((experiment >= dom_crit) & (prediction >= i_crit))
    TN <-  sum((experiment < nondom_crit) & (prediction < i_crit))
    FP <-  sum((experiment < nondom_crit) & (prediction >= i_crit))
    FN <-  sum((experiment >= dom_crit) & (prediction < i_crit))
    df_TP_FP <- rbind(df_TP_FP, c(TP/(TP + FN), FP/(FP + TN)))
  }
  TPrate <- df_TP_FP[,1]
  FPrate <- df_TP_FP[,2]
  auc <- -sum((TPrate[-1] + TPrate[-length(TPrate)])/2 * diff(FPrate))
  
  list(auc, df_TP_FP)
}

f_train_full <- function(a, b, foreign, self, binding){
  colSums(foreign[(a+1):51, ])*colSums(self[(b+1):51, ])/binding
}

f_train_sf <- function(a, b, foreign, self, binding){
  colSums(foreign[(a+1):51, ])*colSums(self[(b+1):51, ])
}

f_train_sb <- function(a, b, foreign, self, binding){
  colSums(self[(b+1):51, ])/binding 
}

f_train_fb <- function(a, b, foreign, self, binding){
  colSums(foreign[(a+1):51, ])/binding 
}

f_train_s <- function(a, b, foreign, self, binding){
  colSums(self[(b+1):51, ]) 
}

f_train_f <- function(a, b, foreign, self, binding){
  colSums(foreign[(a+1):51, ]) 
}

test_acc <- function(training_indicies, test_indicies, model){
  
  if (model == "full"){
    corr_a_training <- sapply(a_seq, function(a){
      sapply(b_seq, function(b){
        cor(x = all_response_percent[training_indicies], y = f_train_full(a , b, all_foreignsimilarity_topN[, training_indicies], all_selfsimilarity_topN[, training_indicies], all_binding[training_indicies]), method = "spearman")
      })
    })
  }else if (model == "sf"){
    corr_a_training <- sapply(a_seq, function(a){
      sapply(b_seq, function(b){
        cor(x = all_response_percent[training_indicies], y = f_train_sf(a , b, all_foreignsimilarity_topN[, training_indicies], all_selfsimilarity_topN[, training_indicies], all_binding[training_indicies]), method = "spearman")
      })
    })
  }else if (model == "sb"){
    corr_a_training <- sapply(a_seq, function(a){
      sapply(b_seq, function(b){
        cor(x = all_response_percent[training_indicies], y = f_train_sb(a , b, all_foreignsimilarity_topN[, training_indicies], all_selfsimilarity_topN[, training_indicies], all_binding[training_indicies]), method = "spearman")
      })
    })
  }else if (model == "fb"){
    corr_a_training <- sapply(a_seq, function(a){
      sapply(b_seq, function(b){
        cor(x = all_response_percent[training_indicies], y = f_train_fb(a , b, all_foreignsimilarity_topN[, training_indicies], all_selfsimilarity_topN[, training_indicies], all_binding[training_indicies]), method = "spearman")
      })
    })
  }else if (model == "s"){
    corr_a_training <- sapply(a_seq, function(a){
      sapply(b_seq, function(b){
        cor(x = all_response_percent[training_indicies], y = f_train_s(a , b, all_foreignsimilarity_topN[, training_indicies], all_selfsimilarity_topN[, training_indicies], all_binding[training_indicies]), method = "spearman")
      })
    })
  }else if (model == "f"){
    corr_a_training <- sapply(a_seq, function(a){
      sapply(b_seq, function(b){
        cor(x = all_response_percent[training_indicies], y = f_train_f(a , b, all_foreignsimilarity_topN[, training_indicies], all_selfsimilarity_topN[, training_indicies], all_binding[training_indicies]), method = "spearman")
      })
    })
  }
  
  if (model != "b"){
    ind_max <- which(corr_a_training == max(corr_a_training, na.rm = TRUE), arr.ind = TRUE)[1,]
    b_best <- b_seq[ind_max[1]]
    a_best <- a_seq[ind_max[2]]
  }
  
  
  if (model == "full"){
    test_mag <- f_train_full(a_best, b_best,  all_foreignsimilarity_topN[, test_indicies], all_selfsimilarity_topN[, test_indicies], all_binding[test_indicies])
  }else if (model == "sf"){
    test_mag <- f_train_sf(a_best, b_best,  all_foreignsimilarity_topN[, test_indicies], all_selfsimilarity_topN[, test_indicies], all_binding[test_indicies])
  }else if (model == "sb"){
    test_mag <- f_train_sb(a_best, b_best,  all_foreignsimilarity_topN[, test_indicies], all_selfsimilarity_topN[, test_indicies], all_binding[test_indicies])
  }else if (model == "fb"){
    test_mag <- f_train_fb(a_best, b_best,  all_foreignsimilarity_topN[, test_indicies], all_selfsimilarity_topN[, test_indicies], all_binding[test_indicies])
  }else if (model == "s"){
    test_mag <- f_train_s(a_best, b_best,  all_foreignsimilarity_topN[, test_indicies], all_selfsimilarity_topN[, test_indicies], all_binding[test_indicies])
  }else if (model == "f"){
    test_mag <- f_train_f(a_best, b_best,  all_foreignsimilarity_topN[, test_indicies], all_selfsimilarity_topN[, test_indicies], all_binding[test_indicies])
  }else if (model == "b"){
    test_mag <- 1/all_binding[test_indicies]
  }
  
  data_auc_roc <- get_auc(all_response_percent[test_indicies], test_mag)
  auc <- data_auc_roc[[1]]
  
  
  roc <- data_auc_roc[[2]]
  tpr_new <- approx(x = roc[,2],  y = roc[,1], xout = fpr, method = "linear")$y
  roc_new <- cbind(tpr_new, fpr)
  
  corr_a_test <- cor(x = all_response_percent[test_indicies], y = test_mag, method = "spearman")
  if (model != "b"){
    res_final <- list(corr_a_test, auc, roc_new, a_best, b_best)
  }else{
    res_final <- list(corr_a_test, auc, roc_new)
  }
  res_final
}

dom_crit <- 0.25 # the criteria for labeling positive and negative data
nondom_crit <- 0.25

all_binding <- peptides_MHC$binding
all_response_percent <- peptides_MHC$response
all_selfsimilarity_topN <- selfsimilarity_topN
all_foreignsimilarity_topN <- foreignsimilarity_topN

a_seq <- seq(from = 1, to = 25, by = 1)
b_seq <- seq(from = 1, to = 37, by = 1)

nfold <- 10 # split the data into nfold to perform 10 fold cross validation

corrs <- numeric()
aucs <- numeric()
tprs_full <- numeric()
as <- numeric()
bs <- numeric()
fpr <- seq(0, 1, 0.01) 
nrepeat <- 20
for (irepeat in 1:nrepeat){
  folds <- cut(seq(1,length(all_response_percent)),breaks=nfold,labels=FALSE) # cut will separate the folds as evenly as possible
  folds <- sample(folds)
  for(i in 1:nfold){
    test_indicies <- which(folds==i,arr.ind=TRUE)
    training_indicies <-setdiff(1:length(all_response_percent), test_indicies)
    
    df_acc_full <- test_acc(training_indicies, test_indicies, "full")
    df_acc_sf <- test_acc(training_indicies, test_indicies, "sf")
    df_acc_sb <- test_acc(training_indicies, test_indicies, "sb")
    df_acc_fb <- test_acc(training_indicies, test_indicies, "fb")
    df_acc_s <- test_acc(training_indicies, test_indicies, "s")
    df_acc_f <- test_acc(training_indicies, test_indicies, "f")
    df_acc_b <- test_acc(training_indicies, test_indicies, "b")
    
    corr_test_full <- df_acc_full[[1]]
    corr_test_sf <- df_acc_sf[[1]]
    corr_test_sb <- df_acc_sb[[1]]
    corr_test_fb <- df_acc_fb[[1]]
    corr_test_s <- df_acc_s[[1]]
    corr_test_f <- df_acc_f[[1]]
    corr_test_b <- df_acc_b[[1]]
    
    auc_full <- df_acc_full[[2]]
    auc_sf <- df_acc_sf[[2]]
    auc_sb <- df_acc_sb[[2]]
    auc_fb <- df_acc_fb[[2]]
    auc_s <- df_acc_s[[2]]
    auc_f <- df_acc_f[[2]]
    auc_b <- df_acc_b[[2]]
    
    tpr_full <- df_acc_full[[3]][,1]
    
    a_best_full <- df_acc_full[[4]]
    b_best_full <- df_acc_full[[5]]
    
    corrs <- rbind(corrs, c(corr_test_full, corr_test_sf, corr_test_sb, corr_test_fb, corr_test_s, corr_test_f, corr_test_b))
    aucs <- rbind(aucs, c(auc_full, auc_sf, auc_sb, auc_fb, auc_s, auc_f, auc_b))
    tprs_full <- cbind(tprs_full, tpr_full)
    as <- c(as, a_best_full)
    bs <- c(bs, b_best_full)
  }
}

colnames(corrs) <- c("full", "self_foreign", "self_binding", "foreign_binding", "self", "foreign", "binding")
corrs_mean <- matrix(colMeans(corrs), nrow = 1)
colnames(corrs_mean) <-  c("full", "self_foreign", "self_binding", "foreign_binding", "self", "foreign", "binding")

write.table(x = corrs, row.names = FALSE, quote = FALSE, file = here("model_train_test/chronic/res_allmodels","corrs.txt"))
write.table(x = corrs_mean, row.names = FALSE, quote = FALSE, file = here("model_train_test/chronic/res_allmodels","corrs_mean.txt"))

colnames(aucs) <- c("full", "self_foreign", "self_binding", "foreign_binding", "self", "foreign", "binding")
aucs_mean <- matrix(colMeans(aucs), nrow = 1)
colnames(aucs_mean) <- c("full", "self_foreign", "self_binding", "foreign_binding", "self", "foreign", "binding")

write.table(x = aucs, row.names = FALSE, quote = FALSE, file = here("model_train_test/chronic", paste("res_allmodels/aucs", dom_crit, ".txt", sep = "")))
write.table(x = aucs_mean, row.names = FALSE, quote = FALSE, file = here("model_train_test/chronic", paste("res_allmodels/aucs_mean", dom_crit, ".txt", sep = "")))

write.table(x = tprs_full, row.names = FALSE, col.names = FALSE, quote = FALSE, file = here("model_train_test/chronic", paste("res_allmodels/tprs_full", dom_crit, ".txt", sep = "")))

write.table(data.frame(fpr), row.names = FALSE, col.names = FALSE, file = here("model_train_test/chronic", paste("res_allmodels/fpr_full", dom_crit, ".txt", sep = "")))

write.table(x = cbind(as, bs), row.names = FALSE, col.names = FALSE, quote = FALSE, file = here("model_train_test/chronic", paste("res_allmodels/ab_full", dom_crit, ".txt", sep = "")))

