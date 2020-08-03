# this file is used to run the computational tool to train and test our model and compare it to models
# with one or two of the three terms missing
library(tidyverse)
library(seqinr)
library(here)
library(readxl)
Acute_data_file <-here("model_train_test/acute/acute_HLA_pep_response.xlsx")
peptides_MHC <- read_xlsx(Acute_data_file)
colnames(peptides_MHC) <- c("HLA", "peptide_name", "peptide", "response")

binders_MHC_binding <- c() # store the binding data for all the peptides
# run netMHCpan to get the binding of peptide and MHC
for (i in 1:nrow(peptides_MHC)){
  pep <- peptides_MHC$peptide[i]
  HLAs <- paste("HLA-", peptides_MHC$HLA[i], sep = "")
  
  file_fasta <- here("model_train_test/acute", paste(pep, ".fasta", sep = ""))
  write.fasta(sequences = pep, names = 1, file.out = file_fasta)
  
  path_netmhcpan <- "/Users/anggao/Documents/netMHCpan-4.0/netMHCpan" # this should be replaced by the path of your own installation of netMHCpan4.0
  
  source(here("model_train_test/acute/binding_prediction.R"))
}
peptides_MHC$binding <- binders_MHC_binding

# align the peptides with self and foreign peptides
source(here("model_train_test/acute/get_similarity.R"))


source(here("model_train_test/acute/train_test_all_models.R"))

source(here("model_train_test/acute/plot_auc_corr_comparison.R"))

source(here("model_train_test/acute/plot_roc.R"))

