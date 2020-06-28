# this file is used to run the computational tool to predict the immunogenicity
# of Ahmed peptides
library(tidyverse)
library(seqinr)
library(here)

Ahmed_data_file <-here("examples", "dataAhmed", "Ahmed_peptides_MHC.txt")
peptides_MHC <- read.table(Ahmed_data_file)
colnames(peptides_MHC) <- c("HLA", "peptide")

for (i in 1:nrow(peptides_MHC)){
  pep <- peptides_MHC$peptide[i]
  HLAs <- peptides_MHC$HLA[i]
  
  file_fasta <- here("examples", "dataAhmed", paste(pep, ".fasta", sep = ""))
  write.fasta(sequences = pep, names = 1, file.out = file_fasta)
  
  path_netmhcpan <- "/Users/anggao/Documents/netMHCpan-4.0/netMHCpan" # this should be replaced by the path of your own installation of netMHCpan4.0
  
  final_output_file <- here("examples", "dataAhmed", paste(pep, ".csv", sep = ""))
  source(here("src", "binding_prediction.R"))
  source(here("src", "get_similarity.R"))
  source(here("src", "predict_amp.R"))
}
