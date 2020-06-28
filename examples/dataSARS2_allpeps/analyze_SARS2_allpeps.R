# this is to apply the computational tool to all the possible 8-11mers of SARS-Cov-2
library(tidyverse)
library(seqinr)
library(here)

path_netmhcpan <- "/Users/anggao/Documents/netMHCpan-4.0/netMHCpan" # this should be replaced by the path of your own installation of netMHCpan4.0

HLAs <- read.table(here("examples", "dataSARS2_allpeps", "HLAs.txt"), stringsAsFactors = FALSE)$V1 # a list of 38 HLAs that can cover most population in the world
HLAs <- paste(HLAs, collapse = ",")
file_fasta <- here("examples", "dataSARS2_allpeps", "SARS2_all_peps.fasta") # the input fasta file that contains the peptides of SARS2
final_output_file <- here("examples", "dataSARS2_allpeps", "SARS2_prediction.csv") # the output file that contains the immunogenicity predicition for the peptides

source(here("src", "binding_prediction.R"))
source(here("src", "get_similarity.R"))
source(here("src", "predict_amp.R"))