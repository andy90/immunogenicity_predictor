
path_netmhcpan <- "/Users/anggao/Documents/netMHCpan-4.0/netMHCpan" # this should be replaced by the path of your own installation of netMHCpan4.0

HLAs <- "HLA-A02:01,HLA-A03:01" # HLA alleles that presents the peptides. this should be replaced by the HLAs you are interested; HLAs should be separated by comma, no space allowed.

file_fasta <- here("input.fasta") # the input fasta file that contains the peptides which one want to predict immunogenicity
final_output_file <- here("predicted_immunogenicity.csv") # the output file that contains the immunogenicity predicition for the peptides

source(here("src", "binding_prediction.R"))
source(here("src", "get_similarity.R"))
source(here("src", "predict_amp.R"))

# run main.R