library(tidyverse)
library(seqinr)
library(here)

file_out <- here("model_train_test/chronic", "output.out")

system(paste(path_netmhcpan, file_fasta, "-a", HLAs, "-l 8,9,10,11,12 >", file_out, sep = " " ))

peps <- toupper(unlist(read.fasta(file_fasta, as.string = TRUE)))

for (pep in peps){
  all_content <-readLines(file_out)
  a_binding <- all_content[grep(paste("\\s+",pep,"\\s+",sep=""), all_content)]
  a_binding_table <- read.csv(textConnection(a_binding), stringsAsFactors = FALSE, header = FALSE, sep="")

  df_pep_binding <- a_binding_table[a_binding_table$V3 == pep, 13][1]

  binders_MHC_binding <- c(binders_MHC_binding, df_pep_binding)
}


