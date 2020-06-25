library(tidyverse)
library(seqinr)
library(here)

file_out <- here("src", "output.out")

system(paste(path_netmhcpan, file_fasta, "-a", HLAs, " >", file_out, sep = " " ))

peps <- toupper(unlist(read.fasta(file_fasta, as.string = TRUE)))

binders_MHC_binding <- data.frame()
for (pep in peps){
  all_content <-readLines(file_out)
  a_binding <- all_content[grep(paste("\\s+",pep,"\\s+",sep=""), all_content)]
  a_binding_table <- read.csv(textConnection(a_binding), stringsAsFactors = FALSE, header = FALSE, sep="")

  df_pep_binding <- a_binding_table[a_binding_table$V3 == pep, c(2,13)]
  df_pep_binding$pep <- pep
  
  binders_MHC_binding <- rbind(binders_MHC_binding, df_pep_binding)
}
colnames(binders_MHC_binding)[c(1,2)] <- c("HLA", "binding")


