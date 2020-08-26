rm(list = ls())
library(tidyverse)
library(seqinr)
library(here)

H_229E <- unlist(read.fasta(here("examples/common_corona/HCov-229E.fasta"), seqtype = "AA", as.string = TRUE, seqonly = TRUE))
H_OC43 <- unlist(read.fasta(here("examples/common_corona/HCov-OC43.fasta"), seqtype = "AA", as.string = TRUE, seqonly = TRUE))
H_HL63 <- unlist(read.fasta(here("examples/common_corona/HCov-HL63.fasta"), seqtype = "AA", as.string = TRUE, seqonly = TRUE))
H_HKU1 <- unlist(read.fasta(here("examples/common_corona/HCov-HKU1.fasta"), seqtype = "AA", as.string = TRUE, seqonly = TRUE))

in_229E <- function(pep){
  grepl(pep, H_229E)
}

peps_in_229E <- function(peps){
  sapply(peps, function(pep){
    s <- sum(in_229E(pep))
    s > 0.5
  })
}

in_OC43 <- function(pep){
  grepl(pep, H_OC43)
}

peps_in_OC43 <- function(peps){
  sapply(peps, function(pep){
    s <- sum(in_OC43(pep))
    s > 0.5
  })
}


in_HL63 <- function(pep){
  grepl(pep, H_HL63)
}

peps_in_HL63 <- function(peps){
  sapply(peps, function(pep){
    s <- sum(in_HL63(pep))
    s > 0.5
  })
}

in_HKU1 <- function(pep){
  grepl(pep, H_HKU1)
}

peps_in_HKU1 <- function(peps){
  sapply(peps, function(pep){
    s <- sum(in_HKU1(pep))
    s > 0.5
  })
}

merge_HLA <- function(df){
  paste(df$HLA, collapse = ";")
}
merge_pep <- function(df){
  paste(df$pep, collapse = "; ")
}

all_SARS2<- read_csv(here("examples/dataSARS2_allpeps/SARS2_prediction.csv"))

all_SARS2_immuno <-
  all_SARS2 %>%
  filter(immunogenic == "YES") %>%
  set_names(c("HLA", "pep", "amp", "immuno"))

length(unique(all_SARS2_immuno$pep))

all_SARS2_immuno_shared <- 
  all_SARS2_immuno %>%
  mutate(in229E = peps_in_229E(pep)) %>%
  mutate(inOC43 = peps_in_OC43(pep)) %>%
  mutate(inHL63 = peps_in_HL63(pep)) %>%
  mutate(inHKU1 = peps_in_HKU1(pep)) %>%
  mutate(total_match = in229E + inOC43 + inHL63 + inHKU1) %>%
  filter(total_match > 0.5)

length(unique(all_SARS2_immuno_shared$pep)) # still 46 common peptides
length(unique(all_SARS2_immuno_shared$HLA)) # now we got 31 HLAs

all_SARS2_immuno_shared_in229E <-
  all_SARS2_immuno_shared %>%
  filter(in229E)

length(unique(all_SARS2_immuno_shared_in229E$pep))  # 8 peptides 
length(unique(all_SARS2_immuno_shared_in229E$HLA))  # 24 HLAs
all_SARS2_immuno_shared_in229E %>%
  group_by(pep) %>%
  nest() %>%
  mutate(HLA_merged = map_chr(data, merge_HLA)) %>%
  select(pep, HLA_merged) %>%
  ungroup() %>%
  write.table(file = here("examples/common_corona/peps_in_229E_formatted.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# 90.21% World coverage, 93.6% US coverage

all_SARS2_immuno_shared_in229E %>%
  group_by(HLA) %>%
  nest() %>%
  mutate(n_pep = map_dbl(data, nrow)) %>%
  mutate(pep_merged = map_chr(data, merge_pep)) %>%
  select(HLA, pep_merged, n_pep) %>%
  ungroup() %>%
  summarise(mean(n_pep)) # 1.67

all_SARS2_immuno_shared_inHKU1 <-
  all_SARS2_immuno_shared %>%
  filter(inHKU1)

length(unique(all_SARS2_immuno_shared_inHKU1$pep))  # 31 peptides 
length(unique(all_SARS2_immuno_shared_inHKU1$HLA))  # 29 HLA
all_SARS2_immuno_shared_inHKU1 %>%
  group_by(pep) %>%
  nest() %>%
  mutate(HLA_merged = map_chr(data, merge_HLA)) %>%
  select(pep, HLA_merged) %>%
  ungroup() %>%
  write.table(file = here("examples/common_corona/peps_in_HKU1_formatted.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# 96.6 World Coverage, 97.8 US coverage

all_SARS2_immuno_shared_inHKU1 %>%
  group_by(HLA) %>%
  nest() %>%
  mutate(n_pep = map_dbl(data, nrow)) %>%
  mutate(pep_merged = map_chr(data, merge_pep)) %>%
  select(HLA, pep_merged, n_pep) %>%
  ungroup() %>%
  summarise(mean(n_pep)) # 5.24

all_SARS2_immuno_shared_inHL63 <-
  all_SARS2_immuno_shared %>%
  filter(inHL63)

length(unique(all_SARS2_immuno_shared_inHL63$pep) ) # 17 peptides 
length(unique(all_SARS2_immuno_shared_inHL63$HLA) ) # 28 
all_SARS2_immuno_shared_inHL63 %>%
  group_by(pep) %>%
  nest() %>%
  mutate(HLA_merged = map_chr(data, merge_HLA)) %>%
  select(pep, HLA_merged) %>%
  ungroup() %>%
  write.table(file = here("examples/common_corona/peps_in_HL63_formatted.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# 95.04 World Coverage, 96.7 US Coverage

all_SARS2_immuno_shared_inHL63 %>%
  group_by(HLA) %>%
  nest() %>%
  mutate(n_pep = map_dbl(data, nrow)) %>%
  mutate(pep_merged = map_chr(data, merge_pep)) %>%
  select(HLA, pep_merged, n_pep) %>%
  ungroup() %>%
  summarise(mean(n_pep)) # 2.86

all_SARS2_immuno_shared_inOC43 <-
  all_SARS2_immuno_shared %>%
  filter(inOC43)

length(unique(all_SARS2_immuno_shared_inOC43$pep))  # 20 peptides 
length(unique(all_SARS2_immuno_shared_inOC43$HLA))  # 25 HLA
all_SARS2_immuno_shared_inOC43 %>%
  group_by(pep) %>%
  nest() %>%
  mutate(HLA_merged = map_chr(data, merge_HLA)) %>%
  select(pep, HLA_merged) %>%
  ungroup() %>%
  write.table(file = here("examples/common_corona/peps_in_OC43_formatted.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# 96.32 World Coverage, 97.5 US Coverage

all_SARS2_immuno_shared_inOC43 %>%
  group_by(HLA) %>%
  nest() %>%
  mutate(n_pep = map_dbl(data, nrow)) %>%
  mutate(pep_merged = map_chr(data, merge_pep)) %>%
  select(HLA, pep_merged, n_pep) %>%
  ungroup() %>%
  summarise(mean(n_pep)) # 2.32
  
all_SARS2_immuno_shared %>%
  select(-c(total_match)) %>%
  write_csv(here("examples/common_corona/shared_with_common_corona.csv"))

all_SARS2_immuno_shared_reformatted <-
  all_SARS2_immuno_shared %>%
  group_by(pep) %>%
  nest() %>%
  mutate(HLA_merged = map_chr(data, merge_HLA)) %>%
  select(pep, HLA_merged) %>%
  ungroup()

write.table(x = all_SARS2_immuno_shared_reformatted, file = here("examples/common_corona/peps_shared_common_corona.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# HLA covers 98.6 world, 99.0 US

all_SARS2_immuno_shared_HLAgrouped <- 
  all_SARS2_immuno_shared %>%
  group_by(HLA) %>%
  nest() %>%
  mutate(n_pep = map_dbl(data, nrow)) %>%
  mutate(pep_merged = map_chr(data, merge_pep)) %>%
  select(HLA, pep_merged, n_pep) %>%
  ungroup() 

mean(all_SARS2_immuno_shared_HLAgrouped$n_pep) # 5.6



  
  
