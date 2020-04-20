amp_thresh <- 7024 # this is the point where the crit is 25%, and TPR and TNR are the same

binders_MHC_binding_self_foreign_amp <-
  binders_MHC_binding_self_foreign %>%
  mutate(amp = self*foreign/binding) %>%
  mutate(immunogenic = if_else(amp > amp_thresh, true = "YES", false = "NO"))

data_final <- 
  binders_MHC_binding_self_foreign_amp %>%
  select(HLA, pep, amp, immunogenic) %>%
  set_names(c("HLA", "peptide", "amplitude", "immunogenic")) %>%
  arrange(desc(amplitude) )

write.csv(data_final, file = here("predicted_immunogenicity.csv"), quote = FALSE, row.names = FALSE)

system(paste("rm",  file_out, sep = " " ))


