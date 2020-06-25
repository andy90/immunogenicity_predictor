library("Biostrings")
data("BLOSUM62")
library(doParallel)


##############
# function used to perform local alignment
##############
blosum_align <- function(pep1, pep2){
  res <- pairwiseAlignment(pep1, pep2, type = "local", substitutionMatrix = BLOSUM62,
                           gapOpening = 11, gapExtension = 1) # see wikipedia for the detail of the Smith-Waterman alignment
  res
}

score_blosum_align <- function(pep1, pep2){
  res <- pairwiseAlignment(pep1, pep2, type = "local", substitutionMatrix = BLOSUM62,
                           gapOpening = 11, gapExtension = 1) # see wikipedia for the detail of the Smith-Waterman alignment
  res@score
}


check_similarity <- function(peps, self_peps){
  sapply(peps, function(pep){
    score_blosum_align(self_peps, pep)
  }) 
}

##################
# check the alignment score for the streeck peptides
#################

similarity_seq <- seq(from = 0, to = 50, by = 1)
get_count <- function(v){
  sapply(similarity_seq, function(x){
    sum(v == x)
  })
}

######################
#get the similarity to self and foreign peptides
######################
self_peps <- read.table(file = here("data", "self_peps.txt"), stringsAsFactors = FALSE)
get_self_similarity <- function(seqs){
  self_similarity_topN <- apply(check_similarity(seqs, self_peps$V1), 2, get_count)
  
  if (ncol(self_similarity_topN) > 1){
    Rself <- colSums(self_similarity_topN[26:51, ] ) # 25 is the threshold for self
  }else{
    Rself <- sum(self_similarity_topN[26:51 ] ) # 25 is the threshold for self
  }
  
}

foreign_peps <- read.table(file = here("data", "foreign_peps_noHIV.txt"), stringsAsFactors = FALSE)
get_foreign_similarity <- function(seqs){
  foreign_similarity_topN <- apply(check_similarity(seqs, foreign_peps$V1), 2, get_count)
  if (ncol(foreign_similarity_topN) > 1){
    Rforeign <- colSums(foreign_similarity_topN[18:51, ] ) # 17 is the threshold for foreign
  }else{
    Rforeign <- sum(foreign_similarity_topN[18:51] ) # 17 is the threshold for foreign
  }
}

##########
# the parallel version of get_self/foreign_similarity
###########
registerDoParallel(cores=2)
get_self_similarity_par <- function(seqs){
  seqs_splitted <- split(seqs ,ceiling(seq_along(seqs )/20))
  seqs_self <- 
    foreach(i = 1:length(seqs_splitted), .combine=c) %dopar% {
      get_self_similarity(seqs_splitted[[i]])  
    }
}

get_foreign_similarity_par <- function(seqs){
  seqs_splitted <- split(seqs ,ceiling(seq_along(seqs )/20))
  seqs_foreign <- 
    foreach(i = 1:length(seqs_splitted), .combine=c) %dopar% {
      get_foreign_similarity(seqs_splitted[[i]])  
    }
}

#################
# get the similarity for sette nCov peptides
#################

binders_MHC_binding_self_foreign <- 
  binders_MHC_binding %>%
  mutate("self" = get_self_similarity_par(pep)) %>%
  mutate("foreign" = get_foreign_similarity_par(pep))
