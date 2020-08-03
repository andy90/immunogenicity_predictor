library("Biostrings")
data("BLOSUM62")


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
foreign_peps <- read.table(file = here("data", "foreign_peps_noHIV.txt"), stringsAsFactors = FALSE)

selfsimilarity_topN <- apply(check_similarity(peptides_MHC$peptide, self_peps$V1), 2, get_count)
foreignsimilarity_topN <- apply(check_similarity(peptides_MHC$peptide, foreign_peps$V1), 2, get_count)

