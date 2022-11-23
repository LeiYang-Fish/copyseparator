#' @title copy_assemble
#'
#' @description Assembles a small number of overlapping DNA sequences into their respective gene copies.
#'
#' @param filename A fasta alignment of a small number of overlapping DNA sequences (results from "copy_separate") covering the entire length of the target gene. Check the alignment carefully before proceeding.
#'
#' @param copy_number An integer (e.g. 2,3, or 4) giving the anticipated number of gene copies. Must be the same value as used for "copy_separate".
#'
#' @param verbose Turn on (verbose=1; default) or turn off (verbose=0) the output.
#'
#' @return A fasta alignment of the anticipated number of full-length gene copies.
#'
#' @importFrom seqinr read.fasta write.fasta
#'
#' @importFrom stringr str_count str_order
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet width
#'
#' @importFrom DECIPHER ConsensusSequence RemoveGaps
#'
#' @importFrom beepr beep
#'
#' @examples
#' \dontrun{
#' copy_assemble("inst/extdata/combined_con.fasta",2,1)
#' }
#'
#' @export copy_assemble
#'

copy_assemble<-function(filename,copy_number, verbose=1)
{
filename_short <- gsub("_copies.*","", gsub("[:.:].*","", filename))

sink(paste0(filename_short, "_assemble_log.txt"), append=FALSE, split=TRUE) # begin to record log
error_log_function <- function() {
  cat(geterrmessage(), file="Error_log.txt", append=T)
}

if (copy_number<=1) stop ("The anticipated copy number must be a number larger than one!")

# Define a function to count the total number of ambiguous sites in an alignment
ambiguity_count <- function(x){
  sum(as.integer(stringr::str_count(as.character(DECIPHER::ConsensusSequence(Biostrings::readDNAStringSet(x, format="fasta",nrec=-1L, skip=0L),
                                    threshold = 0.4,ambiguity = TRUE, minInformation=0.8, noConsensusChar = "N")), c("M", "K", "R", "Y","N","W", "S", "H", "V", "D", "B"))))
  }

Consensus_seq_org <- seqinr::read.fasta(file = filename,seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)
subset_num_org <- as.numeric(gsub("\\D", "", gsub("_cluster_.", "", names(Consensus_seq_org))))
subset_sum <- length(unique(subset_num_org))
cluster_list_org <- lapply(sort(unique(subset_num_org)), function(x) which(subset_num_org==x))

# List missing subsets when exist
if (max(subset_num_org)!=subset_sum) {
  if (verbose) {cat(paste0("Note: Subset ", setdiff(1:max(subset_num_org), subset_num_org), " is missing. Subset numbering will be changed.","\n"),sep="")}
}

# Remove ambiguous region from the end of some sequences
Consensus_seq_clean <- list()
for (i in (1:length(Consensus_seq_org))){
  ambiguity_each <- sum(stringr::str_count(stringr::str_sub(Consensus_seq_org[i],-70,-1), c("M", "K", "R", "Y","N","W", "S", "H", "V", "D", "B")))
  if (ambiguity_each<10){
    Consensus_seq_clean <- append(Consensus_seq_clean, Consensus_seq_org[i])
  } else {
    Consensus_seq_clean <- append(Consensus_seq_clean, stringr::str_sub(Consensus_seq_org[i],1, -70))
  }
}
seqinr::write.fasta(sequences = Consensus_seq_clean, names(Consensus_seq_org),file.out="Consensus_seq_clean.fasta")      # saving while re-naming
Consensus_seq_clean2 <- seqinr::read.fasta(file = "Consensus_seq_clean.fasta",seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)
unlink("*_clean.fasta")

# Re-arrange to make sure that all subsets and clusters are arranged in order
Consensus_seq <- list()
for (i in 1:subset_sum){
  Consensus_seq <- append(Consensus_seq, Consensus_seq_clean2[sort(cluster_list_org[[i]])])
}

#seqinr::write.fasta(sequences = Consensus_seq, names(Consensus_seq),file.out="Consensus_seq_.fasta")

subset_names <- gsub("_cluster_.", "", gsub("_consensus.fasta","", names(Consensus_seq)))
cluster_list <- lapply(unique(subset_names), function(x) which(subset_names==x))
cluster_equal <- which(lapply(1:subset_sum, function (x) length(cluster_list[[x]]))==copy_number)  # subsets that copy_separate found = copy_number clusters
cluster_over <- which(lapply(1:subset_sum, function (x) length(cluster_list[[x]]))>copy_number) # subsets that copy_separate found > copy_number clusters

if (verbose) { cat("*************************************************************************\n")}


# If every subset have copy_number of clusters, no need to reduce the dataset
if (length(cluster_equal)==subset_sum) {Consensus_seq_reduced <- Consensus_seq}

# If a subset has > copy_number of clusters/consensus sequences, keep only the copy_number of different sequences
if (length(cluster_equal)<subset_sum && copy_number<=3){

 second_cp <- character(0)
 for (i in cluster_over){
  for (j in 2:length(cluster_list[[i]])){
    seqinr::write.fasta(sequences = as.list(c(stringr::str_sub(Consensus_seq[cluster_list[[i]][1]][[1]],1,-70), stringr::str_sub(Consensus_seq[cluster_list[[i]][j]][[1]],1,-70))),
                        names(c(Consensus_seq[cluster_list[[i]][1]],Consensus_seq[cluster_list[[i]][j]])), file.out = paste0("xnc_",i,"_compare_1_&_",j,".fasta"))
    ambiguity_among <- ambiguity_count(paste0("xnc_",i,"_compare_1_&_",j,".fasta"))
    if (ambiguity_among>4) break
  }
    cat("Searching the second copy for the subset ", i, "\n")
    cat(paste0("--- The sequence ", j, " may be the second copy\n"))
    cat("Number of ambiguous sites between the picked 2nd copy and the 1st copy:", ambiguity_among, "\n")
    second_cp <- append(second_cp, j)
 }
 unlink(list.files(pattern="xnc_"))

 picked_2nd <- character(0)
 for (k in 1:length(cluster_over)){picked_2nd <- append(picked_2nd, c(cluster_list[cluster_over][[k]][1], cluster_list[cluster_over][[k]][as.numeric(second_cp[k])]))}
 list_reduced <- sort(as.numeric(c(unlist(cluster_list[cluster_equal]),picked_2nd)))


 if (copy_number==3){
  third_cp <- character(0)
  for (i in cluster_over){
    which_2nd <- as.numeric(second_cp[which(cluster_over==i)]) # the 2nd copy (may not be the 2nd sequence of the subset)
    #if (which_2nd==length(cluster_list[[i]])){cat("Failed to find the 3rd copy for the subset ", i, "\n")}
    for (m in (which_2nd+1):length(cluster_list[[i]])){
      # comparing with 1st copy
      seqinr::write.fasta(sequences = as.list(c(stringr::str_sub(Consensus_seq[cluster_list[[i]][1]][[1]],1,-70), stringr::str_sub(Consensus_seq[cluster_list[[i]][m]][[1]],1,-70))),
                          names(c(Consensus_seq[cluster_list[[i]][1]],Consensus_seq[cluster_list[[i]][j]])), file.out = paste0("x12nc_",i,"_compare_1_&_",m,".fasta"))
      ambiguity_among13 <- ambiguity_count(paste0("x12nc_",i,"_compare_1_&_",m,".fasta"))
      # comparing with 2nd copy
      seqinr::write.fasta(sequences = as.list(c(stringr::str_sub(Consensus_seq[cluster_list[[i]][which_2nd]][[1]],1,-70), stringr::str_sub(Consensus_seq[cluster_list[[i]][m]][[1]],1,-70))),
                          names(c(Consensus_seq[cluster_list[[i]][1]],Consensus_seq[cluster_list[[i]][j]])), file.out = paste0("x23nc_",i,"_compare_1_&_",m,".fasta"))
      ambiguity_among23 <- ambiguity_count(paste0("x23nc_",i,"_compare_1_&_",m,".fasta"))
     if (ambiguity_among13>4 && ambiguity_among23>4) break
    }
    cat("Searching the third copy for the subset ", i, "\n")
    cat(paste0("--- The sequence ", m, " may be the third copy\n"))
    cat("Number of ambiguous sites between the picked 3rd copy and the 1st copy:", ambiguity_among13, "\n")
    cat("Number of ambiguous sites between the picked 3rd copy and the 2nd copy:", ambiguity_among23, "\n")
    third_cp <- append(third_cp, m)
  }
  unlink(list.files(pattern="x12nc_"))
  unlink(list.files(pattern="x23nc_"))

  picked_3rd <- character(0)
  for (n in 1:length(cluster_over)) {picked_3rd <- append(picked_3rd, cluster_list[cluster_over][[n]][as.numeric(third_cp[n])])}
  list_reduced <- sort(as.numeric(c(unlist(cluster_list[cluster_equal]),picked_2nd,picked_3rd)))
 }

   Consensus_seq_reduced <- Consensus_seq[list_reduced]
   seqinr::write.fasta(sequences = Consensus_seq_reduced, names(Consensus_seq_reduced),file.out = paste0(filename_short,"_reduced.txt"))
}


# When copy_number >3, it is strongly encouraged to do the assembling manually. Below we pick the copy_number of largest clusters to represent each gene copy.
if (length(cluster_equal)<subset_sum && copy_number>3){
   picked_copies <- unlist(lapply(1:length(cluster_over), function(y) {lapply(1:copy_number, function(x) cluster_list[cluster_over][[y]][x])}))
   list_reduced <- sort(as.numeric(c(unlist(cluster_list[cluster_equal]),picked_copies)))
   Consensus_seq_reduced <- Consensus_seq[list_reduced]
   seqinr::write.fasta(sequences = Consensus_seq_reduced, names(Consensus_seq_reduced),file.out = paste0(filename_short,"_reduced.txt"))
}



# Codes below will assemble copy_number of consensus sequences from each subset into copy_number of full length gene copies

if (verbose) { cat("*************************************************************************\n")}

seq_to_con <- character(0)
for (l in 1:copy_number) {
  for (i in seq(l, (subset_sum-1)*copy_number, copy_number)){
    ambiguity_num <- character(0)
    for (j in 1:copy_number) {
      seqinr::write.fasta(sequences = c(Consensus_seq_reduced[i], Consensus_seq_reduced[i+copy_number+j-l]), names(c(Consensus_seq_reduced[i], Consensus_seq_reduced[i+copy_number])),
                          file.out = paste0("cnx_",i,"_",i+copy_number+j-l,".fasta"))
      ambiguity_num <- append(ambiguity_num, ambiguity_count(paste0("cnx_",i,"_",i+copy_number+j-l,".fasta")))

      if (verbose) { cat(paste0("Matching sequences ",i, " & ", i+copy_number+j-l, "\n"))}
    }
    if (verbose) { cat("Number of ambiguous sites for each match:", ambiguity_num, "\n")}
    if (verbose) { cat(paste0("--- The pair ", which(as.numeric(ambiguity_num)==min(as.numeric(ambiguity_num)))[1], " has fewer ambiguous sites and should be assembled together\n"))}

    seq_to_con <- append(seq_to_con, paste0("A",i,"B_A",i+copy_number+which(as.numeric(ambiguity_num)==min(as.numeric(ambiguity_num)))[1]-l,"B"))
  }
}
seq_to_con1 <- seq_to_con[stringr::str_order(seq_to_con, decreasing = FALSE, na_last = TRUE, locale = "en",numeric = TRUE)] # order numerically

unlink("cnx_*") # Delete all intermediate files whose names begin with "cnx_"

if (verbose) { cat("*************************************************************************\n")}
# To find out which sequences belong to which gene copy
Copy_list <- data.frame()
for (i in 1:copy_number) {
  cp <- c(i)
  for (j in 1:(subset_sum-copy_number)) {
    m <- which(gsub(".*_", "", seq_to_con1[cp[j]]) == gsub("_.*", "", seq_to_con1)) # From A1B_A4B to A4B_A5B, then to A5B_A8B, then to A8B_A10B ...; linking 1, 4, 5, 8, 10 ...
    cp <- append(cp,m)
  }
  cp_all <- as.numeric(gsub("[AB]","",c(gsub("_.*", "", seq_to_con1[cp]), gsub(".*_", "", seq_to_con1[cp][subset_sum-copy_number+1]))))
  Copy_list <- append(Copy_list,as.data.frame(cp_all))
}

seq_list <- character(0)
for (i in 1:copy_number) {
  cat("List of sequences to be assembled for gene copy", i, ": ", Copy_list[[i]], "\n")
  for (j in (1:copy_number)[-i]){
    seq_list <- append(seq_list, intersect(as.character(Copy_list[[i]]),as.character(Copy_list[[j]])))
  }
}

if (length(unique(seq_list))>0) {
  cat("--- Sequences involved in the assembling of multiple gene copies: ", stringr::str_sort(unique(seq_list), numeric=TRUE), "\n")
  cat("Please check your input file carefully and run 'copy_assemble' again or do the assembling manually!\n")
}

all_copies_final <- character(0)
for (i in 1:copy_number) {
  copy_subseqs <- lapply(as.numeric(unlist(Copy_list[i])), function(x) Consensus_seq_reduced[x])
  copy_subseqs_name <- lapply(as.numeric(unlist(Copy_list[i])), function(x) names(Consensus_seq_reduced[x]))
  seqinr::write.fasta(sequences = copy_subseqs, copy_subseqs_name,file.out=paste0(filename_short, "_Copy_",i,"_subseqs.fasta"))
  seqinr::write.fasta(sequences = as.character(DECIPHER::ConsensusSequence(Biostrings::readDNAStringSet(paste0(filename_short, "_Copy_",i,"_subseqs.fasta"), format="fasta",nrec=-1L, skip=0L),
                                  threshold = 0.4,ambiguity = TRUE, minInformation=0.8, noConsensusChar = "N")), paste0("Copy_",i,"_final"),file.out=paste0("Copy_",i,"_final.fasta"))
  all_copies_final <- append(all_copies_final, seqinr::read.fasta (paste0("Copy_",i,"_final.fasta"), seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE))
}
seqinr::write.fasta(sequences = all_copies_final, names(all_copies_final),file.out="All_copies_final.fasta")

## Remove shared gaps in the final gene copies
DNAbin <- Biostrings::readDNAStringSet("All_copies_final.fasta")
seqinr::write.fasta(sequences = lapply(1:length(DNAbin), function(x) paste0(c(as.character(DNAbin[x][[1]]), as.character(rep("-",max(Biostrings::width(DNAbin))-Biostrings::width(DNAbin)[x]))), collapse='')),
                    names = names(DNAbin), file.out = "All_copies_same_length_final.fasta")
All_copies_final <- DECIPHER::RemoveGaps(Biostrings::readDNAStringSet("All_copies_same_length_final.fasta", format="fasta"),removeGaps = "common")
Biostrings::writeXStringSet(All_copies_final, paste0(filename_short, "_assembled_", copy_number, "_gene_copies.txt"))
unlink("*_final.fasta")

if (verbose) { cat("*************************************************************************\n")}
cat("Run finished!\n")
cat("Please check '*_subseqs.fasta' files to see if sequences have been linked correctly. Pay attention to nucleotide overhangs introduced by mistake.\n")

beepr::beep(sound = 1, expr = NULL) # make a sound when run finishes
options("error" = error_log_function)
sink() # turn off log
}

