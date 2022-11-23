#' @title sep_assem
#'
#' @description Separates two or more gene copies from short-read Next-Generation Sequencing data into a small number of overlapping DNA sequences and assemble them into their respective gene copies.
#'
#' @param filename A fasta file contains thousands of short reads that have been mapped to a reference. The reference and reads that are not directly mapped to the reference need to be removed after mapping.
#'
#' @param copy_number An integer (e.g. 2,3, or 4) giving the anticipated number of gene copies in the input file.
#'
#' @param read_length An integer (e.g. 250, or 300) giving the read length of your Next-generation Sequencing data. This method is designed for read length >=250bp.
#'
#' @param overlap An integer describing number of base pairs of overlap between adjacent subsets. More overlap means more subsets. Default 225.
#'
#' @param rare_read A positive integer. During clustering analyses, clusters with less than this number of reads will be ignored. Default 10.
#'
#' @param verbose Turn on (verbose=1; default) or turn off (verbose=0) the output.
#'
#' @return A fasta alignment of the anticipated number of full-length gene copies.
#'
#' @importFrom seqinr read.fasta write.fasta
#'
#' @importFrom stringr str_count str_sort
#'
#' @importFrom kmer otu
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet width
#'
#' @importFrom ape as.character.DNAbin read.FASTA
#'
#' @importFrom DECIPHER ConsensusSequence RemoveGaps
#'
#' @importFrom foreach foreach
#'
#' @importFrom parallel makeCluster
#'
#' @importFrom doParallel registerDoParallel
#' 
#' @importFrom beepr beep
#'
#' @examples
#' \dontrun{
#' sep_assem("inst/extdata/toydata.fasta",2,300,225,10,1,1)
#' }
#'
#' @export sep_assem
#'

sep_assem<-function(filename,copy_number,read_length,overlap=225, rare_read=10, core_number=1, verbose=1)
{
for (filename in list.files(pattern=".fasta")) {

# set parallel computing
my.cluster <- parallel::makeCluster(core_number, type = "PSOCK")

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

sink(paste0(gsub("[:.:].*","", filename), "_log.txt"), append=FALSE, split=TRUE) # begin to record log
error_log_function <- function() {
  cat(geterrmessage(), file="Error_log.txt", append=T)
}

if (copy_number<=1) stop ("The anticipated copy number must be a number larger than one!")

if (read_length<250) warning ("This method is designed for read length 250bp or longer. Short reads can easily result in chimeric sequences.")

NGSreads <- seqinr::read.fasta(file = filename,seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)

total_reads <- length(NGSreads); if (verbose) {cat(paste0("Total number of reads imported = ",total_reads,"\n"))}
alignment_length <- nchar(NGSreads[[1]]); if (verbose) {cat(paste0("Length of the alignment = ",alignment_length,"\n"))}

if (verbose) {cat(paste0("Read length = ",read_length,"\n"))}
if (verbose) {cat(paste0("Overlap between adjacent subsets = ",overlap,"\n"))}
if (overlap>=read_length) stop("Overlap between adjacent subsets must be smaller than the read length!")

begin_number <- seq(1,alignment_length-200, by=read_length-overlap)
end_number <- begin_number+read_length-1
number_of_subsets <- length(begin_number); if (verbose) {cat(paste0("Total number of subsets = ",number_of_subsets,"\n"))}

# define a function for subsetting and downsizing
To_subset <- function(i) {
  subset_original <- lapply(1:total_reads, function(x) {substr(NGSreads[x],i,i+read_length-1)}) # begin to subdivide the big alignment into subsets, each has the length of read_length
  subset_small <- subset_original[which(as.character(lapply(1:total_reads, function(x) {substr(subset_original[x],1,10)}))!="----------")]
  subset_smaller <- subset_small[which(as.character(lapply(1:length(subset_small), function(x) {substr(subset_small[x],200,209)}))!="----------")]
  subset_smallest <- subset_smaller[which(lapply(1:total_reads, function(x) {stringr::str_count(substr(subset_smaller[x],1,200),"A")>25})==TRUE)]
  seqinr::write.fasta(sequences = subset_smallest, names = 1:length(subset_smallest), file.out = paste0("Subset_",which(begin_number==i),"_downsized.fasta"))
}

invisible(foreach::foreach(x = begin_number) %dopar% {
  To_subset(x)
})

Subsets <- stringr::str_sort(list.files(pattern="_downsized.fasta"), numeric = TRUE)

# find empty subsets
empty_subset <- list()
for (i in 1:number_of_subsets) {
  if (file.info(Subsets[i])$size==0) {empty_subset <- append(empty_subset, i)}
}

if (length(empty_subset)>0) {
  if (verbose) {cat(paste0("Note: the following subset has no reads: ", empty_subset,"\n"),sep="")}
}

if (verbose) {cat("Performing clustering analyses for each subset...\n")}

# define a function for clustering analyses on a single subset
OTU_subset <- function (i) {
  Sub_set <- ape::as.character.DNAbin(ape::read.FASTA(file=Subsets[i], type = "DNA"))

  # find the threshold range for OTU to find the major clusters (number=copy_number) for each subset
  for (m in seq(0.3,1, by = 0.1)) {
    Subset_OTU <- kmer::otu(Sub_set, k = 5, threshold = m, method = "central", nstart = 20)
    if (length(unique(Subset_OTU))>=copy_number) {break}
  }
  
  # try different threshold values in the range found above
  if (length(unique(Subset_OTU))>copy_number) {
    for (j in seq(m-0.09,m, by = 0.01)) {
      Subset_OTU <- kmer::otu(Sub_set, k = 5, threshold = j, method = "central", nstart = 20)
      if (length(unique(Subset_OTU))>=copy_number) {break}
    }
  }

  reads_each_cluster <- sapply(unique(Subset_OTU), function(x) length(which(Subset_OTU==x)))
  cluster_DF <- cbind(as.data.frame(unique(Subset_OTU)), as.data.frame(reads_each_cluster))[order(-reads_each_cluster),]
  names(cluster_DF) <- c("cluster_num","cluster_total")

 for  (l in (1:length(cluster_DF$cluster_num))){
   Picked_cluster <- Sub_set[which(Subset_OTU==cluster_DF$cluster_num[l])]

   if (length(Picked_cluster)<=rare_read) break

   seqinr::write.fasta(sequences = Picked_cluster, names = labels(Picked_cluster), file.out = paste0("Subset_",i,"_cluster_",l,".fasta"))
   seqinr::write.fasta(sequences = paste0(c(as.character(rep("-",(read_length-overlap)*(i-1))),as.character(DECIPHER::ConsensusSequence(Biostrings::readDNAStringSet(paste0("Subset_",i,"_cluster_",l,".fasta"),
                       format="fasta",nrec=-1L, skip=0L),threshold = 0.4, ambiguity = TRUE, minInformation=0.6, noConsensusChar = "N")[1])), collapse=''),
                       names = paste0("Subset_",i,"_cluster_",l,"_consensus"), file.out = paste0("Subset_",i,"_cluster_",l,"_consensus.fasta"))
 }
}

# for each of the subset_1, 2, 3...
invisible(foreach::foreach(x = setdiff(1:length(Subsets), empty_subset)) %dopar% {
 OTU_subset(x)
})

if (verbose) {cat("Clustering analyses finished\n")}

# put together all the consensus sequences from all subsets into one file
Consensus_list <- stringr::str_sort(list.files(pattern="_consensus.fasta"), numeric = TRUE)
All_consensus <- lapply(1:length(Consensus_list), function (x) seqinr::read.fasta(file = Consensus_list[x], seqtype = "DNA",
                                                                                  as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE, whole.header = TRUE))
filename_short <- gsub("[:.:].*","", filename) # remove file extensions, e.g. ".fasta", ".txt"

# move subset files into the intermediate files folder
dir.create(paste0(filename_short,"_intermediate_files"))
invisible(file.copy(list.files(pattern="Subset_"), paste0(filename_short,"_intermediate_files")))  # use "invisible" so that output do not show here
unlink("Subset_*")

seqinr::write.fasta(sequences=All_consensus, names=Consensus_list, file.out=paste0(filename_short,"_copies",copy_number,"_overlap",overlap,".txt"))




# Runs "copy_separate" above this line, runs "copy_assemble" below this line.
#####################################################################################################################################################################################################################

# Define a function to count the total number of ambiguous sites in an alignment
ambiguity_count <- function(x){
  sum(as.integer(stringr::str_count(as.character(DECIPHER::ConsensusSequence(Biostrings::readDNAStringSet(x, format="fasta",nrec=-1L, skip=0L),
                                                                             threshold = 0.4,ambiguity = TRUE, minInformation=0.6, noConsensusChar = "N")), c("M","K","R","Y","W","S","H","V","D","B"))))
}

Consensus_seq_org <- seqinr::read.fasta(file = paste0(filename_short,"_copies",copy_number,"_overlap",overlap,".txt"),seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)
subset_num_org <- as.numeric(gsub("\\D", "", gsub("_cluster_.", "", names(Consensus_seq_org))))
subset_sum <- length(unique(subset_num_org))
cluster_list_org <- lapply(sort(unique(subset_num_org)), function(x) which(subset_num_org==x))

# List missing subsets when exist
if (max(subset_num_org)!=subset_sum) {
  if (verbose) {cat(paste0("Note: Subset ", setdiff(1:max(subset_num_org), subset_num_org), " is not present in the file: ",paste0(filename_short,"_copies",copy_number,"_overlap",overlap,".txt"),"\n"),sep="")}
}

# Remove ambiguous region from the end of some sequences
Consensus_seq_clean <- list()
for (i in (1:length(Consensus_seq_org))){
  ambiguity_each <- sum(stringr::str_count(stringr::str_sub(Consensus_seq_org[i],-70,-1), c("M", "K", "R", "Y","N","W", "S", "H", "V", "D", "B")))
  if (ambiguity_each<8){
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
      third_cp <- append(third_cp, m)
    }
    unlink(list.files(pattern="x12nc_"))
    unlink(list.files(pattern="x23nc_"))

    picked_3rd <- character(0)
    for (n in 1:length(cluster_over)) {picked_3rd <- append(picked_3rd, cluster_list[cluster_over][[n]][as.numeric(third_cp[n])])}
    list_reduced <- sort(as.numeric(c(unlist(cluster_list[cluster_equal]),picked_2nd,picked_3rd)))
  }

  Consensus_seq_reduced <- Consensus_seq[list_reduced]
  seqinr::write.fasta(sequences = Consensus_seq_reduced, names(Consensus_seq_reduced),file.out = paste0(filename_short,"_copies",copy_number,"_overlap",overlap,"_reduced.txt"))
}


# When copy_number >3, it is strongly encouraged to do the assembling manually. Below we pick the copy_number of largest clusters to represent each gene copy.
if (length(cluster_equal)<subset_sum && copy_number>3){
  picked_copies <- unlist(lapply(1:length(cluster_over), function(y) {lapply(1:copy_number, function(x) cluster_list[cluster_over][[y]][x])}))
  list_reduced <- sort(as.numeric(c(unlist(cluster_list[cluster_equal]),picked_copies)))
  Consensus_seq_reduced <- Consensus_seq[list_reduced]
  seqinr::write.fasta(sequences = Consensus_seq_reduced, names(Consensus_seq_reduced),file.out = paste0(filename_short,"_copies",copy_number,"_overlap",overlap,"_reduced.txt"))
}


# Code below will assemble copy_number of consensus sequences from each subset into copy_number of full length gene copies

seq_to_con <- character(0)
for (l in 1:copy_number) {
  for (i in seq(l, (subset_sum-1)*copy_number, copy_number)){
    ambiguity_num <- character(0)
    for (j in 1:copy_number) {
      seqinr::write.fasta(sequences = c(Consensus_seq_reduced[i], Consensus_seq_reduced[i+copy_number+j-l]), names(c(Consensus_seq_reduced[i], Consensus_seq_reduced[i+copy_number])),
                          file.out = paste0("cnx_",i,"_",i+copy_number+j-l,".fasta"))
      ambiguity_num <- append(ambiguity_num, ambiguity_count(paste0("cnx_",i,"_",i+copy_number+j-l,".fasta")))
    }
    seq_to_con <- append(seq_to_con, paste0("A",i,"B_A",i+copy_number+which(as.numeric(ambiguity_num)==min(as.numeric(ambiguity_num)))[1]-l,"B"))
  }
}
seq_to_con1 <- seq_to_con[stringr::str_order(seq_to_con, decreasing = FALSE, na_last = TRUE, locale = "en",numeric = TRUE)] # order numerically

unlink("cnx_*") # Delete all intermediate files whose names begin with "cnx_"

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

cat("Lists of sequences in", paste0(filename_short,"_copies",copy_number,"_overlap",overlap,"_reduced.txt"),"(if not exist, then its un-reduced form) to be assembled for different gene copies", "\n")
seq_list <- character(0)
for (i in 1:copy_number) {
  cat(paste0("Gene_copy_", i, ": "), Copy_list[[i]], "\n")
  for (j in (1:copy_number)[-i]){
    seq_list <- append(seq_list, intersect(as.character(Copy_list[[i]]),as.character(Copy_list[[j]])))
  }
}

if (length(unique(seq_list))>0) {
  shared_seq <- stringr::str_sort(unique(seq_list), numeric=TRUE)
  First_prob_subset <- gsub("_cluster.*", "", names(Consensus_seq_reduced)[as.numeric(shared_seq[1])])
  cat("--- Sequences involved in the assembling of multiple gene copies: ", "\n", shared_seq, "\n")
  cat(paste0("Please check and edit sequences of the ", First_prob_subset, " in the file ", paste0(filename_short,"_copies",copy_number,"_overlap",overlap,".txt"), ", then run 'copy_assemble' on the revised file or do the assembling manually.","\n"))
}

all_copies_final <- character(0)
for (i in 1:copy_number) {
  copy_subseqs <- lapply(as.numeric(unlist(Copy_list[i])), function(x) Consensus_seq_reduced[x])
  copy_subseqs_name <- lapply(as.numeric(unlist(Copy_list[i])), function(x) names(Consensus_seq_reduced[x]))
  seqinr::write.fasta(sequences = copy_subseqs, copy_subseqs_name,file.out=paste0(filename_short, "_Copy_",i,"_subseqs.fasta"))
  seqinr::write.fasta(sequences = as.character(DECIPHER::ConsensusSequence(Biostrings::readDNAStringSet(paste0(filename_short, "_Copy_",i,"_subseqs.fasta"), format="fasta",nrec=-1L, skip=0L),
                                                                           threshold = 0.4,ambiguity = TRUE, minInformation=0.6, noConsensusChar = "N")), paste0("Copy_",i,"_final"),file.out=paste0("Copy_",i,"_final.fasta"))
  all_copies_final <- append(all_copies_final, seqinr::read.fasta (paste0("Copy_",i,"_final.fasta"), seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE))
}
seqinr::write.fasta(sequences = all_copies_final, names(all_copies_final),file.out="All_copies_final.fasta")

## Remove shared gaps in the final gene copies
DNAbin <- Biostrings::readDNAStringSet("All_copies_final.fasta")
seqinr::write.fasta(sequences = lapply(1:length(DNAbin), function(x) paste0(c(as.character(DNAbin[x][[1]]), as.character(rep("-",max(nchar(DNAbin))-nchar(DNAbin)[x]))), collapse='')),
                    names = names(DNAbin), file.out = "All_copies_same_length_final.fasta")
All_copies_final <- DECIPHER::RemoveGaps(Biostrings::readDNAStringSet("All_copies_same_length_final.fasta", format="fasta"),removeGaps = "common")
Biostrings::writeXStringSet(All_copies_final, paste0(filename_short, "_assembled_", copy_number, "_gene_copies.txt"))
unlink("*_final.fasta")

if (verbose) { cat("*************************************************************************\n")}
cat("Run finished!\n")
cat("Please check '*_subseqs.fasta' files to see if sequences have been linked correctly. Pay attention to nucleotide overhangs introduced by mistakes.\n")

# Move the input data and all intermediate and resulting files into a single newly created folder
dir.create(paste0("--- Results_", filename_short))
invisible(file.copy(list.files(pattern=paste0(filename_short,"_")), paste0("--- Results_", filename_short), recursive = TRUE))
invisible(file.copy(filename, paste0("--- Results_", filename_short)))
unlink(list.files(pattern=paste0(filename_short,"_")), recursive=TRUE)
unlink(filename)

beepr::beep(sound = 1, expr = NULL) # make a sound when run finishes
options("error" = error_log_function)
sink() # turn off log

}
