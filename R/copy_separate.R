#' @title copy_separate
#'
#' @description Separates two or more gene copies from short-read Next-Generation Sequencing data into a small number of overlapping DNA sequences.
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
#' @return A fasta alignment of a small number of overlapping DNA sequences covering the entire length of the target gene. Gene copies can be assembled by reordering the alignment manually or use the function "copy_assemble".
#'
#' @importFrom seqinr read.fasta write.fasta
#'
#' @importFrom stringr str_count str_sort
#'
#' @importFrom kmer otu
#'
#' @importFrom Biostrings readDNAStringSet
#'
#' @importFrom ape as.character.DNAbin read.FASTA
#'
#' @importFrom DECIPHER ConsensusSequence
#'
#' @importFrom beepr beep
#'
#' @examples
#' \dontrun{
#' copy_separate("inst/extdata/toydata.fasta",2,300,225,10,1)
#' }
#'
#' @export copy_separate
#'

copy_separate<-function(filename,copy_number,read_length,overlap=225, rare_read=10, verbose=1)
{

sink(paste0(gsub("[:.:].*","", filename), "_separate_log.txt"), append=FALSE, split=TRUE) # begin to record log
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

begin_number <- seq(1,alignment_length-200, by=read_length-overlap); if (verbose) {cat("Beginning position of each subset","\n"); cat(begin_number,"\n")}
end_number <- begin_number+read_length-1; if (verbose) {cat("Ending position of each subset\n"); cat(end_number,"\n")}
number_of_subsets <- length(begin_number); if (verbose) {cat(paste0("Total number of subsets = ",number_of_subsets,"\n"))}

empty_subset <- list()
for (i in begin_number) {
  subset_original <- lapply(1:total_reads, function(x) {substr(NGSreads[x],i,i+read_length-1)}) # begin to subdivide the big alignment into subsets, each has the length of read_length
  subset_small <- subset_original[which(as.character(lapply(1:total_reads, function(x) {substr(subset_original[x],1,10)}))!="----------")]
  subset_smaller <- subset_small[which(as.character(lapply(1:length(subset_small), function(x) {substr(subset_small[x],200,209)}))!="----------")]
  subset_smallest <- subset_smaller[which(lapply(1:total_reads, function(x) {stringr::str_count(substr(subset_smaller[x],1,200),"A")>25})==TRUE)]
  if (verbose) {cat(paste0("Number of reads in subset ", which(begin_number==i), " = ",length(subset_smallest),"\n"))}
  if (length(subset_smallest)==0) {empty_subset <- append(empty_subset, which(begin_number==i))}
  seqinr::write.fasta(sequences = subset_smallest, names = 1:length(subset_smallest), file.out = paste0("Subset_",which(begin_number==i),"_downsized.fasta"))
}

if (length(empty_subset)>0) {
  if (verbose) {cat(paste0("Note: the following subset has no reads: ", empty_subset,"\n"),sep="")}
}

Subsets <- stringr::str_sort(list.files(pattern="_downsized.fasta"), numeric = TRUE)

# for each of the subset_1, 2, 3...
for (i in 1:length(Subsets)) {

  if (verbose) {cat("********************************************\n")}

  if (i %in% as.numeric(empty_subset)) next # skip empty subsets

  Sub_set <- ape::as.character.DNAbin(ape::read.FASTA(file=Subsets[i], type = "DNA"))

  if (verbose) {cat(paste0("Clustering analyses for the Subset ", i,"\n"))}

  # find the threshold range for OTU to find the major clusters (number=copy_number) for each subset
  for (m in seq(0.3,1, by = 0.1)) {
    Subset_OTU <- kmer::otu(Sub_set, k = 5, threshold = m, method = "central", nstart = 20)
    if (verbose) {cat(paste0("threshold = ",m),"\n")}
    if (verbose) {cat(unique(Subset_OTU),"\n")}
    if (length(unique(Subset_OTU))>=copy_number) {break}
  }

  # try different threshold values in the range found above
  if (length(unique(Subset_OTU))>copy_number) {
    for (j in seq(m-0.09,m, by = 0.01)) {
     Subset_OTU <- kmer::otu(Sub_set, k = 5, threshold = j, method = "central", nstart = 20)
     if (verbose) {cat(paste0("threshold = ",j),"\n")}
     if (verbose) {cat( unique(Subset_OTU),"\n")}
     if (length(unique(Subset_OTU))>=copy_number) {break}
     }
  }

  reads_each_cluster <- sapply(unique(Subset_OTU), function(x) length(which(Subset_OTU==x)))
  cluster_DF <- cbind(as.data.frame(unique(Subset_OTU)), as.data.frame(reads_each_cluster))[order(-reads_each_cluster),]
  names(cluster_DF) <- c("cluster_num","cluster_total")

  if (verbose) {cat("Number of reads in each cluster\n")}
  if (verbose) {cat(reads_each_cluster,"\n")}

 for  (l in (1:length(cluster_DF$cluster_num))){
   Picked_cluster <- Sub_set[which(Subset_OTU==cluster_DF$cluster_num[l])]

   if (length(Picked_cluster)<rare_read) break

   if (verbose) {cat(paste0("Number of reads in picked cluster ",l, " = ", length(Picked_cluster),"\n"))}

   seqinr::write.fasta(sequences = Picked_cluster, names = labels(Picked_cluster), file.out = paste0("Subset_",i,"_cluster_",l,".fasta"))
   seqinr::write.fasta(sequences = paste0(c(as.character(rep("-",(read_length-overlap)*(i-1))),as.character(DECIPHER::ConsensusSequence(Biostrings::readDNAStringSet(paste0("Subset_",i,"_cluster_",l,".fasta"),
                       format="fasta",nrec=-1L, skip=0L),threshold = 0.4, ambiguity = TRUE, minInformation=0.6, noConsensusChar = "N")[1])), collapse=''),
                       names = paste0("Subset_",i,"_cluster_",l,"_consensus"), file.out = paste0("Subset_",i,"_cluster_",l,"_consensus.fasta"))
 }
}

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
cat("Run finished!\n")
beepr::beep(sound = 1, expr = NULL) # make a sound when run finishes
options("error" = error_log_function)
sink() # turn off log
}

