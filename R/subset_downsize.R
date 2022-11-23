#' @title subset_downsize
#'
#' @description Subdivides the imported read alignment into subsets and then downsizes each subset by deleting those sequences that have too many gaps or missing data.
#'
#' @param filename A fasta file contains thousands of short reads that have been mapped to a reference. The reference and reads that are not directly mapped to the reference need to be removed after mapping.
#'
#' @param read_length An integer (e.g. 250, or 300) giving the read length of your Next-generation Sequencing data. This method is designed for read length >=250bp.
#'
#' @param overlap An integer describing number of base pairs of overlap between adjacent subsets. More overlap means more subsets.
#'
#' @param verbose Turn on (verbose=1; default) or turn off (verbose=0) the output.
#'
#' @return A number of overlapping subsets (before and after downsizing) of the input alignment.
#'
#' @importFrom seqinr read.fasta write.fasta
#'
#' @importFrom stringr str_count str_sort
#'
#' @importFrom beepr beep
#'
#' @examples
#' \dontrun{
#' subset_downsize("inst/extdata/toydata.fasta",300,225,1)
#' }
#'
#' @export subset_downsize
#'

subset_downsize<-function(filename,read_length,overlap, verbose=1)
{
  sink("log.txt", append=FALSE, split=TRUE) # begin to record log
  error_log_function <- function() {
    cat(geterrmessage(), file="Error_log.txt", append=T)
  }

  if (read_length<250) warning ("This method is designed for read length 250bp or longer. Short reads can easily result in chimeric sequences.")

  NGSreads <- seqinr::read.fasta(file = filename,seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)

  total_reads <- length(NGSreads); if (verbose) { cat(paste0("Total number of reads imported = ",total_reads,"\n"))}
  alignment_length <- nchar(NGSreads[[1]]); if (verbose) { cat(paste0("Length of the alignment = ",alignment_length,"\n"))}

  if (verbose) { cat(paste0("Read length = ",read_length,"\n"))}
  if (verbose) { cat(paste0("Overlap between adjacent subsets = ",overlap,"\n"))}
  if (overlap>=read_length) stop("Overlap between adjacent subsets must be smaller than the read length!")

  begin_number <- seq(1,alignment_length-200, by=read_length-overlap); if (verbose) { cat("Beginning position of each subset","\n"); cat(begin_number,"\n")}
  end_number <- begin_number+read_length-1; if (verbose) { cat("Ending position of each subset\n"); cat(end_number,"\n")}
  number_of_subsets <- length(begin_number); if (verbose) { cat(paste0("Total number of subsets = ",number_of_subsets,"\n"))}

  for (i in begin_number) {
    subset_original <- lapply(1:total_reads, function(x) {substr(NGSreads[x],i,i+read_length-1)}) # begin to subdivide the big alignment into subsets, each has the length of read_length
    seqinr::write.fasta(sequences = subset_original, names = 1:length(subset_original), file.out = paste0("Subset_",which(begin_number==i),"_original_size.fasta"))
    subset_small <- subset_original[which(as.character(lapply(1:total_reads, function(x) {substr(subset_original[x],1,10)}))!="----------")]
    subset_smaller <- subset_small[which(as.character(lapply(1:length(subset_small), function(x) {substr(subset_small[x],200,209)}))!="----------")]
    subset_smallest <- subset_smaller[which(lapply(1:total_reads, function(x) {stringr::str_count(substr(subset_smaller[x],1,200),"A")>25})==TRUE)]
    if (verbose) { cat(paste0("Number of reads in downsized subset ", which(begin_number==i), " = ",length(subset_smallest),"\n"))}
    seqinr::write.fasta(sequences = subset_smallest, names = 1:length(subset_smallest), file.out = paste0("Subset_",which(begin_number==i),"_downsized.fasta"))
  }

  cat("Subsetting and downsizing finished!\n")
  beepr::beep(sound = 1, expr = NULL) # make a sound when run finishes
}


