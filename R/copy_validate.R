#' @title copy_validate
#'
#' @description A tool to help identify incorrectly assembled chimeric sequences.
#'
#' @param filename A DNA alignment in fasta format that contains sequences of two or more gene copies (e.g. results from "copy_assemble").
#'
#' @param copy_number An integer (e.g. 2,3, or 4) giving the number of gene copies in the input file.
#'
#' @param read_length An integer (e.g. 250, or 300) giving the read length of your Next-generation Sequencing data.
#'
#' @param verbose Turn on (verbose=1; default) or turn off (verbose=0) the output.
#'
#' @return A histogram in pdf format showing the relationships between the physical distance between neighboring variable sites and read length.
#'
#' @importFrom seqinr read.fasta
#'
#' @importFrom grDevices dev.off pdf
#'
#' @importFrom graphics abline hist text
#'
#' @importFrom beepr beep
#' 
#' @examples 
#' \dontrun{
#' copy_validate("inst/extdata/Final_two_copies.fasta",2,300,1)
#' }
#'
#' @export copy_validate


copy_validate<-function(filename,copy_number,read_length, verbose=1)
{
  sink("copy_validate_log.txt", append=FALSE, split=TRUE) # begin to record log
  error_log_function <- function() {
    cat(geterrmessage(), file="Error_log.txt", append=T) # begin to record error log
  }

  all_copies <- seqinr::read.fasta(filename, seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)

  all_dist <-  character(0)
  for (i in 1:(copy_number)) {
    for (j in i:(copy_number)) {
      Paralog_diff <- which(unlist(strsplit(all_copies[[i]], ""))!=unlist(strsplit(all_copies[[j]], "")))

      dist_var_sites <- sapply(2:length(Paralog_diff), function(x) Paralog_diff[x]-Paralog_diff[x-1])
      all_dist <- append(all_dist,  dist_var_sites)
      if (verbose) { cat(paste0("Copy ", i, " and Copy ", j, " is different at site: ", Paralog_diff,"\n"))}
    }
  }

  dist_final <- sort(as.numeric(all_dist))

  if (max(dist_final)>read_length)
   cat ("Warning! Distance between some neighboring variable sites is larger than the read length (see the plot in pdf), chimeric sequences may have formed during assembling!", "\n")

  text_to_show <- "read_length"

  grDevices::pdf(file="Distance between neighboring variable sites VS. Read length.pdf", width=8, height=8)

  graphics::hist(dist_final,main="Distance between neighboring variable sites VS. Read length", xlab="Base pairs",
       ylim=c(0,length(dist_final)/2),xlim=c(0,read_length+100),breaks="Sturges",col="blue")

  graphics::abline(v=read_length, col="red", lwd=3, lty=2)
  graphics::text(read_length, length(dist_final)/2, text_to_show, pos = 4, offset = 0.5)

  cat('Run finished! Please check results in the file: "Distance between neighboring variable sites VS. read length.pdf"',"\n")

  beepr::beep(sound = 1, expr = NULL) # make a sound when run finishes
  options("error" = error_log_function)
  sink() # turn off log
  invisible(grDevices::dev.off())
}




