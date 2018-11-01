#' Get Candidate Oligomers
#'
#' This function allows you to provide a fastafile with a non-target sequence and a desired oligomer length and gives back all possible oligomers of that length in the sequence.
#' @param nontarget Filepath for the non-target DNA sequence. No Default
#' @param length Desired length of candidate oligomers to return. Defaults to 30.
#' @param outfile Filepath for the output file name, which will contain the candidate oligomers. Defaults to "Candidate_oligos.fasta"
#' @keywords nontarget
#' @keywords length
#' @keywords outfile
#' @export
#' @examples
#' get_candidate_oligos(nontarget="Athaliana_ITS1_new_wholelength.fasta")

get_candidate_oligos <- function(nontarget, length=30, outfile="Candidate_oligos.fasta"){
	
	#First, get the sequence into an R object - here I will read in a fasta file containing one sequence
	fastafile <- readLines(nontarget)

	if(length(fastafile) > 2){
		#Throw error and quit, this only works with one nontarget sequence
	}
	
	# Now get the candidate kmers of the specified length from the nontarget sequence
	fastaoutput <- c()
	for(i in 1:nchar(fastafile[2])-length){
		
		header <- paste(">kmer_",i, sep="")
		kmer <- substring(fastafile[2], i, i+length)
		q <- i + i-1
		fastaoutput[q] <- header
		fastaoutput[q+1] <- kmer
			
	}
	
	# Output the kmers to the output file
	fileConn<-file(outfile)
	writeLines(fastaoutput, fileConn)
	close(fileConn)
	
}