#' Check Oligomers and Potential Primer Pairs
#'
#' This function takes output from get_candidate_oligos() and blast_candidate_oligos(). It checks each candidate for how many times it hit the target organisms, then sorts potential pairs to minimize the number of hits. The pairs in the output table can then be manually checked for ideal primer characteristics, etc.
#' @param candidates Filepath for the candidate oligomers file. No Default.
#' @param blastoutput Filepath for the output from blasting the candidates against the target database. No Default.
#' @param outputgraph Filepath for the output graph, which will show how many times each oligo hits target sequences. Defaults to "Oligos_Plot.pdf".
#' @param outputfile Filepath for the output table, which will check pairs of oligos and sort them by the minimum total number of times they hit target sequences. Defaults to "OligoPairs_Table.txt".
#' @keywords candidates
#' @keywords blastoutput
#' @keywords outputgraph
#' @keywords outputfile
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#' check_candidate_oligos(candidates="Candidate_oligos.fasta", blastoutput="Candidate_vs_UNITE.out")

check_candidate_oligos <- function(candidates, blastoutput, outputgraph="Oligos_Plot.pdf", outputfile="OligoPairs_Table.txt"){

	# Get the query file into a data frame
	kmers <- grep('kmer', readLines(candidates), value=T)
	seqs <- grep('kmer', readLines(candidates), value=T, invert=T)
	kmers <-as.data.frame(kmers)
	seqs <- as.data.frame(seqs)
	kmers <- cbind(kmers,seqs)

	# Grep the blast result and return the number of hits, put in dataframe
	grepreturn <- grep('hits found', readLines(blastoutput), value=T)
	df <- as.data.frame(grepreturn)
	df2 <- separate(df, grepreturn, c("hash", "numhits", "hits", "found"), " ")

	numcandidates <- length(row.names(df2))
	minhits <- min(as.numeric(df2$numhits))
	maxhits <- max(as.numeric(df2$numhits))

	# Now, check all primer pair sets with left primer between 0 to 1/3 of non-target length and right primer between 2/3 to end of non-target length
	#OnesixthL <- floor((1/6) * numcandidates)
	OnethirdL <- floor((1/3) * numcandidates)
	TwothirdL <- floor((2/3) * numcandidates)
	#FiveSixthsL <- floor((5/6) * numcandidates)

	primerpairs <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("Left", "Right", "LeftSeq","RightSeq","HitsLeft", "HitsRight", "HitsTotal"))

	for(j in 1: OnethirdL){
		for(i in TwothirdL: numcandidates){

			newline <- c(j, i, 	as.character(kmers$seqs[which(kmers$kmers==paste(">kmer_",j,sep=""))]),as.character(kmers$seqs[which(kmers$kmers==paste(">kmer_",i,sep=""))]), as.numeric(df2$numhits)[j], as.numeric(df2$numhits)[i], as.numeric(df2$numhits)[j] + as.numeric(df2$numhits)[i])
			primerpairs[nrow(primerpairs) + 1,] <- newline

		}
	}

	# Sort the primer pairs by the least total number of hits and make the output file
	sortpps <- primerpairs[order(as.numeric(primerpairs$HitsTotal)),]
	write.table(sortpps, file=outputfile, sep="\t", row.names=F)

	# Print a graph of the number of hits per candidate oligos
	pdf(file= outputgraph)
	plot(row.names(df2), df2$numhits, xlab = c("Candidate blocking oligo position on non-target sequence"), ylab=c("Number of hits to target database"))
	dev.off()

}
