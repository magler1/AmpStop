#' Check Oligomers and Potential Primer Pairs
#'
#' This function takes output from get_candidate_oligos() and blast_candidate_oligos(). It checks each candidate for how many times they hits the target organisms (in the example fungi in the Unite ITS database), how many times the hits align along >90% of the length of the primer, and how many hits align at the 3' end of the primer (alignments are counted if they end within 3 bases of the end of hte primer. This type of alignment is the most serious since the polymerase binds there). It outputs a graph (in pdf format) with each of the three counts as well as the same data with the primer sequences in a tab-delimited table. Good candidates can be checked subsequently for the ideal primer characteristics, etc.
#' @param candidates Filepath for the candidate oligomers file. No Default.
#' @param blastoutput Filepath for the output from blasting the candidates against the target database. No Default.
#' @param outputgraph Filepath for the output graph, which will show how many times each oligo hits target sequences. Defaults to "Oligos_Plot.pdf".
#' @param outputfile Filepath for the output table, which will check pairs of oligos and sort them by the minimum total number of times they hit target sequences. Defaults to "OligoPairs_Table.txt".
#' @keywords candidates
#' @keywords blastoutput
#' @keywords outputgraph
#' @keywords outputfile
#' @import tidyr
#' @import ggpubr
#' @import ggplot2
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

  ## Grep the blast result and return the total number of hits, put in dataframe
  grepreturn <- grep('hits found', readLines(blastoutput), value=T)
  df <- as.data.frame(grepreturn)
  df2 <- separate(df, grepreturn, c("hash", "numhits", "hits", "found"), " ")

  # Divide the blast results into the different kmers in individual temporary files
  dir.create("./TempFolder")
  blastoutputlines <- readLines(blastoutput)
  nlines <- length(blastoutputlines)

  queryline <- 0
  for(i in 1:nlines){

    if(grepl("Query", blastoutputlines[i])){
      # Make a new file in the temp directory for each kmer
      currentline <- blastoutputlines[i]
      kmernumber <- strsplit(currentline, split=" ")[[1]][3]
      queryline <- 0
      outfilename <- paste("./TempFolder/",kmernumber,".txt", sep="")
      file.create(outfilename)
      next
    }

    if(queryline == 1){
      # Get the headers into the output file
      currentline <- blastoutputlines[i]
      csfields <- strsplit(currentline, ":")[[1]][2]
      fields <- paste(strsplit(csfields, ", ")[[1]], collapse="\t")
      write(fields,file=outfilename)
      queryline = queryline + 1
      next
    }

    if(queryline > 2){
      # After skipping the 3rd line (hits found) add the rest of the lines for the query, except for the commented line bfore the next "Query" line
      currentline <- blastoutputlines[i]
      if(grepl("#", currentline)){
        queryline = queryline + 1
        next
      }
      write(currentline,file=outfilename,append=TRUE)
      queryline = queryline + 1
      next
    }
    queryline = queryline + 1
  }

  numcandidates <- length(row.names(kmers))

  # Now, get the stats for each potential oligo
  allprimers <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("kmer", "Seq","HitsTotal", "Hits >90% aligned length", "Hits alignd at 3' End"))

  for(i in 1:numcandidates){
    sequence <- as.character(kmers$seqs[which(kmers$kmers==paste(">kmer_",i,sep=""))])
    HitsTotal <- as.numeric(df2$numhits)[i]
    maxalignlength <- 0.9 * nchar(sequence)

    filename <- paste("./TempFolder/kmer_",i,".txt", sep="")
    kmertable <- read.table(file=filename, header=T, sep="\t")
    HitsThreePrime <- length(which(kmertable$q..end > (nchar(sequence) - 3)))

    HitsIdenticalPositions <- length(which(kmertable$identical > maxalignlength))

    newline <- c(i, sequence, HitsTotal, HitsIdenticalPositions, HitsThreePrime)
    allprimers[nrow(allprimers) + 1,] <- newline
  }
  allprimers$kmer <- as.numeric(allprimers$kmer)
  write.table(allprimers, file=outputfile, sep="\t", row.names=F)

  # Print a graph with the 3 different hit counts per candidate oligos
  OnethirdL <- floor((1/3) * numcandidates)
  TwothirdL <- floor((2/3) * numcandidates)

  a <- ggplot(allprimers, aes(x=as.numeric(kmer), y=as.numeric(HitsTotal))) + geom_point() + xlab("kmer") + ylab("Total Hits") + theme_classic() + geom_vline(xintercept = OnethirdL, linetype="dashed", color = "dodgerblue", size=1.5) + geom_vline(xintercept = TwothirdL, linetype="dashed", color = "coral2", size=1.5)
  b <- ggplot(allprimers, aes(x=as.numeric(kmer), y=as.numeric(`Hits >90% aligned length`))) + geom_point() + xlab("kmer") + ylab("Hits > 90% primer length")  + theme_classic() + geom_vline(xintercept = OnethirdL, linetype="dashed", color = "dodgerblue", size=1.5) + geom_vline(xintercept = TwothirdL, linetype="dashed", color = "coral2", size=1.5)
  c <- ggplot(allprimers, aes(x=as.numeric(kmer), y=as.numeric(`Hits alignd at 3' End`))) + geom_point()+ xlab("kmer") + ylab("Hits aligned at 3' End")  + theme_classic() + geom_vline(xintercept = OnethirdL, linetype="dashed", color = "dodgerblue", size=1.5) + geom_vline(xintercept = TwothirdL, linetype="dashed", color = "coral2", size=1.5)
  pdf(file=outputgraph, height=10)
  print(ggarrange(a+ rremove("xlab")+rremove("x.text"),b+ rremove("xlab")+rremove("x.text"),c, labels = c("A","B","C"), ncol=1, nrow=3))
  dev.off()

  # Remove the temporary directory:
  unlink("./TempFolder", recursive=TRUE)

}
