#' Blast Candidate Oligomers
#'
#' This function will blast candidate oligomers (from a nontarget organism) against a database of target organisms and returns a parsable format.
#' @param blastn_path Path to blastn. Defaults to "blastn". If your blast copy is somewhere else, please point to it here.
#' @param blastdb Filepath for the blast databse. No Default
#' @param candidates Filepath for the candidate oligomers file. No Default.
#' @param outfile Filepath for the output file name, which will contain the blast results. Defaults to "Candidates_vs_UNITE.out"
#' @keywords nontarget
#' @keywords length
#' @keywords outfile
#' @export
#' @examples
#' blast_candidate_oligos(blastdb="./Unite_s_011217/UNITE_s_011217")

OS <- Sys.info()['sysname']

blast_candidate_oligos <- function(blastn_path = "blastn", blastdb, candidates, outfile="Candidate_vs_UNITE.out"){

	# create blast command function. The parameters are currently not adjustable, but will be printed
  cmd <- paste(blastn_path, " -db ", blastdb ," -query ", candidates ," -out ", outfile," -outfmt \"7 nident std\" -perc_identity 50 -word_size 12 -dust no -evalue 1000 -max_target_seqs 1000 ", sep = "")
  paste("Running Blast command: ",cmd, sep = "")

  # Run the blast command. This is OS-specific
  if(OS=="Darwin" || OS=="Linux"){
	  #Run the command
	  sapply(cmd, system)
  }

  if(OS=="Windows"){
    #Run the command
    shell(cmd)
  }

}
