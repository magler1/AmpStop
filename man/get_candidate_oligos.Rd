% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_Oligo_Candidates_Func.r
\name{get_candidate_oligos}
\alias{get_candidate_oligos}
\title{Get Candidate Oligomers}
\usage{
get_candidate_oligos(
  nontarget,
  length = 30,
  outfile = "Candidate_oligos.fasta"
)
}
\arguments{
\item{nontarget}{Filepath for the non-target DNA sequence. No Default}

\item{length}{Desired length of candidate oligomers to return. Defaults to 30.}

\item{outfile}{Filepath for the output file name, which will contain the candidate oligomers. Defaults to "Candidate_oligos.fasta"}
}
\description{
This function allows you to provide a fastafile with a non-target sequence and a desired oligomer length and gives back all possible oligomers of that length in the sequence.
}
\examples{
get_candidate_oligos(nontarget="Athaliana_ITS1_new_wholelength.fasta")
}
\keyword{length}
\keyword{nontarget}
\keyword{outfile}
