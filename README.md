# AmpStop

# What and why?: 
“AmpStop”, an “R” package to automate design of blocking oligos. It makes it easy to quickly and easily design oligos to block nontarget amplification in amplicon sequencing projects. AmpStop can be used by anyone with R and BLAST+ installed on their computer. It requires as input only the nontarget, amplified region and a target sequence database that is BLAST-formatted (an example can be found in the folder "Testing"). By running three R functions, users can within minutes generate a list of all possible blocking oligos, a figure showing how many times each oligo “hits” target templates (fewer hits indicates higher specificity for non-targets) and a ranked list of the most promising blocking oligo pairs.

# Installation:

```
devtools::install_github("magler1/AmpStop")
library(AmpStop)
```

# Running:

First, download the folder "Testing", then:

```
setwd("Path/to/Testing/")

?get_candidate_oligos
get_candidate_oligos(nontarget="Athaliana_ITS1_new_wholelength.fasta")

?blast_candidate_oligos
blast_candidate_oligos(blastn_path="/Path/to/blastn", blastdb="./Unite_s_011217/UNITE_s_011217", candidates="Candidate_oligos.fasta")

?check_candidate_oligos
check_candidate_oligos(candidates="Candidate_oligos.fasta", blastoutput="Candidate_vs_UNITE.out")
```
