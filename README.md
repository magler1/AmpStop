# AmpStop

# What and why?: 
“AmpStop” is an “R” package to automate design of blocking oligos. It makes it easy to quickly and easily design oligos to block nontarget amplification in amplicon sequencing projects. AmpStop can be used by anyone with R and BLAST+ installed on their computer. It requires as input only the nontarget, amplified region and a target sequence database that is BLAST-formatted (an example can be found in the folder "Testing"). By running three R functions, users can within minutes generate a list of all possible blocking oligos, a figure showing how many times each oligo “hits” target templates (fewer hits indicates higher specificity for non-targets) and sortable table of the blocking oligos pairs. The goal is to find a primer set that specifically amplifies the non-target organism and avoids the target organisms and which is nested inside "universal" primers that amplify target and non-target organisms. By using the blocking oligos in amplicon sequencing library preparation, non-target amplification can be effectively avoided, saving valuable sequencing depth and increasing diversity discovery.

# Installation:
1. AmpStop uses tools from the following R packages, so you will need them installed:
* tidyr
* ggplot2
* ggpubr

2. AmpStop needs a local install of blastn, which is part of NCBI's BLAST+ package. You can find that here: https://www.ncbi.nlm.nih.gov/books/NBK279671/. You will need to know the path to blastn on your local machine.

3. AmpStop installation itself requires "devtools", so you will need that for the following AmpStop installation command:
```
devtools::install_github("magler1/AmpStop")
library(AmpStop)
```

# Running:

## 1. Set the working directory
First, download the folder "Testing", then set the working directory:
```
setwd("Path/to/Testing/")
```
## 2. Make a list of candidate blocking oligos
The first script divides what we refer to as a "non-target" sequence into all possible oligos (kmers) of a given length (the default is 30). The non-target is the sequence of the organism for which you want to block amplification. In this example it is the Arabidopsis thlaiana Col-0 ITS sequence. These will be output into a single fasta-formatted text file called "Candidate_oligos.fasta". The filename can be changed.
```
?get_candidate_oligos   # To get more information about the script and see all options
get_candidate_oligos(nontarget="Athaliana_ITS1_new_wholelength.fasta", length=30, outfile="Candidate_oligos.fasta")
```

## 3. Blast the candidates against target organism sequences
The next script runs a blast search, aligning each of the possible oligos to a database of what we refer to as "target" sequences to find out how often they migth align with one of these sequences. This depends on the locus you are working with, but could be bacterial 16S, etc. In this example we use the UNITE fungal ITS database.

The default BLAST parameters used here allow fairly poor hits with the idea of erring on the side of caution (it is better to find potential matches than not). The default parameters are:
* percid = 25 (can be adjusted in the parameters)
* wordsize = 7 (can be adjusted in the parameters)
* evalue = 100000 (not currently adjustable)

The script will produce a blast text output containing all alignments for each of the oligos. The default filename for the output is "Candidate_vs.UNITE.out". This filename can be changed.

Note: in Windows, blastn_path will looks something like "C:\\path\\to\\blastn.exe"

```
?blast_candidate_oligos
blast_candidate_oligos(blastn_path="/Path/to/blastn", blastdb="./Unite_s_011217/UNITE_s_011217", candidates="Candidate_oligos.fasta", outfile = "Candidates_vs_UNITE.out", percid = 25, wordsize = 7)
```

## 4. Summarize the results of the BLAST search
The final script parses the blast results and outputs a summary of hit results for each candidate oligo. The output is in both graph and tab-delimited table that can be opened in excel. Essentially, for every possible oligo three parameters are calculated. Ideally, a "good" blocking oligo candidate would minimize these values.:
1. The total number of hits of the oligo to the "target" database. 
2. The total number of hits where >90% of the candidate oligo length aligned to a target sequence. These are probably very problematic regions where amplification is quite likely.
3. The total number of hits where the oligo aligned at the 3' end. Even oligos that only partially align to a template can result in amplfication if they anneal at the 3' end because this is where the polymerase attaches to the oligo/template dimer. Thus, alignments at this end are the most serious and are therefore counted separately. Unfortunately, there is no way to estimate this perfectly from a BLAST alignment. The approach used here can detect alignments where at least 7 oligo bases align to a target, with a couple of possible mismatches. We include the count if it aligns up to 3 bases from the end of the oligo (if your candidates are 30bp, an alignment of 7 bases up to base 27 will be counted as a 3' alignment). This approach should be quite conservative and give a good idea of 3' alignment issues.

```
?check_candidate_oligos
check_candidate_oligos(candidates="Candidate_oligos.fasta", blastoutput="Candidate_vs_UNITE.out")
```
### Example graphical output:
![Example Output Figure](https://github.com/magler1/AmpStop/blob/master/Testing/Oligos_Plot.jpg)

Note that the blue and red dotted lines signify 1/3 and 2/3 of the length of the whole non-target sequence. A left and right primer may be best found outside these two lines.

# What do I do with the results?
In the end, the goal is to find a useful primer pair. You will need to pick a forward and reverse primer out of the candidate oligos. Picking a couple of each is probably desireable for testing purposes. *The reverse primer will need to be reverse complimented or else you won't get amplification.* 

If you've designed an ITS blocking oligo set, as in this example, you will need to check that they amplify the non-target gDNA and do not amplify some examples of target gDNA. This can preliminarily be done by amplification and viewing results on a gel. Full testing will probably require checking a real sample with and without blocking oligos via amplicon sequencing to see if there is an effect on the diversity. 
