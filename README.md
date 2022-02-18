# *triple-liftOver*

*triple-liftOver* is a PERL script which takes an input file in PLINK .bim format and converts the genomic coordinate for three consecutive bases at each chromosomal position between genome builds using UCSC’s tool liftOver (https://genome.ucsc.edu/cgi-bin/hgLiftOver). It outputs the new positions in the destination build with a category column indicating whether a variant is found in a Between-Builds Inverted Sequence (BBIS) region. Variants that cannot be lifted over (in the unmapped output of liftOver) or are no longer on the same chromosome are excluded from the output file. For a focal SNV to qualify as an inverted site, the heuristic requires either the succeeding base in the old build to become the preceding base in the new build or the preceding base in the old build to become the succeeding base in the new build. 

## Usage ##
perl tripleliftover_v12.pl --bim input.bim --outprefix outprefix

## Example run ##
perl tripleliftover_v12.pl --bim 1000GP_Phase3_chr10.test.bim --outprefix 1000GP_Phase3_chr10.test

## Output
Upon successful completion of the script, an output file  named "outprefix.invr.txt" will be generated. For the example run, the output file is 1000GP_Phase3_chr10.test.invr.txt. It is a 5 column tab delimited text file which contains all the variants that can be lifted over to the same chromosome by UCSC liftOver. These five columns are
- snpid: original variant name from the bim file
- chr: chromosome
- pos_hg38: position in the destination build (current setting is GRCh38/hg38)
- pos_hg19: position in the source build (current setting is GRCh37/hg19)
- category: either "inverted" (is an inverted site) or "lifted"(is NOT an inverted site or cannot be deduced based on current heuristic) 

## Notes ##

-  Currently the locations for UCSC liftOver program (line 19) and the chain file "hg19ToHg38.over.chain.gz" (line 20) are hardwired in the script. Both the program and different chain files can be downloaded from UCSC (https://genome.ucsc.edu/cgi-bin/hgLiftOver) by following the links under "Command Line Tool".
-  For the input bim file, please make sure the variant ID column is populated and each row has an unique identifier. Otherwise the results for these locations may not be as expected.
-  *triple-liftOver* does not verify alleles against the reference alleles in the source and destination reference genome builds. It is presumed that your input bim file contains alleles on the genome forward strand in the source build and the output flags any SNV site falls within BBIS regions. For these inverted SNVs to have the forward strand alleles in the destination build, their original alleles need to be flipped.
-  Current algorithm treats an indels as a SNV and may flag it as inverted, if you don't want to flip strand for indels, exclude them from your strand flip list
 
 ## Citation ##
