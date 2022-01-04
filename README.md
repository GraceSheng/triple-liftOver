# *triple-liftOver*

*triple-liftOver* is a PERL script which takes an input file in PLINK .bim format, lifts over three consecutive bases at each chromosomal position between human genome reference builds using UCSC’s tool liftOver (https://genome.ucsc.edu/cgi-bin/hgLiftOver) and outputs the new positions in the destination build with a category column indicating whether it is an inverted site. Variants that cannot be lifted over (in the unmapped output of liftOver) or are no longer on the same chromosome are excluded from the output file. For a focal SNP to qualify as an inverted site, it requires either the succeeding base in the old build to become the preceding base in the new build or the preceding base in the old build to become the succeeding base in the new build.

## Usage ##
perl tripleliftover_v12.pl --bim input.bim --outprefix outprefix

## Output
Upon successful completion of the script, an output file  named "outprefix.invr.txt" will be generated. It is a 5 column tab delimited text file which contains all the variants that can be lifted over to the same chromosome by UCSC liftOver. These five columns are
- snpid: original variant name from the bim file
- chr: chromosome
- pos_hg38: position in destination build (current setting is GRCh38/hg38)
- pos_hg19: position in source build (current setting is GRCh37/hg19)
- category: either "inverted" (is an inverted site) or "lifted"(is NOT an inverted site or cannot be deduced based on current heuristic) 

## Notes ##

-  Currently the locations for UCSC liftOver program and the chain file "hg19ToHg38.over.chain.gz" are hardwired in the script.
-  For the input bim file, please make sure the variant ID column is populated and each row has an unique identifier. Otherwise the results for these locations may not be as expected.
-  *triple-liftOver* does not verify alleles against the reference alleles in the source and destination reference genome builds. It is presumed that your input bim file contains alleles on the genome forward strand in the source build and the output flags any SNV site falls within Between-Builds Inverted Sequence (BBIS) regions in the destination build. For these inverted SNVs to have the forward strand alleles in the destination build, their original alleles need to be flipped.
