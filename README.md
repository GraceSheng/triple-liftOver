# triple-liftOver

triple-liftOver is a PERL script which takes an input file in PLINK .bim format, lifts over three consecutive bases at each chromosomal position between genome builds using UCSC’s tool liftOver (https://genome.ucsc.edu/cgi-bin/hgLiftOver) and outputs the new positions in the destination build with a category column indicating whether it is an inverted site. Variants that cannot be lifted over (in the unmapped output of liftOver) or are no longer on the same chromosome are excluded from the output file. For a focal SNP to qualify as an inverted site, it requires either the succeeding base in the old build to become the preceding base in the new build or the preceding base in the old build to become the succeeding base in the new build.

## Usage ##
perl tripleliftover_v12.pl --bim input.bim --outprefix outprefix

## Notes ##

1. Currently the locations for UCSC liftOver program and the chain file "hg19ToHg38.over.chain.gz" are hardwired in the script.
2. triple-liftOver does not verify alleles against the reference allele in the source and destination reference genome builds. It is presumed that your input bim file contains alleles on the genome forward strand in the source build and the output flags any SNV site falls within Between-Builds Inverted Sequence (BBIS) regions. For these inverted SNVs to have the forward strand alleles in the destination build, their original alleles need to be flipped.
