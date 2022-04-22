# *triple-liftOver*
*triple-liftOver* is a utility that takes an input file in PLINK .bim format (https://www.cog-genomics.org/plink/1.9/) and converts the genomic coordinate for three consecutive bases at each chromosomal position between genome builds using UCSC’s tool liftOver (https://genome.ucsc.edu/cgi-bin/hgLiftOver). It outputs the new positions in the destination build with a category column indicating whether a variant is found in a Between-Builds Inverted Sequence (BBIS) region. Variants that cannot be lifted over (in the unmapped output of liftOver) or are no longer on the same chromosome are excluded from the output file. For a focal SNV to qualify as an inverted site, the heuristic requires either the succeeding base in the old build to become the preceding base in the new build or the preceding base in the old build to become the succeeding base in the new build.

## Usage ## 
`perl tripleliftover_v133.pl --bim input.bim --base base_genome_build --target target_genome_build [--download-chainfile] --outprefix outprefix`

`[--download-chainfile]` tag is optional, and it allows an interactive run of the program so liftOver chain files can be downloaded to library/chainfiles/. Once the chain files are downloaded, a 2nd run of the program without this tag would generate the final output file for *triple-liftOver*. Note that downloading the chain files requires internet access. If you are running a job on the compute node without internet access, please follow this 2-step approach. Otherwise, the *triple-liftOver* program can be run in one step without using this tag.

## Output ##
Upon successful execution of the script, an output file named "outprefix.invr.txt" will be generated. For the example run, the output file is 1000GP_Phase3_chr10.test.invr.txt. It is a 5-column tab delimited text file which contains all the variants that can be lifted over to the same chromosome by UCSC liftOver. These five columns are
- snpid: original variant name from the bim file
- chr: chromosome
- pos_[target_genome_build]: position in the destination build (target_genome_build is Hg38 in the example)
- pos_[base_genome_build]: position in the source build (base_genome_build is hg19 in the example)
- category: either "inverted" (is an inverted site) or "lifted" (is NOT an inverted site or cannot be deduced based on *triple-liftOver* heuristic)

## Brief Example of Using *triple-liftOver* ##
To use *triple-liftOver* to correctly convert the genomic coordinates of a PLINK dataset with variants in hg19 genome forward (PLUS) strand alleles to PLINK dataset with coordinates in hg38 genome forward (PLUS) strand alleles:
- Run *triple-liftOver* for your bim file with hg19 as base build and hg38 as target build: 

`perl tripleliftover_v133.pl --bim 1000GP_Phase3_chr10.test.bim --base hg19 --target hg38 --outprefix 1000GP_Phase3_chr10.test`

- Then use the variant name in the output to exclude any variant that is not in the *triple-liftOver* output (these variants either failed the lift over or were lifted to a different chromosome). It can be executed with a command in PLINK like this:

`plink --bfile 1000GP_Phase3_chr10.test --extract variant.name.in.liftOver.output --make-bed --out 1000GP_Phase3_chr10.test.lifted`

- Use the variant name and pos_[Hg38] from the output to update the positions in the hg19 PLINK files. This step generates the hg38 PLINK file with genome forward strand alleles (with the inverted site issue between hg19 and hg38):

`plink --bfile 1000GP_Phase3_chr10.test.lifted --update-map variant.name.hg38.position --make-bed --out 1000GP_Phase3_chr10.test.lifted.hg38`

- Parse out the “inverted” variant names from *triple-liftOver* output, limit them to SNVs only (as it may include indels). Execute the strand flip command in PLINK for these inverted SNVs. This will generate the final hg38 PLINK dataset for imputation (fixes the inverted site issue between hg19 and hg38):

`plink --bfile 1000GP_Phase3_chr10.test.lifted.hg38 --flip inverted.SNV.name --make-bed --out 1000GP_Phase3_chr10.test.lifted.hg38.final`

## Notes ##
- Current repository includes the following files
  - tripleliftover_v133.pl: the latest *triple-liftOver* program
  - main.py: downloader program called by tripleliftover_v133.pl if the chain file between user specified base and target genome builds is not available from library/chainfiles/, it will prompt the user with all the available target genome build chain files for his base build so he can choose which one to download.
  - urlParser.py: URL parsing program which parses the UCSC liftOver chain file download web page and is called by the main.py.
  - library/liftOver: precompiled UCSC liftOver program
  - library/chainfiles/
    - hg18ToHg38.over.chain.gz: hg18 to hg38 liftOver chain file
    - hg19ToHg38.over.chain.gz: hg19 to hg38 liftOver chain file
  - Examples/
    - 1000GP_Phase3_chr10.test.bim : one example PLINK bim file
    - 1000GP_Phase3_chr10.test.invr.txt: output file from running *triple-liftOver* for the above bim file
-  After cloning the repository, please add executable permission to the pre-compiled “liftOver” program which is under library/ 
-  For the input bim file, please make sure the variant ID column is populated and each row has a unique identifier. Otherwise, the results for these locations may not be as expected.
-  The *triple-liftOver* program does not use any allele info from the bim file for its procedure, nor does it verify these alleles against the reference alleles in the source and destination genome reference builds. Its heuristic finds the inverted sites between source and target genome reference builds based on the relative position change between the focal SNV and its two neighboring bases. If your input bim file contains genome forward (PLUS) strand alleles in the source build, you can use this program to identify any SNV site falls within BBIS regions. For these inverted SNVs to have the genome forward strand (PLUS) alleles in the destination build, their original alleles need to be converted to the reverse complement allele (i.e. A to a T, C to a G, and vice versa). Current algorithm treats an indel as a SNV and may flag it as inverted, if you don't want to flip strand for indels, please exclude them from your strand flip list. 

## Citation ##
Please cite Sheng et al. bioRxiv 2022 (doi: 10.1101/2022.02.19.481174) if you use this tool.
