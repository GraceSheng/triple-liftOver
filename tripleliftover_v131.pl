#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use File::Basename;
use File::Spec;
use Getopt::Long;

##perl tripleliftover.pl --bim input.bim --base [base build, e.g. hg19] --target [target build] [--download-chainfile] --outprefix outputprefix
##version 1.31
##script assumes liftover is installed in library
##script assumes relevant chain files are installed in library/chainfiles, otherwise it would ask for downloads
##script will output all SNPs and indicate whether they can be directly lifted over or note that they are inverted. SNPs that did not lift over successfully (in the unmapped output of liftover) or that lifted to a different chr will not be included.
##created by - charleston chiang & grace sheng 11/10/2021
##modified by - charleston 11/14/2021 - output lifted over positions too
##modified by - charleston 11/22/2021 - now check final output before deleting intermediate files
##modified by - charleston 4/20/2022 - interface with Jordan's script to download on the fly chainfiles, fixing relative paths, etc.

my ($bimfile, $useroutprefix, $basebuild, $targetbuild, $download);
GetOptions("bim=s" => \$bimfile,
	   "base=s" => \$basebuild,
	   "target=s" => \$targetbuild,
	   "download-chainfile" => \$download,
	   "outprefix=s" => \$useroutprefix);

die "usage: perl $0 --bim input.bim --base base_genome_build --target target_genome_build [--download-chainfile] --outprefix outprefix\n" if ! $bimfile;

print STDERR "\n\n";
print STDERR "===============================================\n";
print STDERR "triple-liftOver, version 1.3\n\n";
print STDERR "please cite Sheng et al. bioRxiv 2022 (doi: 10.1101/2022.02.19.481174) if you use this tool\n\n";
print STDERR "for the input bim file, make sure the variant ID column is populated and each row has an unique identifier!\n";
print STDERR "make sure python can be invoked from your path!\n";
print STDERR "if you don't have the chainfiles necessary for liftover, you can run the script in interactive mode with --download-chainfile option\n";
print STDERR "===============================================\n\n";

my $liftover = "$Bin/library/liftOver";
my $chainfiledir = "$Bin/library/chainfiles";
($basebuild, $targetbuild) = buildnameconvention($basebuild, $targetbuild);
my $chainfile = $basebuild . "To" . $targetbuild . '.over.chain.gz';
my $rstr = randomstring(6);
my ($outprefix, $outdir, $outsuffix) = fileparse(File::Spec->rel2abs($useroutprefix));

print STDERR "analysis started: ", getcurrenttime(), "\n";

##check if chainfile exists:
print STDERR "checking if required chainfile exists [$basebuild => $targetbuild]\n";
if($download){
    print STDERR "interactive mode only to download chainfiles...\n";
    system("python $Bin/main.py $basebuild");
    die "chainfiles download attempted. If successful, please run this script again to proceed if you do not wish to continue in interactive mode (e.g. you need to submit the job to a compute cluster)\n";
}elsif(! -e "$chainfiledir/$chainfile"){
    print STDERR "chainfile not found in $chainfiledir, attempt to download it...\n";
    system("python $Bin/main.py $basebuild");
}

print STDERR "converting bim file into UCSC bed format...\n";
my $tmpbedfile = "$outdir/$outprefix.$rstr.txt";
convertbim2ucscbedformat($bimfile, $tmpbedfile);

print STDERR "running liftOver...\n";
my $tmpliftedfile = "$outdir/$outprefix.$rstr.lifted.txt";
my $tmpunmappedfile = "$outdir/$outprefix.$rstr.unmapped.txt";
my $command = "$liftover $tmpbedfile $chainfiledir/$chainfile $tmpliftedfile $tmpunmappedfile";
system($command);

print STDERR "building the structure of lifted coordinates...\n";
my %db;
readinliftedstructure($tmpliftedfile, \%db);
#print STDERR scalar(keys %db), "\n";

print STDERR "identifying SNPs in inverted regions...\n";
identifyinvertedpos(\%db, "$outdir/$outprefix.$rstr.step1.invr.txt");

print STDERR "tidying up the output...\n";
tidyup("$outdir/$outprefix.$rstr.step1.invr.txt", "$outdir/$outprefix.$rstr.step2.invr.txt");

##sort output by pos
print STDERR "sorting the output file...\n";
system("sort -k2n -k3n -o $outdir/$outprefix.invr.txt $outdir/$outprefix.$rstr.step2.invr.txt");
if(-s "$outdir/$outprefix.invr.txt"){
    system("rm $outdir/*$rstr*");
    print STDERR "analysis finished: ", getcurrenttime(), "\n";
}




sub tidyup {
    my $tmpinfile = shift @_;
    my $tmpoutfile = shift @_;
    open(IN, $tmpinfile) || die "cannot open $tmpinfile: $!\n";
    my %finaldb;
    while(<IN>){
	chomp;
	my ($chr, $pos0, $pos, $pos1, $id, $cat) = split("\t", $_); #cat has 'inverted' or 'lifted'
	$chr =~ s/chr//;
	$finaldb{$id}[0] = $chr;
	$finaldb{$id}[1] = $pos;
	$finaldb{$id}[2] = $cat;
    }
    close IN;

    ##check to make sure the SNP is still mapped to the same chromosome
    open(OUT, "> $tmpoutfile") || die "cannot open $tmpoutfile: $!\n";
    print OUT join "\t", "snpid", "chr", "pos_[$targetbuild]", "pos_[$basebuild]", "category"; print OUT "\n";    
    open(BIM, $bimfile) || die "cannot open $bimfile: $!\n";
    while(<BIM>){
	chomp;
	my ($chr, $id, $gd, $pos, $a1, $a2) = split('\s+', $_);
	$chr = chrconvert($chr);
	if(defined($finaldb{$id})){
	    if($finaldb{$id}[0] eq $chr){
		print OUT join "\t", $id, $finaldb{$id}[0], $finaldb{$id}[1], $pos, $finaldb{$id}[2];
		print OUT "\n";
	    }
	}
    }
    close BIM;
    close OUT;
}

sub identifyinvertedpos {
    my $d = shift @_;
    my $tmpoutfile = shift @_;
    open(OUT, "> $tmpoutfile") || die "cannot open $tmpoutfile: $!\n";
    foreach my $snp (keys %$d){
	if(!defined($$d{$snp}[2])){ ##values are chr, pos0, pos, pos1
	    next; #the current snp did not liftover, skip
	}elsif($$d{$snp}[1] ne "NA" && $$d{$snp}[1] - $$d{$snp}[2] == 1){ #previous bp inverted
	    print OUT join "\t", @{$$d{$snp}}, $snp, "inverted";
	    print OUT "\n";
	}elsif($$d{$snp}[3] ne "NA" && $$d{$snp}[3] - $$d{$snp}[2] == -1){ #succeeding bp inverted
	    print OUT join "\t", @{$$d{$snp}}, $snp, "inverted";
	    print OUT "\n";
	}else{
	    print OUT join "\t", @{$$d{$snp}}, $snp, "lifted";
	    print OUT "\n";
	}
    }
    close OUT;
}

sub readinliftedstructure {
    my $liftedfile = shift @_;
    my $d = shift @_;
    open(FILE, $liftedfile) || die "cannot open $liftedfile: $!\n";
    while(<FILE>){
	chomp;
	my ($chr, $pos0, $pos, $id) = split("\t", $_);
	my @ids = split("_", $id);
	my $tag = pop @ids;
	my $realid = join "_", @ids;
	if($tag eq 'b'){ #the snp we want
	    $$d{$realid}[0] = $chr; #record its chr; if it's different from chr in the original build, it will be removed later
	    $$d{$realid}[2] = $pos;
	    $$d{$realid}[1] = "NA";
	    $$d{$realid}[3] = "NA";
	}elsif(defined($$d{$realid}) && $chr eq $$d{$realid}[0]){ #preceding and succeeding bp will always come after the main bp of interest
	    if($tag eq 'a'){ #preceding bp
		$$d{$realid}[1] = $pos;
	    }elsif($tag eq 'c'){ #succeeding bp
		$$d{$realid}[3] = $pos;
	    }
	}
    }
    close FILE;
}

sub convertbim2ucscbedformat {
    my $bimfile = shift @_;
    my $outbedfile = shift @_;
    open(BIM, $bimfile) || die "cannot open $bimfile: $!\n";
    open(OUT, "> $outbedfile");
    while(<BIM>){
        chomp;
        my ($chr, $id, $gmap, $pos, $a1, $a2) = split('\s+', $_);
        my $pos0 = $pos - 1;
	$chr = chrconvert($chr);

	#this is the SNP we want
        print OUT "chr$chr\t$pos0\t$pos\t$id" . "_b\n";

	#preceding SNP
	$pos0 --; $pos --;
	print OUT "chr$chr\t$pos0\t$pos\t$id" . "_a\n";

	#succeeding SNP
	$pos0 += 2; $pos += 2;
	print OUT "chr$chr\t$pos0\t$pos\t$id" . "_c\n";
    }
    close BIM;
    close OUT;
}

sub buildnameconvention {
    ##UCSC naming convention, base build is always lower case
    ##target build is capitalize first letter or first and third letter
    my $b = shift @_;
    my $t = shift @_;
    $b = lc($b);
    my @t = split("", $t);
    $t[0] = uc($t[0]);
    if(scalar(@t) > 6){
	$t[3] = uc($t[3]);
    }
    my $nt = join "", @t;
    return($b, $nt);
}

##misc utilities
sub getcurrenttime {
    my $format = shift @_;
    my @f = localtime(time);
    $f[2] = "0" . $f[2] if $f[2] < 10;
    $f[1] = "0" . $f[1] if $f[1] < 10;
    $f[0] = "0" . $f[0] if $f[0] < 10;
    my $time;
    if(defined($format) && $format eq "fn"){ #better to use in filename
        $time = join "", $f[5]+1900, "_", $f[4]+1, "_", $f[3], ".$f[2]-$f[1]-$f[0]";
    }else{
        $time = join "", $f[5]+1900, '/', $f[4]+1, '/', $f[3], ",", " $f[2]:$f[1]:$f[0]";
    }
    return($time);
}

sub randomstring {
    my $len = shift;
    my @chars = ("A".."Z", "a".."z", 1..9, 0);
    my $string;
    $string .= $chars[rand @chars] for 1..$len;
    return($string);
}

sub chrconvert {
    my $c = shift @_;
    if($c eq '23' || $c eq '25'){
	return("X");
    }elsif($c eq '24'){
	return("Y");
    }elsif($c eq '26'){
	return("M");
    }elsif($c eq 'XY'){
	return("X");
    }elsif($c eq 'MT'){
	return("M");
    }else{
	return($c);
    }
}
