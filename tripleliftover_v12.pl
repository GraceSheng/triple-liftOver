#!/usr/bin/perl -w
use strict;
use Getopt::Long;

##perl tripleliftover.pl --bim input.bim --outprefix outputprefix
##version 1.2
##script will output all SNPs and indicate whether they can be directly lifted over or note that they are inverted. SNPs that did not lift over successfully (in the unmapped output of liftover) or that lifted to a different chr will not be included.
##created by - charleston chiang & xin sheng 11/10/2021
##modified by - charleston 11/14/2021 - output lifted over positions too
##modified by - charleston 11/22/2021 - now check final output before deleting intermediate files

my ($bimfile, $outprefix);
GetOptions("bim=s" => \$bimfile,
	   "outprefix=s" => \$outprefix);

die "usage: perl $0 --bim input.bim --outprefix outprefix\n" if ! $bimfile;
warn "for the input bim file, make sure the variant ID column is populated and each row has an unique identifier!\n";

my $liftover = "/project/chia657_28/programs/liftOver/liftOver";
my $chainfile = "/project/chia657_28/programs/liftOver/hg19ToHg38.over.chain.gz";
my $rstr = randomstring(6);

print STDERR "analysis started: ", getcurrenttime(), "\n";
print STDERR "converting bim file into UCSC bed format...\n";
my $tmpbedfile = "$outprefix.$rstr.txt";
convertbim2ucscbedformat($bimfile, $tmpbedfile);

print STDERR "running liftOver...\n";
my $tmpliftedfile = "$outprefix.$rstr.lifted.txt";
my $tmpunmappedfile = "$outprefix.$rstr.unmapped.txt";
my $command = "$liftover $tmpbedfile $chainfile $tmpliftedfile $tmpunmappedfile";
system($command);

print STDERR "building the structure of lifted coordinates...\n";
my %db;
readinliftedstructure($tmpliftedfile, \%db);
#print STDERR scalar(keys %db), "\n";

print STDERR "identifying SNPs in inverted regions...\n";
identifyinvertedpos(\%db, "$outprefix.$rstr.step1.invr.txt");

print STDERR "tidying up the output...\n";
tidyup("$outprefix.$rstr.step1.invr.txt", "$outprefix.$rstr.step2.invr.txt");

##sort output by pos
print STDERR "sorting the output file...\n";
system("sort -k2n -k3n -o $outprefix.invr.txt $outprefix.$rstr.step2.invr.txt");
if(-s "$outprefix.invr.txt"){
    system("rm *$rstr*");
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
    print OUT join "\t", "snpid", "chr", "pos_hg38", "pos_hg19", "category"; print OUT "\n";    
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
