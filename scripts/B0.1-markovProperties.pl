#!/usr/bin/perl

#markovProperties:
#        Calculate bitscores for each kmer of a given length. Frequencies are provided by an input sequence. 

use warnings;
use strict;
use Getopt::Long;

my $k = 2; #size of the k-mer
my $alpha = 'ACGT'; # alphabet to use 
my ($inFile, $help, $verbose); 

&GetOptions( 
    "h|help"              => \$help,
    "v|verbose"           => \$verbose,
    "i|seqfile=s"         => \$inFile, 
    "k|kmer=s"            => \$k, 
    "a|alphabet=s"        => \$alpha
    );

if( $help ) {
    &help();
    exit(1);
}

die "FATAL: [$k] must be greater than 1!" if $k <= 1; 

#Read sequence from file:
#Pad with a pseudosequence that contains all 2mers: 
# AACAGAT
#   CCGCT
#     GGT
#TA    

#AACAGATCCGCTGGTTA    

#my ($seq, $name) = ('AACAGATCCGCTGGTTA','seqName'); #This added pseudocounts for each k-mer 
my ($seq, $name) = ('','seqName');
open(IN,  "< $inFile") || die "FATAL: cannot open [$inFile]\n[$!]";
while(my $fa = <IN>){
    chomp($fa);
     if($fa =~ /^>(\S+)/){
	 $name = $1;  
     }
     else{
	 $seq .= $fa;
     }
}
close(IN);

#Check sequence is legal:
my $checkedSeq = ''; 
my @seq = split(//, $seq); 
foreach my $s (@seq){
    if( isDNA($s)){
	$checkedSeq .= $s;
    }
    else {
	print STDERR "Warning: [$s] is not a legal DNA character!\n";
    }
}
die "FATAL: [$inFile] contains no legal DNA sequence info!" if length($checkedSeq) == 0; 
$checkedSeq =~ tr/a-z/A-Z/;

#printf "%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", log(10), log(1000), log(16), log(2.718282), log(16)/log(2); 

my %kmerFreq = %{extractNmerFreqs($checkedSeq, $k  , $alpha)};
my %jmerFreq = %{extractNmerFreqs($checkedSeq, $k-1, $alpha)};
my $psK = 0.01;        #A pseudocount for each kmer 
my $psJ = sqrt($psK);  #A pseudocount for each k-1mer 

foreach my $kmer (sort keys %kmerFreq) {
    my $jmer1 = substr($kmer, 0, $k-1);
    my $jmer2 = substr($kmer, 1, $k-1);
    
    
    if( defined($kmerFreq{$kmer}) && defined($jmerFreq{$jmer1}) && defined($jmerFreq{$jmer2}) ){
	my $bits = ( log($kmerFreq{$kmer} + $psK) - log($jmerFreq{$jmer1} + $psJ) - log($jmerFreq{$jmer2} + $psJ) )/log(2);
	printf "$kmer\t%0.2f", $bits;
	printf "\t\t$kmer\[%0.2f\]\t$jmer1\[%0.2f\]\t$jmer2\[%0.2f\]\t\t", $kmerFreq{$kmer}, $jmerFreq{$jmer1}, $jmerFreq{$jmer2} if defined($verbose);
	#printf "( log($kmerFreq{$kmer} + $psK) - log($jmerFreq{$jmer1} + $psJ) - log($jmerFreq{$jmer2} + $psJ) )/log(2);" if defined($verbose); #DEBUGGING
	print "\n"; 
    }
}

exit(0);   
######################################################################

#extractNmerFreqs: return a pointer to a hash of nmer counts for an input sequence 
sub extractNmerFreqs {
    my $seq     = shift;
    my $kmerLen = shift;
    my $alpha   = shift;
    my $len = length($seq);
    my @alpha = split(//, $alpha); 
    
    my $kmerFreq = initialiseKmerHash($kmerLen, $alpha);
    my %kmerFreq = %{$kmerFreq};
    
    if($kmerLen==1){
	foreach my $a (@alpha){
	    $kmerFreq{$a}=0;           #initialise hash, make sure "0" freq nucleotides are captured 
	}
    }
    elsif($kmerLen==2){
	foreach my $a (@alpha){
	    foreach my $a2 (@alpha){
		$kmerFreq{$a . $a2}=0; #initialise hash, make sure "0" freq dinucleotides are captured 
	    }
	}
    }
    
    for (my $i = 0; $i < $len+1-$kmerLen; $i++){
	my $subseq = substr($seq, $i, $kmerLen);
	if(defined($kmerFreq{$subseq})){
	    $kmerFreq{$subseq}++;
	}
	else{
	    $kmerFreq{$subseq}=0;
	}
    }
    print "Raw counts:\n" if defined($verbose);
    foreach my $kmer (sort keys %kmerFreq) {
	printf "$kmer\t$kmerFreq{$kmer}\n" if defined($verbose);
	$kmerFreq{$kmer} = $kmerFreq{$kmer}/($len+1-$kmerLen);
    }
	
    return \%kmerFreq; 
}

######################################################################
#initialiseKmerHash: Given a kmer length and alphabet, initialise a hash with all possible kmers of the given length & alphabet 
sub initialiseKmerHash {
    my $kmerLen = shift;
    my $alpha   = shift;
    my %kmerFreq; 

    my @alpha = split( //, $alpha );
    my @kmers = @alpha;
    for (my $i = 0; $i < ($kmerLen-1); $i++){
	my @kmersLonger = ();
	foreach my $k (@kmers){
	    foreach my $a (@alpha){
		push(@kmersLonger, $k . $a);
		#print "[$k]\t[$a]\n";
	    }
	}
	@kmers =  @kmersLonger; #sort()
    }
    
    foreach my $k (@kmers){
	$kmerFreq{$k} = 0;
    }
    return \%kmerFreq;
}

######################################################################
#returns true if input character is a DNA nucleotide (excludes IUPAC codes):
sub isDNA {
    my $a = shift;
    
    if (defined($a)){
        $a =~ tr/a-z/A-Z/;
    }
    
    if (defined($a) && length($a) && ($a =~ /[ACGT]/) ){
        return 1;
    }
    else {
        return 0;
    }
    
}


######################################################################
sub help {
    print STDERR <<EOF;

    markovProperties.pl: compute bitscores for each kmer of a given
	length over an input sequence.  Background probabilities are
	computed using 'k-1'mers.  Pseudocounts are added to each
	frequency to address possible log(0) issues.

      Usage:
	markovProperties.pl -i <seqfile> -k 2 -a 'ACGT'

      Options:
	       -h|--help Show this help.
               -i|--seqfile  <file>           a fasta file 
               -k|--kmer     <num>            kmer value, >1 
               -a|--alphabet <string>         the alphabet to be used, e.g. ACGT

	Returns log2( F_k / (F_{k-1} * F_{k+1}) ) for each kmer. 

	Where F_{k-1} is the frequency of the first k-1mer of each kmer,
	and F_{k+1} is the frequency of the last k-1mer of each kmer.
	
TODO:


    Paul P. Gardner
    
    
EOF
}
