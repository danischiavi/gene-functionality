#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $k = 2; #size of the k-mer 
my ($inFile, $help, $verbose); 

&GetOptions( 
    "h|help"              => \$help,
    "v|verbose"           => \$verbose,
    "i|seqfile=s"         => \$inFile, 
    "k|kmer=s"            => \$k
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

my ($seq, $name) = ('AACAGATCCGCTGGTTA','seqName');
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

my %kmerFreq = %{extractNmerFreqs($checkedSeq, $k  )};
my %jmerFreq = %{extractNmerFreqs($checkedSeq, $k-1)};
foreach my $kmer (sort keys %kmerFreq) {
    my $jmer1 = substr($kmer, 0, $k-1);
    my $jmer2 = substr($kmer, 1, $k-1);

    if( $kmerFreq{$kmer} && $jmerFreq{$jmer1} && $jmerFreq{$jmer2} ){
    	my $bits = ( log($kmerFreq{$kmer}) - log($jmerFreq{$jmer1}) - log($jmerFreq{$jmer2}) )/log(2);
    	printf "%0.2f,", $bits;
    	printf "\t\t$kmer\[%0.2f\]\t$jmer1\[%0.2f\]\t$jmer2\[%0.2f\]\n", $kmerFreq{$kmer}, $jmerFreq{$jmer1}, $jmerFreq{$jmer2} if defined($verbose);
	
    }
}

exit(0);   
######################################################################

#extractNmerFreqs: return a pointer to a hash of nmer counts for an input sequence 
sub extractNmerFreqs {
    my $seq  = shift;
    my $kmerLen = shift; 
    my $len = length($seq);
    my %kmerFreq;
    for (my $i = 0; $i < $len+1-$kmerLen; $i++){
    	my $subseq = substr($seq, $i, $kmerLen);
    	if(defined($kmerFreq{$subseq})){
    	    $kmerFreq{$subseq}++;
    	}
    	else{
    	    $kmerFreq{$subseq}=1;
    	}
    }

    foreach my $kmer (keys %kmerFreq) {
	   $kmerFreq{$kmer} = $kmerFreq{$kmer}/($len+1-$kmerLen);
    }
	
    return \%kmerFreq; 
}



#########################
#calculate_genomic_signature: return a pointer to a hash of genomic signature for input sequence
sub calculate_genomic_signature {
  my $seq = shift;

  # Calculate individual nucleotide frequencies
  my %nucleotide_freq;
  for (my $i = 0; $i < length($seq); $i++) {
    my $nucleotide = substr($seq, $i, 1);
    $nucleotide_freq{$nucleotide}++;
  }
  # Normalize nucleotide frequencies
  foreach my $nucleotide (keys %nucleotide_freq) {
    $nucleotide_freq{$nucleotide} /= length($seq);
  }

  # Calculate dinucleotide frequencies
  my %dinucleotide_freq;
  for (my $i = 0; $i < length($seq) - 1; $i++) {
    my $dinucleotide = substr($seq, $i, 2);
    $dinucleotide_freq{$dinucleotide}++;
  }

  # Calculate genomic signatures
  my %genomic_signatures;
  foreach my $dinucleotide (keys %dinucleotide_freq) {
    my $nuc1 = substr($dinucleotide, 0, 1);
    my $nuc2 = substr($dinucleotide, 1, 1);
    my $signature = $dinucleotide_freq{$dinucleotide} / 
                    ($nucleotide_freq{$nuc1} * $nucleotide_freq{$nuc2});
    $genomic_signatures{$dinucleotide} = $signature;
  }

  return \%genomic_signatures;
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

markovProperties.pl: compute bitscores for each kmer of a given length.

Usage:   markovProperties.pl -i <seqfile> -k 2
Options:       -h|--help                     Show this help.

               -i|--seqfile <file>           a fasta file 
               -k|--kmer    <num>            kmer value, >1 
TODO:

EOF
}
