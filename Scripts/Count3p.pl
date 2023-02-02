#!/usr/local/bin/perl -w

use strict;

#how to run this file:
#cat input.bed|perl Count3p.pl input_genome.genome
#input bed file (BAmtobed)
#input Genome size file


my %genome_hash=();
my $genome_size=0;

open FILE1, $ARGV[0] or die "Can't open genome file!\n";

#create an array for each chromosome - each bp initialized to 0
#genome file format chrom_name<tab>size_of_chrom
while (<FILE1>){
  chomp;
  my ($name,$genome_size)=split(/\t/,$_);
  @{$genome_hash{$name}}=(0) x $genome_size;
}
close FILE1;


while (<STDIN>){
  chomp;
  my ($chrom,$start,$end,$read)=split(/\t/,$_);
  $genome_hash{$chrom}[$end]+=1;
}

foreach my $g (sort keys %genome_hash){
  for (my $x=0; $x <= $#{$genome_hash{$g}};$x++){
    if (defined $genome_hash{$g}[$x]){
      print $g, "\t", $x, "\t",  $genome_hash{$g}[$x], "\n";
    }
    else{
      print $g, "\t", $x, "\t", "0\n";
    }
  }
}
