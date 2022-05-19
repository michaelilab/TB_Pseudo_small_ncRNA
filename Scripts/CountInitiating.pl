#!/usr/local/bin/perl -w

use strict;
#input bed file (BAmtobed)
#input Genome size file
#my $name="lcl|Tb927_02_v4:259325-270494";
#my $genome_size=10258;
#my @array=(0) x $genome_size;

my %genome_hash=();
my $genome_size=0;

open FILE1, $ARGV[0] or die "Can't open control genome file!\n";

while (<FILE1>){
  chomp;
  my ($name,$genome_size)=split(/\t/,$_);
  @{$genome_hash{$name}}=(0) x $genome_size;
}
close FILE1;

while (<STDIN>){
  chomp;
  my ($chrom,$start,$end,$read)=split(/\t/,$_);
  $genome_hash{$chrom}[$start]+=1;
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
