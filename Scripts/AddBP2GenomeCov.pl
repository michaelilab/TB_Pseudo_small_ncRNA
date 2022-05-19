#!/usr/local/bin/perl -w

use strict;

#input fasta
open FILE0, $ARGV[0] or die "Can't open fasta file!\n";
my %seq_hash=();
my $seq;

while (<FILE0>){
  chomp;
  if (/^>(.*?)$/){
#  if (/^>(.*?)\s+/){
    $seq=$1;
    #print $seq, "\n";
    $seq_hash{$seq}=();
  }
  else{
      s/\s+//g;
      s/\t//g;
      $seq_hash{$seq}.=$_;
   }
}
close FILE0;

open FILE1, $ARGV[1] or die "Can't open genomecov file!\n";
while (<FILE1>){
  chomp;
  my ($name,$num,$count)=split(/\t/,$_);
  if (exists $seq_hash{$name}){
      my @seq=split(//,$seq_hash{$name});
      print "$name\t$num\t",$seq[$num-1],"\t$count\n";
  }
  else{
      print STDERR $name, "\n";
  }
}
close FILE1;
