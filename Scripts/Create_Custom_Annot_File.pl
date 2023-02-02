#!/usr/local/bin/perl

my %seq_hash=();
my %mod_table;

#read in rRNA sequence file
open FILE0, $ARGV[0] or die "cannot open rRNA sequence file!\n";
my $seq_name;
while (<FILE0>){
  chomp;
  if (/^>(.*)/){
      $seq_name=$1;
      $seq_hash{$seq_name}="";
  }
  else{
      $seq_hash{$seq_name}.=$_;
  }
}
close FILE0;
open FILE1, $ARGV[1] or die "cannot open rRNA description file!\n";
my @rRNA_array=();
my %conversion=();
my %reverse=();

while (<FILE1>){
  chomp;
  my ($rRNA,$start,$end)=split(/\t/,$_);
  for (my $y=$start; $y <= $end; $y++){
    $rRNA_array[$y]=$rRNA;
    $conversion{$y}=$y-$start+1;
    $reverse{$rRNA}{$y-$start+1}=$y;
  }
}
close FILE1;

open FILE2, $ARGV[2] or die "cannot open rRNA modification table file!\n";
my @modified  = (0) x @rRNA_array;
my @sno_array =  (0) x @rRNA_array;
#open FILE2 - known sites
while (<FILE2>){
  chomp;
  my ($rRNA,$mod,$type,$rRNA_seq,$ref,$sno)=split(/\t/,$_);
  $modified[$reverse{$rRNA}{$mod}]=10;
  $sno_array[$reverse{$rRNA}{$mod}]=$sno;
}
close FILE2;

#Put it all together
print "rRNA\tposition\trel_count\tmodified\tsnoRNA\tBP\n";
foreach $a (keys %seq_hash){
    #print $a, "\n";
    my @temp=();
    @temp= split(//,$seq_hash{$a});
    for (my $x = 1; $x < $#temp;$x++){
	my $bp=$temp[$x-1];
	print $rRNA_array[$x], "\t", $x, "\t", $conversion{$x},"\t", $modified[$x],"\t", $sno_array[$x], "\t", $bp, "\n";

    }
}
