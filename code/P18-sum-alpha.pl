#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
# ------------------------
my $outDir = "../analysis/P18-sum-alpha";
`rm -r $outDir`;
`mkdir $outDir`;
# ------------------------
my $files = `ls ../analysis/P15-clust-tree-cm-goods/*/exp*/alpha-diversity.tsv`;
chomp($files);
my @files = split "\n", $files;
# ------------------------
my %totals  = ();
my %data    = ();
my %allestimators = ();
foreach my $f (@files){
  my @f      = split /\//, $f;
  my $v      = $f[3];
  my $est    = $f[4];
  $est       =~ s/\_vector.qza.exp//g;
  $allestimators{$est} = 1;
  my @header = ();
  open IN, "$f" or die;
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    if (!defined($header[1])){
      @header = ("SampleID", @A);
    }else{
      for my $i (1 .. $#A){
        $data{$v}{$A[0]}{$est} = $A[$i];
      }
    }
  }
  close IN;
}
# print dataset ---------------------
open OUT, ">$outDir/alpha-diversity-per-sample-region.txt" or die "Error: Cannot open $outDir/alpha-diversity-per-sample-region.txt\n";
print OUT "Region\tSampleID";
foreach my $t (sort keys %allestimators){
  print OUT "\t$t";
}
print OUT "\n";
foreach my $v (sort keys %data){
  foreach my $s (sort keys %{$data{$v}}){
    print OUT "$v\t$s";
    foreach my $t (sort keys %allestimators){
      print OUT "\t$data{$v}{$s}{$t}";
    }
    print OUT "\n";
  }
}
close OUT;
