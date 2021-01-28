#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
# ------------------------
my $outDir = "../analysis/P03-summarize-qc";
`rm -r $outDir`;
`mkdir $outDir`;
# ------------------------
my $files = `ls ../analysis/P02-dada2/v*/exported-denoise-stats/chip*/stats.tsv`;
chomp($files);
my @files = split "\n", $files; # now the @files array holds the path to each sample
# ------------------------
my %data        = (); # data holds: {region}{chip}{sample}{measure} = value
my @allmeasures = ();
my %sample2chip = ();
foreach my $f (@files){
  my @f      = split /\//, $f;
  my $v      = $f[3];
  my $chip   = $f[5];
  $chip      = substr($chip,0,5);

  my @header = ();
  open IN, "$f" or die;
  while(<IN>){
    chomp($_);
    $_ =~ s/\r//g;
    my @A = split "\t", $_;
    if (!defined($header[1])){
      @header = @A;
      @allmeasures = @A;
    }else{
      next if ($_ =~ /^#/);
      for my $i (1 .. $#A){
        $data{$v}{$A[0]}{$header[$i]} = $A[$i];
      }
      $sample2chip{$A[0]} = $chip;
    }
  }
  close IN;
}

# print dataset ---------------------
open OUT, ">$outDir/qc-per-sample-region.txt" or die "Error: Cannot open $outDir/qc-per-sample-region.txt\n";
print OUT "Region\tSampleID\tChip";
foreach my $t (@allmeasures[1..$#allmeasures]){
  my $newt = $t;
  $newt =~ s/\-/\./g;
  $newt =~ s/\ /\./g;
  print OUT "\t$newt";
}
print OUT "\n";
foreach my $v (sort keys %data){
  foreach my $s (sort keys %{$data{$v}}){
    print OUT "$v\t$s\t$sample2chip{$s}";
    foreach my $t (@allmeasures[1..$#allmeasures]){
      if (!defined($data{$v}{$s}{$t})){
        print "$v $s $t\n";
      }
      print OUT "\t$data{$v}{$s}{$t}";
    }
    print OUT "\n";
  }
}
close OUT;
