#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
# ------------------------
my $outDir = "../analysis/P07-sum-taxa-db4.0.pl";
`rm -r $outDir`;
`mkdir $outDir`;


my @regions = ("v2", "v3", "v4", "v67", "v8", "v9");
foreach my $v (@regions){
  # make directory for results
  print "$v...\n";
  `mkdir $outDir/$v`;
  # ------------------------
  # specify path to feature biom file
  my $file = "../analysis/P06-convert-tax-class/$v/exported-feature-table/feature-table.biom";
  my $tx   = "../analysis/P06-convert-tax-class/$v/taxonomy.tsv";
  # ------------------------
  # convert to tsv --
  `biom convert -i $file -o $outDir/$v/feature-table.txt --to-tsv`;
  # load into data
  my %data    = ();
  my %totals  = ();
  my @header  = ();
  open IN, "$outDir/$v/feature-table.txt" or die;
  while(<IN>){
    chomp($_);
    next if ($_ =~ /Constructed/);
    my @A = split "\t", $_;
    if (!defined($header[1])){
      @header = @A;
    }else{
      for my $i (1 .. $#A){
        $data{$header[$i]}{$A[0]} += $A[$i];
        $totals{$header[$i]}      += $A[$i];
      }
    }
  }
  close IN;

  # collection names of OTUs from full file
  my %otumap = ();
  open IN, "$tx" or die;
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    next if ($A[1] eq "taxonomy");
    $otumap{$A[0]} = $A[1];
  }

  # create normalized dataset ---------------------
  my %norm = ();
  my %otus = ();
  foreach my $s (sort keys %data){
    foreach my $t (sort keys %{$data{$s}}){
      if (!defined($data{$s}{$t})){
        $data{$s}{$t} = 0;
      }
      $otus{$t}     = 1;
      $norm{$s}{$t} = 100*$data{$s}{$t}/$totals{$s};
    }
  }
  # create normalized dataset at otu level with taxonomy ---------------------
  open OUT, ">$outDir/$v/$v.otu-table.txt" or die "Error: Cannot write to $outDir/$v/$v.otu-table.txt\n";
  print OUT "Region\tSampleID";
  foreach my $t (sort keys %otus){
    print OUT "\t$t";
  }
  print OUT "\n";
  foreach my $s (sort keys %data){
    print OUT "$v\t$s";
    foreach my $t (sort keys %otus){
      print OUT "\t$norm{$s}{$t}";
    }
    print OUT "\n";
  }
  close OUT;

  # now aggregate at phylum class order family genus species levels
  my @levels = qw/kingdom phylum class order family genus species/;
  my @ls     = qw/k p c o f g s/;
  for my $i (1 .. 6){
    print "$v\t$levels[$i]\n";
    my %agg    = ();
    my %txstrs = ();
    # summarize at this level
    foreach my $t (sort keys %otus){
      my @t     = split /\; /, $otumap{$t};
      my $txstr = "Unassigned";
      if (defined($t[$i])){
        $txstr  = join("; ", @t[0..$i]);
      }else{
        $txstr  = $otumap{$t};
        for my $j (($#t+1)..$i){
          $txstr .= "; $ls[$j]\__unassigned";
        }
      }
      $txstrs{$txstr} = 1;
      foreach my $s (sort keys %norm){
        $agg{$s}{$txstr} += $norm{$s}{$t};
      }
    } # end of aggregation
    # zero out ---
    foreach my $s (sort keys %agg){
      foreach my $t (sort keys %txstrs){
        if (!defined($agg{$s}{$t})){
          $agg{$s}{$t} = 0;
        }
      }
    }

    open OUT, ">$outDir/$v/$v.$levels[$i].txt" or die "Error: Cannot write to $outDir/$v/$v.$levels[$i].txt\n";
    print OUT "Region\tSampleID";
    foreach my $t (sort keys %txstrs){
      print OUT "\t$t";
    }
    print OUT "\n";
    foreach my $s (sort keys %agg){
      print OUT "$v\t$s";
      foreach my $t (sort keys %txstrs){
        print OUT "\t$agg{$s}{$t}";
      }
      print OUT "\n";
    }
    close OUT;
  } # end of levels loop
} # end of regions loop

# --------
# given all regions, aggregate into merged results
`mkdir $outDir/final`;
my @flevels = qw/phylum class order family genus species otu-table/;
foreach my $fl (@flevels){
  my %fldata = ();
  my %tx     = ();
  foreach my $v (@regions){
    my @header = ();
    open IN, "$outDir/$v/$v.$fl.txt" or die "$outDir/$v/$v.$fl.txt";
    while(<IN>){
      chomp($_);
      my @A = split "\t", $_;
      if (!defined($header[1])){
        @header = @A;
      }else{
        for my $i (2 .. $#A){
          $fldata{$v}{$A[1]}{$header[$i]} = $A[$i];
          $tx{$header[$i]} = 1;
        }
      }
    }
    close IN;
  } # end of regions

  # zero out ---
  foreach my $v (sort keys %fldata){
    foreach my $s (sort keys %{$fldata{$v}}){
      foreach my $t (sort keys %tx){
        if (!defined($fldata{$v}{$s}{$t})){
          $fldata{$v}{$s}{$t} = 0;
        }
      }
    }
  }
  # print final ---
  open OUT, ">$outDir/final/$fl.txt" or die "Error: Cannot write to $outDir/final/$fl.txt\n";
  print OUT "Region\tSampleID";
  foreach my $t (sort keys %tx){
    print OUT "\t$t";
  }
  print OUT "\n";
  foreach my $v (sort keys %fldata){
    foreach my $s (sort keys %{$fldata{$v}}){
      print OUT "$v\t$s";
      foreach my $t (sort keys %tx){
        print OUT "\t$fldata{$v}{$s}{$t}";
      }
      print OUT "\n";
    }
  }
  close OUT;
} # end of levels
