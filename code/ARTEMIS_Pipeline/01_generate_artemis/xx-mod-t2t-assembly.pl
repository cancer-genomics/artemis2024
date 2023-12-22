#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;


open OUT, ">chrm13.mod.fna" or die;
open IN, "gunzip -c GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz |"  or die;
while(<IN>){
  chomp($_);
  if ($_ =~ /^>/){
   my @A = split " ", $_;
   $A[6] =~ s/\,//g;
   print OUT ">chr$A[6]\n"; 
  }else{
    my $line = uc($_);
    print OUT "$line\n";
  }
}
close IN;
close OUT;
