#!/usr/bin/perl
use strict;
use warnings;

my @files = `ls pali/`;
chomp @files; 

for(my $i=0; $i<scalar(@files); $i++){
  print("mkdir results/$files[$i]_f ; cd results/$files[$i]_f ; ../../ipa -x 1000 -f ../../pali/$files[$i]\n");
  #for(my $j=$i+1; $j<scalar(@files); $j++){
  #  print("mkdir results/$files[$i]_$files[$j] ; cd results/$files[$i]_$files[$j] ; ../../ipa -x 1000 ../../pali/$files[$i] ../../pali/$files[$j]\n");
  #}
}
