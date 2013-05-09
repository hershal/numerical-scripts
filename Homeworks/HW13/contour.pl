#!/usr/bin/perl -w

use warnings;
use strict;

if (!defined($ARGV[0]) || !defined($ARGV[1])) {
  die("usage: <format> <data files...>");
} 

my $format = shift(@ARGV);
  
foreach my $file (@ARGV) {
  open(my $in, "<$file") or die("Could not open $file for reading.");
  my $basename = $file;
  $basename =~ s/(^.*)\..*?$/$1/;
  $basename =~ s/_/-/g;

  my $outFile = $basename.".tmp";
  my $gnuOutputFile = $basename.".".$format;
    
  open(my $out, ">$outFile") or die("Could not open $outFile for writing");

  my $old1;
  while (<$in>) {
    if (defined($old1) && $1 ne $old1) {
      print($out "\n")
    }
    if ($_ =~ m/(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)/g) {
      print($out "$1 $2 $3\n");
    }
    $old1 = $1;
  }

  my $plotcmd = <<PERLEOF;
#!/bin/sh
gnuplot << GNUEOF
set style line 11 lc rgb 'black' lt 1
set style line 12 lc rgb '#808080' lt 0 lw 1
set tics nomirror

set nokey

set terminal $format
set output "$gnuOutputFile"

set style line 1 lc rgb 'black' lt 1 lw 2

set pm3d map

# Matlab-style colorbars
# Source: http://www.gnuplotting.org/matlab-colorbar-with-gnuplot/
set palette defined ( 0 "#000090", 1 "#000fff", 2 "#0090ff",\\
                      3 "#0fffee", 4 "#90ff70", 5 "#ffee00",\\
                      6 "#ff7000", 7 "#ee0000", 8 "#7f0000")

splot '$outFile'

GNUEOF

PERLEOF
  
  system("$plotcmd");
}

printf("done\n");


