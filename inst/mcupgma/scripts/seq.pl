#!/usr/bin/perl

use strict;
use Getopt::Long;

my $seperator = "\n";
my $format;

my ($help,$log_file,$out_file,$verbose);
(my $progname = $0) =~ s'^.*/'';
my $USAGE = "$progname [-help] [-log_file <FILE>] [-out_file <FILE>] [-verbose|-noverbose]

A perl implementation imitating the unix command seq. 


------------------------ Here's the seq man page: ---------------------------------
NAME
       seq - print a sequence of numbers

SYNOPSIS
       seq [OPTION]... LAST
       seq [OPTION]... FIRST LAST
       seq [OPTION]... FIRST INCREMENT LAST

DESCRIPTION
       Print numbers from FIRST to LAST, in steps of INCREMENT.

       -f, --format=FORMAT
              use printf style floating-point FORMAT (default: %g)

       -s, --separator=STRING
              use STRING to separate numbers (default: \\n)

       -w, --equal-width
              equalize width by padding with leading zeroes

       --help display this help and exit

       --version
              output version information and exit

       If  FIRST  or INCREMENT is omitted, it defaults to 1.  That is, an omitted INCREMENT defaults to 1 even when LAST
       is smaller than FIRST.  FIRST, INCREMENT, and LAST are interpreted as floating point values.  INCREMENT  is  usu‐
       ally  positive  if  FIRST  is smaller than LAST, and INCREMENT is usually negative if FIRST is greater than LAST.
       When given, the FORMAT argument must contain exactly one of the printf-style, floating point output  formats  %e,
       %f, %g

";


GetOptions(
           'help'       => \$help,
           'log_file=s' => \$log_file,
           'out_file=s' => \$out_file,
           'verbose!'   => \$verbose,
	   'separator=s' => \$seperator,
	   'format=s'   => \$format
          );
die("$USAGE") if ($help);
if ($log_file){
  open(STDERR,"> $log_file") or die("could not open log file:$!");
  my $oldfh = select(STDERR); $| = 1; select($oldfh);
}
open(LOG,">&STDERR") or die;

if ($out_file) {
  open(STDOUT,"> $out_file") or die ("could not open out file:$!");
}


# ====================================================
my ($start,$inc,$end);
if      (scalar(@ARGV) == 3) {
        $start = shift;
        $inc   = shift;
        $end   = shift;
} elsif (scalar(@ARGV) == 2) {
	$start = shift;
	$end   = shift;
	$inc = 1;
} elsif (scalar(@ARGV) == 1) {
	$end = shift;
	$start = 1;
	$inc = 1;
} else { die("not 1,2 or 3 arguments, see -help");}

unless(defined($format)){
	if (int($inc) == $inc) { # integer steps ==> integer output
		$format = '%d';
	} else {                 # non integer steps ==> %g formatted output
		$format = '%g';
	}
}



for (my $cur = $start; $cur <= $end ; $cur += $inc) {
	print $seperator unless ($cur == $start);
	printf($format,$cur);
}
print "\n";
	
