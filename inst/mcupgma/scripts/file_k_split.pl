#!/usr/bin/perl

##########################################################################
# MC-UPGMA  - Accurate huge scale clustering by Memory Constrained UPGMA #
#             Loewenstein et al. Bioinformatics. 2008 Jul 1;24(13):i41-9.#
#                                                                        #
#  usage: emulates mixed cat and zcat command by filename extension      #
#                                                                        #
#    exec's to zcat if all extensions are .gz, cat if all or not.        #
#    Otherwise (a mixture of .gz files and non gzipped), each input file #
#    is piped to stdout using cat or zcat - this is slower               #
#                                                                        #
# Copyright (C) 2007, Yaniv Loewenstein                                  #
#                School of Computer Science And Engineering              #
#                Hebrew University of Jerusalem                          #
#                                                                        #
#      All Rights Reserved                                               #
#                                                                        #
#      This source code is distributed under the terms of the            #
#      GNU General Public License. See the file LICENSE                  #
#      for details.                                                      #
#                                                                        #
##########################################################################

use strict;
use Getopt::Long;

use IO::File;


my ($help,$log_file,$out_file,$verbose);
(my $progname = $0) =~ s'^.*/'';
my $USAGE = "$progname [-help] [-log_file <FILE>] [-out_file <FILE>] [-verbose|-noverbose]
Split a file to k seperate files, based on line numbers modulo k.

-k <INT>                            - number of output files
-ofile_base|base <STRING>           - filename base name for output files
-ofile_extension|extension <STRING> - filename extension for output files
-[no]ogz - gzip output (and add the .gz filename extension to the specified one).
-[no]filenames - echoes the list of created output files to stdout

";

my $K;
my $fn_base = 'splitted_file';
my $fn_ext;
my $is_gz_output;
my $is_echo_filenames;

GetOptions(
           'help'       => \$help,
           'log_file=s' => \$log_file,
           'out_file=s' => \$out_file,
           'verbose!'   => \$verbose,
	   'K|k=i'        => \$K,
	   'ofile_base|base=s' => \$fn_base,
	   'ofile_extension|extension=s' => \$fn_ext,
	   'ogz!' => \$is_gz_output,
           "echo_filenames|filenames!" => \$is_echo_filenames
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
die("-k not specified") unless(defined($K));
die("k must be positive") unless($K > 0);


my @filehandles;
my @files_opened;

for (my $k = 1 ; $k <= $K ; $k++) {
      my $filename = "$fn_base.$k";
      $filename .= ".$fn_ext" if ($fn_ext);
      $filename .= ".gz" if ($is_gz_output);
      $filehandles[$k] = new IO::File (($is_gz_output)  ? "| gzip -f > $filename" : "> $filename");
      push(@files_opened,$filename) if ($is_echo_filenames);
      my $fh = $filehandles[$k];
      die("couldn't open $filename for write") unless (defined $filehandles[$k]);
}

my $fh;
while(<>) {
	$fh = $filehandles[($.  % $K) or $K]; ## $. is the input line number starting at 1
	die("can't write line:$. to file number :".(($.  % $K) or $K)) unless(defined($fh));
	print $fh $_;
}

# close files loop
for (my $k = 1 ; $k <= $K ; $k++) {
	$filehandles[$k]->close;
}

print "@files_opened\n" if ($is_echo_filenames);

