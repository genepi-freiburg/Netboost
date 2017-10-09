#!/usr/bin/perl
########## START AUTOMATIC CODE ##########################################
use strict;
use Getopt::Long;

my ($help,$log_file,$out_file,$verbose);
(my $progname = $0) =~ s'^.*/'';
my $USAGE = "$progname [-help] [-log_file <FILE>] [-out_file <FILE>] [-verbose|-noverbose] EDGE_FILE_1 .. EDGE_FILE_N
Checks that input edges are in the correct format, and obey input assumption, i.e. 
	1. No self edges
	2. No duplicate edges 
	3. Line formatting: 'cluster_id1 cluster_id2 dist'

This script will naively hold all edges in memory to find duplicate edges, and is useful for debugging on 
buggy small examples, rather than huge data that does not fit in memory. Gzipped files are allowed as input 
arguments".

GetOptions(
           'help'       => \$help,
           'log_file=s' => \$log_file,
           'out_file=s' => \$out_file,
           'verbose!'   => \$verbose
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

############### END AUTOAMTIC CODE ####################################
die("no filenames given") unless scalar(@ARGV);
my %E;
while(my $file = shift) {
	warn("Edges filenames must have the filename extension: .edges.gz") unless ($file =~ /\.edges\.gz$/);

	if ($file =~ /.gz$/){
		print LOG "Reading gzipped file $file using zcat pipe\n";
		open(A,"zcat $file |") or die("error zcat-ing $file:$!");
	} else {
		print LOG "Reading plain file $file\n";
		open(A,"< $file") or die("error openning $file:$!");
	}
	while(<A>) {
		chomp;
		my $line_desc = "$file line $.:'$_'";
		die("not in the right format:$_") unless (/^(\d+)\s+(\d+)\s+([\d\.eE\-\+]+)$/);
		my $i = int($1); # gets rid of leading 0's
		my $j = int($2);
		my $dist = $3;
		warn("illegal zero index in $line_desc") if ($i == 0 || $j ==0);
		warn("self edge in $line_desc") if ($i == $j);
		if ($j < $i) { # swap
			my $temp = $i;
			$i = $j;
			$j = $temp;
		}
		if (exists($E{$i}{$j})) { 
			warn ("Input edge in $line_desc duplicated with $E{$i}{$j}");
			$E{$i}{$j} .= ", $line_desc";
		} else {
			$E{$i}{$j} = $line_desc;
		}
	}
	close(A) or warn("error closing $file");
}
	


