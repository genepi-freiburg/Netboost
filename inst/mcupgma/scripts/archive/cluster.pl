#!/usr/bin/perl

##########################################################################
# MC-UPGMA  - Accurate huge scale clustering by Memory Constrained UPGMA #
#             Loewenstein et al. Bioinformatics. 2008 Jul 1;24(13):i41-9.#
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

my ($help,$log_file,$out_file,$verbose);
(my $progname = $0) =~ s'^.*/'';
(my $this_script_path = $0) =~ s/[^\/]+$//;
$this_script_path = './' unless $this_script_path;

# the name the make executable (e.g. use gmake for FreeBSD, not make)
my $MAKE = (`uname -s` =~ /BSD/) ? 'gmake' : 'make'; 

my $max_singleton_id ;              #= 600;
my $jumpstart;                      # start from this iteration if var defined
my $MAX_N_RETRIES = 0;              # how many times to try running make, when exiting with error (non 0)
my $MAX_CLUSTERING_ITERATIONS = 20; # 
my $NUM_EDGES2LOAD = 5000000;       # default number of edges to load
my $jobs;                           # number of concurent jobs passed to make process
my $isInputSorted;
my $H = 40;                         # number of files to split the edges for external merging to
                                    # higher numbers, reduce the input size per process and allow more
                                    # concurent merging processes
my $sleep_seconds_per_round = 0;
my $num_split_unmodified_edges; # see help for the -k option
my $allow_non_dendrogram;

my $MAX_DIST;    # the dissimilarity cutoff (names \psi in the paper).
my $mosix_flag;  # flags to pass to mosrun when running via the MOSIX grid.
my $inc; # a dir with all the scripts used by this script
my $is_print_filesizes = 1;
my $output_tree_file = "cluster_tree";

my $USAGE = "$progname [-help] [-log_file <FILE>] [-verbose|-noverbose]   *.edges.gz
Runs MC-UPGMA multi-round clustering makefiles.

The clustering is done in iterations. In each iteration the minimal M edges are loaded,
and clustering is done to the point where correctness (i.e. as in O(|E|) memory regular UPGMA) 
can no longer be guaranteed.

Input edges filenames must be have the extension .edges.gz and be gzipped.

Self edges, and duplicate edges (the edges i,j and j,i are also regarded duplicates) are not allowed as input, and 
will lead to undefined behavior (crashes) . It is the user's responsibility to eliminate these cases.

The multi-round clustering process, is based on the makefile mechanism which evaluates dependencies
between input and output files, and on their modification time-stamps. Therefor, modifying file 
time-stamps externally (by the user, not by this program) may invalidate the process, and is not 
recommended.

Use dense cluster IDs to improve performance.

Configurable Parameters:


   -max_distance|x|psi <FLOAT> - Upper bound on all distance in input, used as missing edge's distance. Specification of 
                                this option is mandatory.

   -max_singleton <UINT> - maximum cluster ID in input (>0). Can be any number N s.t. N > n for any input singleton ID n. 
                          N+1 is the first merger ID. Dense 

   -iterations <UINT>     - maximum number of clustering iterations.
 
   -M|heap_size  <UINT>   - maximum number of edges to load in a clustering iteration (dictated by hardware limits)
                            the larger M, the less clustering iterations needed.

                         ** NOTE: **  order of the input files specification is important when this option is used.

   -H|num_hash_buckets <UINT>  - how many output files of \"thin\" edges per input
   
   -jobs <UINT>          - run make with -j for parallelization (default: none). The value of this parameter is passed to 
                          make's '-j' option. This option allows a considerable speed-up by parallelization of the 
			  'external merging' procedures on a multi-processor machine.

   -retries <UINT>       - if make returns with an an error code, retry this number of times (default: $MAX_N_RETRIES)

   -otree|output_tree_file <FILE>         - output tree filename (default: $output_tree_file)

   -split_unmodified_edges|k <UINT>  - When this option is set, the file containing the unmodified cluster edges
				     passed from round to round is split into k parts, each part can then be 
				     be processed in parallel to speed-up the process on a multi-processor machine.
				     This is useful when this file contains the vast majority of the data and 
				     becomes the bottleneck, e.g. when very little clustering is done per round. An 
				     optimal value for k would be the number of concurrent jobs passed to the -jobs
				     option.
Advanced options:

   -allow_non_dendrogram -  Allows merging of provably minimal edge intervals, even if the exact merge score
                            (cluster height in dendrogram) is not known at merge time due to partial knowledge
                            implied by the memory constraint. This option allows further clustering per round, 
			    thus speeding up the whole clustering process considerably, by posing a less strict
                            requirement on the output - now the cluster heights are no longer required to be exact.
  			    Note that this renders the cluster heights (merge scores) in the output tree 
  			    meaningless. It does preserve the correct tree topology of the clustering solution 
			    however, and provides a useful speedup if don't need the dendrogram.

   -[no]sorted          - Assume input is sorted, uses UNIX's `head -M` to prevent reading all files by clustering 
                          code, currently only applies to first iteration, since output of merging is not resorted.
  
   -mosix               - Run external parallel merging code under mosix (currently, only the native processes are 
                          mosrun'ed. not the zcat). This capability is still half-baked, since most IO is done 
			  through pipes, which are inefficeinet using mosix. We need to use mospipe in the *.mk 
			  files to make the whole multi-round process efficient under mosix.

   -scriptsdir <DIR> - Use this directory for auxiliary scripts and makefiles, instead of using the same directory as this script's 
                       location (i.e. $this_script_path)

   -[no]filesizes    - Print the filesizes (in bytes) of input edge files per round (default: toggled) . This may slow 
                       down the round post-processing. It is useful for monitoring the rate in which the cluster similarity 
                       graph shrinks, or for identifying situations in which little or no clustering is done per round (e.g. when debugging).

   -start|t|T <INT>     - jump start from iteration T (don't run make for iterations < T).  Useful when debugging.

   -sleep <UINT>     - sleep this number of seconds after executing a round of clustering. Useful for 
			maintaining distinguishable timestamps for tiny inputs where a round lasts less than a 
			second, but shell operations (e.g. touch or cp) provide only whole seconds timestamp 
			resolutions, which might falsify make's dependencies. Sleeping 1 second between rounds is 
			useful when testing on very small inputs (see definitions.mk for more info, or to set from 
			within the makefiles).
Dependencies:
  clustering.mk - clustering makefile (which uses other supplied makefiles as well).
  wc            - standard UNIX word count executable.
  smart_cat     - Perl script choosing zcat if filename ends with .gz, cat otherwise.
  seq | seq.pl	- Linux executable, and Perl implementation providing a similar interface for other systems like FreeBSD.
  make (gmake)  - the GNU makefile program, invoked by make (Linux, Mac-Darwin) or gmake (FreeBSD).

------------------------------------------------------------------------------------------------------------

If you find this code useful, please cite:

Loewenstein Y, Portugaly E, Fromer M, Linial M.
Efficient algorithms for accurate hierarchical clustering of huge datasets: tackling the entire protein space.
Bioinformatics. 2008 Jul 1;24(13):i41-9.

-------------------------------------------------------------------------------------------------------------

Revisions:

Aug 2008, by Yaniv Loewenstein:
     *  Assume clustering makefiles are at the same directory as this script  (i.e. $this_script_path), 
        unless otherwise specified by the -scriptsdir option, which overrides this assumption, 
        and allows running this script with a different set of makefiles (and other auxiliaries).
     * add the print file-sizes option
     * added -otree option (output tree filename)
     * help section revised
Sept 2008:
     * added -allow_non_dendrogram option
-------------------------------------------------------------------------------------------------------------

";



GetOptions(
           'help'         => \$help,
           'log_file=s'   => \$log_file,
           'out_file=s'   => \$out_file,
           'verbose!'     => \$verbose,
	   'jobs=i'                  => \$jobs,
	   'retries=i'               => \$MAX_N_RETRIES,
	   'max_singelton|max_singleton=i'         => \$max_singleton_id, # support older version typo...
	   'iterations=i'            => \$MAX_CLUSTERING_ITERATIONS,
	   'M|heap_size|size_heap=i' => \$NUM_EDGES2LOAD,
	   'sorted!'                 =>\$isInputSorted,
	   'H|num_hash_buckets=i'    =>\$H,
	   'max_distance|x|psi=f'        =>\$MAX_DIST,
	   'mosix!'                  =>\$mosix_flag,
	   't|T|start=i'             =>\$jumpstart,
           'scriptdir=s'             =>\$inc,
	   'filesizes!' =>\$is_print_filesizes,
	   'otree|output_tree_file|tree=s' =>\$output_tree_file,
	   'split_unmodified_edges|k=i' => \$num_split_unmodified_edges,
	   'allow_non_dendrogram!'  => \$allow_non_dendrogram,
	   'sleep=i' => \$sleep_seconds_per_round
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



### Add the scripts dir to the PATH env variable, and to the make include option (-I) ###########
if ($inc) {
  warn("$inc - doesn't look like a directory") unless (-d $inc);
} else {
  $inc = $this_script_path;
}
$ENV{PATH}.=":$inc";
#################################################################################################




die if ($jobs < 0);
die if ($MAX_N_RETRIES    < 0);
die if ($max_singleton_id < 0);
die if ($MAX_CLUSTERING_ITERATIONS <= 0);
die('H < 1 is not allowed') if ($H < 1); 
die('H > 1000 is not allowed [values are restricted to 1000 to prevent overload on filesystem and CPU]')  if ($H > 1000); 
die("max distance option is required (try `$progname -help`)") unless(defined($MAX_DIST));
die("jump start must be for an iteration >= 1")  if (defined($jumpstart) && $jumpstart <1);
die('please pass a value > 1 , and < 200 , for -k (-split_unmodified_edges) [values are restricted to 200 to prevent overload on filesystem and CPU]')  
  if (defined($num_split_unmodified_edges) && (($num_split_unmodified_edges <= 1) ||($num_split_unmodified_edges >= 200)));
warn("You have specified a non positive maximum distace:$MAX_DIST") if ($MAX_DIST <= 0);
$jobs = "-j $jobs" if ($jobs);

my $SEPERATOR =  '
===================================================================================
=====================  Round %4d  ================================================
===================================================================================
';

# the --load-average is passed to the makefile to prevent clogging of this node, see manual page for GNU make
$mosix_flag = "--load-average=8 'USE_GRID = 1'" if ($mosix_flag);

# loop variables:
my $t = 1;
my $n_retries = 0;
my @input_files = @ARGV;
my $com;

die('no inputs') unless (scalar(@input_files));


unless (defined($max_singleton_id)){
  print LOG "$progname:finding the maximum singleton ID by reading all edges - this might take a while (you can supply it externally)\n";
  open(A,"smart_cat @input_files | ") or die;
  while(<A>){
    split;
    $max_singleton_id = $_[0] if ($_[0] > $max_singleton_id);
    $max_singleton_id = $_[1] if ($_[1] > $max_singleton_id);
  }
  close(A) or die;
  $max_singleton_id = $_[1] if ($_[1] > $max_singleton_id);
  die unless( defined($max_singleton_id));
  print LOG "$progname: maximum singleton ID in input edges:$max_singleton_id\n";
}



#######################
# 
#######################
print LOG "output cluster tree file: $output_tree_file (removing file, and building it from scratch)\n";
if (!$jumpstart) {
	print LOG `rm -fv $output_tree_file`; 
} else {
	warn("the tree file does not exist, but we're in jumpstart mode!") unless (-e $output_tree_file);
        print LOG "not removing existing tree, since jumpstart is set\n";
}


#######################
# initial sizes file ##
#######################
#print LOG "$progname:building input clusters size file, assuming all 1..$max_singleton_id clusters are of size = 1\n";
my $sizes_file ;#= 'isizes';
#warn("overwriting existing sizes file $sizes_file") if (-s $sizes_file);
#my $isizes_com = "seq  --format='%.0f' 1 $max_singleton_id  | awk '{print \$1,1}' >  $sizes_file";
#print LOG "running:$isizes_com\n";
#print LOG `$isizes_com`;
#die('something went wrong when creating the singleton sizes files:'.$!) if ($?);






# main clustering loop (runs rounds of clustering) is aborted if 
#    (1) error and reached max number of retries
#    (2) reached maximum number of rounds (iterations)
#    (3) no edges remain after last round of merging (successful completion).

while(1) {
  printf LOG $SEPERATOR,$t;
  if ($is_print_filesizes) { 
    my $total_size = 0;
    $total_size += (-s $_) foreach (@input_files); 
    print LOG "total size in bytes of the ".scalar(@input_files)." input files is $total_size".
		(("@input_files" =~ /^(\s*\S+\.gz)+\s*$/) ? ' (gzipped)' : '')."\n";
  } 

  if ($n_retries > $MAX_N_RETRIES){
    print LOG "$progname:reached $n_retries retries - aborting\n";
    exit(1);
  }  

  if ($t >= $MAX_CLUSTERING_ITERATIONS){
    print LOG "reached maximum number of iterations = $MAX_CLUSTERING_ITERATIONS - aborting\n";
    last;
  }	

  ##########################################
  # Check if we've finished clustering:
  # if one cluster in $sizes_file
  #     if no edges
  #         end of clustering
  
  ## Count edges - abort if any edges found (do not read the whole huge file with `wc -l`)
  open(EDGES,"smart_cat @input_files | head -1 |") or die;
  my $is_edges_remain = 0;
  while(<EDGES>){
    ++$is_edges_remain; last;
  }
  close(EDGES) or die("can't close edge stram $is_edges_remain");
  my $nclusters;
 if ($sizes_file) { # skip the first
    chomp($nclusters = `wc -l $sizes_file `);
    print LOG "checking the number of alive clusters - `wc -l` says:'$nclusters'\n";
    die unless ($nclusters =~ /^\s*(\d+)\s+\S+\s*$/);
    $nclusters = $1;
  } else { $nclusters = "undefined"; } # were in the first round

  if (!$is_edges_remain){
    if ($nclusters == 1){
      print LOG "$progname: one cluster in $sizes_file and edges files seem empty - finished clustering successfully!\n";
    } else {
      print LOG "$progname: edges files seem empty (altough more than one live cluster remains) - finished clustering successfully. Note that the output is a forest, not a tree since the input contained more than one connected component.";
    } 
    last;
  } else {
    print LOG "$progname: one cluster in $sizes_file (`wc -l` said:'$nclusters') but edges file are not empty - continue clustering\n"; 
    warn("*** Warning *** - Suspiciously, edges remain but there's only one cluster in the cluster sizes file (this might be normal if there are some singleton clusters having edges, which are absent from the cluster sizes file). ") if ($nclusters == 1);
  }

  # the make command line spec:
  $com = "$MAKE -I $inc --warn-undefined-variables ".
    (($isInputSorted && $t == 1) ? "'filt = head -$NUM_EDGES2LOAD' " : ' ').
      " $jobs 'PSI = $MAX_DIST' $mosix_flag 'H = $H' 't = $t' 'M = $NUM_EDGES2LOAD' ".
	"'SCRIPTS_DIR = $inc' ".
	"'TREE = $output_tree_file' ".
	((defined($num_split_unmodified_edges)) ? "'K = $num_split_unmodified_edges' " : ' ').
	((defined($allow_non_dendrogram)) ? "'ALLOW_INEXACT_MERGES := 1' " : ' ').
	"'INPUT_EDGES_FILES =  @input_files' ".
	  (($t > 1 ) ? "'ISIZES_FILE = $sizes_file' " : "'FIRST_MERGE_ID = $max_singleton_id' ").
	  "-f $inc/clustering.mk  2>&1 ";

  if (!defined($jumpstart) || $t >= $jumpstart){
    # run make and print its output to the log
    #---#if ($t > 161){ # debug code
    print LOG "$progname:executing : `$com`\n";
    open(A,"$com |") or die("cant open `$com |`");
    print LOG "make:$_" while(<A>);
    close(A);
    sleep($sleep_seconds_per_round);
    print LOG "finished executing\n";
    #---#}
  }
  if ($?){ # make retuned with non 0 exit code (i.e. error).
    ++$n_retries;
    print LOG "$progname:make failure #$n_retries (executed:`$com`)\n";
  } else {
    print LOG "$progname:returned 0 - looks like success!\n";

    # A list of the up-to-date edges files, connecting only valid clusters, is given in iteration_$t/outputs_listing
    unless (-s "iteration_$t/outputs_listing"){ # the file is empty or non existing, meaning the process wasn't complete successfully 
      print LOG "$progname:alas! the iteration_$t/outputs_listing file is missing! let's retry\n";
      ++$n_retries;
      next; # jump to the beginning of this loop to try again without incrementing the iteration counter $t.
    }	
    @input_files = `cat iteration_$t/outputs_listing `;     die if $?;
    $sizes_file = "iteration_$t/cluster_sizes";
    chomp(@input_files);
    die if (scalar(@input_files) > 10); #this is the number of lines
    @input_files = split(/\s+/,"@input_files"); # now we have seperate file per array element again
    print LOG "$progname:inputs to next stage = @input_files\n";
    $t++;

    #    last;
  }	
}

