#! /usr/bin/perl

use strict;
use warnings;
use diagnostics;

#usage:
#specify in the command line:
#1) tab-separated input file with IDs of proteins belonging to one or more functional clusters in the following format:
#  "cluster	terms-in-cluster	genes-in-cluster	genes-in-cluster+rep	genes-in-cluster-rep"
# where "cluster" is the cluster name,
# "terms-in-cluster" is the list of GO/HPO/GWAS terms/disease associations related to the cluster,
# "genes-in-cluster" is the comma-separated list of genes/protein associated with the cluster,
# "genes-in-cluster+rep" is the comma-separated list of genes/protein associated with the cluster that contain the repeat,
# "genes-in-cluster-rep" is the comma-separated list of genes/protein associated with the cluster that do not contain the repeat

#the scritpt generates n OUTPUT files, one for each cluster listed in the input file, containing the comma-separated list of gene/protein IDs related to the cluster that contain the repeat.

#______________________________________________________________________________________________

open( INPUT, "<$ARGV[0]" );

#______________________________________________________________________________________________

my @elements = ();
my $line     = '';
my $elements = '';

#______________________________________________________________________________________________

my @LINES = <INPUT>;
shift @LINES;

foreach $line (@LINES) {

    chomp $line;
    @elements = split '\t', $line;
    open( OUTPUT, ">output-0-ID-list-$ARGV[0]-$elements[0].txt" );
    print $elements[0], "\t", $elements[3], "\n";
    print OUTPUT $elements[3];
    close OUTPUT;

}

#______________________________________________________________________________________________

close INPUT;

