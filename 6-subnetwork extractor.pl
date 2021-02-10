#! /usr/bin/perl

use strict;
use warnings;
use diagnostics;

#usage:
#specify in the command line:
#1) tab-separated input file whose lines specify protein-protein interactions in this format: "BIOGRID INTERACTION:	109314	118845	physical	rep-rep	HTT	SETD2", where "BIOGRID INTERACTION:" is a fixed string, "109314" is the BioGrid ID of one protein, "118845" is the BioGrid ID of the other protein, "physical" is the type of interaction between the two proteins, and "HTT" and "SETD2" are the gene names related to the two proteins, respectively;
#2) comma-separated input file with a set of IDs of interest;

#the script generates an OUTPUT file containing those lines, i.e. protein-protein interactions, of the first input file in which both proteins are encoded by genes whose IDs are both contained in the input file #2; this output file can be used to generate and analyze quantitatively protein-protein interaction networks using Cytoscape;

#_______________________________________________________________________________________

open( INPUTnet, "<$ARGV[0]" );
open( INPUTids, "<$ARGV[1]" );
open( OUTPUT,   ">output-0-SUBNETWORK-$ARGV[1].txt" );

#_______________________________________________________________________________________

my @LINESnet = ();
my @LINESids = ();
my @elements = ();
my @temp     = ();
my $temp     = '';
my $lineNET  = '';
my $elements = '';

#_______________________________________________________________________________________

@LINESnet = <INPUTnet>;
@temp     = <INPUTids>;
chomp @temp;
@LINESids = split ",", $temp[0];

foreach $lineNET (@LINESnet) {

    chomp $lineNET;
    @elements = split '\t', $lineNET;
    chomp $elements[5];
    chomp $elements[6];

    if ( grep { $_ eq $elements[5] } @LINESids ) {
        if ( grep { $_ eq $elements[6] } @LINESids ) {
            print $lineNET, "\n";
            print OUTPUT $lineNET, "\n";
        }
    }

}

#_______________________________________________________________________________________

close INPUTnet;
close INPUTids;
close OUTPUT;
