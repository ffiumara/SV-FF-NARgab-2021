#! /usr/bin/perl

use strict;
use warnings;
use diagnostics;
use List::MoreUtils qw(uniq);

#usage:
#specify in the command line:
#1) Species name abbreviation (e.g. "HomSap" for Homo sapiens);
#2) repeat type (e.g. "polyQ");
#3) input file with protein IDs (BioGrid IDs, one per line; e.g. "116520");
#4) input file with ID mapping (tab-separated BioGrid and Uniprot IDs, one ID pair per line; e.g. "Q2M2I8 116520");
#5) input file with BioGrid interactome (in BioGrid ".tab2" format);

#the scritpt generates 5 output files:
#1) OUTPUTFILE1: results summary (column headers as in OUTPUTFILE4);
#2) OUTPUTFILE2: list of Biogrid physical protein-protein interactions in which both proteins contain the repeat;
#3) OUTPUTFILE3: list of Biogrid physical protein-protein interactions in which only one protein contains the repeat;
#4) OUTPUTFILE4: results summary with column headers;
#5) OUTPUTFILE5: sum of OUTPUTFILE2 and OUTPUTFILE3;

#_______________________________________________________________________________________________________

open( INPUTFILEids,      "<$ARGV[2]" );
open( INPUTFILEidCONV,   "<$ARGV[3]" );
open( INPUTFILEinteract, "<$ARGV[4]" );

open( OUTPUTFILE1, ">>output-2-interactome-biogrid-TOT-$ARGV[0].txt" );
open( OUTPUTFILE2,
    ">output-2-interactome-biogrid-for-graph-Rep-Rep-$ARGV[0]-$ARGV[1].txt" );
open( OUTPUTFILE3,
    ">output-2-interactome-biogrid-for-graph-Rep-NO-Rep-$ARGV[0]-$ARGV[1].txt"
);
open( OUTPUTFILE4, ">output-2-interactome-biogrid-$ARGV[0]-$ARGV[1].txt" );
open( OUTPUTFILE5,
    ">output-2-interactome-biogrid-for-graph-tot-$ARGV[0]-$ARGV[1].txt" );

#_______________________________________________________________________________________________________

my @elementilinea                  = ();
my @IDgenenames                    = ();
my @IDs                            = ();
my @interactionsRepNORep           = ();
my @interactionsRepRep             = ();
my @LINESidCONV                    = ();
my @LINESids                       = ();
my @LINESinteract                  = ();
my @NOrepNODES                     = ();
my @repNODES                       = ();
my @repNODESincludingNOreprepNODES = ();
my @TOTAL                          = ();

my $condizione           = "";
my $IDgenename           = "";
my $IDsingolo            = "";
my $IDsource             = "";
my $IDtarget             = "";
my $indicator            = "";
my $interaction          = "";
my $interaction2         = "";
my $linea                = "";
my $pernodeinteractionsA = 0;
my $pernodeinteractionsB = 0;
my $RATIO1               = 0;
my $variabile1           = 0;
my $variabile2           = 0;
my $variabile3           = 0;

#_____________________________________________________________________________________________________

@LINESids      = <INPUTFILEids>;
@LINESidCONV   = <INPUTFILEidCONV>;
@LINESinteract = <INPUTFILEinteract>;

foreach $linea (@LINESids) {
    @elementilinea = split( '\|', $linea );
    $IDsingolo = $elementilinea[1];
    push( @IDs, $elementilinea[1] );
}

foreach $linea (@LINESidCONV) {
    @elementilinea = split( '\s', $linea );
    $IDgenename = $elementilinea[1];
    push( @IDgenenames, $elementilinea[1] );
}

$condizione = "physical";

foreach $linea (@LINESinteract) {
    @elementilinea = split( '\t', $linea );
    $IDsource = $elementilinea[3];
    chomp $IDsource;
    $IDtarget = $elementilinea[4];
    chomp $IDtarget;
    chomp $elementilinea[12];
    if ( $elementilinea[12] eq $condizione ) {
        if ( grep { $_ eq $IDsource } @IDgenenames ) {
            unless ( grep { $_ eq $IDsource } @repNODESincludingNOreprepNODES )
            {
                push @repNODESincludingNOreprepNODES, $IDsource;
            }
        }
        if ( grep { $_ eq $IDtarget } @IDgenenames ) {
            unless ( grep { $_ eq $IDtarget } @repNODESincludingNOreprepNODES )
            {
                push @repNODESincludingNOreprepNODES, $IDtarget;
            }
        }
    }
}

#_______________________________________________________________________________________________________

$indicator = "BIOGRID INTERACTION:";

foreach $linea (@LINESinteract) {
    @elementilinea = split( '\t', $linea );
    $IDsource = $elementilinea[3];
    chomp $IDsource;
    $IDtarget = $elementilinea[4];
    print ".", "\n";
    chomp $elementilinea[12];

    $condizione = "physical";

    if ( $elementilinea[12] eq $condizione ) {

        if ( grep { $_ eq $IDsource } @IDgenenames ) {
            if ( grep { $_ eq $IDtarget } @IDgenenames ) {
                $interaction  = join "\t", $indicator, $IDsource, $IDtarget;
                $interaction2 = join "\t", $indicator, $IDtarget, $IDsource;
                unless ( grep { $_ eq $interaction2 } @interactionsRepRep ) {
                    unless ( grep { $_ eq $interaction } @interactionsRepRep ) {
                        push @interactionsRepRep, $interaction;
                        print OUTPUTFILE2 $interaction, "\t",
                          $elementilinea[12], "\t", "rep-rep", "\n";
                        print $interaction, "\t", $elementilinea[12], "\t",
                          "rep-rep", "\n";
                        print OUTPUTFILE5 $interaction, "\t",
                          $elementilinea[12], "\t", "rep-rep", "\n";
                        push @TOTAL, $IDsource;
                        push @TOTAL, $IDtarget;
                    }
                }
            }
        }
        @repNODES = uniq @TOTAL;

        if ( grep { $_ eq $IDsource } @IDgenenames ) {
            unless ( grep { $_ eq $IDtarget } @IDgenenames ) {
                $interaction  = join "\t", $indicator, $IDsource, $IDtarget;
                $interaction2 = join "\t", $indicator, $IDtarget, $IDsource;
                unless ( grep { $_ eq $interaction } @interactionsRepNORep ) {
                    unless ( grep { $_ eq $interaction2 }
                        @interactionsRepNORep )
                    {
                        push @interactionsRepNORep, $interaction;
                        print OUTPUTFILE3 $interaction, "\t",
                          $elementilinea[12], "\t", "rep-NOrep", "\n";
                        print $interaction, "\t", $elementilinea[12], "\t",
                          "rep-NOrep", "\n";
                        print OUTPUTFILE5 $interaction, "\t",
                          $elementilinea[12], "\t", "rep-NOrep", "\n";
                        unless ( grep { $_ eq $IDtarget } @NOrepNODES ) {
                            push @NOrepNODES, $IDtarget;
                        }
                    }
                }
            }
        }

        if ( grep { $_ eq $IDtarget } @IDgenenames ) {
            unless ( grep { $_ eq $IDsource } @IDgenenames ) {
                $interaction  = join "\t", $indicator, $IDsource, $IDtarget;
                $interaction2 = join "\t", $indicator, $IDtarget, $IDsource;
                unless ( grep { $_ eq $interaction } @interactionsRepNORep ) {
                    unless ( grep { $_ eq $interaction2 }
                        @interactionsRepNORep )
                    {
                        push @interactionsRepNORep, $interaction;
                        print OUTPUTFILE3 $interaction, "\t",
                          $elementilinea[12], "\t", "NOrep-rep", "\n";
                        print $interaction, "\t", $elementilinea[12], "\t",
                          "NOrep-rep", "\n";
                        print OUTPUTFILE5 $interaction, "\t",
                          $elementilinea[12], "\t", "NOrep-rep", "\n";
                        unless ( grep { $_ eq $IDsource } @NOrepNODES ) {
                            push @NOrepNODES, $IDsource;
                        }
                    }
                }
            }
        }

    }
}

$variabile1 = scalar @repNODESincludingNOreprepNODES;
$variabile2 = scalar @repNODES;
$variabile3 = scalar @IDs;

unless ( $variabile1 == 0 ) {
    $pernodeinteractionsA =
      scalar @interactionsRepRep / scalar @repNODESincludingNOreprepNODES;
}
else { $pernodeinteractionsA = 0 }
unless ( $variabile2 == 0 ) {
    $pernodeinteractionsB = scalar @interactionsRepRep / scalar @repNODES;
}
else { $pernodeinteractionsB = 0 }
unless ( $variabile3 == 0 ) { $RATIO1 = scalar @IDgenenames / scalar @IDs; }
else                        { $RATIO1 = 0 }

#______________________________________________________________________________________________________

print $ARGV[0], "\t", $ARGV[1], "\t", scalar @IDs, "\t", scalar @IDgenenames,
  "\t", $RATIO1, "\t", scalar @repNODESincludingNOreprepNODES, "\t",
  scalar @repNODES, "\t", scalar @NOrepNODES, "\t", scalar @interactionsRepRep,
  "\t", scalar @interactionsRepNORep, "\t", $pernodeinteractionsA, "\t",
  $pernodeinteractionsB, "\n";

print OUTPUTFILE1 $ARGV[0], "\t", $ARGV[1], "\t", scalar @IDs, "\t",
  scalar @IDgenenames, "\t", $RATIO1, "\t",
  scalar @repNODESincludingNOreprepNODES, "\t", scalar @repNODES, "\t",
  scalar @NOrepNODES, "\t", scalar @interactionsRepRep, "\t",
  scalar @interactionsRepNORep, "\t", $pernodeinteractionsA, "\t",
  $pernodeinteractionsB, "\n";

print OUTPUTFILE4 "species", "\t",
  "repeat",                                                       "\t",
  "number of proteins with repeat in proteome",                   "\t",
  "number of proteins with repeat in proteome mapped to Biogrid", "\t",
  "mapping efficiency",                                           "\t",
"number of proteins with the repeat establishing interactions in the Biogrid interactome",
  "\t",
"number of proteins with the repeat establishing interactions with at least one other protein containing the repeat in the Biogrid interactome",
  "\t",
"number of proteins without the repeat establishing interactions with proteins containing the repeat in the Biogrid interactome",
  "\t",
"number of interactions of proteins with the repeat with proteins with the repeat (rep-rep interactions)",
  "\t",
"number of interactions of proteins with the repeat with proteins without the repeat",
  "\t",
  "per-node rep-rep interactions of proteins containing the repeat", "\t",
"per-node rep-rep interactions of proteins containing the repeat that interact at least with one other protein containing the repeat",
  "\n";

print $ARGV[0], "\t", $ARGV[1], "\t", scalar @IDs, "\t", scalar @IDgenenames,
  "\t", $RATIO1, "\t", scalar @repNODESincludingNOreprepNODES, "\t",
  scalar @repNODES, "\t", scalar @NOrepNODES, "\t", scalar @interactionsRepRep,
  "\t", scalar @interactionsRepNORep, "\t", $pernodeinteractionsA, "\t",
  $pernodeinteractionsB, "\n";

print OUTPUTFILE4 $ARGV[0], "\t", $ARGV[1], "\t", scalar @IDs, "\t",
  scalar @IDgenenames, "\t", $RATIO1, "\t",
  scalar @repNODESincludingNOreprepNODES, "\t", scalar @repNODES, "\t",
  scalar @NOrepNODES, "\t", scalar @interactionsRepRep, "\t",
  scalar @interactionsRepNORep, "\t", $pernodeinteractionsA, "\t",
  $pernodeinteractionsB, "\n";

#_____________________________________________________________________________________________________

close INPUTFILEids;
close INPUTFILEidCONV;
close INPUTFILEinteract;

close OUTPUTFILE1;
close OUTPUTFILE4;
close OUTPUTFILE2;
close OUTPUTFILE3;
close OUTPUTFILE5;

#_____________________________________________________________________________________________________
