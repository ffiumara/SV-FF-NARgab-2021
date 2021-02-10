#! /usr/bin/perl

use strict;
use warnings;
use diagnostics;

#usage:
#specify in the command line:
#1) input file with tab-separated matrix of r correlation coefficients;

#the script generates an OUTPUT file containing a matrix of z coefficients derived from the r coefficient matrix in the input file, after Fisher Z-transformation;

#_______________________________________________________________________________________

open( INPUT,  "<$ARGV[0]" );
open( OUTPUT, ">OUTPUT-0-$ARGV[0]-Z-converted.txt" );

#_______________________________________________________________________________________

my @LINES        = ();
my @elements     = ();
my @counter      = ();
my $z            = 0;
my $totelements  = 0;
my $rassoluto    = 0;
my $r            = 0;
my $counter      = 0;
my $line         = '';
my $intestazione = '';

#_______________________________________________________________________________________

@LINES        = <INPUT>;
$intestazione = shift @LINES;
print OUTPUT $intestazione;

foreach $line (@LINES) {

    chomp $line;
    @elements    = split '\t', $line;
    $totelements = scalar @elements - 1;
    @counter     = ( 1 .. $totelements );
    print $elements[0], "\t";
    print OUTPUT $elements[0], "\t";

    foreach $counter (@counter) {

        $r         = $elements[$counter];
        $rassoluto = abs($r);

        unless ( $r eq "NA" ) {
            if ( $rassoluto != 1 ) {
                $z = 0.5 * ( log( 1 + $r ) - log( 1 - $r ) );
            }
            else { $z = "NA" }
        }
        else { $z = "NA" }

        print $z, "\t";
        print OUTPUT $z, "\t";

    }

    print "\n";
    print OUTPUT "\n";

}

#_______________________________________________________________________________________

close INPUT;
close OUTPUT;
