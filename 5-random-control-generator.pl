#! /usr/bin/perl

use strict;
use warnings;
use diagnostics;

#usage:
#specify in the command line:
#1) comma-separated input file with a set of IDs;
#2) comma-separated input file with a set of IDs to select from;

#the script generates an OUTPUT file containing a comma-separated list of IDs randomly selected from input file #2, equinumerous with the list of IDs in input file #1, for control analyses;

#_________________________________________________________________________________________________________

open( INPUT1, "<$ARGV[0]" );
open( INPUT2, "<$ARGV[1]" );
open( OUTPUT, ">RANDOM-CONTROL-$ARGV[0]" );

#_________________________________________________________________________________________________________

my @totalIDs      = ();
my @UNIQUE        = ();
my @TEMP12        = ();
my @LINES         = ();
my @GENERALidLIST = ();
my $RANDOMid      = '';
my $LINES         = '';
my $GENERALidLIST = '';
my $temp          = 0;
my $targetnumber  = 0;
my $randomX       = 0;
my $counter       = 0;

#_________________________________________________________________________________________________________

@LINES         = <INPUT1>;
@TEMP12        = <INPUT2>;
@GENERALidLIST = split ",", $TEMP12[0];
$temp          = scalar @GENERALidLIST;
@totalIDs      = split ",", $LINES[0];
$targetnumber  = scalar @totalIDs;

until ( $counter == $targetnumber ) {

    $randomX  = int( rand($temp) );
    $RANDOMid = $GENERALidLIST[$randomX];
    chomp $RANDOMid;

    unless ( grep { $_ eq $RANDOMid } @UNIQUE ) {
        push @UNIQUE, $RANDOMid;
        $counter = $counter + 1;
        print $counter, "\t", $RANDOMid, "\n";
        print OUTPUT $RANDOMid . ",";
    }

}

#_________________________________________________________________________________________________________

close INPUT1;
close INPUT2;
close OUTPUT;

#_________________________________________________________________________________________________________
