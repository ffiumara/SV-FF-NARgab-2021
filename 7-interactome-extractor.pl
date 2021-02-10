#! /usr/bin/perl

use strict;
use warnings;
use diagnostics;

#usage:
#specify in the command line:
#1) input file with protein-protein interaction network of interest in the following tab-separated format: "BIOGRID INTERACTION:	109314	118845	physical	#rep-rep	HTT	SETD2", where "BIOGRID INTERACTION:" is a fixed line header, "109314" and "118845" are the Biogrid IDs of the two interacting #proteins, "physical" is the BioGrid interaction type, "rep-rep" is an indicator that the interaction is between two proteins bearing a repeat of interest, #and "HTT" and "SETD2" are the gene names associated with the two interacting proteins;
#2) the desired number of randomly selected "seed" nodes (e.g. 10);
#3) text descriptor of the dataset that will be reported in the output file names;

#the script generates 2 output files:
#1) output file (OUTPUTint) reporting 10 seed nodes found randomly by the script and their n direct interactors;
#2) output file (OUTPUTrnd) reporting  the 10 same seed nodes (see above) found randomly by the script and n (see above) randomly selected proteins;
#3) output file (OUTPUTbasicnodes) reporting the IDs of the 10 seed nodes common to the other two output files;

#_______________________________________________________________________________________

open( INPUT1,           "<$ARGV[0]" );
open( OUTPUTint,        ">OUTPUT-1-INTERACTING-$ARGV[2]-N-$ARGV[1]-$ARGV[0]" );
open( OUTPUTrnd,        ">OUTPUT-1-RANDOM-$ARGV[2]-N-$ARGV[1]-$ARGV[0]" );
open( OUTPUTbasicnodes, ">OUTPUT-1-COMMON-NODES-$ARGV[2]-N-$ARGV[1]-$ARGV[0]" );

#_______________________________________________________________________________________

my @unique       = ();
my @unique2      = ();
my @RNDunique    = ();
my @noINT        = ();
my @LINESnetwork = ();
my @INT          = ();
my @INTunique    = ();
my @elements     = ();

my $temp     = '';
my $proteina = '';
my $line     = '';
my $node     = '';

my $unique  = 0;
my $Ntarget = 0;
my $N       = 0;
my $counter = 0;
my $casuale = 0;

#_______________________________________________________________________________________

@LINESnetwork = <INPUT1>;
$N            = $ARGV[1];
chomp $N;

foreach $line (@LINESnetwork) {

    @elements = split '\t', $line;
    chomp $elements[5];
    chomp $elements[6];
    print $elements[5], "\t", $elements[6], "\n";
    unless ( grep { $_ eq $elements[5] } @unique ) {
        push @unique, $elements[5];
    }
    unless ( grep { $_ eq $elements[6] } @unique ) {
        push @unique, $elements[6];
    }

}

print scalar @unique, "\n";

$counter = 0;

until ( $counter == $ARGV[1] ) {
    $casuale = rand( scalar @unique - 1 );
    $temp    = $unique[$casuale];
    unless ( grep { $_ eq $temp } @unique2 ) {
        push @unique2, $temp;
        print OUTPUTrnd $temp,        ",";
        print OUTPUTbasicnodes $temp, ",";
        $counter = $counter + 1;
    }
}

print scalar @unique2, "\n";

#_______________________________________________________________________________________

foreach $line (@LINESnetwork) {

    @elements = split '\t', $line;
    chomp $elements[5];
    chomp $elements[6];
    if ( grep { $_ eq $elements[5] } @unique2 ) {
        push @INT, $elements[5];
        push @INT, $elements[6];
    }
    if ( grep { $_ eq $elements[6] } @unique2 ) {
        push @INT, $elements[5];
        push @INT, $elements[6];
    }

}

foreach $node (@INT) {

    chomp $node;
    unless ( grep { $_ eq $node } @INTunique ) {
        push @INTunique, $node;
        print OUTPUTint $node, ",";
    }

}

print scalar @INTunique, "\n";

#_______________________________________________________________________________________

foreach $proteina (@unique) {

    unless ( grep { $_ eq $proteina } @INTunique ) {
        push @noINT, $proteina;
    }

}

print scalar @noINT, "\n";
@RNDunique = @unique2;
$Ntarget   = scalar @INTunique;
$counter   = scalar @RNDunique;

until ( $counter == $Ntarget ) {
    $casuale = rand( scalar @noINT - 1 );
    $temp    = $noINT[$casuale];
    unless ( grep { $_ eq $temp } @RNDunique ) {
        push @RNDunique, $temp;
        print OUTPUTrnd $temp, ",";
        $counter = $counter + 1;
    }
}

print scalar @RNDunique, "\n";

#_______________________________________________________________________________________

close INPUT1;
close OUTPUTint;
close OUTPUTrnd;
close OUTPUTbasicnodes;

