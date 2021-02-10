#! /usr/bin/perl

use strict;
use warnings;
use diagnostics;

#usage:
#specify in the command line:
#1) input file with tab-separated matrix of z correlation coefficients;
#2) input file with comma-separated list of gene names related to polyQ proteins of interest;
#3) input file with comma-separated list of gene names related to a control group of randomly selected polyQ proteins, equinumerous to the list in input file 2;

#the script generates 8 output files:
#1) output file reporting all the pairwise combinations of gene/protein names in the input file 2 and the associated z scores, as derived from input file 1;
#2) output file reporting all the pairwise combinations of gene/protein names in the input file 3 and the associated z scores, as derived from input file 1;
#3) output file as output 1 without those pairs for which a z score was not available (NA) in the input file 1;
#4) output file as output 2 without those pairs for which a z score was not available (NA) in the input file 1;
#5) output file with a summary of the analyses reported in output files 1 and 3;
#6) output file with a summary of the analyses reported in output files 2 and 4;
#7) output file as output file 5; the file is not overwritten in case of repeated use of the script with different pairs of input files;
#8) output file as output file 6; the file is not overwritten in case of repeated use of the script with different pairs of input files;

#_______________________________________________________________________________________

chomp $ARGV[1];
chomp $ARGV[2];

open( INPUT1, "<$ARGV[0]" );
open( INPUT2, "<$ARGV[1]" );
open( INPUT3, "<$ARGV[2]" );

open( OUTPUTset1, ">OUTPUT-1a-Z-COEFF-$ARGV[1]" );
open( OUTPUTrnd1, ">OUTPUT-1a-Z-COEFF-$ARGV[2]" );
open( OUTPUTset2, ">OUTPUT-1b-CLEAN-Z-COEFF-$ARGV[1]" );
open( OUTPUTrnd2, ">OUTPUT-1b-CLEAN-Z-COEFF-$ARGV[2]" );
open( OUTPUTset3, ">OUTPUT-1c-SUMMARY-Z-COEFF-$ARGV[1]" );
open( OUTPUTrnd3, ">OUTPUT-1c-SUMMARY-CLEAN-Z-COEFF-$ARGV[2]" );
open( OUTPUTset4, ">>OUTPUT-1d-SUMMARY-Z-COEFF-TOTAL-SET.txt" );
open( OUTPUTrnd4, ">>OUTPUT-1d-SUMMARY-Z-COEFF-TOTAL-RND.txt" );

#_______________________________________________________________________________________

my @SETrnd             = ();
my @SETids             = ();
my @POSITIONStargetIDs = ();
my @POSITIONSsourceIDs = ();
my @positionsSETids    = ();
my @positionsRNDids    = ();
my @listtabIDs         = ();
my @LINEStab           = ();
my @LINESset           = ();
my @LINESrnd           = ();
my @intestazione       = ();
my @elements2          = ();
my @counterTABids      = ();

my $tabID = '';
my $setID = '';
my $rndID = '';

my $ZsetTOTALabs     = 0;
my $ZsetTOTAL        = 0;
my $ZrndTOTALabs     = 0;
my $ZrndTOTAL        = 0;
my $totaltableIDs    = 0;
my $targetposition   = 0;
my $sourceposition   = 0;
my $numberofSETids   = 0;
my $numberofRNDids   = 0;
my $Nsetna           = 0;
my $NsetCLEAN        = 0;
my $Nset             = 0;
my $Nrndna           = 0;
my $NrndCLEAN        = 0;
my $Nrnd             = 0;
my $meanZsetTOTALabs = 0;
my $meanZsetTOTAL    = 0;
my $meanZset         = 0;
my $meanZrndTOTALabs = 0;
my $meanZrndTOTAL    = 0;
my $meanZrnd         = 0;
my $meanabsZset      = 0;
my $meanabsZrnd      = 0;
my $counter          = 0;

#_______________________________________________________________________________________section0

@LINEStab = <INPUT1>;
@LINESset = <INPUT2>;
@LINESrnd = <INPUT3>;

@listtabIDs = split '\t', $LINEStab[0];
shift @listtabIDs;
@SETids = split ",", $LINESset[0];
@SETrnd = split ",", $LINESrnd[0];

print OUTPUTset1 $ARGV[1], "\n";
print OUTPUTset2 $ARGV[1], "\n";
print OUTPUTrnd1 $ARGV[2], "\n";
print OUTPUTrnd2 $ARGV[2], "\n";

foreach $tabID (@listtabIDs) { print "gene name: ", $tabID, "\n"; }

print "\n";

foreach $setID (@SETids) { print "gene name: ", $setID, "\n"; }

print "\n";

foreach $rndID (@SETrnd) { print "gene name: ", $rndID, "\n"; }

#_______________________________________________________________________________________section1

$totaltableIDs = scalar @listtabIDs - 1;
@counterTABids = ( 0 .. $totaltableIDs );

foreach $counter (@counterTABids) {

    #print $counter,"\t",$listtabIDs[$counter],"\n";

    chomp $listtabIDs[$counter];

    foreach $setID (@SETids) {
        if ( $setID eq $listtabIDs[$counter] ) {
            print "FOUND setID: ", $setID, "\t", $listtabIDs[$counter], "\n";
            push @positionsSETids, $counter;
        }
    }

    foreach $rndID (@SETrnd) {
        if ( $rndID eq $listtabIDs[$counter] ) {
            print "FOUND rndID: ", $rndID, "\t", $listtabIDs[$counter], "\n";
            push @positionsRNDids, $counter;
        }
    }

}

$numberofSETids = scalar @positionsSETids;
$numberofRNDids = scalar @positionsRNDids;

#_______________________________________________________________________________________section2

$Nset             = 0;
$NsetCLEAN        = 0;
$Nsetna           = 0;
$ZsetTOTAL        = 0;
$ZsetTOTALabs     = 0;
$meanZsetTOTAL    = 0;
$meanZsetTOTALabs = 0;

@intestazione = split '\t', $LINEStab[0];
@POSITIONSsourceIDs = @positionsSETids;
shift @positionsSETids;
@POSITIONStargetIDs = @positionsSETids;

foreach $sourceposition (@POSITIONSsourceIDs) {

    @elements2 = split '\t', $LINEStab[ $sourceposition + 1 ];

    foreach $targetposition (@POSITIONStargetIDs) {

        print "SET: ", "\t", $intestazione[ $sourceposition + 1 ], "\t",
          $intestazione[ $targetposition + 1 ], "\t",
          $elements2[ $targetposition + 1 ],    "\n";
        print OUTPUTset1 $intestazione[ $sourceposition + 1 ], "\t",
          $intestazione[ $targetposition + 1 ], "\t",
          $elements2[ $targetposition + 1 ],    "\n";
        $Nset = $Nset + 1;

        if ( $elements2[ $targetposition + 1 ] ne "NA" ) {
            print OUTPUTset2 $intestazione[ $sourceposition + 1 ], "\t",
              $intestazione[ $targetposition + 1 ], "\t",
              $elements2[ $targetposition + 1 ],    "\n";
            $NsetCLEAN = $NsetCLEAN + 1;
            $ZsetTOTAL = $ZsetTOTAL + $elements2[ $targetposition + 1 ];
            $ZsetTOTALabs =
              $ZsetTOTALabs + abs( $elements2[ $targetposition + 1 ] );
        }
        else { $Nsetna = $Nsetna + 1; }

    }

    shift @POSITIONStargetIDs;

}

if ( $NsetCLEAN > 0 ) {
    $meanZset    = $ZsetTOTAL / $NsetCLEAN;
    $meanabsZset = $ZsetTOTALabs / $NsetCLEAN;
}
else { $meanZset = "NA"; $meanabsZset = "NA"; }

#_______________________________________________________________________________________section3

$Nrnd             = 0;
$NrndCLEAN        = 0;
$Nrndna           = 0;
$ZrndTOTAL        = 0;
$ZrndTOTALabs     = 0;
$meanZrndTOTAL    = 0;
$meanZrndTOTALabs = 0;

@intestazione = split '\t', $LINEStab[0];
@POSITIONSsourceIDs = @positionsRNDids;
shift @positionsRNDids;
@POSITIONStargetIDs = @positionsRNDids;

foreach $sourceposition (@POSITIONSsourceIDs) {

    @elements2 = split '\t', $LINEStab[ $sourceposition + 1 ];

    foreach $targetposition (@POSITIONStargetIDs) {

        print "RND: ", "\t", $intestazione[ $sourceposition + 1 ], "\t",
          $intestazione[ $targetposition + 1 ], "\t",
          $elements2[ $targetposition + 1 ],    "\n";

        print OUTPUTrnd1 $intestazione[ $sourceposition + 1 ], "\t",
          $intestazione[ $targetposition + 1 ], "\t",
          $elements2[ $targetposition + 1 ],    "\n";
        $Nrnd = $Nrnd + 1;

        if ( $elements2[ $targetposition + 1 ] ne "NA" ) {
            print OUTPUTrnd2 $intestazione[ $sourceposition + 1 ], "\t",
              $intestazione[ $targetposition + 1 ], "\t",
              $elements2[ $targetposition + 1 ],    "\n";
            $NrndCLEAN = $NrndCLEAN + 1;
            $ZrndTOTAL = $ZrndTOTAL + $elements2[ $targetposition + 1 ];
            $ZrndTOTALabs =
              $ZrndTOTALabs + abs( $elements2[ $targetposition + 1 ] );
        }
        else { $Nrndna = $Nrndna + 1; }

    }

    shift @POSITIONStargetIDs;

}

if ( $NrndCLEAN > 0 ) {
    $meanZrnd    = $ZrndTOTAL / $NrndCLEAN;
    $meanabsZrnd = $ZrndTOTALabs / $NrndCLEAN;
}
else { $meanZrnd = "NA"; $meanabsZrnd = "NA"; }

#_______________________________________________________________________________________section4

print "ID set", "\t", "n", "\t", "total pairwise combinations", "\t",
  "pairwise combinations with z", "\t", "pairwise combinations without z",
  "\t", "total z", "\t", "total absolute z", "\t", "mean z", "\t",
  "mean absolute z", "\n";

print $ARGV[1], "\t", $numberofSETids, "\t", $Nset, "\t", $NsetCLEAN, "\t",
  $Nsetna, "\t", $ZsetTOTAL, "\t", $ZsetTOTALabs, "\t", $meanZset, "\t",
  $meanabsZset, "\n";

print OUTPUTset3 "ID set", "\t", "n", "\t", "total pairwise combinations",
  "\t", "pairwise combinations with z", "\t",
  "pairwise combinations without z", "\t", "total z", "\t", "total absolute z",
  "\t", "mean z", "\t", "mean absolute z", "\n";

print OUTPUTset3 $ARGV[1], "\t", $numberofSETids, "\t", $Nset, "\t",
  $NsetCLEAN, "\t", $Nsetna, "\t", $ZsetTOTAL, "\t", $ZsetTOTALabs, "\t",
  $meanZset, "\t", $meanabsZset, "\n";

print OUTPUTset4 "ID set", "\t", "n", "\t", "total pairwise combinations",
  "\t", "pairwise combinations with z", "\t",
  "pairwise combinations without z", "\t", "total z", "\t", "total absolute z",
  "\t", "mean z", "\t", "mean absolute z", "\n";

print OUTPUTset4 $ARGV[1], "\t", $numberofSETids, "\t", $Nset, "\t",
  $NsetCLEAN, "\t", $Nsetna, "\t", $ZsetTOTAL, "\t", $ZsetTOTALabs, "\t",
  $meanZset, "\t", $meanabsZset, "\n";

#_______________________________________________________________________________________section5

print "ID set", "\t", "n", "\t", "total pairwise combinations", "\t",
  "pairwise combinations with z", "\t", "pairwise combinations without z",
  "\t", "total z", "\t", "total absolute z", "\t", "mean z", "\t",
  "mean absolute z", "\n";

print $ARGV[2], "\t", $numberofRNDids, "\t", $Nrnd, "\t", $NrndCLEAN, "\t",
  $Nrndna, "\t", $ZrndTOTAL, "\t", $ZrndTOTALabs, "\t", $meanZrnd, "\t",
  $meanabsZrnd, "\n";

print OUTPUTrnd3 "ID set", "\t", "n", "\t", "total pairwise combinations",
  "\t", "pairwise combinations with z", "\t",
  "pairwise combinations without z", "\t", "total z", "\t", "total absolute z",
  "\t", "mean z", "\t", "mean absolute z", "\n";

print OUTPUTrnd3 $ARGV[2], "\t", $numberofRNDids, "\t", $Nrnd, "\t",
  $NrndCLEAN, "\t", $Nrndna, "\t", $ZrndTOTAL, "\t", $ZrndTOTALabs, "\t",
  $meanZrnd, "\t", $meanabsZrnd, "\n";

print OUTPUTrnd4 "ID set", "\t", "n", "\t", "total pairwise combinations",
  "\t", "pairwise combinations with z", "\t",
  "pairwise combinations without z", "\t", "total z", "\t", "total absolute z",
  "\t", "mean z", "\t", "mean absolute z", "\n";

print OUTPUTrnd4 $ARGV[2], "\t", $numberofRNDids, "\t", $Nrnd, "\t",
  $NrndCLEAN, "\t", $Nrndna, "\t", $ZrndTOTAL, "\t", $ZrndTOTALabs, "\t",
  $meanZrnd, "\t", $meanabsZrnd, "\n";

#_______________________________________________________________________________________

close INPUT1;
close INPUT2;
close INPUT3;

close OUTPUTset1;
close OUTPUTrnd1;
close OUTPUTset2;
close OUTPUTrnd2;
close OUTPUTset3;
close OUTPUTrnd3;
close OUTPUTset4;
close OUTPUTrnd4;

#_______________________________________________________________________________________

