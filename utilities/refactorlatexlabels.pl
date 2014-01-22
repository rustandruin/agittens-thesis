#!/usr/bin/env perl

use File::Basename qw(fileparse);
use File::Spec::Functions qw(catfile);

if ($#ARGV != 1) {
	print "usage: refactorlatexlabels.pl file.tex labels\n";
	exit;
}

($texinfname, $path, $suffix) = fileparse($ARGV[0], "\.[^.]*");

if (lc($suffix) ne ".tex") {
	print "usage: refactorlatexlabels.pl file.tex labels\n";
	exit;
}

$outfname = catfile($path, "$texinfname$suffix.refactored");

open TEXFIN, "<$ARGV[0]" or die "Cannot open tex file: $!";
open LABELSIN, "<$ARGV[1]" or die "Cannot open label file: $!";
open TEXFOUT, ">$outfname" or die "Cannot open output tex file: $!";

@labels = <LABELSIN>;

foreach $pair (@labels) {
	chomp($pair);
	@labelpair = split('\t', $pair);
	$labelsubs{$labelpair[0]} = $labelpair[1];
	$labelmatches{$labelpair[0]} = 0;
}

while(<TEXFIN>) {
	foreach $inlabel (keys(%labelsubs)) {
		if (m/\Q$inlabel\E/) {
			$labelmatches{$inlabel} =
	$labelmatches{$inlabel} +  ( s/\Q$inlabel\E/$labelsubs{$inlabel}/g );
		}
	}
	print TEXFOUT $_;
}

close(TEXFIN);
close(LABELSIN);
close(TEXFOUT);

print "\n";
foreach $inlabel (keys(%labelsubs)) {
	print "Replaced $inlabel with $labelsubs{$inlabel} $labelmatches{$inlabel} times\n";
}

