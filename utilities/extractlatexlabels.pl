#!/usr/bin/env perl

use File::Basename qw(fileparse);
use File::Spec::Functions qw(catfile);

if ($#ARGV != 0) {
	print "usage: extractlatexlabels.pl file.tex\n";
	exit;
}

($texinfname, $path, $suffix) = fileparse($ARGV[0], "\.[^.]*");

if (lc($suffix) ne ".tex") {
	print "usage: extractlatexlabels.pl file.tex\n";
	exit;
}

$outfname = catfile($path, "$texinfname$suffix.labels");

open TEXFIN, "<$ARGV[0]" or die "Cannot open tex file: $!";
open LABELSOUT, ">$outfname" or die "Cannot open label file: $!";

while(<TEXFIN>) {
	if (m/\\ref{(.*)
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

