#!/usr/bin/perl -w
open (IN, $ARGV[0]) || die "Can't open $ARGV[0]";
open (OUT, ">", "$ARGV[0].csv");
print OUT "t;selfs;non-selfs;Bcells;antibodies;Tcells;interleukins;bone marrows\n";
while (<IN>) {
	if ($_ =~ /t=\s+(\S+)\s+nW=\s+(\S+)\s+nR=\s+(\S+)\s+nB=\s+(\S+)\s+nAb=\s+(\S+)\s+nTh=\s+(\S+)\s+nIL=\s+(\S+)\s+nM=\s+(\S+)/ ) {
		#print $_;
		print OUT "$1;$2;$3;$4;$5;$6;$7;$8\n";
	}
}
close IN;
close OUT;