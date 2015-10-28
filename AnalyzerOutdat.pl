
use strict;
use warnings;

my $filename="";
my $won=0;

open OUT,">","$ARGV[0]-AnalyzerOutdat-$ARGV[1].txt";
print "$ARGV[0]-AnalyzerOutdat-$ARGV[1].txt";
open(IN1,"$ARGV[0]");
	while (<IN1>) {
		chomp;
		$filename=$_;
		open(IN2,"$filename");
			$won=0;
			while (<IN2>) {
				if ($_ =~ /t=\s+(\d+)\s+nW=\s+\d+\s+nR=\s+(\d+)/ ) {
					if (($won== 0) && ($1 > 3000) && ($2<=$ARGV[1])) {
						print OUT "$filename;$1\n";
						$won=1;
					}
				}
			}
			if ($won == 0) {
				print OUT "$filename;9999\n";
			}
		close IN2;
	}
close IN1;
close OUT;

