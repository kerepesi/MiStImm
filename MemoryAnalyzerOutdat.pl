:#!C:\strawberry\perl\bin\perl.exe

use strict;
use warnings;

my $filename="";
my $won1=0;
my $won2=0;
my $won1_time=0;
my $won2_time=0;

open(IN1,"$ARGV[0]");
	while (<IN1>) {
		chomp;
		$filename=$_;
		open(IN2,"$filename");
			$won1=0;
			$won2=0;
			while (<IN2>) {
				if ($_ =~ /t=\s+(\d+)\s+nW=\s+\d+\s+nR=\s+(\d+)/ ) {
					if (($won1== 0) && ($1 > 3000) && ($1<=3150) && ($2<=50)) {
						$won1=1;
						$won1_time=$1;
					}elsif (($won2== 0) && ($1 > 3150) && ($2<=50)) {
						$won2=1;
						$won2_time=$1;
					}
				}
			}
			if ($won1) {
				my $t1=$won1_time-3000;
				print"$filename;$t1;";
			}else {
				print"$filename;9999;";
			}
			if ($won2) {
				my $t2=$won2_time-3150;
				print"$t2\n";
			}else {
				print"9999\n";
			} 
		close IN2;
	}
close IN1;

