
use strict;
use warnings;

my $num_win=0;
my $win_time=0;
my $sum_win_time=0;
my $avg_win_time=0;

open OUT, ">", "$ARGV[0]-CountStat.txt";
print "$ARGV[0]-CountStat.txt";
open(IN,"$ARGV[0]");
	while (<IN>) {
		if ($_ =~ /\S+;(\d+)/) {
		$win_time=$1;
			if ($win_time < 9999 ) {
				$num_win++;
				$sum_win_time=$sum_win_time+$win_time;
			}
		}
	}
close IN;


if ($num_win >0) {
	$avg_win_time=$sum_win_time/$num_win;
	print OUT "Number of wins = $num_win\n";
	print OUT "Average win time: $avg_win_time\n";
}
close OUT;