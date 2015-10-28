
use strict;
use warnings;

open IN, '>', "randomnumbers.txt";
for (my $i=1;$i<=500;$i++) {
	my $rnd=int(rand(20000));
	system ("immune3_dec17_2014 $rnd");
	print IN "$rnd\n";
}
close IN;
