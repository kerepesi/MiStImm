
use strict;
use warnings;

open IN, '>', "randomnumbers.txt";
for (my $i=1;$i<=20;$i++) {
	my $rnd=int(rand(20000));
	system ("MiStImm_oct18_2017 $rnd");
	print IN "$rnd\n";
}
close IN;
