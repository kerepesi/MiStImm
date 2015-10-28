#!/usr/bin/perl -w
#0: inputfile
#1: start
#2: end
#3: step
#Typical running command: perl Out2Vizu.pl out2-1400092876 0 3500 100 area > ki
my %sum_reprod=();
my %num_reprod_type=();
my $num_reprod_type=0;
my $step=0;
print(qq~<html>
    <head>
    <script type="text/javascript" src="https://www.google.com/jsapi"></script>
    <script type="text/javascript">
      google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = google.visualization.arrayToDataTable([
~);
open (IN, $ARGV[0]) || die "Can't open $ARGV[0]";
while (<IN>) { 
        if ($_ =~ /^\s+(\D+reprod)\D+(\d+\.\d+)/ ) {
            push(@reprod_type,$1);
            push(@reprod_time,$2); 
            if (! exists $num_reprod_type{$1}) {
                $num_reprod_type{$1}=1;
                $num_reprod_type+=1;
            }
        }
}
close IN;
if (! $#reprod_type == $#reprod_time ) {
    die "$reprod_type not equal $reprod_time ";
}
print"['Time steps'";
foreach my $key (sort keys %num_reprod_type) {
    print", '$key'";
}
print"],\n";
my $j=0;
for($i = $ARGV[1]; $i < $ARGV[2]; $i+=$ARGV[3]) {
    $step=(2*$i+$ARGV[3])/2;
    print("['$step'");
    while ($reprod_time[$j] <= $i) {
        if ($reprod_time[$j] >= $ARGV[1]) {
            $sum_reprod{$reprod_type[$j]}{$step}+=1;
        }
        $j+=1;
    }
    foreach my $key (sort keys %num_reprod_type) {
        if (exists $sum_reprod{$key}{$step}) {
                print(", $sum_reprod{$key}{$step}");
        } else {
            print", 0";
        }
    }
    print("],\n")
}
if ($ARGV[4] eq "line") {
 print(qq~        ]);
        var options = {
          title: 'Number of the cell reproductions'
        };

        var chart = new google.visualization.LineChart(document.getElementById('chart_div'));
        chart.draw(data, options);
      }
        </script>
    </head>
    <body>
    <div id="chart_div" style="width: 1200px; height: 700px;"></div>
    </body>
    </html>
 ~);
} elsif ($ARGV[4] eq "area") {
 print(qq~ ]);
        var options = {
          title: 'Number of the cell reproductions',
          vAxis: {title: 'Accumulated'},
          isStacked: true
        };

        var chart = new google.visualization.SteppedAreaChart(document.getElementById('chart_div'));
        chart.draw(data, options);
      }
    </script>
    </head>
    <body>
    <div id="chart_div" style="width: 1200px; height: 700px;"></div>
    </body>
    </html>
 ~);
}
