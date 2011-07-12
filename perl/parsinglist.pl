my $s=shift @ARGV;
use Text::Balanced qw (extract_bracketed);

$s = "12,x[1,1],[123],x[2+3,2]-y[1,2],sin(f,x),sdsd";
   
print "$s\n\n";


my @cols;

my %counter;
$counter{"{}"}=0;
$counter{"[]"}=0;
$counter{"()"}=0;

print $start=0;
for (my $i=0;$i<length($s);$i++) {
  $char = substr($s,$i,1);
  for my $key (keys %counter) { 
    $counter{$key}++ if ($char eq substr($key,0,1));
    $counter{$key}-- if ($char eq substr($key,1,1));
  }
  print $char;
  print ":";
  print $counter{"{}"};  print $counter{"[]"};  print $counter{"()"};print "\n";
  if ($char eq ',') {
    my $nest=0;
    $nest = $nest || $_ for (values %counter);
    if ($nest==0) {
      push @cols, substr($s,$start,$i-$start);
      $start = $i+1;
    }
  }
}
push @cols, substr($s,$start,length($s)-$start);
print "$s\n\n";
print ":$_:\n" for @cols;

