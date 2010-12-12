use warnings;

print "Content-type: image/svg+xml";
 
 use SVG::Graph;
use Data::Dumper;
  use SVG::Graph::Data;
  use SVG::Graph::Data::Datum;

  #create a new SVG document to plot in...
  my $graph = SVG::Graph->new(width=>1200,height=>600,margin=>30);

  #and create a frame to hold the data/glyphs

my $frame0 = $graph->add_frame;
  my $frame = $frame0->add_frame;


open IN ,"tail -10 my.log |";

my $prev=0;
my $first=0;

my @data;

while (<IN>) {
/: ([\d.]+)/;
push @data, SVG::Graph::Data::Datum->new(x=>($1-$first)/60,y => 3.6/($1-$prev)) if $first;
$first=$1 unless $first;

$prev=$1;
}


my $data = SVG::Graph::Data->new(data => \@data);

  #put the xy data into the frame
  $frame->add_data($data);

$frame->add_glyph('axis',x_fractional_ticks=>4,y_fractional_ticks=>6,'stroke'=>'black','stroke-width'=>2);


  #print the graphic
  print $graph->draw;


close IN;
