open IN ,"tail -10 my.log |";

my $prev=0;
my $first=0;
print <<XML;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

<svg width="1200px" height="600px" version="1.1"
viewBox="0 -4 10 4"
xmlns="http://www.w3.org/2000/svg">
<g transform="scale(5,-1)">
<path d="M
XML

while (<IN>) {
/: ([\d.]+)/;
printf("%0.4f %0.4f ", ($1-$first)/60, 3.6/($1-$prev)) if $first;
$first=$1 unless $first;

$prev=$1;
}

print <<XML;
" style="fill:white;stroke:blue;stroke-width:0.02;"/>
<line y1="0" x1="0%" x2="100%" y2="0" style="stroke:red;stroke-width:0.01"/>
<line y1="1" x1="0%" x2="100%" y2="1" style="stroke:red;stroke-width:0.01"/>
<line y1="2" x1="0%" x2="100%" y2="2" style="stroke:red;stroke-width:0.01"/>
<line y1="3" x1="0%" x2="100%" y2="3" style="stroke:red;stroke-width:0.01"/>
<line y1="4" x1="0%" x2="100%" y2="4" style="stroke:red;stroke-width:0.01"/>
</g>
</svg> 
XML

close IN;
