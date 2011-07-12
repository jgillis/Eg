use Text::Balanced qw (extract_bracketed);

   
my $s="[1,2[5,6],3],[4,5,6],[7,[5,6],9],[4,[3],6]";


while (my $row = extract_bracketed( $s, '[]' )) {
  $s = substr($s,1);
  print ":$row:\n";
}

