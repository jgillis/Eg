use Regexp::Common qw /balanced/;

my $_='abc (fo + /* 12 / (5+30) -t^2 + (18+g)^2)^4 - sqrt(x_ff0^2+y^2) +sin(4+3)^2';
  while (s/((\w*)$RE{balanced}{-parens=>'()'})\^(\d)/'pow(' . ($2?$1:substr($1,1,-1)) . ',' . $3 . ')'/ge) {};
  while (s/sqrt($RE{balanced}{-parens=>'()'})/'pow(' . substr($1,1,-1) . ',0.5)'/ge) {};
  while (s/(\w+)\^(\d)/power($1,$2)/g) {};
print;
