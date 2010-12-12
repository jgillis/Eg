use Sys::Syslog;  
my $name="ikke";
openlog($name, "ndelay", "local0");

syslog("info|local0", "testje 123....blahblah error hier en daar ");
