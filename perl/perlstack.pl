#!/usr/bin/perl -w
# server using IO::Socket
#---------------------
# http://docstore.mik.ua/orelly/perl/advprog/ch12_03.htm
use strict;
use IO::Socket;
$SIG{CHLD} = sub {wait ()};

print "Perlstack\n";
my $sock = new IO::Socket::INET(
LocalHost => '192.168.2.13',
LocalPort => 5743,
Proto => 'tcp',
Listen => SOMAXCONN,
Reuse => 1);

$sock or die "no socket :$!";
$sock->autoflush();


my $size=10;
my @buffer=(1..$size);

my($new_sock, $c_addr, $buf);

while (($new_sock, $c_addr) = $sock->accept())
{
	my $pid = fork();
    	die "Cannot fork: $!" unless defined($pid);
	if ($pid == 0) {
		my ($client_port, $c_ip) =sockaddr_in($c_addr);
		my $client_ipnum = inet_ntoa($c_ip);
		my $client_host =gethostbyaddr($c_ip, AF_INET);
		print "got a connection from: [$client_ipnum] \n";
		if ($client_ipnum ne "192.168.2.13") {
			print $new_sock join("\n",@buffer);
			print $new_sock "\n";
			exit(0);
		} else {
			while (defined ($buf = <$new_sock>))
			{	
				chomp $buf;
				shift @buffer;
				push @buffer, $buf;
				print join("\n",@buffer);
				print "\n";
			}
			exit(0);
		}
	}

} 
