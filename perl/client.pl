#!/usr/bin/perl -w
# cliIO.pl
# a simple client using IO:Socket
#----------------


# notes: always print newlines
use strict;
use IO::Socket;

#my $host = 'localhost';
my $host = '192.168.2.13';
my $port = 5743;
my $sock = new IO::Socket::INET( PeerAddr => $host, PeerPort => $port, Proto => 'tcp');
$sock or die "no socket :$!";
$sock->autoflush();

while (my $buf=<$sock>) {
	print $buf;
}


close $sock;

