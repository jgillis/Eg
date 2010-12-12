
    $FIFO = 'apipe';

    while (1) {
        unless (-p $FIFO) {
            unlink $FIFO;
            require POSIX;
            POSIX::mkfifo($FIFO, 0700)
                or die "can't mkfifo $FIFO: $!";
        }

        # next line blocks until there's a reader
        open (FIFO, "> $FIFO") || die "can't write $FIFO: $!";
        print FIFO "John Smith (smith\@host.org)\n", `fortune -s`;
        close FIFO;
        sleep 2;    # to avoid dup signals
    }
