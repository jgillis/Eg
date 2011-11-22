#!/bin/bash
./server 5001 &
echo "Hello world" | ./client localhost 5001
kill %1
