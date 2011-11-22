#!/bin/bash
python server.py & # Start the server
python client.py foobar
python client.py baz
kill %1 # Kill the server
