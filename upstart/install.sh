#!/bin/bash


sudo cp helloworld.conf /etc/init
# A symbolic link will not suffice


# If this is the very first run, reboot your PC or use this:
sudo initctl reload-configuration
# If you are just updating helloworld.conf, it's done now. upstart will have seen the change and applied it.
