#!/bin/bash

# Get and extract
wget http://bitbucket.org/gavinband/bgen/get/master.tar.gz
tar -zxvf master.tar.gz

# Configure and build
cd gavinband-bgen-d1f03a2c308a
./waf configure
./waf

# bgenix will be installed in build/apps/bgenix

# If you have ~/bin in your $PATH then add bgenix to ~/bin, e.g.
cd ~/bin
ln -s ~/programs/gavinband-bgen-d1f03a2c308a/build/apps/bgenix

