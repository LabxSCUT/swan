#!/bin/bash

#automated script to run SWAN test
#the log of tests will be in all.log

if [ ! -e "swan_test.tgz" ]; then wget http://meta.usc.edu/softs/swan/swan_test.tgz; tar -zxvf swan_test.tgz; fi
($SWAN_BIN/paired.sh all && $SWAN_BIN/single.sh all) >all.log 2>&1
