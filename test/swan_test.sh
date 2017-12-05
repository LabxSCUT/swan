#!/bin/bash

#automated script to run SWAN test
#the log of tests will be in all.log

if [ ! -e "swan_test.tgz" ]; then wget https://s3-us-west-2.amazonaws.com/lixiabucket/swan_test.tgz; fi
tar -zxvf swan_test.tgz
(paired.sh all && single.sh all) >all.log 2>&1
