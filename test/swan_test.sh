#!/bin/bash

#automated script to run SWAN test
#the log of tests will be in all.log
#please provide $SWAN_BIN in $1

if [ ! -e "swan_test.tgz" ]; then wget https://s3-us-west-2.amazonaws.com/lixiabucket/swan_test.tgz; fi
tar -zxvf swan_test.tgz
SWAN_BIN=$1
echo SWAN_BIN=$SWAN_BIN
echo "./paired.sh $SWAN_BIN all >all.paired.log 2>&1"
./paired.sh $SWAN_BIN all >all.paired.log 2>&1 
echo "./single.sh $SWAN_BIN all >all.single.log 2>&1"
./single.sh $SWAN_BIN all >all.single.log 2>&1
