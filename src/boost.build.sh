#!/bin/sh
cd boost_tree
./bootstrap.sh
./b2 --with-regex
cd ..
