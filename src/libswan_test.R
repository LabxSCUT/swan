#!/usr/bin/env Rscript
suppressMessages(library(swan))
.Call( "swan_unit_test", PACKAGE = "swan" )
