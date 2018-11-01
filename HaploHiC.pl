#!/usr/bin/perl -w
use strict;
use HaploHiC::HaploHiC;
use HaploHiC::RunFunc;

my ($VERSION, $DATE, $AUTHOR, $EMAIL);

#----- version --------
$VERSION = "0.05";
$DATE = '2018-11-01';
#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

HaploHiC::RunFunc->options_alert_and_run( argv_Aref => \@ARGV );