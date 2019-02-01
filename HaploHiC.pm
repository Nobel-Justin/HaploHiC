package HaploHiC::HaploHiC;

use strict;
use warnings;

# basic
use HaploHiC::LoadOn;
# functions
use HaploHiC::Extensions::JuicerDB;
use HaploHiC::Extensions::ShellForJuicer;
use HaploHiC::Extensions::JuicerDump;
use HaploHiC::Extensions::FRAG2Gene;
use HaploHiC::PhasedHiC::HaploDivideReads;
use HaploHiC::PhasedHiC::adjustDumpContacts;

#--- 
1; ## tell the perl script the successful access of this module.
