package HaploHiC::HaploHiC;

use strict;
use warnings;

# basic
use HaploHiC::LoadOn;
# functions
use HaploHiC::Extensions::ShellForJuicer;
use HaploHiC::Extensions::JuicerDumpBP;
use HaploHiC::Extensions::JuicerDumpFRAG;
use HaploHiC::Extensions::FRAG2Gene;
use HaploHiC::PhasedHiC::HaploDivideReads;

#--- 
1; ## tell the perl script the successful access of this module.
