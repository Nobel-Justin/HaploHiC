package HaploHiC::Check;

use strict;
use warnings;
use BioFuse::Util::Sys qw/ check_java_version /;
use HaploHiC::LoadOn;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
	          check
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::Check';
#----- version --------
$VERSION = "0.02";
$DATE = '2019-01-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
						check
					 /;

#--- check para and database
sub check{

	my $second_cmd = $V_Href->{command};

    # version check, this part consumes lots of vmem: use '-Xmx1m', ~1.5G; OR >9G.
    if(    $second_cmd eq 'run_juicer'
        || $second_cmd eq 'juicerDump'
      ){
        check_java_version(javaPath => $V_Href->{java}, minVer => $V_Href->{java_minimum_version});
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
