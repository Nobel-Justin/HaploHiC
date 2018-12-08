package HaploHiC::LoadOn;

use strict;
use warnings;
use FindBin qw/ $RealBin /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              $V_Href
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::LoadOn';
#----- version --------
$VERSION = "0.50";
$DATE = '2018-11-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#------ variants -------
our $V_Href = {

    #++++++++++++++#
    # Main Version #
    #++++++++++++++#

    MainName => 'HaploHiC.pl',
    Version => '0.14',
    Date => '2018-11-14',
    AUTHOR => 'Wenlong Jia',
    EMAIL => 'wenlongkxm@gmail.com',

    # functions
    RealBin => $RealBin,
    func => {},
    command => undef,
    argv_Aref => undef,
    run_mode => 0,

    ## manual
    HELP => 0,
    HELP_INFO => {},

    # for debug to keep some key intermediate folders/files
    in_debug => 0,

    #++++++++++#
    # settings #
    #++++++++++#

    # global setting
    ref_version => undef,
    # setting about dump func
    dump_allowNormMtd => {
                            'NONE'  => 1,
                            'VC'    => 1,
                          'VC_SQRT' => 1,
                            'KR'    => 1
                         },
    dump_allowBinSize => {
                            BP =>   {
                                        '2.5MB' => 2500000,
                                        '1MB'   => 1000000,
                                        '500KB' => 500000,
                                        '250KB' => 250000,
                                        '100KB' => 100000,
                                        '50KB'  => 50000,
                                        '25KB'  => 25000,
                                        '10KB'  => 10000,
                                        '5KB'   => 5000
                                    },
                            FRAG => {
                                        '500'   => 500,
                                        '200'   => 200,
                                        '100'   => 100,
                                        '50'    => 50,
                                        '20'    => 20,
                                        '5'     => 5,
                                        '2'     => 2,
                                        '1'     => 1
                                    }
                        },

    #++++++++++#
    # variants #
    #++++++++++#

    # intermediated containers
    ## for dump (JuicerDumpBP and JuicerDumpFRAG)
    ### ChrThings -> { $chr => { chr => $chr, len => $length, turn => $. } }
    ChrThings => {},

    #++++++++++#
    # software #
    #++++++++++#

    # only allows xx.xx, donot write third number
    java => 'java',
    java_minimum_version => 1.8
};

#--------- functions in this pm --------#
my @functoion_list = qw/
                        load_variants_dict
                     /;

#--- load variant dict Href ---
sub load_variants_dict{
    return $V_Href;
}

#--- 
1; ## tell the perl script the successful access of this module.
