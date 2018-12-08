package HaploHiC::Extensions::JuicerDumpFRAG;

use strict;
use warnings;
use Getopt::Long;
use File::Spec::Functions qw/ catfile /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use HaploHiC::LoadOn;
use HaploHiC::Extensions::JuicerDumpBP qw/ check_juicer_files make_workspace load_chr_Things /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              JuicerDumpFRAGWrok
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::Extensions::JuicerDumpFRAG';
#----- version --------
$VERSION = "0.04";
$DATE = '2018-11-14';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        JuicerDumpFRAGWrok
                        get_merged_dump_frag
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} dump_FRAG <[Options]>
     
     Options:

       # run juicer_tools dump func in FRAG mode and merge all chromosomes to one final file #

       # Inputs and Outputs #
        -juicer [s]  path of juicer folder, which has sub-folder, e.g., 'CPU'. <required>
        -topdir [s]  topdir stores juicer results. <required>
                      this func creates '$V_Href->{workspace_name}-<postfix>' folder under topdir.
        -db_dir [s]  database folder made by 'juicer_db' function. <required>
        -ref_v  [s]  version of reference genome, see contents under 'db_dir'. <required>

       # Options #
        -hic_fl [s]  prefix of '.hic' file: 'inter_30' or 'inter'. [inter_30]
        -bin    [s]  bin size. [1]
                      Note: FRAG mode allows: 500, 200, 100, 50, 20, 5, 2, 1
        -ctype  [s]  'observed' or 'oe'. [observed]
        -norm   [s]  kind of normalization, must be one of NONE/VC/VC_SQRT/KR. ['NONE']
                      notes from juicer_tools dump usage:
                       1) VC is vanilla coverage,
                       2) VC_SQRT is square root of vanilla coverage;
                       3) KR is Knight-Ruiz or Balanced normalization.
        -gzip        gzip compress the contact txt files. [disabled]

        -h|help      Display this help info.

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
    $V_Href->{ $_->[0] } = $_->[1] for
        map {
            if( !exists $V_Href->{$_->[0]} ){
                ( $_ );
            }
            else{
                warn_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
            }
        }
        (
            # input/output
            [ juicer_dir => undef ],
            [ JuicerTopdir => undef ],
            [ db_dir => undef ],

            # options
            [ hic_type => 'inter_30' ],
            [ bin_size => '1' ],
            [ count_type => 'observed' ],
            [ norm_method => 'NONE' ],
            [ gzip_output => 0 ],

            # settings
            [ dump_resolUnitMode => 'FRAG' ],

            # intermediate variants
            [ hic_file => undef ],
            [ chrLenFile => undef ],
            [ juicer_tool_Jar => undef ],
            [ workspace_name => 'extract_matrix' ],
            [ workspace => undef ],
            [ MergedFile => 'InterIntraChr.merged.txt' ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['JuicerTopdir'],
                                  ['db_dir'],
                                  ['juicer_dir']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-juicer:s" => \$V_Href->{juicer_dir},
        "-topdir:s" => \$V_Href->{JuicerTopdir},
        "-db_dir:s" => \$V_Href->{db_dir},
        "-ref_v:s"  => \$V_Href->{ref_version},
        # options
        "-hic_fl:s" => \$V_Href->{hic_type},
        "-bin:s"    => \$V_Href->{bin_size},
        "-ctype:s"  => \$V_Href->{count_type},
        "-norm:s"   => \$V_Href->{norm_method},
        "-gzip"     => \$V_Href->{gzip_output},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !defined $V_Href->{juicer_dir} || !-d $V_Href->{juicer_dir}
             || !defined $V_Href->{JuicerTopdir} || !-d $V_Href->{JuicerTopdir}
             || !defined $V_Href->{db_dir} || !-d $V_Href->{db_dir}
             || !defined $V_Href->{ref_version}
             || !exists $V_Href->{dump_allowBinSize}->{$V_Href->{dump_resolUnitMode}}->{$V_Href->{bin_size}}
             || !exists $V_Href->{dump_allowNormMtd}->{$V_Href->{norm_method}}
             || ( $V_Href->{hic_type} ne 'inter' && $V_Href->{hic_type} ne 'inter_30' )
             || ( $V_Href->{count_type} ne 'observed' && $V_Href->{count_type} ne 'oe' )
            );
}

#--- run juicer_tools dump func in FRAG mode and merge all chromosomes ---
sub JuicerDumpFRAGWrok{
    check_juicer_files;

    make_workspace;

    load_chr_Things;

    &get_merged_dump_frag;
}

#--- get intar-chr and inter-chr and merge ---
sub get_merged_dump_frag{

    my @chr_Href = sort { $a->{turn} <=> $b->{turn} } values %{$V_Href->{ChrThings}};

    for my $i ( 0 .. $#chr_Href ){
        my $i_chr = $chr_Href[$i]->{chr};
        my $i_chr_workspace = catfile( $V_Href->{workspace}, "$i_chr" );
        `rm -rf $i_chr_workspace && mkdir -p $i_chr_workspace`;
        for my $j ( $i .. $#chr_Href ){
            my $j_chr = $chr_Href[$j]->{chr};
            my $ij_FRAG_file = catfile( $i_chr_workspace, "$i_chr.vs.$j_chr.FRAG.txt" );
            my $command  = "($V_Href->{java} -jar $V_Href->{juicer_tool_Jar} dump $V_Href->{count_type} $V_Href->{norm_method} $V_Href->{hic_file} $i_chr $j_chr FRAG $V_Href->{bin_size} $ij_FRAG_file)";
               $command .= " && (gzip $ij_FRAG_file)" if( $V_Href->{gzip_output} );
            &trible_run_for_success( $command, "$i_chr-vs-$j_chr-dump-FRAG" );
            # recode
            $chr_Href[$i]->{frag}->{$j_chr} = $V_Href->{gzip_output} ? "$ij_FRAG_file.gz" : $ij_FRAG_file;
        }
        # inform
        stout_and_sterr "[INFO]\t".`date`
                             ."\tdump FRAG deals with $i_chr OK\n";
    }

    # merge to one
    $V_Href->{MergedFile} = catfile( $V_Href->{workspace}, $V_Href->{MergedFile} );
    $V_Href->{MergedFile} .= '.gz' if( $V_Href->{gzip_output} );
    open (MERGE, Try_GZ_Write($V_Href->{MergedFile})) || die "cannot write merge file: $!\n";
    for my $i ( 0 .. $#chr_Href ){
        my $i_chr = $chr_Href[$i]->{chr};
        for my $j ( 0 .. $#chr_Href ){
            my $j_chr = $chr_Href[$j]->{chr};
            next unless exists($chr_Href[$i]->{frag}->{$j_chr});
            my $ij_FRAG_file = $chr_Href[$i]->{frag}->{$j_chr};
            open (IJCHR, Try_GZ_Read($ij_FRAG_file)) || die "fail to read ijchr file: $!\n";
            while(<IJCHR>){
                my ($i_idx, $j_idx, $contacts) = (split);
                print MERGE join("\t", $i_chr, $i_idx, $j_chr, $j_idx, $contacts)."\n";
            }
            close IJCHR;
        }
    }
    close MERGE;
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tmerge dump FRAG results OK\n"
                         ."\t$V_Href->{MergedFile}\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
