package HaploHiC::Extensions::JuicerDumpBP;

use strict;
use warnings;
use Getopt::Long;
use File::Spec::Functions qw/ catfile /;
use BioFuse::Util::Sys qw/ file_exist trible_run_for_success /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read /;
use HaploHiC::LoadOn;
use HaploHiC::GetPath qw/ GetPath /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              JuicerDumpBPWrok
              make_workspace
              load_chr_Things
              check_juicer_files
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::Extensions::JuicerDumpBP';
#----- version --------
$VERSION = "0.06";
$DATE = '2018-08-05';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        JuicerDumpBPWrok
                        merge_inter_intra
                        convert_intraChr_PosToBinIdx
                        get_intraChr
                        get_interChr
                        make_workspace
                        load_chr_Things
                        check_juicer_files
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} dump_BP <[Options]>
     
     Options:

       # run juicer_tools dump func in BP mode and merge all chromosomes to one final file #

       # Inputs and Outputs #
        -juicer [s]  path of juicer folder, which has sub-folder, e.g., 'CPU'. <required>
        -topdir [s]  topdir stores juicer results. <required>
                      this func creates '$V_Href->{workspace_name}-<postfix>' folder under topdir.
        -db_dir [s]  database folder made by 'juicer_db' function. <required>

       # Options #
        -ref_v  [s]  version of reference genome. [hg19]
                      Note: this depends on the contents under 'db_dir'.
        -hic_fl [s]  prefix of '.hic' file: 'inter_30' or 'inter'. [inter_30]
        -bin    [s]  bin size. [1MB for BP]
                      Note: BP mode allows: 2.5MB, 1MB, 500KB, 250KB, 100KB, 50KB, 25KB, 10KB, 5KB
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
            [ bin_size => '1MB' ],
            [ count_type => 'observed' ],
            [ norm_method => 'NONE' ],
            [ gzip_output => 0 ],

            # settings
            [ dump_resolUnitMode => 'BP' ],

            # intermediate variants
            [ hic_file => undef ],
            [ chrLenFile => undef ],
            [ juicer_tool_Jar => undef ],
            [ workspace_name => 'extract_matrix' ],
            [ workspace => undef ],
            [ chrBinRangeInJuicer => 'chrBinRangeInJuicer.txt' ],
            [ ALLtoALLFile => 'ALLtoALL.txt' ],
            [ intraChr_dir => 'intraChr' ],
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
        # options
        "-ref_v:s"  => \$V_Href->{ref_version},
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
             || ( !defined $V_Href->{juicer_dir} || !-d $V_Href->{juicer_dir} )
             || ( !defined $V_Href->{JuicerTopdir} || !-d $V_Href->{JuicerTopdir} )
             || ( !defined $V_Href->{db_dir} || !-d $V_Href->{db_dir} )
             || !exists( $V_Href->{dump_allowBinSize}->{$V_Href->{dump_resolUnitMode}}->{$V_Href->{bin_size}} )
             || !exists( $V_Href->{dump_allowNormMtd}->{$V_Href->{norm_method}} )
             || ( $V_Href->{hic_type} ne 'inter' && $V_Href->{hic_type} ne 'inter_30' )
             || ( $V_Href->{count_type} ne 'observed' && $V_Href->{count_type} ne 'oe' )
            );
}

#--- run juicer_tools dump func in BP mode and merge all chromosomes ---
sub JuicerDumpBPWrok{
    &check_juicer_files;

    &make_workspace;

    &load_chr_Things;

    &get_interChr;

    &get_intraChr;

    &convert_intraChr_PosToBinIdx;

    &merge_inter_intra;
}

#--- check files ---
## also applied by JuicerDumpFRAG.pm
sub check_juicer_files{
    $V_Href->{hic_file} = GetPath(filekey => 'hic_file');
    $V_Href->{chrLenFile} = GetPath(filekey => 'chrLenFile');
    $V_Href->{juicer_tool_Jar} = GetPath(filekey => 'juicer_tool_Jar');
    # check existence
    if(    !-e $V_Href->{hic_file}
        || !-e $V_Href->{chrLenFile}
        || !-e $V_Href->{juicer_tool_Jar}
    ){
        warn "Please make sure the existence of files below.\n";
        warn "  $V_Href->{hic_file}\n";
        warn "  $V_Href->{chrLenFile}\n";
        warn "  $V_Href->{juicer_tool_Jar}\n";
        exit;
    }

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tfiles checking DONE\n";
}

#--- prepare workspace
## also applied by JuicerDumpFRAG.pm
sub make_workspace{
    $V_Href->{workspace} = GetPath(filekey => 'JuicerAlignedDir');
    if( !-w $V_Href->{workspace} ){
        warn "no write authority in topdir/aligned fold:\n";
        warn "  $V_Href->{workspace}\n";
        exit;
    }

    # output folder
    $V_Href->{workspace_name} = join('-', $V_Href->{workspace_name}, $V_Href->{hic_type}, $V_Href->{count_type}, $V_Href->{norm_method}, $V_Href->{dump_resolUnitMode}.$V_Href->{bin_size} );
    $V_Href->{workspace} = catfile( $V_Href->{JuicerTopdir}, 'aligned', $V_Href->{workspace_name} );
    `rm -rf $V_Href->{workspace} && mkdir -p $V_Href->{workspace}`;

    # for dump BP
    if( defined $V_Href->{chrBinRangeInJuicer} ){
        $V_Href->{chrBinRangeInJuicer} = catfile( $V_Href->{workspace}, $V_Href->{chrBinRangeInJuicer} );
    }

    # use real number of bin_size
    $V_Href->{bin_size} = $V_Href->{dump_allowBinSize}->{$V_Href->{dump_resolUnitMode}}->{ $V_Href->{bin_size} };

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tcreate workspace OK\n"
                         ."\t$V_Href->{workspace}\n";
}

#--- read chromosome info ---
## also applied by JuicerDumpFRAG.pm
sub load_chr_Things{
    # the chrLenFile chr-order is very important, we assume that it is identical to the juicer-order
    open (CHRLEN, $V_Href->{chrLenFile}) || die "fail read chrLenFile: $!\n";
    while(<CHRLEN>){
        next if(/^#/);
        my ($chr, $length) = (split)[0, 1];
        $V_Href->{ChrThings}->{$chr} = { chr=>$chr, len=>$length, turn=>$. };
    }
    close CHRLEN;
}

#--- get inter-chr
sub get_interChr{
    $V_Href->{ALLtoALLFile} = catfile( $V_Href->{workspace}, $V_Href->{ALLtoALLFile} );
    my $command  = "($V_Href->{java} -jar $V_Href->{juicer_tool_Jar} dump $V_Href->{count_type} $V_Href->{norm_method} $V_Href->{hic_file} ALL ALL BP $V_Href->{bin_size} $V_Href->{ALLtoALLFile})";
       $command .= " && (gzip $V_Href->{ALLtoALLFile})" if( $V_Href->{gzip_output} );
    &trible_run_for_success( $command, 'ALLtoALL-interChr' );

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tget ALL inter-chr matrix-count OK\n";
}

#--- get intra-chr
sub get_intraChr{
    $V_Href->{intraChr_dir} = catfile( $V_Href->{workspace}, $V_Href->{intraChr_dir} );
    `rm -rf $V_Href->{intraChr_dir} && mkdir -p $V_Href->{intraChr_dir}`;

    my @chr_Href = sort { $a->{turn} <=> $b->{turn} } values %{$V_Href->{ChrThings}};
    for my $Href ( @chr_Href ){
        my $chr = $Href->{chr};
        my $ctChrFile = catfile( $V_Href->{intraChr_dir}, "$chr.txt" );
        $Href->{PosCt_file} = $ctChrFile;
        my $command  = "($V_Href->{java} -jar $V_Href->{juicer_tool_Jar} dump $V_Href->{count_type} $V_Href->{norm_method} $V_Href->{hic_file} $chr $chr BP $V_Href->{bin_size} $ctChrFile)";
           $command .= " && (gzip $ctChrFile)" if( $V_Href->{gzip_output} );
        &trible_run_for_success( $command, "$chr-intraChr" );
        # inform
        stout_and_sterr "[INFO]\t".`date`
                         ."\tget intra-$chr matrix-count OK\n";
    }
}

#--- convert the intraChr from Pos to BinIdx
sub convert_intraChr_PosToBinIdx{

    open (CHRBR, ">$V_Href->{chrBinRangeInJuicer}") || die "fail write chr bin-range file.\n";

    my $preChrBinSum = 0;
    my @chr_Href = sort { $a->{turn} <=> $b->{turn} } values %{$V_Href->{ChrThings}};
    for my $Href ( @chr_Href ){
        my $chr = $Href->{chr};
        my $chrlen = $Href->{len};
        my $ctChrFile = $Href->{PosCt_file};
        $ctChrFile .= '.gz' if( $V_Href->{gzip_output} );
        my $ctChrFileBinIdx = $Href->{PosCt_file}.".byBinIdx";
        $Href->{BinCt_file} = $ctChrFileBinIdx;
        open (BINCT, ">$ctChrFileBinIdx") || die "fail write bin-index count file of $chr\n";
        open (POSCT, Try_GZ_Read($ctChrFile)) || die "fail read count file of $chr\n";
        while(<POSCT>){
            my ($i_BinStP, $j_BinStP, $ij_Count) = (split)[0,1,2];
            if(    $i_BinStP % $V_Href->{bin_size} != 0
                || $j_BinStP % $V_Href->{bin_size} != 0
            ){
                warn "Not integer times of BinSize ($V_Href->{bin_size}):\n";
                warn "  i_BinStP: $i_BinStP; j_BinStP=$j_BinStP\n";
                warn "  $ctChrFile\n";
                exit;
            }
            # Pos to BinIdx
            my $i_BinIdx = $i_BinStP / $V_Href->{bin_size} + $preChrBinSum;
            my $j_BinIdx = $j_BinStP / $V_Href->{bin_size} + $preChrBinSum;
            print BINCT join("\t", $i_BinIdx, $j_BinIdx, $ij_Count)."\n";
        }
        close POSCT;
        close BINCT;
        # compress?
        &trible_run_for_success( "gzip $ctChrFileBinIdx", "$chr-binIdx" ) if( $V_Href->{gzip_output} );
        # update
        print CHRBR "$chr\t$preChrBinSum\t";
        my $thisChrBinCount = int( $chrlen / $V_Href->{bin_size} );
        $thisChrBinCount++ if( $chrlen % $V_Href->{bin_size} != 0 );
        $preChrBinSum += $thisChrBinCount;
        print CHRBR ($preChrBinSum-1)."\n";
        # inform
        stout_and_sterr "[INFO]\t".`date`
                         ."\t$chr finishes, its BinCount is $thisChrBinCount, preChrBinSum updates to $preChrBinSum\n";
    }

    close CHRBR;
}

#--- convert interChr and intraChr result
sub merge_inter_intra{
    $V_Href->{MergedFile} = catfile( $V_Href->{workspace}, $V_Href->{MergedFile} );
    # compress?
    my @chr_Href = sort { $a->{turn} <=> $b->{turn} } values %{$V_Href->{ChrThings}};
    if( $V_Href->{gzip_output} ){
        `cp $V_Href->{ALLtoALLFile}.gz $V_Href->{MergedFile}.gz`;
        `zcat $_->{BinCt_file}.gz | gzip -c >> $V_Href->{MergedFile}.gz` for @chr_Href;
    }
    else{
        `cp $V_Href->{ALLtoALLFile} $V_Href->{MergedFile}`;
        `cat $_->{BinCt_file} >> $V_Href->{MergedFile}` for @chr_Href;
    }
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tmerge inter-chr and intra-chr OK:\n"
                         ."\t$V_Href->{MergedFile}\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
