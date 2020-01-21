package HaploHiC::Extensions::JuicerDump;

use strict;
use warnings;
use Getopt::Long;
use File::Spec::Functions qw/ catfile /;
use BioFuse::Util::Sys qw/ file_exist trible_run_for_success /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use HaploHiC::LoadOn;
use HaploHiC::GetPath qw/ GetPath /;
use HaploHiC::PhasedHiC::dumpContacts qw/ get_dumpHeader load_enzyme_site_list write_dumpBinLog get_contacts_idx /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
                JuicerDumpWrok
                check_juicer_files
                make_workspace
                load_chr_Things
                getPairChrDump
                mergePairChrDump
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::Extensions::JuicerDump';
#----- version --------
$VERSION = "0.11";
$DATE = '2019-03-23';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        JuicerDumpWrok
                        check_juicer_files
                        make_workspace
                        load_chr_Things
                        getPairChrDump
                        mergePairChrDump
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} juicerDump <[Options]>
     
     Options:

       # run juicer_tools dump func and merge all chromosomes to one final file #

       # Inputs and Outputs #
        -juicer [s]  path of juicer folder, which has sub-folder, e.g., 'CPU'. <required>
        -topdir [s]  topdir stores juicer results. <required>
                        this func creates '$V_Href->{workspace_name}-<postfix>' folder under topdir.
        -db_dir [s]  database folder made by 'juicer_db' function. <required>
        -ref_v  [s]  version of reference genome, see contents under 'db_dir'. <required>

       # Options #
        -hic_fl [s]  prefix of '.hic' file: 'inter_30' or 'inter'. [inter_30]
        -dpmode [s]  mode of dump, 'BP' or 'FRAG'. [BP]
        -enzyme [s]  enzyme type. <required in FRAG mode>
        -bin    [s]  bin size. [1MB]
                       1)   BP mode: 2.5MB, 1MB, 500KB, 250KB, 100KB, 50KB, 25KB, 10KB, 5KB
                       2) FRAG mode: 500, 200, 100, 50, 20, 5, 2, 1
        -ctype  [s]  'observed' or 'oe'. [observed]
        -norm   [s]  kind of normalization, must be one of NONE/VC/VC_SQRT/KR. ['NONE']
                       Notes from juicer_tools dump usage:
                        1) VC is vanilla coverage,
                        2) VC_SQRT is square root of vanilla coverage;
                        3) KR is Knight-Ruiz or Balanced normalization.

        -h|help      Display this help info.

     Version:
        $VERSION on $DATE

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
            [ count_type => 'observed' ],
            [ norm_method => 'NONE' ],
            [ dumpMode => 'BP' ],
            [ dumpBinSize => '1MB' ],
            [ enzyme_type => undef ],

            # output
            [ dumpHeader => undef ],
            [ dumpBinLog => 'chrBinRange.txt' ],
            [ dumpMerged => 'InterIntraChr.merged.txt.gz' ],

            # intermediate variants
            [ hic_file => undef ],
            [ chrLenFile => undef ],
            [ juicer_tool_Jar => undef ],
            [ enzyme_site => undef ],
            [ chr2enzymePos => {} ],
            [ workspace_name => 'extract_matrix' ],
            [ workspace => undef ],

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
        "-dpmode:s" => \$V_Href->{dumpMode},
        "-bin:s"    => \$V_Href->{dumpBinSize},
        "-ctype:s"  => \$V_Href->{count_type},
        "-norm:s"   => \$V_Href->{norm_method},
        "-enzyme:s" => \$V_Href->{enzyme_type},
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
             || !exists $V_Href->{dump_allowBinSize}->{$V_Href->{dumpMode}}->{$V_Href->{dumpBinSize}}
             || ( $V_Href->{dumpMode} eq 'FRAG' && !defined $V_Href->{enzyme_type} )
             || !exists $V_Href->{dump_allowNormMtd}->{$V_Href->{norm_method}}
             || ( $V_Href->{hic_type} ne 'inter' && $V_Href->{hic_type} ne 'inter_30' )
             || ( $V_Href->{count_type} ne 'observed' && $V_Href->{count_type} ne 'oe' )
            );
}

#--- run juicer_tools dump func and merge all chromosomes ---
sub JuicerDumpWrok{
    &check_juicer_files;

    &make_workspace;

    &load_chr_Things;

    &getPairChrDump;

    get_dumpHeader;
    load_enzyme_site_list if($V_Href->{dumpMode} eq 'FRAG');
    write_dumpBinLog;

    &mergePairChrDump;
}

#--- check files ---
## also applied by JuicerDumpFRAG.pm
sub check_juicer_files{
    $V_Href->{hic_file} = GetPath(filekey => 'hic_file');
    $V_Href->{chrLenFile} = GetPath(filekey => 'chrLenFile');
    $V_Href->{juicer_tool_Jar} = GetPath(filekey => 'juicer_tool_Jar');
    $V_Href->{enzyme_site} = GetPath(filekey => 'enzyme_site') if $V_Href->{dumpMode} eq 'FRAG';
    # check existence
    if(    !file_exist(filePath => $V_Href->{hic_file})
        || !file_exist(filePath => $V_Href->{chrLenFile})
        || !file_exist(filePath => $V_Href->{juicer_tool_Jar})
        || ($V_Href->{dumpMode} eq 'FRAG' && !file_exist(filePath => $V_Href->{enzyme_site}))
    ){
        warn "Please make sure the existence of files below.\n";
        warn "  $V_Href->{hic_file}\n";
        warn "  $V_Href->{chrLenFile}\n";
        warn "  $V_Href->{juicer_tool_Jar}\n";
        warn "  $V_Href->{enzyme_site}\n" if $V_Href->{dumpMode} eq 'FRAG';
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
    $V_Href->{workspace_name} = join('-', $V_Href->{workspace_name}, $V_Href->{hic_type}, $V_Href->{count_type}, $V_Href->{norm_method}, $V_Href->{dumpMode}.$V_Href->{dumpBinSize} );
    $V_Href->{workspace} = catfile( $V_Href->{JuicerTopdir}, 'aligned', $V_Href->{workspace_name} );
    `rm -rf $V_Href->{workspace} && mkdir -p $V_Href->{workspace}`;

    # dump outputs
    $V_Href->{dumpBinLog} = catfile( $V_Href->{workspace}, $V_Href->{dumpBinLog} );
    $V_Href->{dumpMerged} = catfile( $V_Href->{workspace}, $V_Href->{dumpMerged} );

    # use real number of dumpBinSize
    $V_Href->{dumpBinSize} = $V_Href->{dump_allowBinSize}->{$V_Href->{dumpMode}}->{ $V_Href->{dumpBinSize} };

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tcreate workspace OK\n"
                         ."\t$V_Href->{workspace}\n";
}

#--- read chromosome info ---
## also applied in haplo_div
sub load_chr_Things{
    # the chrLenFile chr-order is very important, same with refseg order in reference genome
    # juicer applies this order in its outputs
    open (CHRLEN, $V_Href->{chrLenFile}) || die "fail read chrLenFile: $!\n";
    while(<CHRLEN>){
        next if(/^#/);
        my ($chr, $length) = (split)[0, 1];
        # avoid comma in chr-name
        if($chr =~ /,/){
            warn_and_exit "<ERROR>\tchromosome $chr has comma (',') in its name, not allowed.\n";
        }
        # record
        $V_Href->{ChrThings}->{$chr} = { chr=>$chr, len=>$length, turn=>$. };
        push @{$V_Href->{sortedChr}}, $chr;
    }
    close CHRLEN;
}

#--- get contact of paired chr ---
sub getPairChrDump{
    # juicer_tools.jar func 'dump' follows ref-seg order, i.e., the '@SQ' in SAM header
    # no matter what i_chr,j_chr or j_chr,i_chr it receive, always follow the ref-seg order
    my @chr_Href = sort { $a->{turn} <=> $b->{turn} } values %{$V_Href->{ChrThings}};

    for my $i ( 0 .. $#chr_Href ){
        my $i_chr = $chr_Href[$i]->{chr};
        my $i_chr_workspace = catfile($V_Href->{workspace}, $i_chr);
        `rm -rf $i_chr_workspace && mkdir -p $i_chr_workspace`;
        for my $j ( $i .. $#chr_Href ){
            my $j_chr = $chr_Href[$j]->{chr};
            my $ij_dump_file = catfile($i_chr_workspace, "$i_chr.vs.$j_chr.$V_Href->{dumpMode}.txt");
            # dump outputs from juicer_tools.jar func 'dump'
            my $command  = "($V_Href->{java} -jar $V_Href->{juicer_tool_Jar} dump $V_Href->{count_type} $V_Href->{norm_method} $V_Href->{hic_file} $i_chr $j_chr $V_Href->{dumpMode} $V_Href->{dumpBinSize} $ij_dump_file.tmp)";
            &trible_run_for_success( $command, "$i_chr-vs-$j_chr-dump-$V_Href->{dumpMode}" );
            # add chr
            open (NDUMP, Try_GZ_Write("$ij_dump_file.gz")) || die "fail write $i_chr.vs.$j_chr.$V_Href->{dumpMode}.txt.gz: $!\n";
            # header
            if($V_Href->{dumpMode} eq 'BP'){
                print NDUMP join("\t", '#chr_i', 'chrBinStPos_i', 'chr_j', 'chrBinStPos_j', "contacts\n");
            }
            else{
                print NDUMP join("\t", '#chr_i', 'chrFragIdx_i',  'chr_j', 'chrFragIdx_j',  "contacts\n");
            }
            open (ODUMP, Try_GZ_Read("$ij_dump_file.tmp")) || die "fail read $i_chr.vs.$j_chr.$V_Href->{dumpMode}.txt.tmp: $!\n";
            while(<ODUMP>){
                my ($i_chrBin, $j_chrBin, $contact) = (split);
                print NDUMP join("\t", $i_chr, $i_chrBin, $j_chr, $j_chrBin, $contact) . "\n";
            }
            close ODUMP;
            close NDUMP;
            # delete tmp
            `rm -f $ij_dump_file.tmp`;
            # record
            push @{$chr_Href[$i]->{dump}}, [$j_chr, "$ij_dump_file.gz"];
        }
        # inform
        stout_and_sterr "[INFO]\t".`date`
                             ."\tdump($V_Href->{dumpMode}) of $i_chr OK\n";
    }
}

#--- merge pairChr dump results and add whole-genome bin-index ---
sub mergePairChrDump{
    open (CNTC, Try_GZ_Write($V_Href->{dumpMerged})) || die "cannot write merged dump file: $!\n";
    # header
    print CNTC $V_Href->{dumpHeader};
    print CNTC join("\t", '#chr_i', 'chrBinIdx_i', "wgBinIdx_i", 'chr_j', 'chrBinIdx_j', "wgBinIdx_j", "contacts\n");
    # load dump file in refseg order
    my @chr_Href = sort { $a->{turn} <=> $b->{turn} } values %{$V_Href->{ChrThings}};
    for my $chr_Href ( @chr_Href ){
        my $i_chr = $chr_Href->{chr};
        my $chr_i_preChrBinSum = $chr_Href->{preChrBinSum};
        next unless exists $chr_Href->{dump};
        for my $Aref (@{$chr_Href->{dump}}){
            my ($j_chr, $ij_dump_file) = @$Aref;
            my $chr_j_preChrBinSum = $V_Href->{ChrThings}->{$j_chr}->{preChrBinSum};
            open (IJDUMP, Try_GZ_Read($ij_dump_file)) || die "fail read $i_chr.vs.$j_chr.$V_Href->{dumpMode}.txt.gz: $!\n";
            while(<IJDUMP>){
                next if /^#/;
                my ($i_chrBin, $j_chrBin, $contact) = (split)[1,3,4];
                # if it's BP mode, the x_chrBin is the bin-StPos, transform to binIdx
                if($V_Href->{dumpMode} eq 'BP'){
                    $i_chrBin = &get_contacts_idx(chr => $i_chr, pos => $i_chrBin);
                    $j_chrBin = &get_contacts_idx(chr => $j_chr, pos => $j_chrBin);
                }
                # output
                print CNTC join("\t", $i_chr,
                                      $i_chrBin,
                                      $i_chrBin + $chr_i_preChrBinSum,
                                      $j_chr,
                                      $j_chrBin,
                                      $j_chrBin + $chr_j_preChrBinSum,
                                      $contact)."\n";
            }
            close IJDUMP;
        }
        # inform
        stout_and_sterr "[INFO]\t".`date`
                             ."\tmerge $i_chr related chrPairs dump OK.\n";
    }
    close CNTC;
}

#--- 
1; ## tell the perl script the successful access of this module.
