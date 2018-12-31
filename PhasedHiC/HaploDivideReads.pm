package HaploHiC::PhasedHiC::HaploDivideReads;

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use List::Util qw/ max /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist trible_run_for_success /;
use BioFuse::BioInfo::Objects::Bam_OB;
use HaploHiC::LoadOn;
use HaploHiC::GetPath qw/ GetPath /;
use HaploHiC::Extensions::JuicerDumpBP qw/ load_chr_Things /;
use HaploHiC::PhasedHiC::phasedMutWork qw/ load_phased_VCF release_phaseMut_OB /;
use HaploHiC::PhasedHiC::splitPairBam qw/ divide_pairBam /;
use HaploHiC::PhasedHiC::dEndUkHapConfirm qw/ dEndUK_get_HapLink /;
use HaploHiC::PhasedHiC::sEndSoloHapConfirm qw/ sEndhx_get_HapLink /;
use HaploHiC::PhasedHiC::mergeHaploReads qw/ merge_haplo_reads /;
use HaploHiC::PhasedHiC::dumpContacts qw/ dump_contacts /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              DivideHiCreadsToHaplotypes
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::HaploDivideReads';
#----- version --------
$VERSION = "0.15";
$DATE = '2018-12-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        DivideHiCreadsToHaplotypes
                        check_files
                        prepare
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} haplo_div <[Options]>

     Options:

       # assign Hi-C reads to paternal and maternal haplotypes via Phased-Het-Mutation (PhaHetMut).  #
       # allocate PhaHetMut-uncovered Hi-C PE reads (PhaUcPE) based on Local Hap-Divide-Ratio (LHDR) #

       # Inputs and Outputs #
        -phsvcf  [s]  phased sorted VCF file. <required>
        -outdir  [s]  folder to store results and temporary files. <required>
        -bam     [s]  paired SE bam files. <required>
                       Note: 1) PE Hi-C reads are always aligned seperatedly, e.g., BWA mem.
                             2) inputs as: -bam a_R1.bam,a_R2.bam -bam b_R1.bam,b_R2.bam
                             3) can applied for multiple runs: run-a, run-b, .., run-n

       # Tools and Database #
        -samt    [s]  SAMtools. <required>
        -db_dir  [s]  database folder made by 'juicer_db' function. <required>
        -ref_v   [s]  version of reference genome, see contents under 'db_dir'. <required>

       # General Options #
        -enzyme  [s]  enzyme type. <required>
        -hapct   [i]  N-ploid of the sample. [2]
                       Note: default 2 is for homo sapiens.
        -use_indel    use phased InDels. [disabled]
        -min_idd [i]  minimum distance to accept adjacent indels. [5]
        -fork    [i]  to run operations with N forks in parallel. [1]
                       Note: more forks consume more memory.
        -st_step [i]  start wrokflow from certain step. [1]
        -ed_step [i]  stop workflow at certain step. [$V_Href->{totalStepNO}]
                       Note: step-NO. list:  1: splitBam;   2: sEndUKtoHap;  3: dEndUKtoHap; 
                                             4: readsMerge; 5: dumpContacts

       # Options of step NO.1 #
        -use_sp       use PE-reads having supplementary alignments, suggested to enable. [disabled]
        -use_sd       use PE-reads having secondary alignments. [disabled]
        -use_mm       use multiple mapped reads. [disabled]
        -min_mq  [i]  the minimum mapping quality. [20]
        -min_re  [i]  the minimum distance to bilateral reads-edges to judge one het-allele. [5]
        -min_bq  [i]  the minimum base quality to accept mut-allele on one reads. [20]
        -min_ad  [i]  the minimum distance to nearby alteration on reads to accept one het-allele. [3]
        -qual    [i]  base quality offset. [33]
                       Note: 1) set 33 for Sanger and Illumina 1.8+.
                             2) set 64 for Solexa, Illumina 1.3+ and 1.5+.
        -max_rd  [s]  the maximum distance to determine close alignments. [1E3]
        -use_caln     accept PE-reads whose two ends are close alignments ('-max_rd'). [disabled]
                       Note: default to put them in 'invalid PE'.

       # Options of step NO.2-3 #
        -skipddp      skip de-duplication of phased reads. [disabled]
                       ** Note: this option ('-skipddp') is also effective in step NO.5.
        -ucrfut  [s]  the flanking region to calculate LHDR for PhaUcPE. [1E4]
                       Note: 1) the minimum value allowed is 5E3.
                             2) it's unilateral size, and apply it bilaterally.
        -ucrpha  [i]  at most use such amount of flanking PhaHetMut to calculate LHDR for PhaUcPE. [5]
                       Note: 1) set 0 to disable, i.e., only '-ucrfut' works.
        -mpwr    [f]  ratio of '-ucrfut' to set as window size to store phased contacts. [0.1]
                       Note: 1) the less this option is set, the more memory and cpu-time is consumed.
                             2) available interval: (0, 0.5].
        -slbmpf  [s]  only deal with bam files with provided prefix. [all]
                       Note: 1) still load haplotypes contacts got from all bam files.
                             2) prefix's regex: /^(PREFIX)[\\._]R[12].*\\.bam\$/
                             3) effective in only run at step NO.2 or NO.3.

       # Options of step NO.5 #
        -dpmode  [s]  mode of dump, 'BP' or 'FRAG'. [BP]
        -dpbin   [s]  bin size of contacts dump. [1MB]
                       Note: 1) BP   mode allows: 2.5MB, 1MB, 500KB, 250KB, 100KB, 50KB, 25KB, 10KB, 5KB
                             2) FRAG mode allows: 500, 200, 100, 50, 20, 5, 2, 1
        -dphpcmb [s]  select any haplo combination matched string set by this option. default to use all.
                       Note: e.g., use 'Intra' to get intra-haplotype contacts of all chromosomes.
                                   use 'Inter' to get inter-haplotype contacts of all chromosomes.
                                   use 'h2Intra' to get only intra-h2 contacts of all chromosomes.

        -h|help       Display this help info.

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
        # -sampid  [s]  sample ID. <required>
        # -use_sgum    use pe-reads having one unmapped end. [disabled]
        # -ucrftm  [s]  allow to extend the '-ucrfut' time by time till to this times. [10]
        #                Note: try next extend only when current flanking region has no phased contacts.
                             # 2) option '-ucrftm' also works on this setting.
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
            # [ sampleID => undef ],
            [ phased_vcf => undef ],
            [ outdir => undef ],
            [ PairBamList => [] ],

            # software and database
            [ samtools => undef ],
            [ db_dir => undef ],
            [ GenomeRefFa => undef ],
            [ enzyme_type => undef ],
            [ enzyme_site => undef ],

            # options
            ## general
            [ haploCount => 2 ],
            [ forkNum => 1 ],
            # 1: splitBam;   2: sEndUKtoHap; 3: dEndUKtoHap;
            # 4: readsMerge; 5: dumpContacts
            [ stepToStart => 1 ],
            [ stepToStop => 5 ],
            [ totalStepNO => 5 ],
            [ SelectBamPref => [] ], # step=2, change to Href soon
            ## phased-Mut
            [ phaMutPosWinSize => 1E3 ],
            [ use_InDel => 0 ],
            [ min_InDelDist => 5 ],
            [ closeInDelSameHap => 0 ], # close-InDel must be on the same haplo (0: N, 1: Y)
            ## reads selection
            [ min_mapQ => 20 ],
            [ use_spmap => 0 ],
            [ use_sdmap => 0 ],
            [ use_multmap => 0 ],
            [ maxCloseAlignDist => 1E3 ],
            [ skipCloseAlignFromInvalid => 0 ],
            # [ use_single_ump => 0 ], # might use supplementary alignments in future
            ## allele on reads
            [ min_baseQ => 20 ],
            [ min_distToRedge => 5 ],
            [ min_distToNearbyAlt => 3 ],
            [ baseQ_offset => 33 ],
            ## unknown reads operation
            [ max_readLength => 150 ],
            [ mapPosWinRatio => 0.1 ],
            [ mapPosWinSize => undef ], # reset as 'UKreadsFlankRegUnit' * 'mapPosWinRatio'
            [ phasePEdetails => {} ], # record PE-map info and for de-dup
            [ phasePEcontact => {} ], # just record counts from 'phasePEdetails' hash
            [ UKreadsFlankRegUnit => 1E4 ],
            [ UKreadsFlankRegUnitMaxTimes => 10 ], # 2**time
            [ UKreadsMaxPhasedHetMut => 5 ],
            [ SkipDeDupPhasedReads => 0 ],
            ## dump contacts
            [ dumpMode => 'BP' ],
            [ dumpBinSize => '1MB' ],
            [ dumpPEdetails => {} ], # record PE-map info and for de-dup, similar to 'phasePEdetails'
            [ dumpPEcontact => {} ], # just record counts from 'dumpPEdetails' hash, similar to 'phasePEcontact'
            [ dumpSubDir => 'dumpContacts' ],
            [ dumpOutput => undef ], # output
            [ dumpBinLog => undef ], # output
            [ dumpHapComb => undef ], # select HapComb (h[x]Intra and hInter) to dump contacts
            ### FRAG specific
            [ chr2enzymePos => {} ],

            # intermediate variants
            [ allHapComb => [] ],
            [ intraHapComb => [] ],
            [ chrLenFile => undef ],
            [ phasedMut_report => undef ],
            [ VCFmutStat => {
                                i01_skipChr => 0,
                                i02_unPhaMut => 0,
                                i03_homPhaMut => 0,
                                i04_skipINDEL => 0,
                                i05_iniPhaMut => 0, # will delete
                                i05_hetPhaMut => 0,
                                i07_INDEL_rEdge => 0,
                                i08_INDEL_close => 0
                            } ],
            # [ { R1_bam => $R1_bamOB,
            #     R2_bam => $R2_bamOB,
            #     prefix => $bam_prefix,
            #     no => $no,
            #     workSpace => $workspace_folder,
            #     PEsplit_report => $report_path,
            #     splitBam => {tag1=>$splitBamOB, tag2=>$splitBamOB, ..},
            #     bamToMerge => {h[x]Intra=>[$splitBamOB,..], hInter=>[$splitBamOB,..]},
            #     mergeBam => {h[x]Intra=>$mergeBamOB, hInter=>$mergeBamOB}
            #   } ]
            ## tags: phMut-[ds]End-hx, phMut-[ds]End-hInter, discarded, unknown, invalid
            [ PairBamFiles => [] ],
            [ rOB_AbufferSize => 1E5 ], # 5E5
            [ chrPairSort_peOB_AbufferSize => 2E6 ], # 5E6
            [ peC_ReportUnit => 1E6 ],
            [ PEsplit_report => undef ],
            [ PEsplitStat => {
                                # 'phMut-dEnd-h1' => 0, # see func 'prepare'
                                # 'phMut-sEnd-h1' => 0,
                                'phMut-dEnd-hInter' => 0,
                                'phMut-sEnd-hInter' => 0,
                                'discarded' => 0,
                                'remove_SP' => 0,
                                'unknown' => 0,
                                'invalid-DanglingEnd' => 0,
                                'invalid-SelfCircle' => 0,
                                'invalid-DumpedPairFw' => 0,
                                'invalid-DumpedPairRv' => 0
                            } ],
            # {$chr}->{$PosIdx} = [ $PhasedMut_OB_1, $PhasedMut_OB_1, .. ]
            # check BioFuse::BioInfo::Objects::PhasedMut_OB
            [ PhasedMut => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['outdir'],
                                  ['db_dir'],
                                  ['samtools'],
                                  ['phased_vcf']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        # "-sampid:s" => \$V_Href->{sampleID},
        "-phsvcf:s" => \$V_Href->{phased_vcf},
        "-outdir:s" => \$V_Href->{outdir},
        "-bam:s"    => \@{$V_Href->{PairBamList}},
        # tools and datbase
        "-samt:s"   => \$V_Href->{samtools},
        "-db_dir:s" => \$V_Href->{db_dir},
        "-ref_v:s"  => \$V_Href->{ref_version},
        # options
        ## general
        "-enzyme:s" => \$V_Href->{enzyme_type},
        "-hapct:i"  => \$V_Href->{haploCount},
        "-fork:i"   => \$V_Href->{forkNum},
        "-st_step:i"=> \$V_Href->{stepToStart},
        "-ed_step:i"=> \$V_Href->{stepToStop},
        "-slbmpf:s" => \@{$V_Href->{SelectBamPref}},
        "-rbfsize:s"=> \$V_Href->{rOB_AbufferSize}, # hidden option
        ## phased-Mut
        "-use_indel"=> \$V_Href->{use_InDel},
        "-min_idd:i"=> \$V_Href->{min_InDelDist},
        ## reads selection
        "-min_mq:i" => \$V_Href->{min_mapQ},
        "-use_sp"   => \$V_Href->{use_spmap},
        "-use_sd"   => \$V_Href->{use_sdmap},
        "-use_mm"   => \$V_Href->{use_multmap},
        "-use_caln" => \$V_Href->{skipCloseAlignFromInvalid},
        # "-use_sgum"   => \$V_Href->{use_single_ump},
        "-max_rd:s" => \$V_Href->{maxCloseAlignDist},
        ## allele on reads
        "-min_bq:i" => \$V_Href->{min_baseQ},
        "-min_re:i" => \$V_Href->{min_distToRedge},
        "-min_ad:i" => \$V_Href->{min_distToNearbyAlt},
        "-qual:i"   => \$V_Href->{baseQ_offset},
        ## unknown reads operation
        "-ucrfut:s" => \$V_Href->{UKreadsFlankRegUnit},
        # "-ucrftm:i" => \$V_Href->{UKreadsFlankRegUnitMaxTimes},
        "-ucrpha:i" => \$V_Href->{UKreadsMaxPhasedHetMut},
        "-mpwr:f"   => \$V_Href->{mapPosWinRatio},
        "-skipddp"  => \$V_Href->{SkipDeDupPhasedReads},
        ## dump contacts
        "-dpmode:s" => \$V_Href->{dumpMode},
        "-dpbin:s"  => \$V_Href->{dumpBinSize},
        "-dphpcmb:s"=> \$V_Href->{dumpHapComb},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{phased_vcf})
             || !defined $V_Href->{outdir}   || !-d $V_Href->{outdir}
             || !defined $V_Href->{samtools} || !-e $V_Href->{samtools}
             || !defined $V_Href->{db_dir}   || !-d $V_Href->{db_dir}
             || !defined $V_Href->{ref_version}
             || scalar(@{$V_Href->{PairBamList}}) == 0
             || !defined $V_Href->{enzyme_type}
             || $V_Href->{haploCount} < 2
             || ( $V_Href->{baseQ_offset} != 33 && $V_Href->{baseQ_offset} != 64 )
             || $V_Href->{UKreadsFlankRegUnit} < 5E3
             # || $V_Href->{UKreadsFlankRegUnitMaxTimes} < 1
             || $V_Href->{UKreadsMaxPhasedHetMut} < 0
             || $V_Href->{mapPosWinRatio} <= 0 || $V_Href->{mapPosWinRatio} > 0.5
             || $V_Href->{stepToStart} < 1 || $V_Href->{stepToStart} > $V_Href->{totalStepNO}
             || $V_Href->{stepToStop}  < 1 || $V_Href->{stepToStop}  > $V_Href->{totalStepNO}
             || $V_Href->{stepToStop}  < $V_Href->{stepToStart}
             || (    scalar(@{$V_Href->{SelectBamPref}}) != 0
                  && !(    $V_Href->{stepToStart} == $V_Href->{stepToStop} # run at only one step
                        && $V_Href->{stepToStart} =~ /^[23]$/ # run at step NO.2 or NO.3
                      )
                )
             || $V_Href->{dumpMode} !~ /^(BP|FRAG)$/
             || !exists $V_Href->{dump_allowBinSize}->{$V_Href->{dumpMode}}->{uc($V_Href->{dumpBinSize})}
            );
}

#--- run juicer_tools dump func in FRAG mode and merge all chromosomes ---
sub DivideHiCreadsToHaplotypes{

    # basic
    &check_files;
    &prepare;
    load_chr_Things;

    # phased Het-mutation
    load_phased_VCF;

    # PE to categories
    # [ds]End-h[x], [ds]End-hInter, unkown, discard
    divide_pairBam;

    # release memory
    release_phaseMut_OB;

    # one-side confirmed contacts
    ## sEnd-h[x]
    sEndhx_get_HapLink;

    # non-side confirmed contacts
    ## unkown
    dEndUK_get_HapLink;

    # merge reads of each haplotype
    merge_haplo_reads;

    # dump contacts
    dump_contacts;
}

#--- check files ---
sub check_files{
    # chr len file
    $V_Href->{chrLenFile}  = GetPath(filekey => 'chrLenFile');
    # ref fasta file
    $V_Href->{GenomeRefFa} = GetPath(filekey => 'GenomeRefFa');
    # enzyme site list
    $V_Href->{enzyme_site} = GetPath(filekey => 'enzyme_site');
    # check existence
    if(    !file_exist(filePath => $V_Href->{chrLenFile})
        || !file_exist(filePath => $V_Href->{GenomeRefFa})
        || !file_exist(filePath => $V_Href->{enzyme_site})
    ){
        warn_and_exit "<ERROR>\tPlease make sure the existence of files below.\n"
                            ."\t$V_Href->{chrLenFile}\n"
                            ."\t$V_Href->{GenomeRefFa}\n"
                            ."\t$V_Href->{enzyme_site}\n";
    }

    # PairBamFiles
    my %preBamFiles;
    for my $PairBamList (@{$V_Href->{PairBamList}}){
        if ($PairBamList =~ /^([^,]+R1[^,]*\.bam),([^,]+R2[^,]*\.bam)$/){
            my ($R1_bam_path, $R2_bam_path) = ($1, $2);
            if(    !file_exist(filePath=>$R1_bam_path)
                || !file_exist(filePath=>$R2_bam_path)
            ){
                warn_and_exit "<ERROR>\tbam file does not existfrom input:\n"
                                    ."\t$PairBamList\n";
            }
            if(    exists $preBamFiles{$R1_bam_path}
                || exists $preBamFiles{$R2_bam_path}
            ){
                warn_and_exit "<ERROR>\tencounter same bam file again from input:\n"
                                    ."\t$PairBamList\n";
            }
            # prepare output prefix
            ( my $R1_bam_prefix = basename($R1_bam_path) ) =~ s/[\.\_]R1.*\.bam$//;
            ( my $R2_bam_prefix = basename($R2_bam_path) ) =~ s/[\.\_]R2.*\.bam$//;
            if(    $R1_bam_prefix ne $R2_bam_prefix
                || length($R1_bam_prefix) == 0
                || exists $preBamFiles{$R1_bam_prefix}
            ){
                warn_and_exit "<ERROR>\trequires valid prefix of bam files from input:\n"
                                    ."\t$PairBamList\n";
            }
            # bam object
            my $R1_bam = BioFuse::BioInfo::Objects::Bam_OB->new(filepath => $R1_bam_path, tag => $R1_bam_prefix);
            my $R2_bam = BioFuse::BioInfo::Objects::Bam_OB->new(filepath => $R2_bam_path, tag => $R2_bam_prefix);
            # records
            my $i = scalar @{$V_Href->{PairBamFiles}};
            push @{$V_Href->{PairBamFiles}}, {R1_bam => $R1_bam, R2_bam => $R2_bam, prefix => $R1_bam_prefix, no => $i};
            $preBamFiles{$_} = 1 for ( $R1_bam_path, $R2_bam_path, $R1_bam_prefix );
        }
        else{
            warn_and_exit "<ERROR>\tcannot recognize R1/R2.bam files from input:\n"
                                ."\t$PairBamList\n";
        }
    }

    # selected bam to process in step NO.2 or NO.3
    my %temp;
    for my $bam_prefix ( @{$V_Href->{SelectBamPref}} ){
        unless( exists $preBamFiles{$bam_prefix} ){
            warn_and_exit "<ERROR>\tcannot recognize selected bam-prefix ('-slbmpf' option):\n"
                                ."\t$bam_prefix\n";
        }
        $temp{$bam_prefix} = 1;
    }
    $V_Href->{SelectBamPref} = \%temp;

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tfiles checking DONE\n";
}

#--- prepare some workflow settings ---
sub prepare{
    # window size to record phased contacts
    if(    $V_Href->{stepToStart} <= 3
        && $V_Href->{stepToStop}  >= 2
    ){
        # set mapPosWinSize
        $V_Href->{mapPosWinSize} = int(max($V_Href->{UKreadsFlankRegUnit} * $V_Href->{mapPosWinRatio} / 100, 1)) * 100;
        ## check
        if($V_Href->{mapPosWinSize} < 200){
            warn_and_exit "<ERROR>\tcurrent window size ($V_Href->{mapPosWinSize}) to record phased contacts is less than 200.\n"
                                ."\tplease set option '-mpwr' as larger value.\n";
        }
        ## inform
        stout_and_sterr "[INFO]\t".`date`
                             ."\twindow size to record phased contacts is set as $V_Href->{mapPosWinSize}.\n";
    }

    # prepare haplotype combination: hx,hx
    for my $hap_i (1 .. $V_Href->{haploCount}){
        push @{$V_Href->{intraHapComb}}, "h$hap_i,h$hap_i";
        for my $hap_j (1 .. $V_Href->{haploCount}){
            push @{$V_Href->{allHapComb}}, "h$hap_i,h$hap_j";
        }
    }

    # prepare PEsplitStat (dEnd and sEnd) of step s01
    for my $hap_i (1 .. $V_Href->{haploCount}){
        $V_Href->{PEsplitStat}->{"phMut-dEnd-h$hap_i"} = 0;
        $V_Href->{PEsplitStat}->{"phMut-sEnd-h$hap_i"} = 0;
    }

    # inform selected bam
    if(    (   $V_Href->{stepToStart} == 2
            || $V_Href->{stepToStop}  == 3
           )
        && scalar(keys %{$V_Href->{SelectBamPref}}) != 0
    ){
        stout_and_sterr "[INFO]\t".`date`
                             ."\tSelect bam files with below prefixes to operate.\n";
        stout_and_sterr       "\t$_\n" for sort keys %{$V_Href->{SelectBamPref}};
    }

    # delete possible previous results
    if( $V_Href->{stepToStart} <= 1 ){
        `rm -rf $V_Href->{outdir}/*`;
    }
    else{
        if( $V_Href->{stepToStart} <= 2 ){
            if( scalar(keys %{$V_Href->{SelectBamPref}}) == 0 ){
                `rm -rf $V_Href->{outdir}/*-workspace/*phMut-sEnd-h[0-9]*.h*Int*.bam`;
            }
            else{
                `rm -rf $V_Href->{outdir}/*-workspace/$_.phMut-sEnd-h[0-9]*.h*Int*.bam` for keys %{$V_Href->{SelectBamPref}};
            }
        }
        if( $V_Href->{stepToStart} <= 3 ){
            if(    scalar(keys %{$V_Href->{SelectBamPref}}) == 0
                || $V_Href->{stepToStart} != 3
            ){
                `rm -rf $V_Href->{outdir}/*-workspace/*unknown.h*Int*.bam`;
            }
            else{
                `rm -rf $V_Href->{outdir}/*-workspace/$_.unknown.h*Int*.bam` for keys %{$V_Href->{SelectBamPref}};
            }
        }
        if( $V_Href->{stepToStart} <= 4 ){
            `rm -rf $V_Href->{outdir}/*.merge.h*.bam`;
        }
        # if only step NO.5,
        # it might use different setting,
        # or else will overwirte former results.
        if( $V_Href->{stepToStart} < 5 ){
            `rm -rf $V_Href->{outdir}/$V_Href->{dumpSubDir}`;
        }
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
