package HaploHiC::PhasedHiC::HaploDivideReads;

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use List::Util qw/ max /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist trible_run_for_success /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use HaploHiC::LoadOn;
use HaploHiC::GetPath qw/ GetPath /;
use HaploHiC::Extensions::JuicerDump qw/ load_chr_Things /;
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
$VERSION = "0.34";
$DATE = '2021-08-14';

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
                        delete_prev_results
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} haplo_div <[Options]>

     Options:

       # assign Hi-C reads to paternal and maternal haplotypes via Phased-Het-Mutation (PhaHetMut).  #
       # allocate haplotype-unknown Hi-C PE reads (HaploUKPE) based on Local Hap-Divide-Ratio (LHDR) #

       # Inputs and Outputs #
        -phsvcf  [s]  phased sorted VCF file. <required>
        -outdir  [s]  folder to store results and temporary files. <required>
        -bamlist [s]  file path list of paired SE bam or sortN PE bam. <required>
                       Note: 1) PE Hi-C reads are always aligned seperatedly.
                             2) list multiple runs, one row one run.
                             3) instance run_c: run_c_R1.bam,run_c_R2.bam
                             4) instance run_k: run_k.sortN.bam

       # Tools and Database #
        -samt    [s]  SAMtools. <required>
        -db_dir  [s]  database folder made by 'juicer_db' function. <required>
        -ref_v   [s]  version of reference genome, see contents under 'db_dir'. <required>

       # General Options #
        -enzyme  [s]  enzyme type. <required>
        -use_indel    use phased InDels. [disabled]
        -min_idd [i]  minimum distance to accept adjacent indels. [5]
        -st_step [i]  start wrokflow from certain step. [1]
        -ed_step [i]  stop workflow at certain step. [$V_Href->{totalStepNO}]
                       Note: step-NO. list:  1: splitBam;   2: sEndUKtoHap;  3: dEndUKtoHap; 
                                             4: readsMerge; 5: dumpContacts
        -fork    [i]  to run operations with N forks in parallel. [5]
                       Note: 1) more forks, more memory consumed.
                             2) step NO.1/4: maximum is the number of Hi-C PE bam files input.
                                step NO.2/3: maximum is the chromosome count to deal.
                             3) if set a large number, it's better to run NO.1 step alone first.

       # Options of step NO.1 #
        -flt_pe       discard PE-reads once either read has unqualitified alignment. [disabled]
                       Note: default to discard unqualitified alignment only.
        -use_sp       use PE-reads having supplementary alignments, suggested to enable. [disabled]
        -use_sd       use PE-reads having secondary alignments. [disabled]
        -use_mm       use multiple mapped reads. [disabled]
        -min_mq  [i]  the minimum mapping quality. [30]
        -min_re  [i]  the minimum distance to bilateral reads-edges to judge one het-allele. [5]
        -min_bq  [i]  the minimum base quality to accept mut-allele on one reads. [20]
        -min_ad  [i]  the minimum distance to nearby alteration on reads to accept one het-allele. [3]
        -qual    [i]  base quality offset. [33]
                       Note: '33' for Sanger and Illumina 1.8+; '64' for Solexa, Illumina 1.3+ and 1.5+.
        -max_rd  [s]  the maximum distance to determine close alignments. [1E3]
        -use_caln     accept PE-reads whose two ends are close alignments ('-max_rd'). [disabled]
                       Note: default to put them in 'invalid PE'.

       # Options of step NO.2-3 #
        -ucrfut  [s]  the unit size of unilateral local flanking region, minimum: 5E3. [1E4]
        -ucrfmx  [s]  the maximum size of unilateral local flanking region when extends. [1E7]
        -mpwr    [f]  ratio of '-ucrfut' to set as window size to store phased contacts. (0, [0.5]]
                       Note: the less '-mpwr' is, the more memory and cpu-time consumed.
        -min_ct  [i]  the minimum contact count from phased pairs to confirm phased local regions. >=[5]

       # Options of step NO.5 #
        -dump    [s]  how to dump. [BP:1MB]
                       Note: 1) format is DumpMode[:BinSize:HapComb], e.g., BP, FRAG:1, BP:500KB:h2Intra
                             2) can be used multiple times, e.g., -dump BP:500KB -dump FRAG:1 -dump BP
                             3) DumpMode: BP or FRAG.
                             4) BP   BinSize: 2.5MB, [1MB], 500KB, 250KB, 100KB, 50KB, 25KB, 10KB, 5KB
                             5) FRAG BinSize: 500, 200, 100, 50, 20, 5, 2, [1]
                             6) HapComb is to select any haplo combination matched this string. default all.
                                e.g., use 'Intra' to get intra-haplotype contacts of all chromosomes.
                                      use 'Inter' to get inter-haplotype contacts of all chromosomes.
                                      use 'h2Intra' to get only intra-h2 contacts of all chromosomes.

        -h|help       Display this help info.

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
        # -hapct   [i]  N-ploid of the sample. [2]
        #                Note: default 2 is for homo sapiens.
        # -skipddp      skip de-duplication of phased reads. [disabled]
        #                ** Note: this option ('-skipddp') is also effective in step NO.5.
        # -ucrpha  [i]  at most use such amount of flanking PhaHetMut to calculate LHDR for HaploUKPE. [5]
        #                Note: 1) set 0 to disable, i.e., only '-ucrfut' works.
        # -slbmpf  [s]  only deal with bam files with provided prefix. [all]
        #                Note: 1) still load haplotypes contacts got from all bam files.
        #                      2) prefix's regex: /^(PREFIX)[\\._]R[12].*\\.bam\$/
        #                      3) effective in only run at step NO.2 or NO.3.
        # -sampid  [s]  sample ID. <required>
        # -use_sgum    use pe-reads having one unmapped end. [disabled]
        # -ucrftm  [s]  allow to extend the '-ucrfut' time by time till to this times. [10]
        #                Note: try next extend only when current flanking region has no phased contacts.
        #                      2) option '-ucrftm' also works on this setting.
        # -bam     [s]  paired SE bam files. <required>
        #                Note: 1) PE Hi-C reads are always aligned seperatedly, e.g., BWA mem.
        #                      2) inputs as: -bam a_R1.bam,a_R2.bam -bam b_R1.bam,b_R2.bam
        #                      3) can applied for multiple runs: run-a, run-b, .., run-n
        # -add_r   [f]  uniform addition ratio to each haplotype combination in phased local region. [0]<=0.1
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
            [ PairBamList => undef ],
            [ SourceBam => [] ],

            # software and database
            [ samtools => undef ],
            [ db_dir => undef ],
            [ GenomeRefFa => undef ],
            [ enzyme_type => undef ],
            [ enzyme_site => undef ],

            # options
            ## general
            [ haploCount => 2 ],
            [ forkNum => 5 ], # 1st,4th: NumOfPairedBam; 2nd,3rd: chrCount*2 (Intra+Inter).
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
            [ skipS01Report => 0 ], # skip over-writing when use '-slbmpf' in s01
            ## reads selection
            [ mustMapChr_i => {Af => [], Hf => {}} ],
            [ mustMapChr_j => {Af => [], Hf => {}} ],
            [ min_mapQ => 30 ],
            [ use_spmap => 0 ],
            [ use_sdmap => 0 ],
            [ use_multmap => 0 ],
            [ filter_pe => 0 ],
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
            [ mapPosWinRatio => 0.5 ],
            [ mapPosWinSize => undef ], # reset as 'UKreadsFlankRegUnit' * 'mapPosWinRatio'
            [ phasePEdetails => {} ], # record PE-map info and for de-dup
            [ phasePEcontact => {} ], # just record counts from 'phasePEdetails' hash
            [ UKreadsFlankRegUnit => 1E4 ],
            [ UKreadsFlankRegMax => 1E7 ], # 10mb
            [ UKreadsFlankRegUnitMaxTimes => undef ],
            [ UKreadsFlankRegUnitMidTimes => undef ],
            [ UKreadsFlankRegUnitIniTimes => undef ],
            # [ UKreadsFlankReg => [] ],
            # [ UKreadsMaxPhasedHetMut => 0 ], # deprecated
            [ SkipDeDupPhasedReads => 0 ],
            [ hapCombMinLinkForPhaReg => 5 ],
            [ uniformAddRatioForHapComb => 0 ],
            ## dump contacts
            [ dump => ['BP:1MB'] ], # integrated option, mode:size:hapComb
            [ dumpMode => undef ],
            [ dumpBinSize => undef ],
            [ dumpHapComb => undef ], # select HapComb (h[x]Intra and hInter) to dump contacts
            [ norm_method => 'NONE' ],
            [ dumpPEdetails => {} ], # record PE-map info and for de-dup, similar to 'phasePEdetails'
            [ dumpPEcontact => {} ], # just record counts from 'dumpPEdetails' hash, similar to 'phasePEcontact'
            [ dumpSubDir => 'dumpContacts' ],
            [ dumpHeader => undef ], # output
            [ dumpFilePrefix => undef ], # output prefix
            [ dumpOutput => undef ], # output
            [ dumpBinLog => undef ], # output
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
            # check BioFuse::BioInfo::Objects::Allele::PhasedMut_OB
            [ PhasedMut => {} ],
            # {flankSize}->{PhasedLinkCount} = $number
            [ LocRegPhased => {} ],

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
        "-bamlist:s"=> \$V_Href->{PairBamList},
        "-bam:s"    => \@{$V_Href->{SourceBam}}, # hidden option
        # tools and datbase
        "-samt:s"   => \$V_Href->{samtools},
        "-db_dir:s" => \$V_Href->{db_dir},
        "-ref_v:s"  => \$V_Href->{ref_version},
        # options
        ## general
        "-enzyme:s" => \$V_Href->{enzyme_type},
        "-hapct:i"  => \$V_Href->{haploCount}, # hidden option
        "-fork:i"   => \$V_Href->{forkNum},
        "-st_step:i"=> \$V_Href->{stepToStart},
        "-ed_step:i"=> \$V_Href->{stepToStop},
        "-slbmpf:s" => \@{$V_Href->{SelectBamPref}}, # hidden option
        "-rbfsize:s"=> \$V_Href->{rOB_AbufferSize}, # hidden option
        ## phased-Mut
        "-use_indel"=> \$V_Href->{use_InDel},
        "-min_idd:i"=> \$V_Href->{min_InDelDist},
        "-sks1rep"  => \$V_Href->{skipS01Report}, # hidden option
        ## reads selection
        "-chr_i:s"  => \@{$V_Href->{mustMapChr_i}->{Af}}, # hidden option
        "-chr_j:s"  => \@{$V_Href->{mustMapChr_j}->{Af}}, # hidden option
        "-min_mq:i" => \$V_Href->{min_mapQ},
        "-use_sp"   => \$V_Href->{use_spmap},
        "-use_sd"   => \$V_Href->{use_sdmap},
        "-use_mm"   => \$V_Href->{use_multmap},
        "-flt_pe"   => \$V_Href->{filter_pe},
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
        "-ucrfmx:s" => \$V_Href->{UKreadsFlankRegMax},
        # "-ucrpha:i" => \$V_Href->{UKreadsMaxPhasedHetMut}, # hidden option
        "-mpwr:f"   => \$V_Href->{mapPosWinRatio},
        "-skipddp"  => \$V_Href->{SkipDeDupPhasedReads}, # hidden option
        "-min_ct:i" => \$V_Href->{hapCombMinLinkForPhaReg},
        "-add_r:f"  => \$V_Href->{uniformAddRatioForHapComb}, # hidden option
        ## dump contacts
        "-dump:s"   => \@{$V_Href->{dump}},
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
             || ( !file_exist(filePath=>$V_Href->{PairBamList}) && scalar(@{$V_Href->{SourceBam}}) == 0 )
             || !defined $V_Href->{enzyme_type}
             || $V_Href->{haploCount} < 2
             || ( $V_Href->{baseQ_offset} != 33 && $V_Href->{baseQ_offset} != 64 )
             || $V_Href->{UKreadsFlankRegUnit} < 5E3
             || $V_Href->{UKreadsFlankRegMax} < $V_Href->{UKreadsFlankRegUnit}
             # || $V_Href->{UKreadsMaxPhasedHetMut} < 0
             || $V_Href->{hapCombMinLinkForPhaReg} < 5
             || $V_Href->{uniformAddRatioForHapComb} > 0.1
             || $V_Href->{mapPosWinRatio} <= 0 || $V_Href->{mapPosWinRatio} > 0.5
             || $V_Href->{stepToStart} < 1 || $V_Href->{stepToStart} > $V_Href->{totalStepNO}
             || $V_Href->{stepToStop}  < 1 || $V_Href->{stepToStop}  > $V_Href->{totalStepNO}
             || $V_Href->{stepToStop}  < $V_Href->{stepToStart}
             || (    scalar(@{$V_Href->{SelectBamPref}}) != 0
                  && !(    $V_Href->{stepToStart} == $V_Href->{stepToStop} # run at only one step
                        && $V_Href->{stepToStart} <= 3 # run at step NO.1, NO.2, or NO.3
                      )
                )
            );
}

#--- assign Hi-C pair with haplotype-unknown reads to intra/inter-haplotype ---
sub DivideHiCreadsToHaplotypes{

    # basic
    &check_files;
    &prepare;
    load_chr_Things;

    # PE to categories
    # [ds]End-h[x], [ds]End-hInter, unkown, discard
    divide_pairBam;

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
    ## load list of pair-se-bam if provide
    if(file_exist(filePath=>$V_Href->{PairBamList})){
      open (BAMLIST, Try_GZ_Read($V_Href->{PairBamList})) || die "fail reading PairBamList: $!\n";
      while(<BAMLIST>){
        next if /^#/;
        chomp;
        push @{$V_Href->{SourceBam}}, $_;
      }
      close BAMLIST;
    }
    ## record pair-se-bam
    my %preBamFiles;
    for my $SourceBam (@{$V_Href->{SourceBam}}){
        if ($SourceBam =~ /^([^,]+R1[^,]*\.bam),([^,]+R2[^,]*\.bam)$/){
            my ($R1_bam_path, $R2_bam_path) = ($1, $2);
            # check existence
            if(    !file_exist(filePath=>$R1_bam_path)
                || !file_exist(filePath=>$R2_bam_path)
            ){
                warn_and_exit "<ERROR>\tbam file does not exist:\n"
                                    ."\t$SourceBam\n";
            }
            # check duplicated input
            if(    exists $preBamFiles{$R1_bam_path}
                || exists $preBamFiles{$R2_bam_path}
            ){
                warn_and_exit "<ERROR>\tencounter same bam file again from input:\n"
                                    ."\t$SourceBam\n";
            }
            # prepare output prefix
            ( my $R1_bam_prefix = basename($R1_bam_path) ) =~ s/[\.\_]R1.*\.bam$//;
            ( my $R2_bam_prefix = basename($R2_bam_path) ) =~ s/[\.\_]R2.*\.bam$//;
            if(    $R1_bam_prefix ne $R2_bam_prefix
                || length($R1_bam_prefix) == 0
                || exists $preBamFiles{$R1_bam_prefix}
            ){
                warn_and_exit "<ERROR>\trequires valid prefix of bam files from input:\n"
                                    ."\t$SourceBam\n";
            }
            # bam object
            my $R1_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $R1_bam_path, tag => $R1_bam_prefix);
            my $R2_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $R2_bam_path, tag => $R2_bam_prefix);
            # records
            my $i = scalar @{$V_Href->{PairBamFiles}};
            push @{$V_Href->{PairBamFiles}}, {R1_bam => $R1_bam, R2_bam => $R2_bam, prefix => $R1_bam_prefix, no => $i};
            $preBamFiles{$_} = 1 for ( $R1_bam_path, $R2_bam_path, $R1_bam_prefix );
        }
        elsif(   basename($SourceBam) =~ /[\.\-_][R_]1[\.\-_]/
              || basename($SourceBam) =~ /[\.\-_][R_]2[\.\-_]/
        ){  # potentially avoid that one line has only one SE bam
            warn_and_exit "<ERROR>\tcannot recognize R1/R2.bam files from input:\n"
                                ."\t$SourceBam\n";
        }
        else{
            # check existence
            if(!file_exist(filePath=>$SourceBam)){
                warn_and_exit "<ERROR>\tbam file does not exist:\n"
                                    ."\t$SourceBam\n";
            }
            # check duplicated input
            if(exists $preBamFiles{$SourceBam}){
                warn_and_exit "<ERROR>\tencounter same bam file again:\n"
                                    ."\t$SourceBam\n";
            }
            # prepare output prefix
            ( my $bam_prefix = basename($SourceBam) ) =~ s/\.bam$//;
            if(    length($bam_prefix) == 0
                || exists $preBamFiles{$bam_prefix}
            ){
                warn_and_exit "<ERROR>\trequires valid prefix of bam files:\n"
                                    ."\t$SourceBam\n";
            }
            # bam object
            my $PE_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $SourceBam, tag => $bam_prefix);
            # records
            my $i = scalar @{$V_Href->{PairBamFiles}};
            push @{$V_Href->{PairBamFiles}}, {PE_bam => $PE_bam, prefix => $bam_prefix, no => $i};
            $preBamFiles{$_} = 1 for ( $SourceBam, $bam_prefix );
        }
    }

    # selected bam to process in step NO.1, NO.2, or NO.3
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

    # chr selection
    for my $chr_ij (qw/ mustMapChr_i mustMapChr_j /){
        if(scalar @{$V_Href->{$chr_ij}->{Af}}){
            $V_Href->{$chr_ij}->{Hf}->{$_} = 1 for @{$V_Href->{$chr_ij}->{Af}};
        }
        else{
            delete $V_Href->{$chr_ij};
        }
    }

    # prepare local region size
    $V_Href->{UKreadsFlankRegUnitMaxTimes} = POSIX::ceil($V_Href->{UKreadsFlankRegMax} / $V_Href->{UKreadsFlankRegUnit});
    $V_Href->{UKreadsFlankRegUnitMidTimes} = int((1 + $V_Href->{UKreadsFlankRegUnitMaxTimes}) / 2);
    $V_Href->{UKreadsFlankRegUnitIniTimes} = $V_Href->{UKreadsFlankRegUnitMidTimes};

    # prepare PEsplitStat (dEnd and sEnd) of step s01
    for my $hap_i (1 .. $V_Href->{haploCount}){
        $V_Href->{PEsplitStat}->{"phMut-dEnd-h$hap_i"} = 0;
        $V_Href->{PEsplitStat}->{"phMut-sEnd-h$hap_i"} = 0;
    }

    # inform selected bam
    if(    $V_Href->{stepToStop} <= 3
        && scalar(keys %{$V_Href->{SelectBamPref}}) != 0
    ){
        stout_and_sterr "[INFO]\t".`date`
                             ."\tSelect bam files with below prefixes to operate.\n";
        stout_and_sterr       "\t$_\n" for sort keys %{$V_Href->{SelectBamPref}};
    }

    # dump options
    shift @{$V_Href->{dump}} if @{$V_Href->{dump}} > 1; # user has input, so discard default
    for my $i (0 .. scalar(@{$V_Href->{dump}})-1){
        my @dumpOpt = split /:/, uc($V_Href->{dump}->[$i]);
        # check dump mode
        if($dumpOpt[0] !~ /^(BP|FRAG)$/i){
            warn_and_exit "<ERROR>\tDumpMode allows 'BP' or 'FRAG'.\n";
        }
        # check dump size
        if(defined $dumpOpt[1]){
            if(!exists $V_Href->{dump_allowBinSize}->{$dumpOpt[0]}->{$dumpOpt[1]}){
                warn_and_exit "<ERROR>\tPlease use required dump BinSize.\n";
            }
        }
        else{ # default
            $dumpOpt[1] = $dumpOpt[0] eq 'BP' ? '1MB' : '1';
        }
        # check haplo combination
        if(defined $dumpOpt[2]){
            $dumpOpt[2] = lc($dumpOpt[2]);
            $dumpOpt[2] =~  s/int/Int/;
        }
        else{ # default
            $dumpOpt[2] = 'all';
        }
        # over-write
        $V_Href->{dump}->[$i] = \@dumpOpt;
    }
    # inform
    if($V_Href->{stepToStop} >= 5){
        stout_and_sterr "[INFO]\t".`date`
                             ."\tdump work details:\n";
        stout_and_sterr       "\t".join(':', @$_)."\n" for @{$V_Href->{dump}};
    }

    # delete possible previous results
    &delete_prev_results;
}

# delete possible previous results
sub delete_prev_results{
    my $SelectBamCount = scalar(keys %{$V_Href->{SelectBamPref}});
    if( $V_Href->{stepToStart} <= 1 ){
        if($SelectBamCount == 0){
            `rm -rf $V_Href->{outdir}/*`;
        }
        else{
            `rm -rf $V_Href->{outdir}/$_-workspace` for keys %{$V_Href->{SelectBamPref}};
        }
    }
    if( $V_Href->{stepToStart} <= 2 ){
        if($SelectBamCount == 0){
            `rm -rf $V_Href->{outdir}/*-workspace/*phMut-sEnd-h[0-9]*.[hs][0-9It]*`; # bam, statOf...
            `rm -rf $V_Href->{outdir}/*-workspace/*phMut-sEnd-h[0-9]*.*chrPair*`; # chrPair dir, list
        }
        else{
            `rm -rf $V_Href->{outdir}/*-workspace/$_.phMut-sEnd-h[0-9]*.[hs][0-9It]*` for keys %{$V_Href->{SelectBamPref}};
            `rm -rf $V_Href->{outdir}/*-workspace/$_.phMut-sEnd-h[0-9]*.*chrPair*` for keys %{$V_Href->{SelectBamPref}};
        }
    }
    if( $V_Href->{stepToStart} <= 3 ){
        if(    $SelectBamCount == 0
            || $V_Href->{stepToStart} != 3
        ){
            `rm -rf $V_Href->{outdir}/*-workspace/*unknown.[hs][0-9It]*`; # bam, statOf...
            `rm -rf $V_Href->{outdir}/*-workspace/*unknown.*chrPair*`; # chrPair dir, list
        }
        else{
            `rm -rf $V_Href->{outdir}/*-workspace/$_.unknown.[hs][0-9It]*` for keys %{$V_Href->{SelectBamPref}};
            `rm -rf $V_Href->{outdir}/*-workspace/$_.unknown.*chrPair*` for keys %{$V_Href->{SelectBamPref}};
        }
    }
    if( $V_Href->{stepToStart} <= 4 ){
        `rm -rf $V_Href->{outdir}/*.merge.h*.bam`;
        `rm -rf $V_Href->{outdir}/*.statOfPhasedLocReg.gz`;
    }
    # if only step NO.5,
    # it might use different setting,
    # or else will overwirte former results.
    if( $V_Href->{stepToStart} < 5 ){
        `rm -rf $V_Href->{outdir}/$V_Href->{dumpSubDir}`;
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
