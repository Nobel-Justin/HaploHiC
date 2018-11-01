package HaploHiC::PhasedHiC::splitPairBam;

use strict;
use warnings;
use List::Util qw/ max min sum first /;
use Data::Dumper;
use Parallel::ForkManager;
use File::Spec::Functions qw/ catfile /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::Quality qw/ baseQ_char2score /;
use BioFuse::BioInfo::Objects::HicReads_OB;
use BioFuse::BioInfo::Objects::HicPairEnd_OB;
use HaploHiC::LoadOn;
use HaploHiC::PhasedHiC::phasedMutWork qw/ PosToPhasedMut /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              divide_pairBam
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::splitPairBam';
#----- version --------
$VERSION = "0.12";
$DATE = '2018-10-27';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        divide_pairBam
                        prepareSplitBamName
                        prepareSplitBamFH
                        PEreads_to_splitBam
                        judge_reads_alignment
                        judge_reads_by_phMut
                        judge_on_PE
                        selectOneFromMultiHap
                     /;

#--- load pair bams and divide PE-reads ---
sub divide_pairBam{

    # variants
    my $rOB_AbufSize = $V_Href->{rOB_AbufferSize};
    my $peC_ReptUnit = $V_Href->{peC_ReportUnit};
    my $forkNum = min( $V_Href->{forkNum}, scalar(@{$V_Href->{PairBamFiles}}) );

    # prepare splitBam name
    &prepareSplitBamName( pairBamHref => $_ ) for @{$V_Href->{PairBamFiles}};

    # start from step after current step
    return if $V_Href->{stepToStart} > 1;

    # fork manager
    my $pm;
    my $fork_DO = ( $forkNum > 1 );
    if( $fork_DO ){ $pm = new Parallel::ForkManager($forkNum) }
    # read paired bam files and split
    for my $bam_i ( 0 .. scalar(@{$V_Href->{PairBamFiles}})-1 ){
        # fork job starts
        if( $fork_DO ){ $pm->start and next; }
        # prepare split-bams
        my $pairBamHref = $V_Href->{PairBamFiles}->[$bam_i];
        # open splitBamFileHandle
        &prepareSplitBamFH( pairBamHref => $pairBamHref );
        # read source pair BAMs
        ## FLAG: 0x100(sd), 0x400(d), 0x800(sp)
        ## WEDO: -F 0x400(d), as 'sd' and 'sp' alignments is important in HiC
        my %SourBamFH;
        open ($SourBamFH{1},"$V_Href->{samtools} view -F 0x400 $pairBamHref->{R1_bam} |") || die"fail read R1_bam: $!\n";
        open ($SourBamFH{2},"$V_Href->{samtools} view -F 0x400 $pairBamHref->{R2_bam} |") || die"fail read R2_bam: $!\n";
        my (%reads_OB) = ( 1=>[], 2=>[] );
        my $pe_Idx = 0;
        while(1){
            # read bam line
            my %bufAvailSize = map { ( $_, $rOB_AbufSize - scalar(@{$reads_OB{$_}}) ) } (1,2);
            for my $rEnd (1,2){
                my $fh = $SourBamFH{$rEnd};
                while(    $bufAvailSize{$rEnd}--
                       && ( my $bamLine = <$fh> )
                ){
                    push @{$reads_OB{$rEnd}}, BioFuse::BioInfo::Objects::HicReads_OB->new(ReadsLineText => $bamLine, reads_end => $rEnd, _rc_optfd => 1, _rc_rdstr => 1);
                }
            }
            # both bams finish reading?
            my $loop_stop_bool = ( $bufAvailSize{1} >= 0 && $bufAvailSize{2} >= 0 );
            # get last PE-id
            my %LastPid = map { ( ($loop_stop_bool ? '__NA__' : $reads_OB{$_}->[-1]->get_pid), 1 ) } (1,2);
            # get shared PE-id
            my %R1_Pid  = map { ($_,1) } grep { !exists $LastPid{$_}                       } map { $_->get_pid } @{$reads_OB{1}};
            my %Sha_Pid = map { ($_,1) } grep { !exists $LastPid{$_} && exists $R1_Pid{$_} } map { $_->get_pid } @{$reads_OB{2}};
            # load PE-object
            my %pe_OB;
            my %pid2idx;
            for my $rEnd (1,2){
                my $shift_time = 0;
                for my $reads_OB ( @{$reads_OB{$rEnd}} ){
                    my $pid = $reads_OB->get_pid;
                    if(    !exists $Sha_Pid{$pid}
                        ||  exists $LastPid{$pid}
                    ){
                        last;
                    }
                    $shift_time ++;
                    # pe object
                    if( !exists $pid2idx{$pid} ){
                        $pe_Idx++;
                        $pid2idx{$pid} = $pe_Idx;
                        $pe_OB{$pe_Idx} = BioFuse::BioInfo::Objects::HicPairEnd_OB->new( pe_Idx => "$bam_i-$pe_Idx" );
                        # inform
                        stout_and_sterr "[INFO]\t".`date`
                                             ."\tload $pe_Idx PE-reads from '$pairBamHref->{prefix}' bam files.\n" if( $pe_Idx % $peC_ReptUnit == 0 );
                    }
                    # load reads_OB to pe_OB
                    $pe_OB{ $pid2idx{$pid} }->load_reads_OB( reads_OB => $reads_OB );
                }
                # release space
                # @{$reads_OB{$rEnd}} = @{$reads_OB{$rEnd}}[ $shift_time .. scalar(@{$reads_OB{$rEnd}})-1 ];
                splice(@{$reads_OB{$rEnd}}, 0, $shift_time);
            }
            undef %pid2idx; # sweep
            # output PE reads to split bam files
            for my $peIdxToDeal (sort {$a<=>$b} keys %pe_OB){
                &PEreads_to_splitBam(pe_OB => $pe_OB{$peIdxToDeal}, splitBamFH_Href => $pairBamHref->{splitBamFH});
                $pe_OB{$peIdxToDeal} = undef; # release memory
            }
            # stop when both bam files finish reading
            last if( $loop_stop_bool );
        }
        # close file-handles
        close $_ for values %SourBamFH;
        close $_ for values %{$pairBamHref->{splitBamFH}};
        $pairBamHref->{splitBamFH} = {}; # empty splitBamFH key
        # inform
        stout_and_sterr "[INFO]\t".`date`
                             ."\ttotally, load $pe_Idx PE-reads from '$pairBamHref->{prefix}' bam files, and split them OK.\n";
        # write report (1st part)
        open (PESRT, Try_GZ_Write($pairBamHref->{PEsplit_report})) || die "fail write PEtoHaplotype.report: $!\n";
        print PESRT "$_:\t$V_Href->{PEsplitStat}->{$_}\n" for grep !/^unknown_to_/, sort keys %{$V_Href->{PEsplitStat}};
        close PESRT;
        # fork job finishes
        if( $fork_DO ){ $pm->finish; }
    }
    # collect fork jobs
    if( $fork_DO ){ $pm->wait_all_children; }

    # stop at current step
    exit(0) if $V_Href->{stepToStop} == 1;
}

#--- prepare split-bam (file-name) ---
# tags:
## t1, phMut-dEnd-hx:     PE-reads, two ends support same single-haplo.
## t2, phMut-dEnd-hInter: PE-reads, two ends support diff haplotypes (Note: no matter hetMut(s) on one end or two ends).
## t3, phMut-sEnd-hx:     PE-reads, one end  supports single-haplo, another end is unknown. (Note: might assign to hInter)
## t4, phMut-sEnd-hInter: PE-reads, one end  supports  multi-haplo, another end is unknown. (Note: must have supplementary alignment)
## t5, discarded: reads filtered out due to um/mm/sd/sp (depends on options), or at least one end ONLY mapped to chr not needed
## t6, unknown: reads qualified but covering none haplotype.  (Note: assign to hx/hInter latter)
## t7, invalid: invalid Hi-C PE (close alignment)
sub prepareSplitBamName{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};

    ## t1,t3: <haplo> bam of phMut (single-haplo) covered reads
    my @tags = map {("phMut-dEnd-h$_" , "phMut-sEnd-h$_")} (1 .. $V_Href->{haploCount});
    ## t2,t4: <haplo> bam of phMut (multi-haplo) covered reads
    push @tags, "phMut-dEnd-hInter", "phMut-sEnd-hInter";
    ## t5: bam of discarded reads
    push @tags, "discarded";
    ## t6: bam of haplo-unknown reads
    push @tags, "unknown";
    ## t7: invalid Hi-C PE
    push @tags, "invalid";
    # workspace
    $pairBamHref->{workSpace} = catfile( $V_Href->{outdir}, $pairBamHref->{prefix}.'-workspace' );
    (-d $pairBamHref->{workSpace}) || `mkdir -p $pairBamHref->{workSpace}`;
    # split-bam files
    $pairBamHref->{splitBam}->{$_} = catfile( $pairBamHref->{workSpace}, $pairBamHref->{prefix}.".$_.bam" ) for @tags;
    # report
    $pairBamHref->{PEsplit_report} = catfile( $V_Href->{outdir}, $pairBamHref->{prefix}.'.PEtoHaplotype.report' );
    # prepare bam to final merge
    push @{$pairBamHref->{bamToMerge}->{"merge.h${_}Intra"}}, $pairBamHref->{splitBam}->{"phMut-dEnd-h$_"} for (1 .. $V_Href->{haploCount});
    push @{$pairBamHref->{bamToMerge}->{"merge.hInter"}}, $pairBamHref->{splitBam}->{$_} for ("phMut-dEnd-hInter", "phMut-sEnd-hInter");
}

#--- start writing to split-bam (file-w-handle, SAM-header) ---
sub prepareSplitBamFH{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};

    # open split-bam file handles
    for my $tag (sort keys %{$pairBamHref->{splitBam}}){
        open ($pairBamHref->{splitBamFH}->{$tag}, "| $V_Href->{samtools} view -b -o $pairBamHref->{splitBam}->{$tag}")
          || die "fail write $pairBamHref->{prefix} split.bam ($tag): $!\n";
    }
    # write SAM-header
    ## copy original SAM-header
    open (BAMH,"$V_Href->{samtools} view -H $pairBamHref->{R1_bam} |") || die"fail read header of R1_bam: $!\n";
    while( my $Hline = <BAMH> ){
        next unless( $Hline =~ /^\@SQ/ || $Hline =~ /^\@RG/ );
        print {$_} $Hline for values %{$pairBamHref->{splitBamFH}};
    }
    close BAMH;
    ## more in PG info
    my $PG_info  = "\@PG\tID:$MODULE_NAME\tPN:haplo_div\tVN:$VERSION\tCL:$V_Href->{MainName} haplo_div";
       # $PG_info .= " -sampid $V_Href->{sampleID}";
       $PG_info .= " -phsvcf $V_Href->{phased_vcf}";
       $PG_info .= " -outdir $V_Href->{outdir}";
       $PG_info .= " -bam $_" for @{$V_Href->{PairBamList}};
       $PG_info .= " -samt $V_Href->{samtools}";
       $PG_info .= " -db_dir $V_Href->{db_dir}";
       $PG_info .= " -ref_v $V_Href->{ref_version}";
       $PG_info .= " -hapct $V_Href->{haploCount}";
       $PG_info .= " -use_indel" if($V_Href->{use_InDel});
       $PG_info .= " -use_sp" if($V_Href->{use_spmap});
       $PG_info .= " -use_sd" if($V_Href->{use_sdmap});
       $PG_info .= " -use_mm" if($V_Href->{use_multmap});
       $PG_info .= " -min_mq $V_Href->{min_mapQ}";
       $PG_info .= " -min_bq $V_Href->{min_baseQ}";
       $PG_info .= " -min_re $V_Href->{min_distToRedge}";
       $PG_info .= " -qual $V_Href->{baseQ_offset}";
    print {$_} "$PG_info\n" for values %{$pairBamHref->{splitBamFH}};
}

#--- distribute PE-reads to different categories ---
# judgement on each reads' alignment:
## UM: UnMapped
## MM: Multiple-Mapped
## SD: Secondary alignment
## SP: Supplementary alignment
## LM: Low Mapping quality
## UK: haplotype Un-Known
## LB: Low Base quality
## RE: short Reads-Edge
## h1: haplotype-NO.1
## h2: haplotype-NO.2
sub PEreads_to_splitBam{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB = $parm{pe_OB};
    my $splitBamFH_Href = $parm{splitBamFH_Href};

    # make judgement on each reads-end
    my %hap2phMut;
    &judge_on_rEnd( pe_OB => $pe_OB, hap2phMutHref => \%hap2phMut );

    # judge PE-reads
    my %Match_haplo;
    my ($SplitBam_tag, $Split_marker) = &judge_on_PE( pe_OB => $pe_OB, MatchHap_Href => \%Match_haplo );

    # add OWN 'XH:Z:' haplo-id to 'optfd' of each reads_OB
    if( $SplitBam_tag =~ /^phMut-/ ){
        $pe_OB->addHapIDtoReadsOptfd( reads_end => $_ ) for (1,2);
    }
    # write PE to split-bam files
    my $peSAM_Aref = $pe_OB->printSAM;
    ## 'XS:Z:' denotes marker of 'S'plit reason
    $_ .= "\tXS:Z:$Split_marker" for @$peSAM_Aref;
    ## add phMut supporting haplotype
    if( $SplitBam_tag =~ /^phMut-/ ){ # 'dEnd-h' and 'sEnd-h'
        my @all_hapInfo;
        my @Match_haplo = sort keys %Match_haplo;
        for my $hapID ( @Match_haplo ){
            my @hapInfo;
            for my $phMutNO (sort {$a<=>$b} keys %{$hap2phMut{$hapID}}){
                my $phMut_OB = $hap2phMut{$hapID}->{$phMutNO};
                # phasedMut info
                push @hapInfo, join(',', $phMut_OB->get_chr, $phMut_OB->get_pos,
                                         @{ $phMut_OB->get_hapInfo(haplo_ID=>$hapID) } );
            }
            push @all_hapInfo, "$hapID(".join(';',@hapInfo).")";
        }
        # 'XM:Z:' denotes the phased'M'ut used for haplotype detection
        $_ .= "\tXM:Z:".join('|',@all_hapInfo) for @$peSAM_Aref;
    }
    ## output
    print {$splitBamFH_Href->{$SplitBam_tag}} "$_\n" for @$peSAM_Aref;
    ## stat
    $V_Href->{PEsplitStat}->{$SplitBam_tag} ++;
}

#--- judgement on on each reads-end ---
sub judge_on_rEnd{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB = $parm{pe_OB};
    my $hap2phMutHref = $parm{hap2phMutHref};

    # discard abnormal supplementary alignment
    if( $pe_OB->discardAbnormalSP ){
        $V_Href->{PEsplitStat}->{remove_SP} ++;
    }

    # judgement on each end's each alignment
    for my $rEnd (1,2){
        my $reads_OB_Aref = $pe_OB->get_reads_OB( reads_end => $rEnd );
        for my $reads_OB ( @$reads_OB_Aref ){ # first comes the best alignment???
            # qulifying on alignment
            my $alignJudge = &judge_reads_alignment( reads_OB => $reads_OB );
            # record
            $reads_OB->load_AlignJudge( alignJudge => $alignJudge );
            # skip bad scenarios (UM/SD/SP/MM/LM) #
            if(     $alignJudge =~ /UM,/ # UnMapped
                ||  $alignJudge =~ /LM,/ # Low MapQ
                || ($alignJudge =~ /SD,/ && !$V_Href->{use_sdmap}) # Secondary alignment, if donot want
                || ($alignJudge =~ /SP,/ && !$V_Href->{use_spmap}) # Supplementary alignment, if donot want
                || ($alignJudge =~ /MM,/ && !$V_Href->{use_multmap}) # Multiple-Mapped, if donot want
            ){
                last;
            }
            # if good, check haplotype covered of this alignment
            ## get phased Mut(s) this reads covers
            my $phMut_OB_Aref = PosToPhasedMut( chr => $reads_OB->get_mseg, pos => $reads_OB->get_mpos, ext => $reads_OB->get_mRefLen - 1 );
            # no phased-mut covered #
            if( scalar(@$phMut_OB_Aref) == 0 ){
                $reads_OB->load_AlignJudge( alignJudge => 'UK,' ); # haplotype Un-Known
                next;
            }
            # get phased Mut(s) alleles from this reads
            my $phMut_covSign = 0;
            for my $phMut_OB (@$phMut_OB_Aref){
                # see whether reads_OB match against phMut_OB
                my $jdA = &judge_reads_by_phMut( reads_OB => $reads_OB, phMut_OB => $phMut_OB );
                if( $jdA->[0] ne 'MISS' ){
                    $reads_OB->load_AlignJudge( alignJudge => $jdA->[0] );
                    $phMut_covSign = 1;
                    # link phMut_OB to haplotype-ID
                    if( $jdA->[0] =~ /^h\d/ ){
                        (my $hapID = $jdA->[0]) =~ s/,$//;
                        $hap2phMutHref->{$hapID}->{ $phMut_OB->get_phMutNO } = $phMut_OB;
                        # record haplo+alleleInfo to this reads_OB
                        $reads_OB->load_SuppHaplo( hapID => $hapID, allele_OB => $jdA->[1] );
                    }
                }
            }
            # if all phMut(s) not covered #
            $reads_OB->load_AlignJudge( alignJudge => 'UK,' ) unless $phMut_covSign;
        }
    }
}

#--- qulifying read's alignment ---
# UM, SD, SP, MM, LM
sub judge_reads_alignment{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $reads_OB = $parm{reads_OB};

    my $alignJudge = ''; # good

    if( $reads_OB->is_unmap ){ # UnMapped
        $alignJudge .= 'UM,';
    }
    else{
        # Secondary alignment
        if( $reads_OB->is_2ndmap ){
            $alignJudge .= 'SD,';
        }
        # Supplementary alignment
        if( $reads_OB->is_suppmap ){
            $alignJudge .= 'SP,';
        }
        # Multiple-Mapped
        if( $reads_OB->is_mltmap ){
            $alignJudge .= 'MM,';
        }
        # Low Mapping quality
        if( $reads_OB->get_mapQ < $V_Href->{min_mapQ} ){
            $alignJudge .= 'LM,';
        }
    }
    # return
    return $alignJudge;
}

#--- see whether reads_OB match against phMut_OB ---
sub judge_reads_by_phMut{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $reads_OB = $parm{reads_OB};
    my $phMut_OB = $parm{phMut_OB};

    # get allele info of this position from this reads_OB
    my @parm = ( chr => $phMut_OB->get_chr, pos => $phMut_OB->get_pos, refBase => $phMut_OB->get_refAllele, no_warn => 1 );
    push @parm, ( preferINDEL => 1 ) if( $phMut_OB->has_indel );
    my $allele_OB = $reads_OB->get_pos_allele( @parm );

    # not cover this phMut pos #
    if( $allele_OB->get_type eq 'miss' ){ 
        return ['MISS']; # reads_OB miss the phMut
    }
    #----------------------------#
    # low base quality filtering #
    #----------------------------#
    # get mean base quality
    my $alleleQual = $allele_OB->get_qual( Q_offset => $V_Href->{baseQ_offset} );
    # filter
    if( $alleleQual < $V_Href->{min_baseQ} ){
        return ['LB,']; # Low Base quality
    }
    #---------------------------#
    # reads-edge-dist filtering #
    #---------------------------#
    if(    $allele_OB->get_rEdgeDist5 < $V_Href->{min_distToRedge}
        || $allele_OB->get_rEdgeDist3 < $V_Href->{min_distToRedge}
        || $allele_OB->get_rEdgeDist3 < $phMut_OB->get_judgeReadEdgeDist
    ){
        return ['RE,']; # short Reads-Edge
    }
    #---------------------------#
    # close alteration distance #
    #---------------------------#
    if(    $allele_OB->get_nAltDist5 < $V_Href->{min_distToNearbyAlt}
        || $allele_OB->get_nAltDist3 < $V_Href->{min_distToNearbyAlt}
    ){
        return ['AD,']; # Alteration Distance
    }
    #-----------------#
    # which haplotype #
    # might NOT match #
    #-----------------#
    my $haplo_ID = $phMut_OB->allele2haploID( allele_OB => $allele_OB );
    if( defined $haplo_ID ){
        return [ "$haplo_ID,", $allele_OB ];
    }
    else{
        return [ 'NH,' ]; # Non-Halptype matched
    }
}

#--- judgement on pe-reads ---
sub judge_on_PE{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB = $parm{pe_OB};
    my $MatchHap_Href = $parm{MatchHap_Href};

    my $rEndJdgHref = {};
    $rEndJdgHref->{$_}->{J} = $pe_OB->get_rEndWholeAlignJudge(rEnd=>$_) for (1,2);
    $rEndJdgHref->{$_}->{H} = $pe_OB->get_rEndWholeSuppHaplo(rEnd=>$_)  for (1,2);

    my $SplitBam_tag = undef;
    my $Split_marker = undef;

    #----#
    # t5 #
    #----#
    ## t5, has unmap, one or both
    if(    $rEndJdgHref->{1}->{J} =~ /UM/
        || $rEndJdgHref->{2}->{J} =~ /UM/
    ){
        $SplitBam_tag = 'discarded';
        $Split_marker = 'Unmap';
    }
    ## t5, has multiple mapped, one or both
    elsif(   !$V_Href->{use_multmap}
          && (   $rEndJdgHref->{1}->{J} =~ /MM/
              || $rEndJdgHref->{2}->{J} =~ /MM/
             )
    ){
        $SplitBam_tag = 'discarded';
        $Split_marker = 'MultiMap';
    }
    ## t5, has Low Mapping quality, one or both
    elsif(   $rEndJdgHref->{1}->{J} =~ /LM/
          || $rEndJdgHref->{2}->{J} =~ /LM/
    ){
        $SplitBam_tag = 'discarded';
        $Split_marker = 'LowMapQ';
    }
    ## t5, has secondary alignment, one or both
    elsif(   !$V_Href->{use_sdmap}
          && (   $rEndJdgHref->{1}->{J} =~ /SD/
              || $rEndJdgHref->{2}->{J} =~ /SD/
             )
    ){
        $SplitBam_tag = 'discarded';
        $Split_marker = 'SecondAlign';
    }
    ## t5, has supplementary alignment, one or both
    elsif(   !$V_Href->{use_spmap}
          && (   $rEndJdgHref->{1}->{J} =~ /SP/
              || $rEndJdgHref->{2}->{J} =~ /SP/
             )
    ){
        $SplitBam_tag = 'discarded';
        $Split_marker = 'SuppleAlign';
    }
    #-----------------#
    # shortcut return #
    #-----------------#
    return ($SplitBam_tag, $Split_marker) if defined $SplitBam_tag;
    ## t5, test whether has end mapped to not-need chr
    for my $rEnd (1,2){
        my $cnt_Hf = $pe_OB->test_need_RefSeg(reads_end => $rEnd, refSegHref => $V_Href->{ChrThings});
        if($cnt_Hf->{noNeed}){ # has not-need chr mapped
            if(!$cnt_Hf->{isNeed}){ # doesn't have need chr mapped
                $SplitBam_tag = 'discarded';
                $Split_marker = "OnlyNoNeedChr(r$rEnd)";
                last;
            }
            else{ # else, only keep alignment of need chr
                $pe_OB->onlyKeep_need_RefSeg(reads_end => $rEnd, refSegHref => $V_Href->{ChrThings});
            }
        }
    }
    #-----------------#
    # shortcut return #
    #-----------------#
    return ($SplitBam_tag, $Split_marker) if defined $SplitBam_tag;
    #---------------#
    # t7, t1-t4, t6 #
    #---------------#
    ## t7, invalid Hi-C PE, close alignment
    if( $pe_OB->isInValidPair( minSplitReadGap => $V_Href->{min_splitReadGap} ) ){
        $SplitBam_tag = 'invalid';
        $Split_marker = 'invalidPE';
    }
    ## t1-t4, matches one or more haplotype(s)
    elsif(   $rEndJdgHref->{1}->{J} =~ /h\d+/
          || $rEndJdgHref->{2}->{J} =~ /h\d+/
    ){
        # only keep hap-info
        $rEndJdgHref->{$_}->{J} =~ s/[^h,\d][^,\d],//g for (1,2);
        # end-NO <-> hapID
        my %reads2HapID = map { ($_,{map {($_,1)} split /,/, $rEndJdgHref->{$_}->{J}}) } (1,2);
        my %readsHapAmt = map { ($_,scalar(keys %{$reads2HapID{$_}})) } (1,2);
        my %rSoloHapID  = ( 1 => [ grep !exists $reads2HapID{2}->{$_}, keys %{$reads2HapID{1}} ],
                            2 => [ grep !exists $reads2HapID{1}->{$_}, keys %{$reads2HapID{2}} ]  );
        # reset if empty
        $rEndJdgHref->{$_}->{J} ||= 'none' for (1,2);
        # for good output
        $rEndJdgHref->{$_}->{J} =~ s/,$//  for (1,2);

        # Category  R1_Hap   R2_hap   C_Multi  C_Sum  Solo_R1  Solo_R2
        # t1        h1,[uk]  h1,[uk]  1        2      0        0
        # t1        h2,[uk]  h2,[uk]  1        2      0        0
        # t2        h1,h2    h1,[uk]  2        3      1        0
        # t2        h1,h2    h2,[uk]  2        3      1        0
        # t2        h1,h2    h1,h2    4        4      0        0
        # t2        h1,[uk]  h1,h2    2        3      0        1
        # t2        h2,[uk]  h1,h2    2        3      0        1
        # t2        h1,[uk]  h2,[uk]  1        2      1        1
        # t2        h2,[uk]  h1,[uk]  1        2      1        1
        # t3        h1,[uk]  uk       0        1      1        0
        # t3        h2,[uk]  uk       0        1      1        0
        # t3        uk       h1,[uk]  0        1      0        1
        # t3        uk       h2,[uk]  0        1      0        1
        # t4        h1,h2    uk       0        2      2        0
        # t4        uk       h1,h2    0        2      0        2

        ## t1, phMut-dEnd-hx, 2 cases
        ## OR, phMut-sEnd-hx
        if(    $readsHapAmt{1} * $readsHapAmt{2} == 1 # R1 and R2 both have only one hap, and
            && scalar(@{$rSoloHapID{1}}) == 0         # they share same hap
        ){
            my $hapID = (split /,/, $rEndJdgHref->{1}->{J})[0];
            $MatchHap_Href->{$hapID} = 1;
            # check whether at least one end has 'UK'
            # if no,  phMut-dEnd-hx
            # if yes, phMut-dEnd-hx or phMut-sEnd-hx, based on close-alignment
            my $changeTOsEndBool = $pe_OB->dEndSameHapJudge( minSplitReadGap => $V_Href->{min_splitReadGap} );
            if( ! $changeTOsEndBool ){
                $SplitBam_tag = "phMut-dEnd-$hapID";
                $Split_marker = "dEndSoloHap;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
            }
            else{
                $SplitBam_tag = "phMut-sEnd-$hapID";
                $rEndJdgHref->{$_}->{J} =~ s/\b$hapID\b/$hapID(C)/g for (1,2); # add '(C)' to mark 'Close aligned'
                $Split_marker = "dEndSoloHapClose;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
            }
        }
        ## t3, phMut-sEnd-hx, 4 cases
        elsif( $readsHapAmt{1} + $readsHapAmt{2} == 1 ){ # either of R1 and R2 has only one hap
            my $hapID = $rSoloHapID{1}->[0] || $rSoloHapID{2}->[0];
            $SplitBam_tag = "phMut-sEnd-$hapID";
            $Split_marker = "sEndSoloHap;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
            $MatchHap_Href->{$hapID} = 1;
        }
        ## t4, phMut-sEnd-hInter, 2 cases
        ## OR, phMut-sEnd-hx
        elsif(    $readsHapAmt{1} * $readsHapAmt{2} == 0 # either of R1 and R2 is unknown, and
               && $readsHapAmt{1} + $readsHapAmt{2} >  1 # another one has multi-hap
        ){
            # the end that has multi-hap should meet:
            ##  1) have two alignments (i.e., SP) that interChr or large-distance 'heuristic';
            ##  2) one of two alignments is close to the other end's alignment
            ## *3) might have enzyme in its sequence (not in use).
            my ($mhEnd,$uhEnd) = $readsHapAmt{1} ? (1,2) : (2,1); # mh: multi-hap; uh: unkown-hap
            my $doSelectOneHapBool = $pe_OB->sEndInterHapJudge( mhEnd => $mhEnd, uhEnd => $uhEnd, minSplitReadGap => $V_Href->{min_splitReadGap} );
            if( ! $doSelectOneHapBool ){
                $SplitBam_tag = 'phMut-sEnd-hInter';
                $Split_marker = "sEndMultHap;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
                $MatchHap_Href->{$_} = 1 for keys %{$reads2HapID{1}};
                $MatchHap_Href->{$_} = 1 for keys %{$reads2HapID{2}};
            }
            else{
                my $selectHapID = &selectOneFromMultiHap( hapAlleleHf_Aref => [$rEndJdgHref->{$mhEnd}->{H}] );
                $SplitBam_tag = "phMut-sEnd-$selectHapID";
                $rEndJdgHref->{$mhEnd}->{J} =~ s/\b$selectHapID\b/$selectHapID(S)/g; # add '(S)' to mark 'Selected hapID'
                $Split_marker = "sEndSoloHapSelect;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
                $MatchHap_Href->{$selectHapID} = 1;
                # dislink other haplo(s) with all reads_OBs of this pe_OB
                $_->onlyKeep_SuppHaplo( hapID => $selectHapID ) for @{ $pe_OB->get_reads_OB( reads_end => $mhEnd ) };
            }
        }
        ## t2, phMut-dEnd-hInter, 7 cases remaining
        ## OR, phMut-xEnd-hx
        else{
            # two PE-ends should meet:
            ## 1) alignments support same haplo  should be close;
            ## 2) alignments support diff haplos should be not close;
            my $doSelectOneHapBool = $pe_OB->dEndInterHapJudge( minSplitReadGap => $V_Href->{min_splitReadGap} );
            if( ! $doSelectOneHapBool ){
                $SplitBam_tag = 'phMut-dEnd-hInter';
                $Split_marker = "dEndMultHap;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
                $MatchHap_Href->{$_} = 1 for keys %{$reads2HapID{1}};
                $MatchHap_Href->{$_} = 1 for keys %{$reads2HapID{2}};
            }
            else{
                my $selectHapID = &selectOneFromMultiHap( hapAlleleHf_Aref => [$rEndJdgHref->{1}->{H},$rEndJdgHref->{2}->{H}] );
                my $newEndType = (    exists $rEndJdgHref->{1}->{H}->{$selectHapID}
                                   && exists $rEndJdgHref->{2}->{H}->{$selectHapID}
                                 ) ? 'dEnd' : 'sEnd';
                $SplitBam_tag = "phMut-$newEndType-$selectHapID";
                $rEndJdgHref->{$_}->{J} =~ s/\b$selectHapID\b/$selectHapID(S)/g for (1,2); # add '(S)' to mark 'Selected hapID'
                $Split_marker = "${newEndType}SoloHapSelect;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
                $MatchHap_Href->{$selectHapID} = 1;
                # dislink other haplo(s) with all reads_OBs of this pe_OB
                $_->onlyKeep_SuppHaplo( hapID => $selectHapID ) for @{ $pe_OB->get_reads_OB( reads_end => 1 ) };
                $_->onlyKeep_SuppHaplo( hapID => $selectHapID ) for @{ $pe_OB->get_reads_OB( reads_end => 2 ) };
            }
        }
    }
    ## t6, unknown
    else{
        $SplitBam_tag = 'unknown';
        $Split_marker = 'unknownHap';
    }

    return ($SplitBam_tag, $Split_marker);
}

#--- compare multi-haplo alleles to select one ---
sub selectOneFromMultiHap{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $hapAlleleHf_Aref = $parm{hapAlleleHf_Aref};

    # calculate score of each haolo-allele
    ## score = NORM_dist_score + NORM_qual_score
    ## NORM_dist_score = sqrt( rEdgeDist_score * aNearDist_score ) / r_mReadLen
    ## rEdgeDist_score = sqrt( rEdgeDist5 * rEdgeDist3 );
    ## aNearDist_score = sqrt( nAltDist5 * nAltDist3 );
    ## NORM_qual_score = get_qual / 40
    ## PS. base-quality range approximates to 40
    my %score;
    for my $hapAllele_Href ( @$hapAlleleHf_Aref ){
        for my $hapID ( keys %$hapAllele_Href ){
            for my $allele_OB ( @{ $hapAllele_Href->{$hapID} } ){
                # dist as score
                my $score = sqrt(   sqrt( $allele_OB->get_rEdgeDist5 * $allele_OB->get_rEdgeDist3 )
                                  * sqrt( $allele_OB->get_nAltDist5  * $allele_OB->get_nAltDist3  )
                                ) / $allele_OB->get_readsOB->get_mReadLen
                            + $allele_OB->get_qual / 40;
                push @{$score{$hapID}}, $score;
            }
        }
    }

    # compare count
    my %cnt2hap;
    push @{$cnt2hap{scalar(@{$score{$_}})}}, $_ for keys %score;
    my ($maxCnt) = max(keys %cnt2hap);
    if( scalar(@{$cnt2hap{$maxCnt}}) == 1 ){
        return $cnt2hap{$maxCnt}->[0];
    }

    # compare score
    my $maxMeanScore = 0;
    my $selectHapID;
    for my $hapID ( @{$cnt2hap{$maxCnt}} ){
        my $meanScore = int( sum(@{$score{$hapID}}) / $maxCnt );
        if( $meanScore > $maxMeanScore ){
            $maxMeanScore = $meanScore;
            $selectHapID = $hapID;
        }
    }
    return $selectHapID;
}

#--- 
1; ## tell the perl script the successful access of this module.
