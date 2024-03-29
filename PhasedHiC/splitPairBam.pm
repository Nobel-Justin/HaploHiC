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
use BioFuse::Util::Sort qw/ sortByStrAndSubNum /;
use BioFuse::BioInfo::Quality qw/ baseQ_char2score /;
use BioFuse::BioInfo::Objects::SeqData::HicReads_OB;
use BioFuse::BioInfo::Objects::SeqData::HicPairEnd_OB;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use HaploHiC::LoadOn;
use HaploHiC::PhasedHiC::dumpContacts qw/ load_enzyme_site_list /;
use HaploHiC::PhasedHiC::phasedMutWork qw/ load_phased_VCF PosToPhasedMut /; # release_phaseMut_OB

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              divide_pairBam
              getTODOpairBamHrefArray
              forkSetting
              write_peOB_to_chrPairBam
              findContactAnchor
              getChrPairTagAf
              chrPairBamToSortSplitBam
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::splitPairBam';
#----- version --------
$VERSION = "0.32";
$DATE = '2021-08-13';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        divide_pairBam
                        prepareSplitBamObj
                        getTODOpairBamHrefArray
                        forkSetting
                        sourceToSplitBam
                        startWriteSplitBam
                        loadrOBfromSourceBam
                        rOB_to_peOB
                        PEreads_to_splitBam
                        judge_on_rEnd
                        judge_reads_alignment
                        judge_reads_by_phMut
                        judge_on_PE
                        selectOneFromMultiHap
                        sortSplitBamByPEcontact
                        write_peOB_to_chrPairBam
                        getIdxOfrOBforChrPair
                        findContactAnchor
                        chrPairBamToSortSplitBam
                        getChrPairTagAf
                        write_sort_peOB_to_newSplitBam
                     /;

#--- load pair bams and divide PE-reads ---
sub divide_pairBam{

    # prepare splitBam name
    &prepareSplitBamObj(pairBamHref => $_) for @{$V_Href->{PairBamFiles}};

    # start from step after current step
    return if $V_Href->{stepToStart} > 1;

    # all bams or selected
    my @TODOpairBamHref = &getTODOpairBamHrefArray;
    # fork manager
    my ($pm, $fork_DO) = &forkSetting;
    # read paired bam files and split
    for my $i ( 0 .. $#TODOpairBamHref ){
        my $pairBamHref = $TODOpairBamHref[$i];
        # fork job starts
        if($fork_DO){ $pm->start($pairBamHref->{prefix}) and next }
        eval{
            # phased Het-mutation
            $V_Href->{skipS01Report} = 1 if $i != 0; # only the first pairBam produce the mutation report
            load_phased_VCF;
            # load enzyme sites list for 'invalid PE'
            load_enzyme_site_list;
            # distribute Hi-C PE-reads from source pairBam to different splitBam
            &sourceToSplitBam(pairBamHref => $pairBamHref);
            # sort splitBam, sEnd and unknown
            &sortSplitBamByPEcontact(pairBamHref => $pairBamHref);
        };
        if($@){
            if($fork_DO){ warn $@; $pm->finish(1) }
            else{ warn_and_exit $@; }
        }
        # fork job finishes
        if($fork_DO){ $pm->finish(0) }
    }
    # collect fork jobs
    if($fork_DO){ $pm->wait_all_children }

    # stop at current step
    exit(0) if $V_Href->{stepToStop} == 1;
}

#--- prepare split-bam objects ---
# tags:
## t1, phMut-dEnd-hx:     PE-reads, two ends support same single-haplo.
## t2, phMut-dEnd-hInter: PE-reads, two ends support diff haplotypes (Note: no matter hetMut(s) on one end or two ends).
## t3, phMut-sEnd-hx:     PE-reads, one end  supports single-haplo, another end is unknown. (Note: might assign to hInter)
## t4, phMut-sEnd-hInter: PE-reads, one end  supports  multi-haplo, another end is unknown. (Note: must have supplementary alignment)
## t5, discarded: reads filtered out due to um/mm/sd/sp (depends on options), or at least one end ONLY mapped to chr not needed
## t6, unknown: reads qualified but covering none haplotype.  (Note: assign to hx/hInter later)
## t7, invalid: invalid Hi-C PE (same fragment or [close alignment])
sub prepareSplitBamObj{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};

    my $pairBamPrefix = $pairBamHref->{prefix};
    # HaploHiC workspace for this pairBam
    $pairBamHref->{workSpace} = catfile( $V_Href->{outdir}, $pairBamPrefix.'-workspace' );
    (-d $pairBamHref->{workSpace}) || `mkdir -p $pairBamHref->{workSpace}`;
    # HaploHiC report for this pairBam
    $pairBamHref->{PEsplit_report} = catfile( $V_Href->{outdir}, $pairBamPrefix.'.PEtoHaplotype.report' );

    # series of split-bam files with different tags
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
    # record split-bam objects
    for my $tag (@tags){
        my $splitBamPath = catfile($pairBamHref->{workSpace}, $pairBamPrefix.".$tag.bam");
        $pairBamHref->{splitBam}->{$tag} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $splitBamPath, tag => "$pairBamPrefix $tag");
    }

    # prepare bam for final merge
    ## here includes t1,t2,t4
    push @{$pairBamHref->{bamToMerge}->{"merge.h${_}Intra"}}, $pairBamHref->{splitBam}->{"phMut-dEnd-h$_"} for (1 .. $V_Href->{haploCount});
    push @{$pairBamHref->{bamToMerge}->{"merge.hInter"}}, $pairBamHref->{splitBam}->{$_} for ("phMut-dEnd-hInter", "phMut-sEnd-hInter");
}

#--- return Href-array of all bams or selected ---
sub getTODOpairBamHrefArray{
    # all bams or selected
    return   scalar(keys %{$V_Href->{SelectBamPref}}) == 0
           ? @{$V_Href->{PairBamFiles}}
           : grep exists $V_Href->{SelectBamPref}->{$_->{prefix}}, @{$V_Href->{PairBamFiles}};
}

#--- fork manager setting ---
sub forkSetting{
    my $forkNum = min( $V_Href->{forkNum}, scalar(@{$V_Href->{PairBamFiles}}), scalar(keys %{$V_Href->{SelectBamPref}}) || 1E5 );
    if($forkNum > 1){
        my $pm = new Parallel::ForkManager($forkNum);
        $pm->run_on_finish(
            sub {
                my ($pid, $exit_code, $bamPrefix, $signal) = @_;
                warn_and_exit "<ERROR>\t$bamPrefix child-process($pid) failed.\n" if $exit_code || $signal;
            }
        );
        return ($pm, 1);
    }
    else{
        return (undef, 0);
    }
}

#--- distribute reads to splitBam from source pairBam ---
sub sourceToSplitBam{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};

    # open splitBam's FileHandle and write header
    &startWriteSplitBam(pairBamHref => $pairBamHref);

    # read source pair BAMs
    my (%rOBbuff) = ( 1=>[], 2=>[] );
    my $peIdx = 0;
    my $peIdxInform = $V_Href->{peC_ReportUnit};
    ## FLAG: 0x100(sd), 0x400(d), 0x800(sp)
    ## WEDO: -F 0x400(d), as 'sd' and 'sp' alignments is important in HiC
    my %SourBamFH;
    if(exists $pairBamHref->{R1_bam}){
        $SourBamFH{1} = $pairBamHref->{R1_bam}->start_read(samtools => $V_Href->{samtools}, viewOpt => '-F 0x400');
        $SourBamFH{2} = $pairBamHref->{R2_bam}->start_read(samtools => $V_Href->{samtools}, viewOpt => '-F 0x400');
    }
    else{
        $SourBamFH{1} = $pairBamHref->{PE_bam}->start_read(samtools => $V_Href->{samtools}, viewOpt => '-F 0x400 -f 0x40');
        $SourBamFH{2} = $pairBamHref->{PE_bam}->start_read(samtools => $V_Href->{samtools}, viewOpt => '-F 0x400 -f 0x80');
    }
    while(1){
        # load rOB from source bam
        my $fileEnd = &loadrOBfromSourceBam(SourBamFH_Hf => \%SourBamFH, rBuff_Hf => \%rOBbuff);
        # load PE-object
        my %pe_OB;
        &rOB_to_peOB(rBuff_Hf => \%rOBbuff, peOB_Hf => \%pe_OB,
                     fileEnd => $fileEnd, pairBamHref => $pairBamHref,
                     peIdx_Sf => \$peIdx, peIdxInform_Sf => \$peIdxInform);
        # output PE reads to split bam files
        &PEreads_to_splitBam(pe_OB => $pe_OB{$_}, pairBamHref => $pairBamHref) for sort {$a<=>$b} keys %pe_OB;
        # stop when both bam files finish reading
        last if $fileEnd;
    }
    # close file-handles
    close $_ for values %SourBamFH;
    $_->stop_write for values %{$pairBamHref->{splitBam}};
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\ttotally, load $peIdx PE-reads from '$pairBamHref->{prefix}' bam files, and split them OK.\n";

    # write report (1st part)
    open (PESRT, Try_GZ_Write($pairBamHref->{PEsplit_report})) || die "fail write PEtoHaplotype.report: $!\n";
    print PESRT "$_:\t$V_Href->{PEsplitStat}->{$_}\n" for sort keys %{$V_Href->{PEsplitStat}};
    close PESRT;
}

#--- start writing split-bam (file-w-handle, SAM-header) ---
sub startWriteSplitBam{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};

    # prepare SAM-header
    my $bam = $pairBamHref->{R1_bam} || $pairBamHref->{PE_bam};
    my $HeadStr = join('', grep {/^\@SQ/ || /^\@RG/} @{ $bam->header_Af(samtools => $V_Href->{samtools}) });
    ## more in PG info
    my $PG_info  = "\@PG\tID:$MODULE_NAME\tPN:haplo_div\tVN:$VERSION\tCL:$V_Href->{MainName} haplo_div";
       # $PG_info .= " -sampid $V_Href->{sampleID}";
       $PG_info .= " -phsvcf $V_Href->{phased_vcf}";
       $PG_info .= " -outdir $V_Href->{outdir}";
       $PG_info .= " -bam $_" for @{$V_Href->{SourceBam}};
       $PG_info .= " -samt $V_Href->{samtools}";
       $PG_info .= " -db_dir $V_Href->{db_dir}";
       $PG_info .= " -ref_v $V_Href->{ref_version}";
       $PG_info .= " -hapct $V_Href->{haploCount}";
       $PG_info .= " -use_indel" if($V_Href->{use_InDel});
       $PG_info .= " -chr_i ".join(' -chr_i ', @{$V_Href->{mustMapChr_i}->{Af}}) if exists $V_Href->{mustMapChr_i};
       $PG_info .= " -chr_j ".join(' -chr_j ', @{$V_Href->{mustMapChr_j}->{Af}}) if exists $V_Href->{mustMapChr_j};
       $PG_info .= " -flt_pe" if($V_Href->{filter_pe});
       $PG_info .= " -use_sp" if($V_Href->{use_spmap});
       $PG_info .= " -use_sd" if($V_Href->{use_sdmap});
       $PG_info .= " -use_mm" if($V_Href->{use_multmap});
       $PG_info .= " -min_mq $V_Href->{min_mapQ}";
       $PG_info .= " -min_bq $V_Href->{min_baseQ}";
       $PG_info .= " -min_re $V_Href->{min_distToRedge}";
       $PG_info .= " -qual $V_Href->{baseQ_offset}";
    $HeadStr .= "$PG_info\n";

    # open split-bam file handles
    for my $splitBam (values %{$pairBamHref->{splitBam}}){
        $splitBam->start_write(samtools => $V_Href->{samtools});
        $splitBam->write(content => $HeadStr);
    }
}

#--- load reads OB to buff from source bam pair ---
# return file-end bool
sub loadrOBfromSourceBam{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $SourBamFH_Hf = $parm{SourBamFH_Hf};
    my $rBuff_Hf = $parm{rBuff_Hf};

    # read bam line
    my %bufAvailSize = map { ( $_, $V_Href->{rOB_AbufferSize} - scalar(@{$rBuff_Hf->{$_}}) ) } (1,2);
    for my $rEnd (1,2){
        my $fh = $SourBamFH_Hf->{$rEnd};
        while(    $bufAvailSize{$rEnd}--
               && ( my $bamLine = <$fh> )
        ){
            push @{$rBuff_Hf->{$rEnd}}, BioFuse::BioInfo::Objects::SeqData::HicReads_OB->new(ReadsLineText => $bamLine, reads_end => $rEnd, _rc_optfd => 1, _rc_rdstr => 1);
        }
    }
    # both bams finish reading?
    return ( $bufAvailSize{1} >= 0 && $bufAvailSize{2} >= 0 );
}

#--- move reads_OB from buff to pe_OB ---
sub rOB_to_peOB{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};
    my $rBuff_Hf = $parm{rBuff_Hf};
    my $peOB_Hf = $parm{peOB_Hf};
    my $fileEnd = $parm{fileEnd};
    my $peIdx_Sf = $parm{peIdx_Sf};
    my $peIdxInform_Sf = $parm{peIdxInform_Sf};

    my $bam_no = $pairBamHref->{no};
    # get last PE-id
    my %LastPid = map { ( ($fileEnd ? '__NA__' : $rBuff_Hf->{$_}->[-1]->pid), 1 ) } (1,2);
    # get shared PE-id
    my %R1_Pid  = map { ($_,1) } grep { !exists $LastPid{$_}                       } map { $_->pid } @{$rBuff_Hf->{1}};
    my %Sha_Pid = map { ($_,1) } grep { !exists $LastPid{$_} && exists $R1_Pid{$_} } map { $_->pid } @{$rBuff_Hf->{2}};
    undef %R1_Pid; # sweep
    # load PE-object
    my %pid2idx;
    for my $rEnd (1,2){
        my $shift_time = 0;
        for my $rOB ( @{$rBuff_Hf->{$rEnd}} ){
            my $pid = $rOB->pid;
            if(    !exists $Sha_Pid{$pid}
                ||  exists $LastPid{$pid}
            ){
                last;
            }
            $shift_time ++;
            # pe object
            if( !exists $pid2idx{$pid} ){
                $$peIdx_Sf++;
                $pid2idx{$pid} = $$peIdx_Sf;
                $peOB_Hf->{$$peIdx_Sf} = BioFuse::BioInfo::Objects::SeqData::HicPairEnd_OB->new( peIdx => "$bam_no-$$peIdx_Sf" );
            }
            # load rOB to pe_OB
            $peOB_Hf->{$pid2idx{$pid}}->load_reads_OB( reads_OB => $rOB );
        }
        # release space
        splice(@{$rBuff_Hf->{$rEnd}}, 0, $shift_time);
    }
    # inform
    if($$peIdx_Sf >= $$peIdxInform_Sf){
        stout_and_sterr "[INFO]\t".`date`
                             ."\tload $$peIdxInform_Sf PE-reads from '$pairBamHref->{prefix}' bam files.\n";
        $$peIdxInform_Sf += $V_Href->{peC_ReportUnit};
    }
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
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB = $parm{pe_OB};
    my $pairBamHref = $parm{pairBamHref};

    # skip this pe if not match required chr-[pair]
    ## this is preliminary filtration
    my @refSegHref_para;
    if(    exists $V_Href->{mustMapChr_i}
        || exists $V_Href->{mustMapChr_j}
    ){
        push @refSegHref_para, (a_refSegHref => $V_Href->{mustMapChr_i}->{Hf}) if exists $V_Href->{mustMapChr_i};
        push @refSegHref_para, (b_refSegHref => $V_Href->{mustMapChr_j}->{Hf}) if exists $V_Href->{mustMapChr_j};
        return unless $pe_OB->test_pair_RefSeg(@refSegHref_para);
    }

    # make judgement on each reads-end
    my %hap2phMut;
    &judge_on_rEnd(pe_OB => $pe_OB, hap2phMutHref => \%hap2phMut);

    # judge PE-reads
    my %Match_haplo;
    my ($SplitBam_tag, $Split_marker) = &judge_on_PE(pe_OB => $pe_OB, MatchHap_Href => \%Match_haplo);

    # skip this pe if Hi-C links not match required chr-[pair]
    ## this is thorough filtration
    if(    exists $V_Href->{mustMapChr_i}
        || exists $V_Href->{mustMapChr_j}
    ){
        return unless $pe_OB->testLinkRefSeg(@refSegHref_para, chrSortHref => $V_Href->{ChrThings}, chrSortKey  => 'turn');
    }

    # add OWN 'XH:Z:' haplo-id to 'optfd' of each reads_OB
    if( $SplitBam_tag =~ /^phMut-/ ){
        $pe_OB->addHapIDtoReadsOptfd(reads_end => $_) for (1,2);
    }

    # write PE to split-bam files
    my $peSAM_Aref = $pe_OB->printSAM(keep_all=>1);
    ## 'XS:Z:' denotes marker of 'S'plit reason
    $_ .= "\tXS:Z:$Split_marker" for @$peSAM_Aref;
    ## add phMut supporting haplotype
    if( $SplitBam_tag =~ /^phMut-/ ){ # 'dEnd-h' and 'sEnd-h'
        my @all_hapInfo;
        for my $hapID ( sort keys %Match_haplo ){
            my @hapInfo;
            for my $phMutNO (sort {$a<=>$b} keys %{$hap2phMut{$hapID}}){
                my $phMut_OB = $hap2phMut{$hapID}->{$phMutNO};
                # phasedMut info
                push @hapInfo, join(',', $phMut_OB->get_chr, $phMut_OB->get_pos, @{$phMut_OB->get_hapInfo(haplo_ID=>$hapID)});
            }
            push @all_hapInfo, "$hapID(".join(';',@hapInfo).")";
        }
        # 'XM:Z:' denotes the phased'M'ut used for haplotype detection
        $_ .= "\tXM:Z:".join('|',@all_hapInfo) for @$peSAM_Aref;
    }
    ## output
    $pairBamHref->{splitBam}->{$SplitBam_tag}->write(content => join("\n",@$peSAM_Aref)."\n");
    ## stat
    $SplitBam_tag = $Split_marker if $SplitBam_tag eq 'invalid';
    $V_Href->{PEsplitStat}->{$SplitBam_tag} ++;
}

#--- judgement on on each reads-end ---
sub judge_on_rEnd{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB = $parm{pe_OB};
    my $hap2phMutHref = $parm{hap2phMutHref};

    # discard abnormal supplementary alignment
    if( $pe_OB->discardAbnormalSP ){
        $V_Href->{PEsplitStat}->{remove_SP} ++;
    }

    # judgement on each end's each alignment
    for my $rEnd (1,2){
        my $reads_OB_Aref = $pe_OB->rOB_Af(reads_end => $rEnd);
        for my $reads_OB ( @$reads_OB_Aref ){
            # qulifying on alignment
            my $alignJudge = &judge_reads_alignment(reads_OB => $reads_OB);
            # record
            $reads_OB->load_AlignJudge(alignJudge => $alignJudge);
            # skip bad scenarios (UM/SD/SP/MM/LM) #
            if(     $alignJudge =~ /UM,/ # UnMapped
                ||  $alignJudge =~ /LM,/ # Low MapQ
                || ($alignJudge =~ /SD,/ && !$V_Href->{use_sdmap}) # Secondary alignment, if donot want
                || ($alignJudge =~ /SP,/ && !$V_Href->{use_spmap}) # Supplementary alignment, if donot want
                || ($alignJudge =~ /MM,/ && !$V_Href->{use_multmap}) # Multiple-Mapped, if donot want
            ){
                last if $V_Href->{filter_pe};
                next;
            }
            # if good, check haplotype covered of this alignment
            ## get phased Mut(s) this reads covers
            my $phMut_OB_Aref = PosToPhasedMut(chr => $reads_OB->mseg, pos => $reads_OB->mpos, ext => $reads_OB->mRefLen - 1);
            # no phased-mut covered #
            if( scalar(@$phMut_OB_Aref) == 0 ){
                $reads_OB->load_AlignJudge( alignJudge => 'UK,' ); # haplotype Un-Known
                next;
            }
            # get phased Mut(s) alleles from this reads
            my $phMut_covSign = 0;
            for my $phMut_OB (@$phMut_OB_Aref){
                # see whether reads_OB match against phMut_OB
                my $jdA = &judge_reads_by_phMut(reads_OB => $reads_OB, phMut_OB => $phMut_OB);
                if( $jdA->[0] ne 'MISS' ){
                    $reads_OB->load_AlignJudge(alignJudge => $jdA->[0]);
                    $phMut_covSign = 1;
                    # link phMut_OB to haplotype-ID
                    if( $jdA->[0] =~ /^h\d/ ){
                        (my $hapID = $jdA->[0]) =~ s/,$//;
                        $hap2phMutHref->{$hapID}->{ $phMut_OB->get_phMutNO } = $phMut_OB;
                        # record haplo+alleleInfo to this reads_OB
                        $reads_OB->load_SuppHaplo(hapID => $hapID, allele_OB => $jdA->[1]);
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
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
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
        if( $reads_OB->mapQ < $V_Href->{min_mapQ} ){
            $alignJudge .= 'LM,';
        }
    }
    # return
    return $alignJudge;
}

#--- see whether reads_OB match against phMut_OB ---
sub judge_reads_by_phMut{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $reads_OB = $parm{reads_OB};
    my $phMut_OB = $parm{phMut_OB};

    # get allele info of this position from this reads_OB
    my @parm = (chr => $phMut_OB->get_chr, pos => $phMut_OB->get_pos, refBase => $phMut_OB->get_refAllele, no_warn => 1);
    push @parm, (preferINDEL => 1 ) if( $phMut_OB->has_indel);
    my $allele_OB = $reads_OB->get_pos_allele(@parm);

    # not cover this phMut pos #
    if( $allele_OB->get_type eq 'miss' ){ 
        return ['MISS']; # reads_OB miss the phMut
    }
    #----------------------------#
    # low base quality filtering #
    #----------------------------#
    # get mean base quality
    my $alleleQual = $allele_OB->get_qual(Q_offset => $V_Href->{baseQ_offset});
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
    my $haplo_ID = $phMut_OB->allele2haploID(allele_OB => $allele_OB);
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
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
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
    if(    $rEndJdgHref->{1}->{J} =~ /UM,/
        || $rEndJdgHref->{2}->{J} =~ /UM,/
    ){
        $SplitBam_tag = 'discarded';
        $Split_marker = 'Unmap';
    }
    ## t5, has supplementary alignment, one or both
    if(   !$V_Href->{use_spmap}
       && (   $rEndJdgHref->{1}->{J} =~ /SP,/
           || $rEndJdgHref->{2}->{J} =~ /SP,/
          )
    ){
        # no needs to make primary-align
        if(     $V_Href->{filter_pe}
            || !$pe_OB->tryDiscardAlign(type => 'SP')
        ){
            $SplitBam_tag = 'discarded';
            $Split_marker = 'SuppleAlign';
        }
    }
    ## t5, has secondary alignment, one or both
    if(   !$V_Href->{use_sdmap}
       && (   $rEndJdgHref->{1}->{J} =~ /SD,/
           || $rEndJdgHref->{2}->{J} =~ /SD,/
          )
    ){
        # no needs to make primary-align
        if(     $V_Href->{filter_pe}
            || !$pe_OB->tryDiscardAlign(type => 'SD')
        ){
            $SplitBam_tag = 'discarded';
            $Split_marker = 'SecondAlign';
        }
    }
    ## t5, has multiple mapped, one or both
    if(   !$V_Href->{use_multmap}
       && (   $rEndJdgHref->{1}->{J} =~ /MM,/
           || $rEndJdgHref->{2}->{J} =~ /MM,/
          )
    ){
        # see whether R1 and R2 have remained alignments after discarding MM-align
        # if not both have, go to 't5'
        # else, do discard MM-align, and then make primary-align if lacks
        if(     $V_Href->{filter_pe}
            || !$pe_OB->tryDiscardAlign(type => 'MM')
        ){
            $SplitBam_tag = 'discarded';
            $Split_marker = 'MultiMap';
        }
        else{
            $pe_OB->makePrimeAlignment;
        }
    }
    ## t5, has Low Mapping quality, one or both
    if(   $rEndJdgHref->{1}->{J} =~ /LM,/
       || $rEndJdgHref->{2}->{J} =~ /LM,/
    ){
        # see whether R1 and R2 have remained alignments after discarding LM-align
        # if not both have, go to 't5'
        # else, do discard LM-align, and then make primary-align if lacks
        if(     $V_Href->{filter_pe}
            || !$pe_OB->tryDiscardAlign(type => 'LM', min_mapQ => $V_Href->{min_mapQ})
        ){
            $SplitBam_tag = 'discarded';
            $Split_marker = 'LowMapQ';
        }
        else{
            $pe_OB->makePrimeAlignment;
        }
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
    ## t7, invalid Hi-C PE
    my $ivdPEjudge = $pe_OB->isInValidPair(chr2enzymePosHf => $V_Href->{chr2enzymePos},
                                           maxCloseAlignDist => $V_Href->{maxCloseAlignDist},
                                           skipCloseAlign  => $V_Href->{skipCloseAlignFromInvalid});
    if($ivdPEjudge){
        $SplitBam_tag = 'invalid';
        $Split_marker = "invalid-$ivdPEjudge";
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
            my $dEndKeepBool = $pe_OB->dEndSameHapJudge(maxCloseAlignDist => $V_Href->{maxCloseAlignDist});
            if($dEndKeepBool){
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
        ## OR, phMut-dEnd-hx
        elsif( $readsHapAmt{1} + $readsHapAmt{2} == 1 ){ # either of R1 and R2 has only one hap
            my $hapID = $rSoloHapID{1}->[0] || $rSoloHapID{2}->[0];
            $MatchHap_Href->{$hapID} = 1;
            # the has-solo-hap rEnd should meet:
            ## 1)     cannot have two alignments supporting the same hapID
            ## 2) OR, the two alignments supporting the same hapID are not close aligned
            my $rEndHasHap = $readsHapAmt{1} ? 1 : 2;
            my $sEndKeepBool = $pe_OB->sEndSoloHapJudge(shEnd => $rEndHasHap, maxCloseAlignDist => $V_Href->{maxCloseAlignDist});
            if($sEndKeepBool){
                $SplitBam_tag = "phMut-sEnd-$hapID";
                $Split_marker = "sEndSoloHap;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
            }
            else{
                $SplitBam_tag = "phMut-dEnd-$hapID";
                $Split_marker = "sEndSoloHapFarSame;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
                # make ukEnd loads this hapID
                my $rEndNonHap = $readsHapAmt{1} ? 2 : 1;
                $_->load_SuppHaplo(hapID => $hapID, allele_OB => 'NULL') for @{$pe_OB->rOB_Af(reads_end => $rEndNonHap)};
            }
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
            my $doSelectOneHapBool = $pe_OB->sEndInterHapJudge(mhEnd => $mhEnd, uhEnd => $uhEnd, maxCloseAlignDist => $V_Href->{maxCloseAlignDist});
            if( ! $doSelectOneHapBool ){
                $SplitBam_tag = 'phMut-sEnd-hInter';
                $Split_marker = "sEndMultHap;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
                $MatchHap_Href->{$_} = 1 for keys %{$reads2HapID{1}};
                $MatchHap_Href->{$_} = 1 for keys %{$reads2HapID{2}};
            }
            else{
                my $selectHapID = &selectOneFromMultiHap(hapAlleleHf_Aref => [$rEndJdgHref->{$mhEnd}->{H}]);
                $SplitBam_tag = "phMut-sEnd-$selectHapID";
                $rEndJdgHref->{$mhEnd}->{J} =~ s/\b$selectHapID\b/$selectHapID(S)/g; # add '(S)' to mark 'Selected hapID'
                $Split_marker = "sEndSoloHapSelect;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
                $MatchHap_Href->{$selectHapID} = 1;
                # dislink other haplo(s) with all reads_OBs of this pe_OB
                $_->onlyKeep_SuppHaplo( hapID => $selectHapID ) for @{ $pe_OB->rOB_Af( reads_end => $mhEnd ) };
            }
        }
        ## t2, phMut-dEnd-hInter, 7 cases remaining
        ## OR, phMut-xEnd-hx
        else{
            # two PE-ends should meet:
            ## 1) alignments support same haplo  should be close;
            ## 2) alignments support diff haplos should be not close;
            my $doSelectOneHapBool = $pe_OB->dEndInterHapJudge(maxCloseAlignDist => $V_Href->{maxCloseAlignDist});
            if( ! $doSelectOneHapBool ){
                $SplitBam_tag = 'phMut-dEnd-hInter';
                $Split_marker = "dEndMultHap;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
                $MatchHap_Href->{$_} = 1 for keys %{$reads2HapID{1}};
                $MatchHap_Href->{$_} = 1 for keys %{$reads2HapID{2}};
            }
            else{
                my $selectHapID = &selectOneFromMultiHap(hapAlleleHf_Aref => [$rEndJdgHref->{1}->{H},$rEndJdgHref->{2}->{H}]);
                my $newEndType = (    exists $rEndJdgHref->{1}->{H}->{$selectHapID}
                                   && exists $rEndJdgHref->{2}->{H}->{$selectHapID}
                                 ) ? 'dEnd' : 'sEnd';
                $SplitBam_tag = "phMut-$newEndType-$selectHapID";
                $rEndJdgHref->{$_}->{J} =~ s/\b$selectHapID\b/$selectHapID(S)/g for (1,2); # add '(S)' to mark 'Selected hapID'
                $Split_marker = "${newEndType}SoloHapSelect;R1:$rEndJdgHref->{1}->{J};R2:$rEndJdgHref->{2}->{J}";
                $MatchHap_Href->{$selectHapID} = 1;
                # dislink other haplo(s) with all reads_OBs of this pe_OB
                $_->onlyKeep_SuppHaplo( hapID => $selectHapID ) for @{ $pe_OB->rOB_Af( reads_end => 1 ) };
                $_->onlyKeep_SuppHaplo( hapID => $selectHapID ) for @{ $pe_OB->rOB_Af( reads_end => 2 ) };
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
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $hapAlleleHf_Aref = $parm{hapAlleleHf_Aref};

    # calculate score of each haolo-allele
    my %score;
    for my $hapAllele_Href ( @$hapAlleleHf_Aref ){
        for my $hapID ( keys %$hapAllele_Href ){
            for my $allele_OB ( @{ $hapAllele_Href->{$hapID} } ){
                push @{$score{$hapID}}, $allele_OB->get_score(Q_offset => $V_Href->{baseQ_offset});
            }
        }
    }

    # compare count
    my %cnt2hap;
    push @{$cnt2hap{scalar(@{$score{$_}})}}, $_ for keys %score;
    my ($maxCnt) = max(keys %cnt2hap);
    if(scalar(@{$cnt2hap{$maxCnt}}) == 1){
        return $cnt2hap{$maxCnt}->[0];
    }

    # compare score
    my $maxMeanScore = -1;
    my $selectHapID;
    for my $hapID ( @{$cnt2hap{$maxCnt}} ){
        my $meanScore = sum(@{$score{$hapID}}) / $maxCnt;
        if($meanScore > $maxMeanScore){
            $maxMeanScore = $meanScore;
            $selectHapID = $hapID;
        }
    }
    return $selectHapID;
}

#--- sort split.bam for next local region haplo-division ---
## t3(phMut-sEnd-hx) and t6(unknown)
sub sortSplitBamByPEcontact{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};

    ## release memory
    ## no matter in main-process or child-process
    $V_Href->{PhasedMut} = undef;
    ## reload in HaploHiC::PhasedHiC::dumpContacts if need
    $V_Href->{chr2enzymePos} = undef;

    # series of split-bam files with different tags, all contribute to contacts
    ## t1,t3: <haplo> bam of phMut (single-haplo) covered reads
    my @tags = map {("phMut-dEnd-h$_" , "phMut-sEnd-h$_")} (1 .. $V_Href->{haploCount});
    ## t2,t4: <haplo> bam of phMut (multi-haplo) covered reads
    push @tags, "phMut-dEnd-hInter", "phMut-sEnd-hInter";
    ## t6: bam of haplo-unknown reads
    push @tags, "unknown";

    for my $tag (@tags){
        my $splitBam = $pairBamHref->{splitBam}->{$tag};
        # temp workspace
        my $chrPairDir = $splitBam->filepath . '-chrPairDir';
        `rm -rf $chrPairDir` if -d $chrPairDir;
        `mkdir -p $chrPairDir`;
        # read $tag split bam, and write peOB to chr-pair bams
        my %chrPairBam;
        my @subrtOpt = (subrtRef => \&write_peOB_to_chrPairBam,
                        subrtParmAref => [chrPairBamHf => \%chrPairBam, splitBam => $splitBam, chrPairDir => $chrPairDir]);
        $splitBam->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC', quiet => 1, deal_peOB_pool => 1, @subrtOpt);
        # close chr-pair bams
        $_->stop_write for values %chrPairBam;
        # load chr-pair bam and sort corrdinates to splitBam
        &chrPairBamToSortSplitBam(chrPairBamHf => \%chrPairBam, splitBam => $splitBam);
        # sweep
        `rm -rf $chrPairDir`;
        # inform
        my $mark = $splitBam->tag;
        stout_and_sterr "[INFO]\t".`date`
                             ."\tsort $mark bam OK.\n";
    }
}

#--- put given peOB to its chr-pair bam ---
sub write_peOB_to_chrPairBam{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};
    my $chrPairBamHf = $parm{chrPairBamHf};
    my $splitBam = $parm{splitBam};
    my $chrPairDir = $parm{chrPairDir};
    my $sorted = $parm{sorted} || 0;

    my $tag = $splitBam->tag;
    my $useUKbool = ($tag !~ /-dEnd-h/ && $tag !~ /-hInter/); # sEnd-hx, unknown
    my $last_chrPairTag = undef;
    for my $pe_OB (@$pe_OB_poolAf){
        # arranged rOB, turn: R1:5p, R1:3p, R2:3p, R2:5p
        my @rOB = grep $useUKbool || $_->printSAM !~ /\sXH:Z:UK\s/, @{ $pe_OB->arranged_rOB_Af };
        # idx for chrPairTag, sorted by chr(ref-dict),pos(ascending)
        my @idx = &getIdxOfrOBforChrPair(tag => $tag, rOB_arngAf => \@rOB);
        my $chrPairTag = join(',', $rOB[$idx[0]]->mseg, $rOB[$idx[1]]->mseg); # chromosome name barely matches ','
        # check chr-pair bam
        unless(exists $chrPairBamHf->{$chrPairTag}){
            my $chrPairFileTag = join('-', $rOB[$idx[0]]->mseg, $rOB[$idx[1]]->mseg);
            $chrPairBamHf->{$chrPairTag} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => catfile($chrPairDir, "$chrPairFileTag.bam"));
            $chrPairBamHf->{$chrPairTag}->start_write(samtools => $V_Href->{samtools});
            $chrPairBamHf->{$chrPairTag}->write(content => $_) for @{ $splitBam->header_Af(samtools => $V_Href->{samtools}) };
            # deal with last chr-pair if sorted
            $chrPairBamHf->{$last_chrPairTag}->stop_write if $sorted && defined $last_chrPairTag;
        }
        # write to chr-pair bam
        $chrPairBamHf->{$chrPairTag}->write(content => join("\n",@{$pe_OB->printSAM(keep_all=>1)})."\n");
        # record
        $last_chrPairTag = $chrPairTag;
    }
}

#--- prepare index of rOB for chrPairTag ---
# take the bilateral idx (i.e., the first and last) AMSP
# the turn follows the chr(ref-dict),pos(ascending)
sub getIdxOfrOBforChrPair{
    # options
    #shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $tag = $parm{tag};
    my $rOB_arngAf = $parm{rOB_arngAf};

    my @idx = (0, -1);
    if($tag =~ /sEnd-h\d/){ # sEnd-hx
        # recover SuppHaplo attribute
        $_->recover_SuppHaploAttr for @$rOB_arngAf;
        # make sure hap and non-hap
        @idx = &findContactAnchor(rOB_arngAf => $rOB_arngAf);
    }
    # sort by chr+pos for chrPairTag (used for file name)
    return sort {
               sortByStrAndSubNum(
                   Str_a => $rOB_arngAf->[$a]->mseg, Num_a => $rOB_arngAf->[$a]->mpos,
                   Str_b => $rOB_arngAf->[$b]->mseg, Num_b => $rOB_arngAf->[$b]->mpos,
                   SortHref => $V_Href->{ChrThings},
                   SortKey  => 'turn'
               )
           } @idx;
}

#--- find the index in rOB_sortArray of hasHap and nonHap reads_OB for sEnd-hx ---
sub findContactAnchor{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $rOB_arngAf = $parm{rOB_arngAf};

    my (@hasHapIdx, @nonHapIdx);
    for my $i (0 .. scalar(@$rOB_arngAf)-1){
        if($rOB_arngAf->[$i]->has_SuppHaplo){
            push @hasHapIdx, $i;
        }
        else{
            push @nonHapIdx, $i;
        }
    }
    # find proper index of hasHap_rOB and nonHap_rOB
    if( $hasHapIdx[-1] < $nonHapIdx[0] ){
        return ($hasHapIdx[0], $nonHapIdx[-1]);
    }
    elsif( $nonHapIdx[-1] < $hasHapIdx[0] ){
        return ($hasHapIdx[-1], $nonHapIdx[0]);
    }
    else{
        my $hasHaprOB_F = $rOB_arngAf->[$hasHapIdx[0]];
        my $hasHaprOB_L = $rOB_arngAf->[$hasHapIdx[-1]];
        my $nonHaprOB_F = $rOB_arngAf->[$nonHapIdx[0]];
        my $nonHaprOB_L = $rOB_arngAf->[$nonHapIdx[-1]];
        my $hFnLsameChr = ($hasHaprOB_F->mseg eq $nonHaprOB_L->mseg);
        my $hLnFsameChr = ($hasHaprOB_L->mseg eq $nonHaprOB_F->mseg);
        # same chr
        if( $hFnLsameChr && $hLnFsameChr ){
            # larger pos gap
            if(   abs($hasHaprOB_F->mpos - $nonHaprOB_L->mpos)
                > abs($hasHaprOB_L->mpos - $nonHaprOB_F->mpos)
            ){
                return ($hasHapIdx[0], $nonHapIdx[-1]);
            }
            else{
                return ($hasHapIdx[-1], $nonHapIdx[0]);
            }
        }
        # diff chr (solo)
        elsif( $hFnLsameChr && !$hLnFsameChr ){
            return ($hasHapIdx[-1], $nonHapIdx[0]);
        }
        elsif( !$hFnLsameChr && $hLnFsameChr ){
            return ($hasHapIdx[0], $nonHapIdx[-1]);
        }
        # diff chr (both)
        else{
            # longer mapped length
            if(   $hasHaprOB_F->mReadLen + $nonHaprOB_L->mReadLen
                > $hasHaprOB_L->mReadLen + $nonHaprOB_F->mReadLen
            ){
                return ($hasHapIdx[0], $nonHapIdx[-1]);
            }
            else{
                return ($hasHapIdx[-1], $nonHapIdx[0]);
            }
        }
        # if( $nonHaprOB_F->is_suppmap ){ # avoid supp-map nonHap alignment, deprecated
        #     return ($hasHapIdx[0], $nonHapIdx[-1]);
        # }
        # else{
        #     return ($hasHapIdx[-1], $nonHapIdx[0]);
        # }
    }
}

#--- sort chrPair bams and write to new split bam ---
sub chrPairBamToSortSplitBam{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $splitBam = $parm{splitBam};
    my $chrPairBamHf = $parm{chrPairBamHf};

    # start to write new split bam
    my $SAMheadAf = $splitBam->header_Af(samtools => $V_Href->{samtools});
    $splitBam->start_write(samtools => $V_Href->{samtools});
    $splitBam->write(content => $_) for @$SAMheadAf;
    # sort each chrPair bam and write to new split bam
    my $AbufferSize = $V_Href->{chrPairSort_peOB_AbufferSize};
    my @subrtOpt = (subrtRef => \&write_sort_peOB_to_newSplitBam, subrtParmAref => [splitBam => $splitBam]);
    my $chrPairTagAf = &getChrPairTagAf;
    ## each chrPair
    my $tag = $splitBam->tag;
    for my $chrPairTag (@$chrPairTagAf){
        next unless exists $chrPairBamHf->{$chrPairTag};
        $chrPairBamHf->{$chrPairTag}->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC',
                                                      mark => "$tag $chrPairTag", quiet => 1,
                                                      peOB_AbufferSize => $AbufferSize, deal_peOB_pool => 1, @subrtOpt);
    }
    # close new split bam
    $splitBam->stop_write;
}

#--- get ref of chrPair array ---
sub getChrPairTagAf{
    my @chrPairTag;
    ## intra-chr
    push @chrPairTag, "$_,$_" for @{$V_Href->{sortedChr}};
    ## inter-chr
    my $chrCount = scalar @{$V_Href->{sortedChr}};
    for my $i (0 .. $chrCount-1){
        my $chr_a = $V_Href->{sortedChr}->[$i];
        for my $j ($i+1 .. $chrCount-1){
            my $chr_b = $V_Href->{sortedChr}->[$j];
            push @chrPairTag, "$chr_a,$chr_b";
        }
    }
    # return array-ref
    return \@chrPairTag;
}

#--- sort pe-OB (pool) and write new split bam ---
# sort by pos ONLY in this pe_pool (related to 'chrPairSort_peOB_AbufferSize')
sub write_sort_peOB_to_newSplitBam{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};
    my $splitBam = $parm{splitBam};

    my $tag = $splitBam->tag;
    my $useUKbool = ($tag !~ /-dEnd-h/ && $tag !~ /-hInter/); # sEnd-hx, unknown
    # map idx to sort_rOB
    my @peCoord;
    for my $i (0 .. scalar(@$pe_OB_poolAf)-1){
        # arranged rOB, turn: R1:5p, R1:3p, R2:3p, R2:5p
        my @rOB = grep $useUKbool || $_->printSAM !~ /\sXH:Z:UK\s/, @{ $pe_OB_poolAf->[$i]->arranged_rOB_Af };
        # idx for chrPairTag, sorted by chr(ref-dict),pos(ascending)
        my @idx = &getIdxOfrOBforChrPair(tag => $tag, rOB_arngAf => \@rOB);
        push @peCoord, [$i, $rOB[$idx[0]]->mpos, $rOB[$idx[1]]->mpos];
    }
    # sort by pos and write to new split bam
    $splitBam->write(content => join("\n",@{$pe_OB_poolAf->[$_->[0]]->printSAM(keep_all=>1)})."\n") for sort {$a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]} @peCoord;
}

#--- 
1; ## tell the perl script the successful access of this module.
