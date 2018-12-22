package HaploHiC::PhasedHiC::phasedPEtoContact;

use strict;
use warnings;
use List::Util qw/ max min sum first /;
use Data::Dumper;
use Parallel::ForkManager;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Index qw/ Pos2Idx /;
use BioFuse::Util::Sort qw/ sortByStrAndSubNum /;
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
              phasePE_to_contactCount
              smartBam_PEread
              load_phasedPE_contacts
              phasePE_contacts_to_count
              get_rOBpair_HapLinkCount
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::phasedPEtoContact';
#----- version --------
$VERSION = "0.07";
$DATE = '2018-11-16';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        phasePE_to_contactCount
                        smartBam_PEread
                        load_phasedPE_contacts
                        mPosToWinIdx
                        phasePE_contacts_to_count
                        get_rOBpair_HapLinkCount
                     /;

#--- load phased PE from given bam files ---
## load PE's contact => [de-dup PE] => get contact count
sub phasePE_to_contactCount{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $tagToBamHref = $parm{tagToBamHref};

    # load reads from each bam
    for my $tag ( sort keys %$tagToBamHref ){
        for my $Aref ( @{$tagToBamHref->{$tag}} ){
            my ($bamPrefix, $hapSplitBam) = @$Aref;
            # read phased bam
            ## FLAG: 0x100(sd), 0x400(d), 0x800(sp)
            ## WEDO: -F 0x400(d), as 'sd' and 'sp' alignments is important in HiC
            &smartBam_PEread(bam => $hapSplitBam, mark => "'$bamPrefix' $tag", viewOpt => '-F 0x400',
                             subrtRef => \&load_phasedPE_contacts,
                             subrtParmAref => [idxFunc=>\&mPosToWinIdx]);
        }
        # contacts to count
        ## here, do de-dup phased reads (optional)
        &phasePE_contacts_to_count(tag => $tag);
    }
}

#--- read smartPE BAM and do something on each pe_OB ---
sub smartBam_PEread{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $bam = $parm{bam};
    my $mark = $parm{mark} || '';
    my $viewOpt = $parm{viewOpt} || '';
    my $subrtRef = $parm{subrtRef};
    my $subrtParmAref = $parm{subrtParmAref} || [];
    my $peOB_AbufferSize = $parm{peOB_AbufferSize} || 1E5;

    # check existence
    warn_and_exit "<ERROR>\tCannot find $mark bam: $bam\n" unless file_exist(filePath => $bam);
    # read bam
    open (BAM,"$V_Href->{samtools} view $viewOpt $bam |") || die "fail read $mark bam: $!\n";
    my $pe_Count = 0;
    my $pe_CountInform = $V_Href->{peC_ReportUnit};
    my $last_pid = '';
    my $last_peOB = undef;
    my @peOB_pool = ();
    while(<BAM>){
        my $reads_OB = BioFuse::BioInfo::Objects::HicReads_OB->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1 );
        if( $reads_OB->get_pid ne $last_pid ){
            # deal last pe_OB
            $pe_Count++;
            push @peOB_pool, $last_peOB if( defined $last_peOB );
            if( $pe_Count % $peOB_AbufferSize == 0 ){
                &{$subrtRef}(pe_OB => $_, @$subrtParmAref) for @peOB_pool;
                @peOB_pool = ();
                # inform
                if($pe_Count >= $pe_CountInform){
                    stout_and_sterr "[INFO]\t".`date`
                                         ."\tload $pe_CountInform PE-reads from $mark bam.\n";
                    $pe_CountInform += $V_Href->{peC_ReportUnit};
                }
            }
            # create new pe_OB
            $last_peOB = BioFuse::BioInfo::Objects::HicPairEnd_OB->new;
            # update
            $last_pid = $reads_OB->get_pid;
        }
        # load this reads_OB to pe_OB
        $last_peOB->load_reads_OB(reads_OB => $reads_OB);
    }
    close BAM;
    # deal the last one
    if( defined $last_peOB ){
        push @peOB_pool, $last_peOB;
        &{$subrtRef}(pe_OB => $_, @$subrtParmAref) for @peOB_pool;
        @peOB_pool = ();
    }
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tTotally, load $pe_Count PE-reads from $mark bam.\n";
}

#--- load contacts of phased PE-reads ---
sub load_phasedPE_contacts{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB = $parm{pe_OB};
    my $idxFunc = $parm{idxFunc};
    # my $hapSort = $parm{hapSort} || 0; # currently, deprecated

    # get chr-pos ascending sorted all mapped reads_OB
    my $rOB_sortAref = $pe_OB->get_sorted_reads_OB(rEndAref => [1,2], onlyMap => 1,
                                                   chrSortHref => $V_Href->{ChrThings},
                                                   chrSortKey  => 'turn');
    # recover SuppHaplo attribute
    $_->recover_SuppHaploAttr for @$rOB_sortAref;
    # extract haplotype confirmed alignments
    my @hasHaprOB = grep $_->has_SuppHaplo, @$rOB_sortAref;
    # double check
    if( scalar(@hasHaprOB) < 2 ){
        warn_and_exit "<ERROR>\tPhased PE-reads doesn't have at least two haplotype confirmed alignments.\n".Dumper($pe_OB);
    }
    # skip hapComb set by flanking region lacking phased contacts, see func 'get_rOBpair_HapLinkCount'
    if(    $hasHaprOB[ 0]->is_fromUnPhasedRegRand
        || $hasHaprOB[-1]->is_fromUnPhasedRegRand
    ){
        return; # do not take this PE into contact counting
    }
    # as sorted, so take first[0] and last[-1] one as index of contacted region
    my (%chr, %pIdx, %hapID);
    $chr{a}   = $hasHaprOB[ 0]->get_mseg;
    $chr{b}   = $hasHaprOB[-1]->get_mseg;
    return if(!exists $V_Href->{ChrThings}->{$chr{a}} || !exists $V_Href->{ChrThings}->{$chr{b}});
    $pIdx{a}  = &{$idxFunc}(chr => $chr{a}, pos => $hasHaprOB[ 0]->get_mpos);
    $pIdx{b}  = &{$idxFunc}(chr => $chr{b}, pos => $hasHaprOB[-1]->get_mpos);
    $hapID{a} = $hasHaprOB[ 0]->get_SuppHaploStr;
    $hapID{b} = $hasHaprOB[-1]->get_SuppHaploStr;
    # hapID sort
    my ($k1,$k2) = qw/ a b /;
    # ($k1,$k2) = sort {$hapID{$a} cmp $hapID{$b}} ($k1,$k2) if ($hapSort && $hapID{a} ne $hapID{b}); # deprecated
    # get PE-Info string
    my $peInfoStr = join(';', map {( join(',', $_->get_mseg, $_->get_mpos) )} @$rOB_sortAref);
    # record, 'peInfoStr' as key, 'accumulated count' as value
    $V_Href->{phasePEdetails}->{$chr{$k1}}->{$chr{$k2}}->{$pIdx{$k1}}->{$pIdx{$k2}}->{"$hapID{$k1},$hapID{$k2}"}->{$peInfoStr} ++;
}

#--- get window index of mapping position ---
sub mPosToWinIdx{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    return Pos2Idx(pos => $parm{pos}, winSize => $V_Href->{mapPosWinSize});
}

#--- convert PE-contacts to count ---
## note here is to accumulate counts to 'phasePEcontact' container
sub phasePE_contacts_to_count{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $tag = $parm{tag};

    for my $chr_a (sort keys %{$V_Href->{phasePEdetails}}){
        for my $chr_b (sort keys %{$V_Href->{phasePEdetails}->{$chr_a}}){
            for my $pIdx_a (sort {$a<=>$b} keys %{$V_Href->{phasePEdetails}->{$chr_a}->{$chr_b}}){
                for my $pIdx_b (sort {$a<=>$b} keys %{$V_Href->{phasePEdetails}->{$chr_a}->{$chr_b}->{$pIdx_a}}){
                    my $hapCombHf = $V_Href->{phasePEdetails}->{$chr_a}->{$chr_b}->{$pIdx_a}->{$pIdx_b};
                    for my $hapComb (sort keys %$hapCombHf){
                        if( $V_Href->{SkipDeDupPhasedReads} ){
                            # count all PE
                            $V_Href->{phasePEcontact}->{$chr_a}->{$chr_b}->{$pIdx_a}->{$pIdx_b}->{$hapComb} += $_ for values %{$hapCombHf->{$hapComb}};
                        }
                        else{
                            # count all PE de-redundancy by 'peInfo'
                            $V_Href->{phasePEcontact}->{$chr_a}->{$chr_b}->{$pIdx_a}->{$pIdx_b}->{$hapComb} += scalar( keys %{$hapCombHf->{$hapComb}} );
                        }
                    }
                }
                # periodically release memory
                $V_Href->{phasePEdetails}->{$chr_a}->{$chr_b}->{$pIdx_a} = {};
            }
        }
        # periodically release memory
        $V_Href->{phasePEdetails}->{$chr_a} = {};
    }
    # finally release
    $V_Href->{phasePEdetails} = {};
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tDe-dup $tag phased PE-reads OK.\n" unless $V_Href->{SkipDeDupPhasedReads};
}

#--- return haplo link count of given reads_OB pair ---
sub get_rOBpair_HapLinkCount{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $skipSort = $parm{skipSort} || 0;
    my $hapRegex = $parm{hapRegex} || undef;

    my %HapLinkCount;
    # chr,pos ascending sort
    my %rOB;
    if( $skipSort ){
        ($rOB{a}, $rOB{b}) = ($parm{rOB_a}, $parm{rOB_b});
    }
    else{
        ($rOB{a}, $rOB{b}) = sort {
                                sortByStrAndSubNum(
                                    Str_a => $a->get_mseg, Num_a => $a->get_mpos,
                                    Str_b => $b->get_mseg, Num_b => $b->get_mpos,
                                    SortHref => $V_Href->{ChrThings},
                                    SortKey  => 'turn'
                                )
                             } ($parm{rOB_a}, $parm{rOB_b});
    }
    # get region pair
    my %mSeg = map { ($_, $rOB{$_}->get_mseg) } keys %rOB;
    my %mPos = map { ($_, $rOB{$_}->get_mpos) } keys %rOB;
    my %mExp = map { ($_, $mPos{$_}+$rOB{$_}->get_mRefLen-1) } keys %rOB;
    # first check
    if(    ! exists $V_Href->{phasePEcontact}->{$mSeg{a}}
        || ! exists $V_Href->{phasePEcontact}->{$mSeg{a}}->{$mSeg{b}}
    ){
        return (\%HapLinkCount, "LackChrLink($mSeg{a},$mSeg{b})");
    }
    # iteratively find sufficient phased Het-Mut in local flank region
    my $chrHapLinkHref = $V_Href->{phasePEcontact}->{$mSeg{a}}->{$mSeg{b}};
    my %fReg = map { ($_, {pos=>{},mct=>{}}) } keys %rOB;
    for my $time ( 0 .. $V_Href->{UKreadsFlankRegUnitMaxTimes} ){
        my $ext_ratio  = 2 ** $time;
        my $FlankSize  = $ext_ratio * $V_Href->{UKreadsFlankRegUnit};
        my $PhaHetMutC = $ext_ratio * $V_Href->{UKreadsMaxPhasedHetMut};
        for my $Ach (sort keys %rOB){
            # extract phased Het-Mut in 5/3 prime flanking region
            my $phMutOB_p5_Aref = PosToPhasedMut( chr => $mSeg{$Ach}, pos => $mPos{$Ach}, ext => -1 * $FlankSize  );
            my $phMutOB_p3_Aref = PosToPhasedMut( chr => $mSeg{$Ach}, pos => $mExp{$Ach}, ext => $FlankSize );
            # adjust phased Het-Mut count, if set
            my $phMutOB_p5_C = scalar(@$phMutOB_p5_Aref);
            my $phMutOB_p3_C = scalar(@$phMutOB_p3_Aref);
            if( $PhaHetMutC ){
                if( $phMutOB_p5_C > $PhaHetMutC ){
                    @{$phMutOB_p5_Aref} = @{$phMutOB_p5_Aref}[ -1*$PhaHetMutC .. -1 ];
                    $phMutOB_p5_C = $PhaHetMutC;
                }
                if( $phMutOB_p3_C > $PhaHetMutC ){
                    @{$phMutOB_p3_Aref} = @{$phMutOB_p3_Aref}[ 0 .. $PhaHetMutC-1 ];
                    $phMutOB_p3_C = $PhaHetMutC;
                }
            }
            # record
            $fReg{$Ach}{mct}{p5} = $phMutOB_p5_C;
            $fReg{$Ach}{mct}{p3} = $phMutOB_p3_C;
            $fReg{$Ach}{pos}{p5} = ( ( $PhaHetMutC && $phMutOB_p5_C == $PhaHetMutC ) ? $phMutOB_p5_Aref->[ 0]->get_pos : $mPos{$Ach}-$FlankSize );
            $fReg{$Ach}{pos}{p3} = ( ( $PhaHetMutC && $phMutOB_p3_C == $PhaHetMutC ) ? $phMutOB_p3_Aref->[-1]->get_pos : $mExp{$Ach}+$FlankSize );
        }
        # extract the haplo link of a/b paired flank regions
        my %pIdx;
        %{$pIdx{p5}} = map { ( $_, Pos2Idx(pos => $fReg{$_}{pos}{p5}, winSize => $V_Href->{mapPosWinSize}) ) } keys %fReg;
        %{$pIdx{p3}} = map { ( $_, Pos2Idx(pos => $fReg{$_}{pos}{p3}, winSize => $V_Href->{mapPosWinSize}) ) } keys %fReg;
        for my $pIdx_a ( $pIdx{p5}{a} .. $pIdx{p3}{a} ){
            next unless exists $chrHapLinkHref->{$pIdx_a};
            for my $pIdx_b ( $pIdx{p5}{b} .. $pIdx{p3}{b} ){
                next unless exists $chrHapLinkHref->{$pIdx_a}->{$pIdx_b};
                my $pIdxHapLinkHref = $chrHapLinkHref->{$pIdx_a}->{$pIdx_b};
                $HapLinkCount{$_} += $pIdxHapLinkHref->{$_} for keys %$pIdxHapLinkHref;
            }
        }
        # only keep pre-set hapID combination
        if(defined $hapRegex){
            delete $HapLinkCount{$_} for grep !/$hapRegex/, keys %HapLinkCount;
        }
        # find HapLink or last time
        my $hapCombCnt = scalar(keys %HapLinkCount);
        if(    $hapCombCnt
            || $time == $V_Href->{UKreadsFlankRegUnitMaxTimes}
        ){
            my $mark = $hapCombCnt ? 'RegionPhased' : 'RegionNotPhased';
            $mark .= ";($ext_ratio"
                    .",$mSeg{$_}:$fReg{$_}{pos}{p5}-$fReg{$_}{pos}{p5}"
                    .",MCT5p:$fReg{$_}{mct}{p5}"
                    .",MCT3p:$fReg{$_}{mct}{p3})"
                    for sort keys %fReg;
            return (\%HapLinkCount, $mark);
        }
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
