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

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              phasePE_to_contactCount
              load_phasedPE_contacts
              phasePE_contacts_to_count
              get_rOBpair_HapLinkCount
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::phasedPEtoContact';
#----- version --------
$VERSION = "0.15";
$DATE = '2019-03-10';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        phasePE_to_contactCount
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
    for my $tag (sort keys %$tagToBamHref){
        for my $hapSplitBam (@{$tagToBamHref->{$tag}}){
            my $mark = $hapSplitBam->get_tag;
            # read phased bam
            my @lastChrPair = ('__NA__', '__NA__', $mark); # takes $mark by the way
            my @subrtOpt = (subrtRef => \&load_phasedPE_contacts, subrtParmAref => [idxFunc => \&mPosToWinIdx, lastChrPairAf => \@lastChrPair, onlyPha => 1]);
            $hapSplitBam->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC', deal_peOB_pool => 1, @subrtOpt);
            # contacts to count for last chr-pair, release memory
            ## here, do de-dup phased reads (in single run) (optional)
            &phasePE_contacts_to_count(tag => "$mark " . join(',', @lastChrPair[0,1])) if $lastChrPair[0] ne '__NA__';
        }
    }
}

#--- load contacts of phased PE-reads ---
sub load_phasedPE_contacts{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};
    my $idxFunc = $parm{idxFunc};
    my $lastChrPairAf = $parm{lastChrPairAf};
    my $onlyPha = $parm{onlyPha} || 0; # only use phased Hi-C pair

    for my $pe_OB (@$pe_OB_poolAf){
        # get chr-pos ascending sorted all mapped reads_OB
        my $rOB_sortAref = $pe_OB->get_sorted_reads_OB(chrSortHref => $V_Href->{ChrThings}, chrSortKey  => 'turn');
        # recover SuppHaplo attribute
        $_->recover_SuppHaploAttr for @$rOB_sortAref;
        # extract haplotype confirmed alignments
        my @hasHaprOB = grep $_->has_SuppHaplo, @$rOB_sortAref;
        # double check
        if( scalar(@hasHaprOB) < 2 ){
            warn_and_exit "<ERROR>\tPhased PE-reads doesn't have at least two haplotype confirmed alignments.\n".Dumper($pe_OB);
        }
        # skip hapComb set by flanking region lacking phased contacts, see func 'get_rOBpair_HapLinkCount'
        if(    $onlyPha
            && (    $hasHaprOB[ 0]->is_fromUnPhasedRegRand
                ||  $hasHaprOB[-1]->is_fromUnPhasedRegRand
               )
        ){
            next; # do not take this PE into contact counting
        }
        # as sorted, so take first[0] and last[-1] one as index of contacted region
        my (%chr, %pIdx, %hapID);
        $chr{a}   = $hasHaprOB[ 0]->get_mseg;
        $chr{b}   = $hasHaprOB[-1]->get_mseg;
        next if(!exists $V_Href->{ChrThings}->{$chr{a}} || !exists $V_Href->{ChrThings}->{$chr{b}});
        # if new chr-pair, do contacts_to_count on former chr-pair
        if(    $chr{a} ne $lastChrPairAf->[0]
            || $chr{b} ne $lastChrPairAf->[1]
        ){
            # contacts to count for last chr-pair, release memory
            ## here, do de-dup phased reads (in single run) (optional)
            my $tag = "$lastChrPairAf->[2] " . join(',', @$lastChrPairAf[0,1]);
            &phasePE_contacts_to_count(tag => $tag) if $lastChrPairAf->[0] ne '__NA__';
            # update
            $lastChrPairAf->[0] = $chr{a};
            $lastChrPairAf->[1] = $chr{b};
        }
        # go on recording
        $pIdx{a}  = &{$idxFunc}(chr => $chr{a}, pos => $hasHaprOB[ 0]->get_mpos);
        $pIdx{b}  = &{$idxFunc}(chr => $chr{b}, pos => $hasHaprOB[-1]->get_mpos);
        $hapID{a} = $hasHaprOB[ 0]->get_SuppHaploStr;
        $hapID{b} = $hasHaprOB[-1]->get_SuppHaploStr;
        # hapID might have multiple haplotype, then make them really inter-haplotype
        if($hapID{a} =~ /,/ || $hapID{b} =~ /,/){
            ## i=initiative; p=passive
            my ($i,$p) = $hapID{b} =~ /,/ ? qw/ a b / : qw/ b a /;
            $hapID{$i} =~ s/,.+//; # arbitrary, greedy
            $hapID{$p} =  first {$_ ne $hapID{$i}} split /,/, $hapID{$p};
        }
        # get PE-Info string
        my $peInfoStr = join(';', map {( join(',', $_->get_mseg, $_->get_mpos) )} @$rOB_sortAref);
        # record, 'peInfoStr' as key, 'accumulated count' as value
        $V_Href->{phasePEdetails}->{$chr{a}}->{$chr{b}}->{$pIdx{a}}->{$pIdx{b}}->{"$hapID{a},$hapID{b}"}->{$peInfoStr} ++;
    }
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

    my $deDupCount = 0;
    for my $chr_a (sort keys %{$V_Href->{phasePEdetails}}){
        for my $chr_b (sort keys %{$V_Href->{phasePEdetails}->{$chr_a}}){
            for my $pIdx_a (sort {$a<=>$b} keys %{$V_Href->{phasePEdetails}->{$chr_a}->{$chr_b}}){
                for my $pIdx_b (sort {$a<=>$b} keys %{$V_Href->{phasePEdetails}->{$chr_a}->{$chr_b}->{$pIdx_a}}){
                    my $hapCombHf = $V_Href->{phasePEdetails}->{$chr_a}->{$chr_b}->{$pIdx_a}->{$pIdx_b};
                    for my $hapComb (sort keys %$hapCombHf){
                        if( $V_Href->{SkipDeDupPhasedReads} ){
                            # count all PE
                            my $count = sum( values  %{$hapCombHf->{$hapComb}} );
                            $deDupCount += $count;
                            $V_Href->{phasePEcontact}->{$chr_a}->{$chr_b}->{$pIdx_a}->{$pIdx_b}->{$hapComb} += $count;
                        }
                        else{
                            # count all PE de-redundancy by 'peInfo'
                            my $count = scalar( keys %{$hapCombHf->{$hapComb}} );
                            $deDupCount += $count;
                            $V_Href->{phasePEcontact}->{$chr_a}->{$chr_b}->{$pIdx_a}->{$pIdx_b}->{$hapComb} += $count;
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
    my $action = $V_Href->{SkipDeDupPhasedReads} ? 'Load' : 'De-dup';
    stout_and_sterr "[INFO]\t".`date`
                         ."\t$action $tag phased PE-reads and get $deDupCount contacts OK.\n" if $deDupCount >= 10;
}

#--- return haplo link count of given reads_OB pair ---
sub get_rOBpair_HapLinkCount{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $skipSort = $parm{skipSort} || 0;
    my $hapRegex = $parm{hapRegex} || undef;
    my $HapLinkHf = $parm{HapLinkHf};

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
    # first check chr-pair
    if(    ! exists $V_Href->{phasePEcontact}->{$mSeg{a}}
        || ! exists $V_Href->{phasePEcontact}->{$mSeg{a}}->{$mSeg{b}}
    ){
        $HapLinkHf->{stat}->{LackChrLink} ++;
        return [\%HapLinkCount, "LackChrLink($mSeg{a},$mSeg{b})", 'rd', 0];
    }
    # then, pos-pair
    my %mPos = map { ($_, $rOB{$_}->get_mpos) } keys %rOB;
    my %mExp = map { ($_, $mPos{$_}+$rOB{$_}->get_mRefLen-1) } keys %rOB;
    # try quick find from HapLinkC_pool
    my %winIdx = map {($_, Pos2Idx(pos => $mPos{$_}, winSize => $V_Href->{UKreadsFlankRegUnit}))} qw/a b/;
    my %winTag = map {($_, "$mSeg{$_},w$winIdx{$_}")} qw/a b/;
    if(exists $HapLinkHf->{link}->{$winTag{a}}){
        if(exists $HapLinkHf->{link}->{$winTag{a}}->{$winTag{b}}){
            $HapLinkHf->{stat}->{QuickFind} ++;
            my $LocRegInfoAf = $HapLinkHf->{link}->{$winTag{a}}->{$winTag{b}};
            # phased local-region stat
            if($LocRegInfoAf->[1] =~ /^RegionPhased;\(fSize:(\d+),/){
                my $HapLinkCountSum = sum(values %{$LocRegInfoAf->[0]});
                my $HapLinkDetails = join(';', map {"$_:$LocRegInfoAf->[0]->{$_}"} sort keys %{$LocRegInfoAf->[0]});
                $V_Href->{LocRegPhased}->{"$mSeg{a},$mSeg{b}"}->{ ($1*2) }->{$HapLinkCountSum}->{$HapLinkDetails} ++;
            }
            # return quickFind
            return $LocRegInfoAf;
        }
    }
    else{ # not find, must be next a-side bin-Idx
        %{$HapLinkHf->{link}} = ();
    }
    # bisection to find proper size of local flank region to provide sufficient phased-contacts to impute
    my $chrHapLinkHref = $V_Href->{phasePEcontact}->{$mSeg{a}}->{$mSeg{b}};
    my %fReg = map { ($_, {pos=>{}}) } keys %rOB;
    my @FlankSizeTime = (1, $V_Href->{UKreadsFlankRegUnitMaxTimes});
    my $time = undef;
    while(1){
        # use UnitTimes conceived from previously applied
        $time = defined($time) ? int(sum(@FlankSizeTime)/2) : $V_Href->{UKreadsFlankRegUnitIniTimes};
        my $FlankSize;
        my @judgement;
        for my $shift (-1, 0){
            $FlankSize = $V_Href->{UKreadsFlankRegUnit} * ($time + $shift);
            for my $s (sort keys %rOB){
                # deal chromosome ends of 5/3 prime flanking region
                my $mSegLen = $V_Href->{ChrThings}->{$mSeg{$s}}->{len};
                my $FlankSize_p5 = max($mPos{$s}-$FlankSize, 1       ) - $mPos{$s}; # negative value
                my $FlankSize_p3 = min($mExp{$s}+$FlankSize, $mSegLen) - $mExp{$s}; # positive value
                # record
                $fReg{$s}{pos}{p5} = $mPos{$s} + $FlankSize_p5;
                $fReg{$s}{pos}{p3} = $mExp{$s} + $FlankSize_p3;
            }
            # extract the haplo link of a/b paired flank regions
            my %pIdx;
            %{$pIdx{p5}} = map { ( $_, Pos2Idx(pos => $fReg{$_}{pos}{p5}, winSize => $V_Href->{mapPosWinSize}) ) } keys %fReg;
            %{$pIdx{p3}} = map { ( $_, Pos2Idx(pos => $fReg{$_}{pos}{p3}, winSize => $V_Href->{mapPosWinSize}) ) } keys %fReg;
            %HapLinkCount = (); # sweep
            my $quickCheck = 0;
            for my $pIdx_a ( $pIdx{p5}{a} .. $pIdx{p3}{a} ){
                next unless exists $chrHapLinkHref->{$pIdx_a};
                # intra-chr pIdx are sorted in ascending order
                my $pIdx_b_p5 = $pIdx{p5}{b};
                $pIdx_b_p5 = $pIdx_a if $mSeg{a} eq $mSeg{b} && $pIdx_a > $pIdx_b_p5;
                for my $pIdx_b ( $pIdx_b_p5 .. $pIdx{p3}{b} ){
                    next unless exists $chrHapLinkHref->{$pIdx_a}->{$pIdx_b};
                    my $pIdxHapLinkHref = $chrHapLinkHref->{$pIdx_a}->{$pIdx_b};
                    # only keep pre-set hapID combination
                    ## if-else is faster than incorporate '!defined $hapRegex || ' in grep syntax
                    if(defined $hapRegex){
                        $HapLinkCount{$_} += $pIdxHapLinkHref->{$_} for grep /$hapRegex/, keys %$pIdxHapLinkHref;
                    }
                    else{
                        $HapLinkCount{$_} += $pIdxHapLinkHref->{$_} for keys %$pIdxHapLinkHref;
                    }
                }
                # enable quick check
                $quickCheck ||= (scalar(keys %HapLinkCount) != 0) if $shift == -1 && !$quickCheck;
                # do quick check to exits loop
                last if $quickCheck && sum(values %HapLinkCount) >= $V_Href->{hapCombMinLinkForPhaReg};
            }
            # satisfy the minimum link count
            if(    scalar(keys %HapLinkCount) != 0
                && sum(values  %HapLinkCount) >= $V_Href->{hapCombMinLinkForPhaReg}
            ){
                push @judgement, 1;
                last if $shift == -1;
            }
            else{
                push @judgement, 0;
            }
        }
        # go on search?
        my $return = 0;
        my $phased = 0;
        if(    $judgement[0] == 0 # 0-1
            && $judgement[1] == 1
        ){
            $return = 1;
            $phased = 1;
        }
        else{
            if($judgement[0] == 1){ # 1-[1], first half part
                $FlankSizeTime[1] = $time - 1;
            }
            else{ # 0-0
                $FlankSizeTime[0] = $time + 1;
            }
            # check region reliability
            if($FlankSizeTime[0] > $FlankSizeTime[1]){
                $return = 1;
                $phased = 1 if sum(@judgement) > 0;
            }
        }
        # find enough HapLink or last time
        if($return){
            $FlankSize ||= $V_Href->{mapPosWinSize}; # zero flank size uses mappos window
            my $mark = $phased ? 'RegionPhased' : 'RegionNotPhased';
            $mark .= ";(fSize:$FlankSize"
                    .",$mSeg{$_}:$fReg{$_}{pos}{p5}-$fReg{$_}{pos}{p3}"
                    .",mPosWinIdx:$winIdx{$_}"
                    .')'
                    for sort keys %fReg;
            # record
            my $LocRegInfoAf = [\%HapLinkCount, $mark, ($phased?'ph':'rd'), 0];
            $HapLinkHf->{stat}->{Calculate} ++;
            $HapLinkHf->{link}->{$winTag{a}}->{$winTag{b}} = $LocRegInfoAf;
            # 1) phased local-region stat
            # 2) update the UnitTimes for next operation
            if($phased){
                my $HapLinkCountSum = sum(values %HapLinkCount);
                my $HapLinkDetails = join(';', map {"$_:$HapLinkCount{$_}"} sort keys %HapLinkCount);
                $V_Href->{LocRegPhased}->{"$mSeg{a},$mSeg{b}"}->{ ($FlankSize*2) }->{$HapLinkCountSum}->{$HapLinkDetails} ++;
                $V_Href->{UKreadsFlankRegUnitIniTimes} = max( min( $time*2, $V_Href->{UKreadsFlankRegUnitMidTimes} ), 1 );
            }
            else{
                $V_Href->{UKreadsFlankRegUnitIniTimes} = $time;
            }
            # return
            return $LocRegInfoAf;
        }
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
