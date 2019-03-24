package HaploHiC::PhasedHiC::sEndSoloHapConfirm;

use strict;
use warnings;
use List::Util qw/ max min sum first /;
use Data::Dumper;
use Parallel::ForkManager;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Index qw/ Pos2Idx /;
use HaploHiC::LoadOn;
use HaploHiC::PhasedHiC::splitPairBam qw/ forkSetting getTODOpairBamHrefArray write_peOB_to_chrPairBam getChrPairTagAf /;
use HaploHiC::PhasedHiC::phasedPEtoContact qw/ phasePE_to_contactCount get_rOBpair_HapLinkCount /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              sEndhx_get_HapLink
              prepareGetHapBamObj
              getChrPairBam
              chrPair_forkSetting
              chrPair_loadContacts
              chrPair_HapConfirmWork
              chrPair_mergeResult
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::sEndSoloHapConfirm';
#----- version --------
$VERSION = "0.16";
$DATE = '2019-03-24';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        sEndhx_get_HapLink
                        getChrPairBam
                        BamToChrPair
                        confirm_sEndU_Hap
                        prepareGetHapBamObj
                        chrPair_forkSetting
                        chrPair_loadContacts
                        chrPair_HapConfirmWork
                        chrPair_prepareGetHapBam
                        chrPair_sEndU_assignHap
                        findContactAnchor
                        chrPair_PEsplitStat
                        chrPair_LocRegPhased
                        chrPair_mergeResult
                        mergeChrPairResult
                        startWriteGetHapBam
                        writeChrPairReadsToGetHapBam
                        loadChrPairPEsplitStat
                        loadChrPairLocRegPhased
                        writeMergeLocRegPhased
                     /;

#--- assign one-side confirmed contacts PE to certain haplo-contacts ---
## sEnd-h[x]
sub sEndhx_get_HapLink{
    # prepare getHapBam name of all source bams
    for my $hapIdx (1 .. $V_Href->{haploCount}){
        my $tag = "phMut-sEnd-h$hapIdx";
        &prepareGetHapBamObj(pairBamHref => $_, tag => $tag) for @{$V_Href->{PairBamFiles}};
    }

    # current step is in pipeline
    if($V_Href->{stepToStart} <= 2){
        # split sEnd-h[x].bam to chrPair.bam
        ## fork at run level
        my %preChr;
        &getChrPairBam(preChrHf => \%preChr, type => 'sEndU');
        # assign haplo to 'UK' end in sEnd-h[x] PE
        ## fork at chrPair level (preChr + chrRel)
        &confirm_sEndU_Hap(preChrHf => \%preChr);
        # stop at current step
        exit(0) if $V_Href->{stepToStop} == 2;
    }
}

#--- prepare getHapBam object from hapSplit.bam ---
sub prepareGetHapBamObj{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};
    my $tag = $parm{tag};

    my $splitBamPath = $pairBamHref->{splitBam}->{$tag}->get_filepath;
    my $splitBamTag = $pairBamHref->{splitBam}->{$tag}->get_tag;
    # possible haplo combination
    my ($hapID) = ($tag =~ /sEnd-(h\d+)/);
    my @hapComb = $hapID ? ("${hapID}Intra") : (map {("h${_}Intra")} (1 .. $V_Href->{haploCount}));
    push @hapComb, 'hInter';
    for my $hapComb (@hapComb){
        # getHapBam object with filepath and tag
        my $getHapTag = "$tag.$hapComb";
        (my $getHapBamPath = $splitBamPath) =~ s/$tag/$getHapTag/;
        (my $getHapBamTag  = $splitBamTag)  =~ s/$tag/$getHapTag/;
        $pairBamHref->{splitBam}->{$getHapTag} = BioFuse::BioInfo::Objects::Bam_OB->new(filepath => $getHapBamPath, tag => $getHapBamTag);
        # prepare bam for final merge
        push @{$pairBamHref->{bamToMerge}->{"merge.$hapComb"}}, $pairBamHref->{splitBam}->{$getHapTag};
    }
}

#--- split bam to chrPair.bam (fork) ---
sub getChrPairBam{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $type = $parm{type} || 'sEndU';
    my $preChrHf = $parm{preChrHf};

    my @tag;
    if($type eq 'sEndU'){
        @tag = map {"phMut-sEnd-h$_"} (1 .. $V_Href->{haploCount});
    }
    elsif($type eq 'dEndU'){
        @tag = ('unknown');
    }

    # all bams or selected
    my @TODOpairBamHref = getTODOpairBamHrefArray;
    # fork manager
    my ($pm, $fork_DO) = forkSetting;
    # load PE from bam to split
    for my $pairBamHref ( @TODOpairBamHref ){
        # fork job starts
        if($fork_DO){ $pm->start($pairBamHref->{prefix}) and next }
        eval{
            # split to chrPair
            &BamToChrPair(pairBamHref => $pairBamHref, tagAf => \@tag);
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

    # load chrPair.bam obj
    for my $pairBamHref ( @TODOpairBamHref ){
        $pairBamHref->{chrPairBam} = {}; # reset
        for my $tag (@tag){
            my $splitBam = $pairBamHref->{splitBam}->{$tag};
            my $chrPairList = $splitBam->get_filepath . '.chrPairList';
            open (CPL, Try_GZ_Read($chrPairList)) || die "fail read chrPairList: $!\n";
            while(<CPL>){
                my ($chrPairTag, $chrPairBamPath) = (split);
                my $chrPairBamTag = "$pairBamHref->{prefix} $tag $chrPairTag";
                $pairBamHref->{chrPairBam}->{$tag}->{$chrPairTag} = BioFuse::BioInfo::Objects::Bam_OB->new(filepath => $chrPairBamPath, tag => $chrPairBamTag);
                my ($preChr) = (split /,/, $chrPairTag)[0];
                $preChrHf->{$preChr} = 1;
            }
            close CPL;
        }
    }
}

#--- split bam to chrPair ---
sub BamToChrPair{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};
    my $tagAf = $parm{tagAf} || [];

    for my $tag (@$tagAf){
        my $splitBam = $pairBamHref->{splitBam}->{$tag};
        # temp workspace
        my $chrPairDir = $splitBam->get_filepath . '-chrPairDir';
        `rm -rf $chrPairDir` if -d $chrPairDir;
        `mkdir -p $chrPairDir`;
        # read $tag split bam, and write peOB to chr-pair bams
        my %chrPairBam;
        my @subrtOpt = (subrtRef => \&write_peOB_to_chrPairBam,
                        subrtParmAref => [chrPairBamHf => \%chrPairBam, splitBam => $splitBam, chrPairDir => $chrPairDir, sorted => 1]);
        $splitBam->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC', quiet => 1, simpleLoad => 1, deal_peOB_pool => 1, @subrtOpt);
        # close chr-pair bams
        $_->stop_write for values %chrPairBam;
        # write chrPair list
        my $chrPairList = $splitBam->get_filepath . '.chrPairList';
        open (CPL, Try_GZ_Write($chrPairList)) || die "fail write chrPairList: $!\n";
        print CPL join("\t", $_, $chrPairBam{$_}->get_filepath)."\n" for sort keys %chrPairBam;
        close CPL;
        # inform
        my $mark = $splitBam->get_tag;
        stout_and_sterr "[INFO]\t".`date`
                             ."\tsplit $mark bam to chrPair bam OK.\n";
    }
}

#--- confirm sEnd-h[x] PE haplotype contact ---
## use multiple forks at preChr level
## here, deal with all phMut-sEnd-h[x],
## as they will be reload to fill the phased-PE-contact for next dEnd-unknown PE hapAssign
sub confirm_sEndU_Hap{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $preChrHf = $parm{preChrHf};

    my @preChrOpt = map {( [$_, 'intraChr'], [$_, 'interChr'] )}
                    grep exists $preChrHf->{$_}, @{$V_Href->{sortedChr}};
    my $pm = &chrPair_forkSetting;
    # each chrPair of sEnd-h[x].bam
    for my $preChrOpt ( @preChrOpt ){
        # fork job starts
        $pm->start(join(';',@$preChrOpt)) and next;
        eval{
            # options
            my ($preChr, $chrRel) = @$preChrOpt;
            # load phased-contact of this chrPair
            &chrPair_loadContacts(preChr => $preChr, chrRel => $chrRel);
            # confirm sEndU haplo of this chrPair
            &chrPair_HapConfirmWork(preChr => $preChr, chrRel => $chrRel);
        };
        if($@){
            warn $@;
            $pm->finish(1);
        }
        # fork job finishes
        $pm->finish(0);
    }
    # collect fork jobs
    $pm->wait_all_children;

    # merge chrPair data
    &chrPair_mergeResult;
}

#--- fork manager setting for parallel chrPair operations ---
sub chrPair_forkSetting{
    my $pm = new Parallel::ForkManager($V_Href->{forkNum});
    $pm->run_on_finish(
        sub {
            my ($pid, $exit_code, $chrSel, $signal) = @_;
            warn_and_exit "<ERROR>\t$chrSel child-process($pid) failed.\n" if $exit_code || $signal;
        }
    );
    return $pm;
}

#--- load phased contacts of given chrPair ---
sub chrPair_loadContacts{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $preChr = $parm{preChr} || undef;
    my $chrRel = $parm{chrRel} || undef;
    my $type = $parm{type} || 'sEndU';

    # load contacts from PE-reads having dual-side confirmed contacts
    ## dEnd-h[x] and [ds]End-hInter
    my $tagToBamHref = {};
    my @tags = map {("phMut-dEnd-h$_")} ( 1 .. $V_Href->{haploCount} ); # t1, but not t3("phMut-sEnd-h$_")
    push @tags, "phMut-dEnd-hInter", "phMut-sEnd-hInter"; # t2,t4
    # for dEndU, load contacts from PE-reads of sEndSoloHapPE getHap.bam
    ## sEnd-h[x].h[x]Intra and sEnd-h[x].hInter
    if($type eq 'dEndU'){
        push @tags, map {("phMut-sEnd-h$_.h${_}Intra", "phMut-sEnd-h$_.hInter")} ( 1 .. $V_Href->{haploCount} );
    }
    # load phased contacts from tags
    for my $tag ( @tags ){
        push @{$tagToBamHref->{$tag}}, $_->{splitBam}->{$tag} for @{$V_Href->{PairBamFiles}};
    }
    # load contacts of selected chrPair
    phasePE_to_contactCount(tagToBamHref => $tagToBamHref, preChr => $preChr, chrRel => $chrRel);
}

#--- confirm haplotype of U-end of given chrPair ---
sub chrPair_HapConfirmWork{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $preChr = $parm{preChr} || undef;
    my $chrRel = $parm{chrRel} || undef;
    my $assignHapSubrtRef = $parm{assignHapSubrtRef} || \&chrPair_sEndU_assignHap; # default is for sEndU, accept dEndU

    # all bams or selected
    my @TODOpairBamHref = getTODOpairBamHrefArray;
    for my $pairBamHref ( @TODOpairBamHref ){
        for my $tag (sort keys %{$pairBamHref->{chrPairBam}}){
            for my $chrPairTag (sort keys %{$pairBamHref->{chrPairBam}->{$tag}}){
                next if defined $preChr && $chrPairTag !~ /$preChr,/;
                next if defined $chrRel && $chrRel eq 'intraChr' && $chrPairTag !~ /^(.+),\1$/;
                next if defined $chrRel && $chrRel eq 'interChr' && $chrPairTag =~ /^(.+),\1$/;
                my $chrPairBam = $pairBamHref->{chrPairBam}->{$tag}->{$chrPairTag};
                # prepare getHapBam
                my %getHapBam;
                &chrPair_prepareGetHapBam(tag => $tag, chrPairBam => $chrPairBam, getHapBamHf => \%getHapBam);
                # read chrPairBam
                $V_Href->{PEsplitStat} = {}; # reset stat
                $V_Href->{LocRegPhased} = {}; # reset stat
                my $HapLinkHf = { link => {}, stat => {Calculate=>0, LackChrLink=>0, QuickFind=>0} };
                my @subrtOpt = (subrtRef => $assignHapSubrtRef,
                                subrtParmAref => [getHapBamHf => \%getHapBam, HapLinkHf => $HapLinkHf]);
                $chrPairBam->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC', deal_peOB_pool => 1, @subrtOpt);
                # close getHapBam file-handle
                $_->stop_write for values %getHapBam;
                # write PEsplitStat of this chrPairBam
                &chrPair_PEsplitStat(tag => $tag, chrPairBam => $chrPairBam);
                # write LocRegPhased of this chrPairBam
                &chrPair_LocRegPhased(chrPairBam => $chrPairBam);
                # inform
                my $mark = $chrPairBam->get_tag;
                my $HapLinkStat = join('; ', map {("$_:$HapLinkHf->{stat}->{$_}")} sort keys %{$HapLinkHf->{stat}});
                stout_and_sterr "[INFO]\t".`date`
                                     ."\tassign haplotype to UK-end of PE-reads from $mark bam OK.\n"
                                     ."\t$mark HapLinkStat: $HapLinkStat\n";
            }
        }
    }
}

#--- start writing getHapBam of given chrPair (file-w-handle, SAM-header) ---
sub chrPair_prepareGetHapBam{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $tag = $parm{tag};
    my $chrPairBam = $parm{chrPairBam};
    my $getHapBamHf = $parm{getHapBamHf};

    my $chrPairBamPath = $chrPairBam->get_filepath;
    my $chrPairBamHead = join('', @{$chrPairBam->get_SAMheader(samtools => $V_Href->{samtools})}) . "\n";
    # possible haplo combination
    my ($hapID) = ($tag =~ /sEnd-(h\d+)/);
    my @hapComb = $hapID ? ("${hapID}Intra") : (map {("h${_}Intra")} (1 .. $V_Href->{haploCount}));
    push @hapComb, 'hInter';
    for my $hapComb (@hapComb){
        # getHapBam object with filepath and tag
        $getHapBamHf->{$hapComb} = BioFuse::BioInfo::Objects::Bam_OB->new(filepath => "$chrPairBamPath.$hapComb.bam", tag => "$tag.$hapComb");
        $getHapBamHf->{$hapComb}->start_write(samtools => $V_Href->{samtools});
        $getHapBamHf->{$hapComb}->write(content => $chrPairBamHead);
    }
}

#--- assign haplotype to sEndh[x] PE UK-end ---
sub chrPair_sEndU_assignHap{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};
    my $getHapBamHf = $parm{getHapBamHf};
    my $HapLinkHf = $parm{HapLinkHf};

    for my $pe_OB (@$pe_OB_poolAf){
        # get chr-pos ascending sorted all mapped reads_OB
        my $rOB_sortAref = $pe_OB->get_sorted_reads_OB(chrSortHref => $V_Href->{ChrThings}, chrSortKey  => 'turn');
        # recover SuppHaplo attribute
        $_->recover_SuppHaploAttr for @$rOB_sortAref;
        # get index of hasHap and nonHap reads_OB
        my ($hasHapIdx, $nonHapIdx) = &findContactAnchor(rOB_sortAref => $rOB_sortAref);
        my $hasHap_rOB = $rOB_sortAref->[$hasHapIdx];
        my $nonHap_rOB = $rOB_sortAref->[$nonHapIdx];
        my $hasHapID   = $hasHap_rOB->get_SuppHaploStr;
        # get haplo link count of paired rOB, <needs to sort again>
        ## only keep pre-set hapID combination
        my $regex = $hasHapIdx < $nonHapIdx ? "^$hasHapID," : ",$hasHapID\$";
        my $LocRegInfoAf = get_rOBpair_HapLinkCount(rOB_a => $hasHap_rOB, rOB_b => $nonHap_rOB, hapRegex => $regex, HapLinkHf => $HapLinkHf);
        my ($HapLinkC_Hf, $mark, $assignMethod, $modBool) = @$LocRegInfoAf;
        # assign HapID to nonHap_rOB
        ## once local region is not phased, reset mark and loads pre-defined HapLink
        my $isIntraChr = $hasHap_rOB->get_mseg eq $nonHap_rOB->get_mseg;
        unless($modBool){
            if($assignMethod eq 'rd'){
                if($isIntraChr){
                    $HapLinkC_Hf->{"$hasHapID,$hasHapID"} = 1;
                    $mark .= ';IntraChrSetHapIntra';
                }
                else{
                    $HapLinkC_Hf->{$_} = 1 for grep /$regex/, @{$V_Href->{allHapComb}};
                    $mark .= ';InterChrSetAllAvibHap';
                }
            }
            else{
                # uniform addition
                if($V_Href->{uniformAddRatioForHapComb}){
                    my $add = sum(values %$HapLinkC_Hf) * $V_Href->{uniformAddRatioForHapComb};
                    $HapLinkC_Hf->{$_} += $add for grep /$regex/, @{$V_Href->{allHapComb}};
                }
            }
            # update
            $LocRegInfoAf->[1] = $mark;
            $LocRegInfoAf->[3] = 1;
        }
        ## select hapComb
        my $assHapComb;
        my @HapComb = sort keys %$HapLinkC_Hf;
        ## if has only single HapComb, just take it
        if(scalar(@HapComb) == 1){
            $assHapComb = $HapComb[0];
            # add mark
            $mark .= ";SoloHapComb;[$assHapComb](C:$HapLinkC_Hf->{$assHapComb})";
        }
        else{ # luck draw
            # prepare haplotype's luck-draw interval
            my $allCount = 0;
            my %hapCombDraw;
            for my $hapComb (@HapComb){
                $allCount += $HapLinkC_Hf->{$hapComb};
                $hapCombDraw{$hapComb} = $allCount;
            }
            # random pick
            my $luck_draw = int(rand($allCount * 1E3)) / 1E3; # not sprintf
            $assHapComb = first { $hapCombDraw{$_} > $luck_draw } @HapComb;
            # add mark
            $mark .= ";[$_](C:$HapLinkC_Hf->{$_},D:$hapCombDraw{$_})" for @HapComb;
            $mark .= ";RD:$luck_draw";
        }
        ## remove the pre-set hasHap-ID
        (my $assHapID = $assHapComb) =~ s/$regex//;
        ## assign
        $nonHap_rOB->load_SuppHaplo(hapID => $assHapID, allele_OB => 'NULL');
        $nonHap_rOB->addHapIDtoOptfd; # update
        $nonHap_rOB->add_str_to_optfd(str => "\tXU:Z:$mark");
        # write PE to getHapBam files
        my $getHapTag = $assHapID eq $hasHapID ? "${assHapID}Intra" : "hInter";
        $getHapBamHf->{$getHapTag}->write(content => join("\n",@{$pe_OB->printSAM(keep_all=>1)})."\n");
        # stat
        my $chrTag = $isIntraChr ? 'IntraChr' : 'InterChr';
        $V_Href->{PEsplitStat}->{"$getHapTag.$chrTag.$assignMethod"} ++;
    }
}

#--- find the index in rOB_sortArray of hasHap and nonHap reads_OB ---
sub findContactAnchor{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $rOB_sortAref = $parm{rOB_sortAref};

    my (@hasHapIdx, @nonHapIdx);
    for my $i (0 .. scalar(@$rOB_sortAref)-1){
        if($rOB_sortAref->[$i]->has_SuppHaplo){
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
        my $hasHaprOB_F = $rOB_sortAref->[$hasHapIdx[0]];
        my $hasHaprOB_L = $rOB_sortAref->[$hasHapIdx[-1]];
        my $nonHaprOB_F = $rOB_sortAref->[$nonHapIdx[0]];
        my $nonHaprOB_L = $rOB_sortAref->[$nonHapIdx[-1]];
        my $hFnLsameChr = ($hasHaprOB_F->get_mseg eq $nonHaprOB_L->get_mseg);
        my $hLnFsameChr = ($hasHaprOB_L->get_mseg eq $nonHaprOB_F->get_mseg);
        # same chr
        if( $hFnLsameChr && $hLnFsameChr ){
            # larger pos gap
            if(   abs($hasHaprOB_F->get_mpos - $nonHaprOB_L->get_mpos)
                > abs($hasHaprOB_L->get_mpos - $nonHaprOB_F->get_mpos)
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
            if(   $hasHaprOB_F->get_mReadLen + $nonHaprOB_L->get_mReadLen
                > $hasHaprOB_L->get_mReadLen + $nonHaprOB_F->get_mReadLen
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

#--- write PEsplitStat of sEndU of given chrPair ---
sub chrPair_PEsplitStat{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $tag = $parm{tag};
    my $chrPairBam = $parm{chrPairBam};

    my $PEsplitStat = $chrPairBam->get_filepath . '.PEsplitStat';
    open (PES, Try_GZ_Write($PEsplitStat)) || die "fail write PEsplitStat: $!\n";
    print PES "$tag.$_\t$V_Href->{PEsplitStat}->{$_}\n" for sort keys %{$V_Href->{PEsplitStat}};
    close PES;
}

#--- write LocRegPhased of sEndU of given chrPair ---
sub chrPair_LocRegPhased{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chrPairBam = $parm{chrPairBam};

    my $LocRegPhased = $chrPairBam->get_filepath . '.statOfPhasedLocReg.gz';
    open (LRP, Try_GZ_Write($LocRegPhased)) || die "fail write LocRegPhased: $!\n";
    for my $ChrPair (sort keys %{$V_Href->{LocRegPhased}}){
        for my $LocRegSize (sort {$a<=>$b} keys %{$V_Href->{LocRegPhased}->{$ChrPair}}){
            for my $LinkCount (sort {$a<=>$b} keys %{$V_Href->{LocRegPhased}->{$ChrPair}->{$LocRegSize}}){
                for my $LinkDetails (sort keys %{$V_Href->{LocRegPhased}->{$ChrPair}->{$LocRegSize}->{$LinkCount}}){
                    print LRP join( "\t",
                                     $ChrPair,
                                     $LocRegSize,
                                     $LinkCount,
                                     $LinkDetails,
                                     $V_Href->{LocRegPhased}->{$ChrPair}->{$LocRegSize}->{$LinkCount}->{$LinkDetails}
                                   ) . "\n";
                }
            }
        }
    }
    close LRP;
}

#--- merge results of all chrPair ---
sub chrPair_mergeResult{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $type = $parm{type} || 'sEndU';

    # all bams or selected
    my @TODOpairBamHref = getTODOpairBamHrefArray;
    # fork manager
    my ($pm, $fork_DO) = forkSetting;
    # merge chrPair results
    for my $pairBamHref ( @TODOpairBamHref ){
        # fork job starts
        if($fork_DO){ $pm->start($pairBamHref->{prefix}) and next }
        eval{
            # merge chrPair related results (bam, PEsplitStat, LocRegPhased)
            &mergeChrPairResult(pairBamHref => $pairBamHref, type => $type);
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
}

#--- merge chrPair related results (bam, PEsplitStat, LocRegPhased) ---
sub mergeChrPairResult{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};
    my $type = $parm{type} || 'sEndU';

    $V_Href->{PEsplitStat} = {}; # reset
    for my $tag (sort keys %{$pairBamHref->{chrPairBam}}){
        $V_Href->{LocRegPhased} = {}; # reset
        # prepare getHapBam and file-handle
        &startWriteGetHapBam(pairBamHref => $pairBamHref, tag => $tag);
        # possible haplo combination
        my ($hapID) = ($tag =~ /sEnd-(h\d+)/);
        my @hapComb = $hapID ? ("${hapID}Intra") : (map {("h${_}Intra")} (1 .. $V_Href->{haploCount}));
        push @hapComb, 'hInter';
        # each hapComb
        my @ToInform;
        for my $hapComb (@hapComb){
            # getHapBam
            my $getHapTag = "$tag.$hapComb";
            my $getHapBam = $pairBamHref->{splitBam}->{$getHapTag};
            push @ToInform, $getHapBam->get_filepath;
            # prepare PEsplitStat
            for my $chrTag (qw/ IntraChr InterChr /){
                for my $assignMethod (qw/ ph rd /){
                    $V_Href->{PEsplitStat}->{"$tag.$hapComb.$chrTag.$assignMethod"} = 0;
                }
            }
            # load each chrPair results
            for my $chrPairTag (@{ &getChrPairTagAf }){
                next unless exists $pairBamHref->{chrPairBam}->{$tag}->{$chrPairTag};
                my $chrPairBam = $pairBamHref->{chrPairBam}->{$tag}->{$chrPairTag};
                # load chrPair hapComb bam and merge to getHapBam
                my $chrPairGetHapBam = BioFuse::BioInfo::Objects::Bam_OB->new(filepath => $chrPairBam->get_filepath . ".$hapComb.bam", tag => "$tag.$hapComb");
                my @subrtOpt = (subrtRef => \&writeChrPairReadsToGetHapBam, subrtParmAref => [getHapBam => $getHapBam]);
                $chrPairGetHapBam->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC', quiet => 1, simpleLoad => 1, deal_peOB_pool => 1, @subrtOpt);
                # load chrPair PEsplitStat, only once, so let's choose 'hInter'
                &loadChrPairPEsplitStat(chrPairBam => $chrPairBam) if $hapComb eq 'hInter';
                # load chrPair LocRegPhased, only once, so let's choose 'hInter'
                &loadChrPairLocRegPhased(chrPairBam => $chrPairBam) if $hapComb eq 'hInter';
            }
            # close merge getHap bam
            $getHapBam->stop_write;
        }
        # write merge stat data of phased local region
        &writeMergeLocRegPhased(pairBamHref => $pairBamHref, tag => $tag);
        # sweep chrPair folder and list
        my $chrPairList = $pairBamHref->{splitBam}->{$tag}->get_filepath . '.chrPairList';
        my $chrPairDir  = $pairBamHref->{splitBam}->{$tag}->get_filepath . '-chrPairDir';
        # `rm -rf $chrPairDir $chrPairList`; # debug
        # inform
        stout_and_sterr "[INFO]\t".`date`
                             ."\tmerge getHap bam files OK.\n"
                             ."\t".join("\n\t", @ToInform)."\n";
    }

    # write PEsplitStat report (2nd part for sEndU, 3rd part for dEndU)
    ## previous contents
    my @Content = split /\s+/, `cat $pairBamHref->{PEsplit_report}`;
    my $regex = $type eq 'dEndU' ? '^unknown.h[I\d]' : '\.[rp][dh]:$';
    my %Content = map {(@Content[$_,$_+1])} grep {$_%2==0 && $Content[$_] !~ /$regex/} (0 .. $#Content);
    ## add contents of this step
    open (PESRT, Try_GZ_Write($pairBamHref->{PEsplit_report})) || die "fail write '$pairBamHref->{prefix}' report: $!\n";
    print PESRT "$Content[$_]\t$Content{$Content[$_]}\n" for grep {$_%2==0 && exists $Content{$Content[$_]}} (0 .. $#Content);
    print PESRT "$_:\t$V_Href->{PEsplitStat}->{$_}\n" for grep {!exists $Content{"$_:"}} sort keys %{$V_Href->{PEsplitStat}};
    close PESRT;
}

#--- start writing getHapBam (file-w-handle, SAM-header) ---
sub startWriteGetHapBam{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};
    my $tag = $parm{tag};

    # prepare SAM-header
    my $splitBam = $pairBamHref->{splitBam}->{$tag};
    my $HeadAf = $splitBam->get_SAMheader(samtools => $V_Href->{samtools});
    ## add options of second part in @PG line
    my $PGi = first {$HeadAf->[$_] =~ /\@PG/} (0 .. scalar(@$HeadAf)-1);
    chomp($HeadAf->[$PGi]);
    $HeadAf->[$PGi] .= " -ucrfut $V_Href->{UKreadsFlankRegUnit}";
    $HeadAf->[$PGi] .= " -ucrfmx $V_Href->{UKreadsFlankRegMax}";
    $HeadAf->[$PGi] .= " -min_ct $V_Href->{hapCombMinLinkForPhaReg}";
    $HeadAf->[$PGi] .= " -add_r $V_Href->{uniformAddRatioForHapComb}" if $V_Href->{uniformAddRatioForHapComb};
    my $HeadStr = join('', @$HeadAf) . "\n";

    # getHap.bam file-handle
    for my $getHapTag (grep /^$tag\./, sort keys %{$pairBamHref->{splitBam}}){
        my $getHapBam = $pairBamHref->{splitBam}->{$getHapTag};
        $getHapBam->start_write(samtools => $V_Href->{samtools});
        $getHapBam->write(content => $HeadStr);
    }
}

#--- write reads from chrPairBam to getHapBam ---
sub writeChrPairReadsToGetHapBam{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};
    my $getHapBam = $parm{getHapBam};
    $getHapBam->write(content => join("\n",@{$_->printSAM(keep_all=>1)})."\n") for @$pe_OB_poolAf;
}

#--- load chrPair PEsplitStat ---
sub loadChrPairPEsplitStat{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chrPairBam = $parm{chrPairBam};

    my $chrPairPEsplitStat = $chrPairBam->get_filepath . '.PEsplitStat';
    open (PES, Try_GZ_Read($chrPairPEsplitStat)) || die "fail read PEsplitStat: $!\n";
    while(<PES>){
        my ($key, $count) = (split)[0,1];
        $V_Href->{PEsplitStat}->{$key} += $count;
    }
    close PES;
}

#--- load chrPair LocRegPhased ---
sub loadChrPairLocRegPhased{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chrPairBam = $parm{chrPairBam};

    my $chrPairLocRegPhased = $chrPairBam->get_filepath . '.statOfPhasedLocReg.gz';
    open (LRP, Try_GZ_Read($chrPairLocRegPhased)) || die "fail read chrPairLocRegPhased: $!\n";
    while(<LRP>){
        my ($ChrPair, $LocRegSize, $LinkCount, $LinkDetails, $Count) = (split)[0..4];
        $V_Href->{LocRegPhased}->{$ChrPair}->{$LocRegSize}->{$LinkCount}->{$LinkDetails} += $Count;
    }
    close LRP;
}

#--- write stat data of phased local region ---
sub writeMergeLocRegPhased{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $tag = $parm{tag};
    my $pairBamHref = $parm{pairBamHref};

    my $splitBam = $pairBamHref->{splitBam}->{$tag};
    (my $LocRegPhased = $splitBam->get_filepath) =~ s/bam$/statOfPhasedLocReg.gz/;
    (my $localRegionUnit = sprintf "%e", $V_Href->{UKreadsFlankRegUnit} * 2) =~ s/\.?0*e\+0?/E/;
    (my $shortTag = $tag) =~ s/^phMut-//;

    open (LRP, Try_GZ_Write($LocRegPhased)) || die "fail to write statOfPhasedLocReg: $!\n";
    print LRP "##LocalRegionUnit: $localRegionUnit\n";
    print LRP join("\t", '#Tag', 'ChrPair', 'LocalRegionSize', 'PhasedContactsCount', 'Details', 'Amount') . "\n";
    for my $ChrPair (sort keys %{$V_Href->{LocRegPhased}}){
        for my $LocRegSize (sort {$a<=>$b} keys %{$V_Href->{LocRegPhased}->{$ChrPair}}){
            for my $LinkCount (sort {$a<=>$b} keys %{$V_Href->{LocRegPhased}->{$ChrPair}->{$LocRegSize}}){
                for my $LinkDetails (sort keys %{$V_Href->{LocRegPhased}->{$ChrPair}->{$LocRegSize}->{$LinkCount}}){
                    print LRP join( "\t",
                                     $shortTag,
                                     $ChrPair,
                                     $LocRegSize,
                                     $LinkCount,
                                     $LinkDetails,
                                     $V_Href->{LocRegPhased}->{$ChrPair}->{$LocRegSize}->{$LinkCount}->{$LinkDetails}
                                   ) . "\n";
                }
            }
        }
    }
    close LRP;
}

#--- 
1; ## tell the perl script the successful access of this module.
