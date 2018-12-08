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
use HaploHiC::PhasedHiC::phasedPEtoContact qw/ phasePE_to_contactCount smartBam_PEread get_rOBpair_HapLinkCount /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              sEndhx_get_HapLink
              getTODOpairBamHrefArray
              prepare_getHapBamName
              prepare_getHapBamFH
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::sEndSoloHapConfirm';
#----- version --------
$VERSION = "0.04";
$DATE = '2018-11-22';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        sEndhx_get_HapLink
                        getTODOpairBamHrefArray
                        confirm_sEndSoloHapPE_HapLink
                        prepare_getHapBamName
                        prepare_getHapBamFH
                        assign_sEndUKend_haplotype
                        findContactAnchor
                     /;

#--- assign one-side confirmed contacts PE to certain haplo-contacts ---
## sEnd-h[x]
sub sEndhx_get_HapLink{
    # load contacts from PE-reads having dual-side confirmed contacts
    ## dEnd-h[x] and [ds]End-hInter
    ## only run when starts before step4 (readsMerge)
    if($V_Href->{stepToStart} < 4){
        my $tagToBamHref = {};
        my @tags = map {("phMut-dEnd-h$_")} ( 1 .. $V_Href->{haploCount} ); # t1, but not t3("phMut-sEnd-h$_")
        push @tags, "phMut-dEnd-hInter", "phMut-sEnd-hInter"; # t2,t4
        for my $tag ( @tags ){
            push @{$tagToBamHref->{$tag}}, [ $_->{prefix}, $_->{splitBam}->{$tag} ] for @{$V_Href->{PairBamFiles}};
        }
        phasePE_to_contactCount( tagToBamHref => $tagToBamHref );
    }

    # assign haplo to 'UK' end in sEnd-h[x] PE
    ## write getHap.bam
    &confirm_sEndSoloHapPE_HapLink;
}

#--- return Href-array of all bams or selected ---
sub getTODOpairBamHrefArray{
    # all bams or selected
    return   scalar(keys %{$V_Href->{SelectBamPref}}) == 0
           ? @{$V_Href->{PairBamFiles}}
           : grep exists $V_Href->{SelectBamPref}->{$_->{prefix}}, @{$V_Href->{PairBamFiles}};
}

#--- confirm sEnd-h[x] PE haplotype contact ---
## use multiple forks
## here, deal with all phMut-sEnd-h[x],
## as they will be reload to fill the phased-PE-contact for next dEnd-unknown PE hapAssign
sub confirm_sEndSoloHapPE_HapLink{
    # prepare getHap.bam name of all source bams
    for my $hapIdx (1 .. $V_Href->{haploCount}){
        my $tag = "phMut-sEnd-h$hapIdx";
        &prepare_getHapBamName(pairBamHref => $_, tag => $tag, hapID => "h$hapIdx", splitBam => $_->{splitBam}->{$tag}) for @{$V_Href->{PairBamFiles}};
    }

    # start from step after current step
    return if $V_Href->{stepToStart} > 2;

    # all bams or selected
    my @TODOpairBamHref = &getTODOpairBamHrefArray;
    # fork manager
    my $pm;
    my $forkNum = min( $V_Href->{forkNum}, scalar(@TODOpairBamHref) );
    my $fork_DO = ( $forkNum > 1 );
    if( $fork_DO ){ $pm = new Parallel::ForkManager($forkNum) }
    # load PE from sEnd-h[x].bam
    for my $pairBamHref ( @TODOpairBamHref ){
        # fork job starts
        if( $fork_DO ){ $pm->start and next; }
        # load reads from each bam
        my $bamPrefix = $pairBamHref->{prefix};
        for my $hapIdx (1 .. $V_Href->{haploCount}){
            my $tag = "phMut-sEnd-h$hapIdx";
            my $mark = "'$bamPrefix' $tag";
            my $hapSplitBam = $pairBamHref->{splitBam}->{$tag};
            # check existence
            warn_and_exit "<ERROR>\tCannot find $tag bam: $hapSplitBam\n" unless file_exist(filePath => $hapSplitBam);
            # prepare getHap.bam and file-handle
            &prepare_getHapBamFH(pairBamHref => $pairBamHref, tag => $tag, splitBam => $hapSplitBam);
            # read hapSplit.BAM
            ## FLAG: 0x100(sd), 0x400(d), 0x800(sp)
            ## WEDO: -F 0x400(d), as 'sd' and 'sp' alignments is important in HiC
            smartBam_PEread(bam => $hapSplitBam, mark => $mark, viewOpt => '-F 0x400',
                            subrtRef => \&assign_sEndUKend_haplotype,
                            subrtParmAref => [tag=>$tag, pairBamHref => $pairBamHref]);
            # close getHap.bam file-handle
            close $_ for values %{$pairBamHref->{splitBamFH}};
            $pairBamHref->{splitBamFH} = {}; # empty splitBamFH key
            # inform
            stout_and_sterr "[INFO]\t".`date`
                                 ."\tassign haplotype to UK-end of PE-reads from $mark bam OK.\n";
        }
        # write report (2nd part)
        ## previous contents
        my @Content = split /\s+/, `cat $pairBamHref->{PEsplit_report}`;
        my %Content = map {(@Content[$_,$_+1])} grep {$_%2==0 && $Content[$_] !~ /\.[rp][dh]:$/} (0 .. $#Content);
        ## add contents of this step
        open (PESRT, Try_GZ_Write($pairBamHref->{PEsplit_report})) || die "fail write '$bamPrefix' report: $!\n";
        print PESRT "$Content[$_]\t$Content{$Content[$_]}\n" for grep {$_%2==0 && exists $Content{$Content[$_]}} (0 .. $#Content);
        print PESRT "$_:\t$V_Href->{PEsplitStat}->{$_}\n" for grep {!exists $Content{"$_:"}} sort keys %{$V_Href->{PEsplitStat}};
        close PESRT;
        # fork job finishes
        if( $fork_DO ){ $pm->finish; }
    }
    # collect fork jobs
    if( $fork_DO ){ $pm->wait_all_children; }

    # stop at current step
    exit(0) if $V_Href->{stepToStop} == 2;
}

#--- prepare getHap.bam from hapSplit.bam ---
## tag and name
sub prepare_getHapBamName{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};
    my $tag = $parm{tag};
    my $hapID = $parm{hapID};
    my $splitBam = $parm{splitBam};

    # possible haplo combination
    my @hapComb = $hapID ? ("${hapID}Intra") : (map {("h${_}Intra")} (1 .. $V_Href->{haploCount}));
    push @hapComb, 'hInter';
    # getHap.bam tag, name and file-handle
    for my $hapComb (@hapComb){
        # set tag and name
        my $getHapTag = "$tag.$hapComb";
        (my $getHapBam = $splitBam) =~ s/$tag/$getHapTag/;
        $pairBamHref->{splitBam}->{$getHapTag} = $getHapBam;
        # prepare bam to final merge
        push @{$pairBamHref->{bamToMerge}->{"merge.$hapComb"}}, $getHapBam;
    }
}

#--- prepare getHap.bam from hapSplit.bam ---
## file-handle and SAM-header
sub prepare_getHapBamFH{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};
    my $tag = $parm{tag};
    my $splitBam = $parm{splitBam};

    # getHap.bam file-handle
    my @ToInform;
    for my $getHapTag (grep /^$tag\./, sort keys %{$pairBamHref->{splitBam}}){
        my $getHapBam = $pairBamHref->{splitBam}->{$getHapTag};
        # file-handle
        open ($pairBamHref->{splitBamFH}->{$getHapTag}, "| $V_Href->{samtools} view -b -o $getHapBam") || die "fail write getHapBam: $!\n";
        push @ToInform, $getHapBam;
    }
    # SAM header
    open (SPLBAM,"$V_Href->{samtools} view -H $splitBam |") || die"fail read splitBam: $!\n";
    while( my $Hline = <SPLBAM> ){
        if( $Hline !~ /\@PG/ ){
            print {$_} $Hline for values %{$pairBamHref->{splitBamFH}};
        }
        else{ # add options of second part
            chomp($Hline);
            $Hline .= " -ucrfut $V_Href->{UKreadsFlankRegUnit}";
            $Hline .= " -ucrftm $V_Href->{UKreadsFlankRegUnitMaxTimes}";
            $Hline .= " -ucrpha $V_Href->{UKreadsMaxPhasedHetMut}";
            print {$_} "$Hline\n" for values %{$pairBamHref->{splitBamFH}};
        }
    }
    close SPLBAM;
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tgetHap bam files are well prepared.\n"
                         ."\t".join("\n\t", @ToInform)."\n";
}

#--- assign haplotype to sEndh[x] PE UK-end ---
sub assign_sEndUKend_haplotype{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB = $parm{pe_OB};
    my $pairBamHref = $parm{pairBamHref};
    my $tag = $parm{tag};

    # get chr-pos ascending sorted all mapped reads_OB
    my $rOB_sortAref = $pe_OB->get_sorted_reads_OB(rEndAref => [1,2], onlyMap => 1,
                                                   chrSortHref => $V_Href->{ChrThings},
                                                   chrSortKey  => 'turn');
    # recover SuppHaplo attribute
    $_->recover_SuppHaploAttr for @$rOB_sortAref;
    # get index of hasHap and nonHap reads_OB
    my ($hasHapIdx, $nonHapIdx) = &findContactAnchor(rOB_sortAref => $rOB_sortAref);
    my $hasHap_rOB = $rOB_sortAref->[$hasHapIdx];
    my $nonHap_rOB = $rOB_sortAref->[$nonHapIdx];
    my $hasHapID   = $hasHap_rOB->get_SuppHaploStr;
    # get haplo link count of paired rOB, <needs to sort again>
    my ($HapLinkC_Href, $mark) = get_rOBpair_HapLinkCount(rOB_a => $hasHap_rOB, rOB_b => $nonHap_rOB);
    # only keep pre-set hapID combination
    my $regex = $hasHapIdx < $nonHapIdx ? "^$hasHapID," : ",$hasHapID\$";
    delete $HapLinkC_Href->{$_} for grep !/$regex/, keys %$HapLinkC_Href;
    ## once empty, reset mark
    my @HapComb = sort keys %$HapLinkC_Href;
    my $HapCombCnt = scalar(@HapComb);
    $mark .= ";LackRequiredHapLink($hasHapID)" unless $HapCombCnt;
    # assign HapID to nonHap_rOB
    my $assHapID;
    my $assignMethod;
    ## if empty, just set as haplo-intra
    if($HapCombCnt == 0){ # even it matches 'RegionPhased' tag, still might lack $hasHapID related contacts
        $assHapID = $hasHapID;
        $mark .= ';SetHapIntra';
        $assignMethod = 'rd';
    }
    else{
        $assignMethod = 'ph';
        my $assHapComb;
        # if has only single HapComb, just take it
        if($HapCombCnt == 1){
            $assHapComb = $HapComb[0];
            # add mark
            $mark .= ";SoloHapComb;[$assHapComb](C:$HapLinkC_Href->{$assHapComb})";
        }
        else{ # luck draw
            # prepare haplotype's luck-draw interval
            my $allCount = 0;
            my %hapCombDraw;
            for my $hapComb (@HapComb){
                $allCount += $HapLinkC_Href->{$hapComb};
                $hapCombDraw{$hapComb} = $allCount;
            }
            # random pick
            my $luck_draw = int(rand($allCount));
            $assHapComb = first { $hapCombDraw{$_} > $luck_draw } sort keys %hapCombDraw;
            # add mark
            $mark .= ";[$_](C:$HapLinkC_Href->{$_},D:$hapCombDraw{$_})" for keys %hapCombDraw;
            $mark .= ";RD:$luck_draw";
        }
        # remove the pre-set hasHap-ID
        ($assHapID = $assHapComb) =~ s/$regex//;
    }
    ## assign
    $nonHap_rOB->load_SuppHaplo(hapID => $assHapID, allele_OB => 'NULL');
    $nonHap_rOB->addHapIDtoOptfd; # update
    $nonHap_rOB->add_str_to_optfd(str => "\tXU:Z:$mark");
    # output
    my $getHapTag = $assHapID eq $hasHapID ? "$tag.${assHapID}Intra" : "$tag.hInter";
    # write PE to getHap.bam files
    print {$pairBamHref->{splitBamFH}->{$getHapTag}} "$_\n" for @{$pe_OB->printSAM};
    # stat
    $V_Href->{PEsplitStat}->{"$getHapTag.$assignMethod"} ++;
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
    my $time = 0;
    FindAnchor: {
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
            if(    $hasHaprOB_F->get_mseg eq $nonHaprOB_L->get_mseg
                && $hasHaprOB_L->get_mseg eq $nonHaprOB_F->get_mseg
            ){ # same chr
                if(   abs($hasHaprOB_F->get_mpos - $nonHaprOB_L->get_mpos)
                    > abs($hasHaprOB_L->get_mpos - $nonHaprOB_F->get_mpos)
                ){ # larger pos gap
                    return ($hasHapIdx[0], $nonHapIdx[-1]);
                }
                else{
                    return ($hasHapIdx[-1], $nonHapIdx[0]);
                }
            }
            else{ # diff chr
                if( $nonHaprOB_F->is_suppmap ){ # avoid supp-map nonHap alignment
                    return ($hasHapIdx[0], $nonHapIdx[-1]);
                }
                else{
                    return ($hasHapIdx[-1], $nonHapIdx[0]);
                }
            }
        }
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
