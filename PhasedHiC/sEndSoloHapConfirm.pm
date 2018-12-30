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
use HaploHiC::PhasedHiC::phasedPEtoContact qw/ phasePE_to_contactCount get_rOBpair_HapLinkCount /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              sEndhx_get_HapLink
              getTODOpairBamHrefArray
              prepareGetHapBamObj
              startWriteGetHapBam
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::sEndSoloHapConfirm';
#----- version --------
$VERSION = "0.05";
$DATE = '2018-12-29';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        sEndhx_get_HapLink
                        getTODOpairBamHrefArray
                        confirm_sEndSoloHapPE_HapLink
                        prepareGetHapBamObj
                        startWriteGetHapBam
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
            push @{$tagToBamHref->{$tag}}, $_->{splitBam}->{$tag} for @{$V_Href->{PairBamFiles}};
        }
        phasePE_to_contactCount(tagToBamHref => $tagToBamHref);
    }

    # assign haplo to 'UK' end in sEnd-h[x] PE
    ## write getHapBam
    &confirm_sEndSoloHapPE_HapLink;
}

#--- confirm sEnd-h[x] PE haplotype contact ---
## use multiple forks
## here, deal with all phMut-sEnd-h[x],
## as they will be reload to fill the phased-PE-contact for next dEnd-unknown PE hapAssign
sub confirm_sEndSoloHapPE_HapLink{
    # prepare getHapBam name of all source bams
    for my $hapIdx (1 .. $V_Href->{haploCount}){
        my $tag = "phMut-sEnd-h$hapIdx";
        &prepareGetHapBamObj(pairBamHref => $_, tag => $tag, hapID => "h$hapIdx") for @{$V_Href->{PairBamFiles}};
    }

    # start from step after current step
    return if $V_Href->{stepToStart} > 2;

    # all bams or selected
    my @TODOpairBamHref = &getTODOpairBamHrefArray;
    # fork manager
    my $pm;
    my $forkNum = min( $V_Href->{forkNum}, scalar(@TODOpairBamHref) );
    my $fork_DO = ( $forkNum > 1 );
    if($fork_DO){ $pm = new Parallel::ForkManager($forkNum) }
    # load PE from sEnd-h[x].bam
    for my $pairBamHref ( @TODOpairBamHref ){
        # fork job starts
        if($fork_DO){ $pm->start and next }
        # load reads from each hapSplit.bam
        for my $hapIdx (1 .. $V_Href->{haploCount}){
            my $tag = "phMut-sEnd-h$hapIdx";
            # prepare PEsplitStat
            for my $chrTag (qw/ IntraChr InterChr /){
                for my $assignMethod (qw/ ph rd /){
                    $V_Href->{PEsplitStat}->{"$tag.hInter.$chrTag.$assignMethod"} = 0;
                    $V_Href->{PEsplitStat}->{"$tag.h${hapIdx}Intra.$chrTag.$assignMethod"} = 0;
                }
            }
            # prepare getHapBam and file-handle
            &startWriteGetHapBam(pairBamHref => $pairBamHref, tag => $tag);
            # read hapSplit.bam
            my $hapSplitBam = $pairBamHref->{splitBam}->{$tag};
            my $HapLinkHf = { link => {}, stat => {Calculate=>0, LackChrLink=>0, QuickFind=>0} };
            my @subrtOpt = (subrtRef => \&assign_sEndUKend_haplotype,
                            subrtParmAref => [tag => $tag, pairBamHref => $pairBamHref, HapLinkHf => $HapLinkHf]);
            $hapSplitBam->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC', @subrtOpt);
            # close getHapBam file-handle
            $_->stop_write for values %{$pairBamHref->{splitBam}};
            # inform
            my $mark = $hapSplitBam->get_tag;
            my $HapLinkStat = join('; ', map {("$_:$HapLinkHf->{stat}->{$_}")} sort keys %{$HapLinkHf->{stat}});
            stout_and_sterr "[INFO]\t".`date`
                                 ."\tassign haplotype to UK-end of PE-reads from $mark bam OK.\n"
                                 ."\tHapLinkStat: $HapLinkStat\n";
        }
        # write report (2nd part)
        ## previous contents
        my @Content = split /\s+/, `cat $pairBamHref->{PEsplit_report}`;
        my %Content = map {(@Content[$_,$_+1])} grep {$_%2==0 && $Content[$_] !~ /\.[rp][dh]:$/} (0 .. $#Content);
        ## add contents of this step
        open (PESRT, Try_GZ_Write($pairBamHref->{PEsplit_report})) || die "fail write '$pairBamHref->{prefix}' report: $!\n";
        print PESRT "$Content[$_]\t$Content{$Content[$_]}\n" for grep {$_%2==0 && exists $Content{$Content[$_]}} (0 .. $#Content);
        print PESRT "$_:\t$V_Href->{PEsplitStat}->{$_}\n" for grep {!exists $Content{"$_:"}} sort keys %{$V_Href->{PEsplitStat}};
        close PESRT;
        # fork job finishes
        if($fork_DO){ $pm->finish }
    }
    # collect fork jobs
    if($fork_DO){ $pm->wait_all_children }

    # stop at current step
    exit(0) if $V_Href->{stepToStop} == 2;
}

#--- prepare getHapBam object from hapSplit.bam ---
sub prepareGetHapBamObj{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};
    my $tag = $parm{tag};
    my $hapID = $parm{hapID};

    my $splitBamPath = $pairBamHref->{splitBam}->{$tag}->get_filepath;
    my $splitBamTag = $pairBamHref->{splitBam}->{$tag}->get_tag;
    # possible haplo combination
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

#--- return Href-array of all bams or selected ---
sub getTODOpairBamHrefArray{
    # all bams or selected
    return   scalar(keys %{$V_Href->{SelectBamPref}}) == 0
           ? @{$V_Href->{PairBamFiles}}
           : grep exists $V_Href->{SelectBamPref}->{$_->{prefix}}, @{$V_Href->{PairBamFiles}};
}

#--- start writing getHapBam (file-w-handle, SAM-header) ---
sub startWriteGetHapBam{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};
    my $tag = $parm{tag};

    # prepare SAM-header
    my $hapSplitBam = $pairBamHref->{splitBam}->{$tag};
    my $HeadAf = $hapSplitBam->get_SAMheader(samtools => $V_Href->{samtools});
    ## add options of second part in @PG line
    my $PGi = first {$HeadAf->[$_] =~ /\@PG/} (0 .. scalar(@$HeadAf)-1);
    chomp($HeadAf->[$PGi]);
    $HeadAf->[$PGi] .= " -ucrfut $V_Href->{UKreadsFlankRegUnit}";
    $HeadAf->[$PGi] .= " -ucrpha $V_Href->{UKreadsMaxPhasedHetMut}\n";
    my $HeadStr = join('', @$HeadAf);

    # getHap.bam file-handle
    my @ToInform;
    for my $getHapTag (grep /^$tag\./, sort keys %{$pairBamHref->{splitBam}}){
        my $getHapBam = $pairBamHref->{splitBam}->{$getHapTag};
        $getHapBam->start_write(samtools => $V_Href->{samtools});
        $getHapBam->write(content => $HeadStr);
        push @ToInform, $getHapBam->get_filepath;
    }

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
    my $tag = $parm{tag};
    my $pairBamHref = $parm{pairBamHref};
    my $HapLinkHf = $parm{HapLinkHf};

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
    my ($HapLinkC_Href, $mark) = get_rOBpair_HapLinkCount(rOB_a => $hasHap_rOB, rOB_b => $nonHap_rOB, hapRegex => $regex, HapLinkHf => $HapLinkHf);
    # assign HapID to nonHap_rOB
    ## once local region is not phased, reset mark and loads pre-defined HapLink
    my $assignMethod;
    my $isIntraChr = $hasHap_rOB->get_mseg eq $nonHap_rOB->get_mseg;
    if($mark !~ /^RegionPhased/){
        $assignMethod = 'rd';
        if($isIntraChr){
            $HapLinkC_Href->{"$hasHapID,$hasHapID"} = 1;
            $mark .= ';IntraChrSetHapIntra';
        }
        else{
            $HapLinkC_Href->{$_} = 1 for grep /$regex/, @{$V_Href->{allHapComb}};
            $mark .= ';InterChrSetAllAvibHap';
        }
    }
    else{
        $assignMethod = 'ph';
    }
    ## select hapComb
    my @HapComb = sort keys %$HapLinkC_Href;
    my $assHapComb;
    ## if has only single HapComb, just take it
    if(scalar(@HapComb) == 1){
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
    ## remove the pre-set hasHap-ID
    (my $assHapID = $assHapComb) =~ s/$regex//;
    ## assign
    $nonHap_rOB->load_SuppHaplo(hapID => $assHapID, allele_OB => 'NULL');
    $nonHap_rOB->addHapIDtoOptfd; # update
    $nonHap_rOB->add_str_to_optfd(str => "\tXU:Z:$mark");
    # write PE to getHapBam files
    my $getHapTag = $assHapID eq $hasHapID ? "$tag.${assHapID}Intra" : "$tag.hInter";
    $pairBamHref->{splitBam}->{$getHapTag}->write(content => join("\n",@{$pe_OB->printSAM})."\n");
    # stat
    my $chrTag = $isIntraChr ? 'IntraChr' : 'InterChr';
    $V_Href->{PEsplitStat}->{"$getHapTag.$chrTag.$assignMethod"} ++;
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
