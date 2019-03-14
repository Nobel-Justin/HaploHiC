package HaploHiC::PhasedHiC::dEndUkHapConfirm;

use strict;
use warnings;
use List::Util qw/ max min sum first /;
use Data::Dumper;
use Parallel::ForkManager;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Write /;
use BioFuse::Util::Index qw/ Pos2Idx /;
use HaploHiC::LoadOn;
use HaploHiC::PhasedHiC::splitPairBam qw/ forkSetting getTODOpairBamHrefArray /;
use HaploHiC::PhasedHiC::sEndSoloHapConfirm qw/ prepareGetHapBamObj startWriteGetHapBam writeStatOfPhasedLocalRegion /;
use HaploHiC::PhasedHiC::phasedPEtoContact qw/ phasePE_to_contactCount get_rOBpair_HapLinkCount /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              dEndUK_get_HapLink
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::dEndUkHapConfirm';
#----- version --------
$VERSION = "0.24";
$DATE = '2019-03-14';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        dEndUK_get_HapLink
                        confirm_dEndUkHapPE_HapLink
                        dEndUkHapConfirmWork
                        assign_dEndUKend_haplotype
                     /;

#--- assign dual-side unknown PE to to certain haplo-contacts ---
## unknown
sub dEndUK_get_HapLink{
    # load contacts from PE-reads of sEndSoloHapPE getHap.bam
    ## sEnd-h[x].h[x]Intra and sEnd-h[x].hInter
    ## only run when starts before step4 (readsMerge)
    if($V_Href->{stepToStart} < 4){
        my $tagToBamHref = {};
        my @tags = map {("phMut-sEnd-h$_.h${_}Intra", "phMut-sEnd-h$_.hInter")} ( 1 .. $V_Href->{haploCount} );
        for my $tag ( @tags ){
            push @{$tagToBamHref->{$tag}}, $_->{splitBam}->{$tag} for @{$V_Href->{PairBamFiles}};
        }
        phasePE_to_contactCount(tagToBamHref => $tagToBamHref);
    }

    # assign haplo to 'UK' end in dEndUK PE
    ## write getHapBam
    &confirm_dEndUkHapPE_HapLink;
}

#--- confirm unknown PE haplotype contact ---
## use multiple forks
sub confirm_dEndUkHapPE_HapLink{
    # prepare getHapBam name of all source bams
    my $tag = 'unknown';
    prepareGetHapBamObj(pairBamHref => $_, tag => $tag) for @{$V_Href->{PairBamFiles}};

    # start from step after current step
    return if $V_Href->{stepToStart} > 3;

    # all bams or selected
    my @TODOpairBamHref = getTODOpairBamHrefArray;
    # fork manager
    my ($pm, $fork_DO) = forkSetting;
    # load PE from unknown.bam
    for my $pairBamHref ( @TODOpairBamHref ){
        # fork job starts
        if($fork_DO){ $pm->start($pairBamHref->{prefix}) and next }
        eval{
            # confrim haplotype of 'U' end
            &dEndUkHapConfirmWork(pairBamHref => $pairBamHref);
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
    exit(0) if $V_Href->{stepToStop} == 3;
}

#--- main work of this module ---
sub dEndUkHapConfirmWork{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};

    my $tag = 'unknown';
    # prepare PEsplitStat
    for my $chrTag (qw/ IntraChr InterChr /){
        for my $assignMethod (qw/ ph rd /){
            $V_Href->{PEsplitStat}->{"$tag.hInter.$chrTag.$assignMethod"} = 0;
            $V_Href->{PEsplitStat}->{"$tag.h${_}Intra.$chrTag.$assignMethod"} = 0 for (1 .. $V_Href->{haploCount});
        }
    }
    # prepare getHapBam and file-handle
    startWriteGetHapBam(pairBamHref => $pairBamHref, tag => $tag);
    # read unknown.bam
    my $hapSplitBam = $pairBamHref->{splitBam}->{$tag};
    my $HapLinkHf = { link => {}, stat => {Calculate=>0, LackChrLink=>0, QuickFind=>0} };
    $V_Href->{LocRegPhased} = {}; # reset
    my @subrtOpt = (subrtRef => \&assign_dEndUKend_haplotype,
                    subrtParmAref => [tag => $tag, pairBamHref => $pairBamHref, HapLinkHf => $HapLinkHf]);
    $hapSplitBam->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC', deal_peOB_pool => 1, @subrtOpt);
    # close getHapBam file-handle
    $_->stop_write for values %{$pairBamHref->{splitBam}};
    # write stat of phased-local-region (size and linkCount)
    writeStatOfPhasedLocalRegion(tag => 'unknown', hapSplitBam => $hapSplitBam);
    # inform
    my $mark = $hapSplitBam->get_tag;
    my $HapLinkStat = join('; ', map {("$_:$HapLinkHf->{stat}->{$_}")} sort keys %{$HapLinkHf->{stat}});
    stout_and_sterr "[INFO]\t".`date`
                         ."\tassign haplotype to UK-end of PE-reads from $mark bam OK.\n"
                         ."\t$mark HapLinkStat: $HapLinkStat\n";
    # write report (3rd part)
    ## previous contents
    my @Content = split /\s+/, `cat $pairBamHref->{PEsplit_report}`;
    my %Content = map {(@Content[$_,$_+1])} grep {$_%2==0 && $Content[$_] !~ /^unknown.h[I\d]/} (0 .. $#Content);
    ## add contents of this step
    open (PESRT, Try_GZ_Write($pairBamHref->{PEsplit_report})) || die "fail write '$pairBamHref->{prefix}' report: $!\n";
    print PESRT "$Content[$_]\t$Content{$Content[$_]}\n" for grep {$_%2==0 && exists $Content{$Content[$_]}} (0 .. $#Content);
    print PESRT "$_:\t$V_Href->{PEsplitStat}->{$_}\n" for grep {!exists $Content{"$_:"}} sort keys %{$V_Href->{PEsplitStat}};
    close PESRT;
}

#--- assign haplotype to dEndUK PE UK-end ---
sub assign_dEndUKend_haplotype{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};
    my $pairBamHref = $parm{pairBamHref};
    my $tag = $parm{tag};
    my $HapLinkHf = $parm{HapLinkHf};

    for my $pe_OB (@$pe_OB_poolAf){
        # get chr-pos ascending sorted all mapped reads_OB
        my $rOB_sortAref = $pe_OB->get_sorted_reads_OB(chrSortHref => $V_Href->{ChrThings}, chrSortKey  => 'turn');
        # recover SuppHaplo attribute
        $_->recover_SuppHaploAttr for @$rOB_sortAref;
        # get haplo link count of paired rOB, <do not need sort again>
        my $LocRegInfoAf = get_rOBpair_HapLinkCount(rOB_a => $rOB_sortAref->[0], rOB_b => $rOB_sortAref->[-1], skipSort => 1, HapLinkHf => $HapLinkHf);
        my ($HapLinkC_Hf, $mark, $assignMethod, $modBool) = @$LocRegInfoAf;
        # assign HapID to both nonHap_rOB
        ## once local region is not phased, reset mark and loads pre-defined HapLink
        my $isIntraChr = $rOB_sortAref->[0]->get_mseg eq $rOB_sortAref->[-1]->get_mseg;
        unless($modBool){
            if($assignMethod eq 'rd'){
                if($isIntraChr){
                    $HapLinkC_Hf->{$_} = 1 for @{$V_Href->{intraHapComb}};
                    $mark .= ';IntraChrSetHapIntra';
                }
                else{
                    $HapLinkC_Hf->{$_} = 1 for @{$V_Href->{allHapComb}};
                    $mark .= ';InterChrSetAllHap';
                }
            }
            else{
                # uniform addition
                if($V_Href->{uniformAddRatioForHapComb}){
                    my $add = sum(values %$HapLinkC_Hf) * $V_Href->{uniformAddRatioForHapComb};
                    $HapLinkC_Hf->{$_} += $add for @{$V_Href->{allHapComb}};
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
        ## assign
        my @assHapComb = split /,/, $assHapComb;
        for my $i (0, -1){
            $rOB_sortAref->[$i]->load_SuppHaplo(hapID => $assHapComb[$i], allele_OB => 'NULL');
            $rOB_sortAref->[$i]->addHapIDtoOptfd; # update
            $rOB_sortAref->[$i]->add_str_to_optfd(str => "\tXU:Z:$mark");
        }
        # output
        my $getHapTag = $assHapComb[0] eq $assHapComb[-1] ? "$tag.$assHapComb[0]Intra" : "$tag.hInter";
        # write PE to getHapBam files
        $pairBamHref->{splitBam}->{$getHapTag}->write(content => join("\n",@{$pe_OB->printSAM(keep_all=>1)})."\n");
        # stat
        my $chrTag = $isIntraChr ? 'IntraChr' : 'InterChr';
        $V_Href->{PEsplitStat}->{"$getHapTag.$chrTag.$assignMethod"} ++;
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
