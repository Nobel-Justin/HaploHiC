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
use HaploHiC::PhasedHiC::phasedMutWork qw/ PosToPhasedMut /;
use HaploHiC::PhasedHiC::sEndSoloHapConfirm qw/ getTODOpairBamHrefArray prepare_getHapBamName prepare_getHapBamFH /;
use HaploHiC::PhasedHiC::phasedPEtoContact qw/ phasePE_to_contactCount smartBam_PEread get_rOBpair_HapLinkCount /;

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
$VERSION = "0.14";
$DATE = '2018-11-22';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        dEndUK_get_HapLink
                        assign_dEndUKend_haplotype
                     /;

#--- assign dual-side unknown PE to to certain haplo-contacts ---
## unknown
## use multiple forks
sub dEndUK_get_HapLink{
    # load contacts from PE-reads of sEndSoloHapPE getHap.bam
    ## sEnd-h[x].h[x]Intra and sEnd-h[x].hInter
    ## only run when starts before step4 (readsMerge)
    if($V_Href->{stepToStart} < 4){
        my $tagToBamHref = {};
        my @tags = map {("phMut-sEnd-h$_.h${_}Intra", "phMut-sEnd-h$_.hInter")} ( 1 .. $V_Href->{haploCount} );
        for my $tag ( @tags ){
            push @{$tagToBamHref->{$tag}}, [ $_->{prefix}, $_->{splitBam}->{$tag} ] for @{$V_Href->{PairBamFiles}};
        }
        phasePE_to_contactCount( tagToBamHref => $tagToBamHref ); 
    }

    # prepare getHap.bam name of all source bams
    my $tag = 'unknown';
    prepare_getHapBamName(pairBamHref => $_, tag => $tag, splitBam => $_->{splitBam}->{$tag}) for @{$V_Href->{PairBamFiles}};

    # start from step after current step
    return if $V_Href->{stepToStart} > 3;

    # all bams or selected
    my @TODOpairBamHref = getTODOpairBamHrefArray;
    # fork manager
    my $pm;
    my $forkNum = min( $V_Href->{forkNum}, scalar(@TODOpairBamHref) );
    my $fork_DO = ( $forkNum > 1 );
    if( $fork_DO ){ $pm = new Parallel::ForkManager($forkNum) }
    # load PE from unknown.bam
    for my $pairBamHref ( @TODOpairBamHref ){
        # fork job starts
        if( $fork_DO ){ $pm->start and next; }
        # load reads from each bam
        my $bamPrefix = $pairBamHref->{prefix};
        my $mark = "'$bamPrefix' $tag";
        my $hapSplitBam = $pairBamHref->{splitBam}->{$tag};
        # check existence
        warn_and_exit "<ERROR>\tCannot find $tag bam: $hapSplitBam\n" unless file_exist(filePath => $hapSplitBam);
        # prepare getHap.bam and file-handle
        prepare_getHapBamFH(pairBamHref => $pairBamHref, tag => $tag, splitBam => $hapSplitBam);
        # read hapSplit.BAM
        ## FLAG: 0x100(sd), 0x400(d), 0x800(sp)
        ## WEDO: -F 0x400(d), as 'sd' and 'sp' alignments is important in HiC
        smartBam_PEread(bam => $hapSplitBam, mark => $mark, viewOpt => '-F 0x400',
                        subrtRef => \&assign_dEndUKend_haplotype,
                        subrtParmAref => [tag=>$tag, pairBamHref => $pairBamHref]);
        # close getHap.bam file-handle
        close $_ for values %{$pairBamHref->{splitBamFH}};
        $pairBamHref->{splitBamFH} = {}; # empty splitBamFH key
        # inform
        stout_and_sterr "[INFO]\t".`date`
                             ."\tassign haplotype to UK-end of PE-reads from $mark bam OK.\n";
        # write report (3rd part)
        ## previous contents
        my @Content = split /\s+/, `cat $pairBamHref->{PEsplit_report}`;
        my %Content = map {(@Content[$_,$_+1])} grep {$_%2==0 && $Content[$_] !~ /^unknown.h[I\d]/} (0 .. $#Content);
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
    exit(0) if $V_Href->{stepToStop} == 3;
}

#--- assign haplotype to dEndUK PE UK-end ---
sub assign_dEndUKend_haplotype{
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
    # get haplo link count of paired rOB, <do not need sort again>
    my ($HapLinkC_Href, $mark) = get_rOBpair_HapLinkCount(rOB_a => $rOB_sortAref->[0], rOB_b => $rOB_sortAref->[-1], skipSort => 1);
    # assign HapID to both nonHap_rOB
    ## once local region is not phased, reset mark and loads pre-defined HapLink
    my $assignMethod;
    my $isIntraChr = $rOB_sortAref->[0]->get_mseg eq $rOB_sortAref->[-1]->get_mseg;
    if($mark !~ /^RegionPhased/){
        $assignMethod = 'rd';
        if($isIntraChr){
            $HapLinkC_Href->{$_} = 1 for @{$V_Href->{intraHapComb}};
            $mark .= ';IntraChrSetHapIntra';
        }
        else{
            $HapLinkC_Href->{$_} = 1 for @{$V_Href->{allHapComb}};
            $mark .= ';InterChrSetAllHap';
        }
    }
    else{
        $assignMethod = 'ph';
    }
    ## select hapComb
    my $assHapComb;
    my @HapComb = sort keys %$HapLinkC_Href;
    # if has only single HapComb, just take it
    if(scalar(@HapComb) == 1){
        $assHapComb = $HapComb[0];
        # add mark
        $mark .= ";SoloHapComb;[$assHapComb](C:$HapLinkC_Href->{$assHapComb})";
    }
    else{ # luck draw
        ## prepare haplotype's luck-draw interval
        my $allCount = 0;
        my %hapCombDraw;
        for my $hapComb (@HapComb){
            $allCount += $HapLinkC_Href->{$hapComb};
            $hapCombDraw{$hapComb} = $allCount;
        }
        ## random pick
        my $luck_draw = int(rand($allCount));
        $assHapComb = first { $hapCombDraw{$_} > $luck_draw } sort keys %hapCombDraw;
        ## add mark
        $mark .= ";[$_](C:$HapLinkC_Href->{$_},D:$hapCombDraw{$_})" for keys %hapCombDraw;
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
    # write PE to getHap.bam files
    print {$pairBamHref->{splitBamFH}->{$getHapTag}} "$_\n" for @{$pe_OB->printSAM};
    # stat
    my $chrTag = $isIntraChr ? 'IntraChr' : 'InterChr';
    $V_Href->{PEsplitStat}->{"$getHapTag.$chrTag.$assignMethod"} ++;
}

#--- 
1; ## tell the perl script the successful access of this module.
