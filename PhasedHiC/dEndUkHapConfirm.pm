package HaploHiC::PhasedHiC::dEndUkHapConfirm;

use strict;
use warnings;
use List::Util qw/ max min sum first /;
use Data::Dumper;
use HaploHiC::LoadOn;
use HaploHiC::PhasedHiC::sEndSoloHapConfirm qw/ prepareGetHapBamObj getChrPairBam chrPair_forkSetting chrPair_loadContacts chrPair_HapConfirmWork chrPair_mergeResult /;
use HaploHiC::PhasedHiC::phasedPEtoContact qw/ get_rOBpair_HapLinkCount /;

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
$VERSION = "0.25";
$DATE = '2019-03-25';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        dEndUK_get_HapLink
                        confirm_dEndU_Hap
                        chrPair_dEndU_assignHap
                     /;

#--- assign dual-side unknown PE to to certain haplo-contacts ---
## unknown
sub dEndUK_get_HapLink{
    # prepare getHapBam name of all source bams
    my $tag = 'unknown';
    prepareGetHapBamObj(pairBamHref => $_, tag => $tag) for @{$V_Href->{PairBamFiles}};

    # current step is in pipeline
    if($V_Href->{stepToStart} <= 3){
        # split dEnd-U.bam to chrPair.bam
        ## fork at run level
        my %preChr;
        getChrPairBam(preChrHf => \%preChr, type => 'dEndU');
        # assign haplo to 'UK' end in dEnd-U PE
        ## fork at chrPair level (preChr + chrRel)
        &confirm_dEndU_Hap(preChrHf => \%preChr);
        # stop at current step
        exit(0) if $V_Href->{stepToStop} == 3;
    }
}

#--- confirm unknown PE haplotype contact ---
## use multiple forks at preChr level
sub confirm_dEndU_Hap{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $preChrHf = $parm{preChrHf};

    my @preChrOpt = grep exists $preChrHf->{$_->[0]}->{$_->[1]}, # filter
                    map {( [$_, 'intraChr'], [$_, 'interChr'] )}
                    grep exists $preChrHf->{$_}, @{$V_Href->{sortedChr}};
    my $pm = chrPair_forkSetting;
    # each chrPair of dEndU.bam
    for my $preChrOpt ( @preChrOpt ){
        # fork job starts
        $pm->start(join(';',@$preChrOpt)) and next;
        eval{
            # options
            my ($preChr, $chrRel) = @$preChrOpt;
            # only one chrPair?
            my @chrPairTag = keys %{$preChrHf->{$preChr}->{$chrRel}};
            my @chrPairOpt = @chrPairTag == 1 ? (chrPair => $chrPairTag[0]) : ();
            # load phased-contact of this chrPair
            chrPair_loadContacts(preChr => $preChr, chrRel => $chrRel, type => 'dEndU', @chrPairOpt);
            # confirm dEndU haplo of this chrPair
            chrPair_HapConfirmWork(preChr => $preChr, chrRel => $chrRel, assignHapSubrtRef => \&chrPair_dEndU_assignHap);
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
    chrPair_mergeResult(type => 'dEndU');
}

#--- assign haplotype to dEndUK PE UK-end ---
sub chrPair_dEndU_assignHap{
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
        my $getHapTag = $assHapComb[0] eq $assHapComb[-1] ? "$assHapComb[0]Intra" : "hInter";
        # write PE to getHapBam files
        $getHapBamHf->{$getHapTag}->write(content => join("\n",@{$pe_OB->printSAM(keep_all=>1)})."\n");
        # stat
        my $chrTag = $isIntraChr ? 'IntraChr' : 'InterChr';
        $V_Href->{PEsplitStat}->{"$getHapTag.$chrTag.$assignMethod"} ++;
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
