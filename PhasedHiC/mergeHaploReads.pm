package HaploHiC::PhasedHiC::mergeHaploReads;

use strict;
use warnings;
use Data::Dumper;
use Parallel::ForkManager;
use File::Spec::Functions qw/ catfile /;
use List::Util qw/ min /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use HaploHiC::LoadOn;
use HaploHiC::PhasedHiC::splitPairBam qw/ forkSetting /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              merge_haplo_reads
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::mergeHaploReads';
#----- version --------
$VERSION = "0.06";
$DATE = '2019-04-06';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        merge_haplo_reads
                        mergeReadsOfEachHapComb
                        mergeStatOfPhasedLocalRegion
                     /;

#--- merge reads of each haplotype ---
## use multiple forks
sub merge_haplo_reads{
    # prepare file name of merged bam
    for my $pairBamHref ( @{$V_Href->{PairBamFiles}} ){
        my $pairBamPrefix = $pairBamHref->{prefix};
        for my $tag (sort keys %{$pairBamHref->{bamToMerge}}){
            my $mergeBamPath = catfile($V_Href->{outdir}, $pairBamPrefix.".$tag.bam");
            $pairBamHref->{mergeBam}->{$tag} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $mergeBamPath, tag => "$pairBamPrefix $tag");
        }
    }

    # start from step after current step
    return if $V_Href->{stepToStart} > 4;

    # fork manager
    my ($pm, $fork_DO) = forkSetting;
    # merge haploComb bam
    for my $pairBamHref ( @{$V_Href->{PairBamFiles}} ){
        # fork job starts
        if($fork_DO){ $pm->start($pairBamHref->{prefix}) and next }
        eval{
            # merge each type of hapComb
            &mergeReadsOfEachHapComb(pairBamHref => $pairBamHref);
            # merge stat data of phased local region
            &mergeStatOfPhasedLocalRegion(pairBamHref => $pairBamHref);
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
    exit(0) if $V_Href->{stepToStop} == 4;
}

#--- merge reads of each hapComb to mergeBam ---
sub mergeReadsOfEachHapComb{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};

    # merge each type of hapComb
    ## h[x]Intra and hInter
    for my $tag (sort keys %{$pairBamHref->{bamToMerge}}){
        # start writing mergedBam
        my $mergeBam = $pairBamHref->{mergeBam}->{$tag};
        $mergeBam->start_write(samtools => $V_Href->{samtools});
        # merge
        my $bamToMergeAf = $pairBamHref->{bamToMerge}->{$tag};
        ## write SAM-header
        ## copy original SAM-header from last bam
        $mergeBam->write(content => $_) for @{ $bamToMergeAf->[-1]->header_Af(samtools => $V_Href->{samtools}) };
        ## copy reads
        for my $bamToMerge (@$bamToMergeAf){
            ## official process!!!
            my $fh = $bamToMerge->start_read(samtools => $V_Href->{samtools});
            $mergeBam->write(content => $_) while(<$fh>);
            close $fh;
        }
        # stop writing mergedBam
        $mergeBam->stop_write;
        # inform
        my $mark = $mergeBam->tag;
        stout_and_sterr "[INFO]\t".`date`
                             ."\tmerge $mark bam OK.\n";
    }
}

#--- merge stat data of phased local region ---
sub mergeStatOfPhasedLocalRegion{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pairBamHref = $parm{pairBamHref};

    my @statFiles = glob catfile($pairBamHref->{workSpace}, '*.statOfPhasedLocReg.gz');
    return if @statFiles == 0;
    # merge
    my $mergeStatFile = catfile($V_Href->{outdir}, $pairBamHref->{prefix}.'.statOfPhasedLocReg.gz');
    ## header
    `gzip -cd $statFiles[0] | grep '^#' | gzip -c > $mergeStatFile`;
    ## content
    `zcat @statFiles | grep -v '^#' | sort -k 2,2 -k 3n,3 -k 4n,4 | gzip -c >> $mergeStatFile`;
}

#--- 
1; ## tell the perl script the successful access of this module.
