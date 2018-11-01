package HaploHiC::PhasedHiC::mergeHaploReads;

use strict;
use warnings;
use Data::Dumper;
use Parallel::ForkManager;
use File::Spec::Functions qw/ catfile /;
use List::Util qw/ min /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use HaploHiC::LoadOn;
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
$VERSION = "0.02";
$DATE = '2018-10-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        merge_haplo_reads
                     /;

#--- merge reads of each haplotype ---
## use multiple forks
sub merge_haplo_reads{
    # prepare file name of merged bam
    for my $pairBamHref ( @{$V_Href->{PairBamFiles}} ){
        for my $tag (sort keys %{$pairBamHref->{bamToMerge}}){
            $pairBamHref->{mergeBam}->{$tag} = catfile($V_Href->{outdir}, $pairBamHref->{prefix}.".$tag.bam");
        }
    }

    # start from step after current step
    return if $V_Href->{stepToStart} > 4;

    # fork manager
    my $pm;
    my $forkNum = min( $V_Href->{forkNum}, scalar(@{$V_Href->{PairBamFiles}}) );
    my $fork_DO = ( $forkNum > 1 );
    if( $fork_DO ){ $pm = new Parallel::ForkManager($forkNum) }
    # merge haploComb bam
    for my $pairBamHref ( @{$V_Href->{PairBamFiles}} ){
        # fork job starts
        if( $fork_DO ){ $pm->start and next; }
        # load reads from each bam
        my $bamPrefix = $pairBamHref->{prefix};
        # merge each type of hapComb
        ## h[x]Intra and hInter
        for my $tag (sort keys %{$pairBamHref->{bamToMerge}}){
            # open file handle of merged bam
            open (MEGBAM, "| $V_Href->{samtools} view -b -o $pairBamHref->{mergeBam}->{$tag}") || die "fail write $bamPrefix merge.bam ($tag): $!\n";
            # merge
            my $bamToMergeAf = $pairBamHref->{bamToMerge}->{$tag};
            ## check existence
            for my $bamToMerge (@$bamToMergeAf){
                warn_and_exit "<ERROR>\tCannot find bam to merge: $bamToMerge\n" unless file_exist(filePath => $bamToMerge);
            }
            ## write SAM-header
            ## copy original SAM-header from last bam
            open (BAMH,"$V_Href->{samtools} view -H $bamToMergeAf->[-1] |") || die"fail read SAM header: $!\n";
            print MEGBAM while(<BAMH>);
            close BAMH;
            ## copy reads
            for my $bamToMerge (@$bamToMergeAf){
                open (BAM,"$V_Href->{samtools} view $bamToMerge |") || die"fail read bam: $!\n";
                print MEGBAM while(<BAM>);
                close BAM;
            }
            # close file handle of merged bam
            close MEGBAM;
            # inform
            stout_and_sterr "[INFO]\t".`date`
                                 ."\tmerge $bamPrefix $tag bam OK.\n";
        }
        # fork job finishes
        if( $fork_DO ){ $pm->finish; }
    }
    # collect fork jobs
    if( $fork_DO ){ $pm->wait_all_children; }

    # stop at current step
    exit(0) if $V_Href->{stepToStop} == 4;
}

#--- 
1; ## tell the perl script the successful access of this module.
