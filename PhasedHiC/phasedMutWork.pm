package HaploHiC::PhasedHiC::phasedMutWork;

use strict;
use warnings;
use List::Util qw/ max min sum first /;
use Data::Dumper;
use HaploHiC::GetPath qw/ GetPath /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;
use BioFuse::Util::String qw/ getStrRepUnit /;
use BioFuse::Util::Index qw/ Pos2Idx FindOverlapIdxRegion /;
use BioFuse::BioInfo::Objects::Allele::PhasedMut_OB;
use HaploHiC::LoadOn;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              load_phased_VCF
              release_phaseMut_OB
              PosToPhasedMut
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::phasedMutWork';
#----- version --------
$VERSION = "0.18";
$DATE = '2025-10-04';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        load_phased_VCF
                        GetPhaseMutEdgeDist
                        release_phaseMut_OB
                        filterClosePhasedInDel
                        PosToPhasedMut
                     /;

#--- load phased mutations from VCF ---
## requires coordinate-sorted in each chr
sub load_phased_VCF{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $sample = $parm{sample} || undef;

    my $sample_colIdx = defined $sample ? undef : 9;
    my $last_chr;
    my $last_pos = -1; # initialize
    open (VCF, Try_GZ_Read($V_Href->{phased_vcf})) || die "fail to read phased VCF: $!\n";
    while(<VCF>){
        if(/^#/){ # header
            if(!/^##/ && defined $sample){ # theme
                my @colHead = split;
                $sample_colIdx = first { $colHead[$_] =~ /^$sample$/i } (0 .. $#colHead);
                unless(defined $sample_colIdx){
                    warn_and_exit "<ERROR>\tcannot find VCF column for required sample ($sample).\n";
                }
            }
            next;
        }
        my ($chr,$pos,$refSeq,$alt,$pass,$info) = (split)[0,1,3,4,6,$sample_colIdx];
        # not chr wanted
        if( !exists $V_Href->{ChrThings}->{$chr} ){
            $V_Href->{VCFmutStat}->{i01_skipChr} ++;
            next;
        }
        # GT:DP:ADALL:AD:GQ:IGT:IPS:PS
        # not phased
        if( $info !~ /^\d\|\d/ ){
            $V_Href->{VCFmutStat}->{i02_unPhaMut} ++;
            next;
        }
        # allele idx
        # my ($allele_idx) = ( $info =~ /^([^:]+):/ );
        my ($allele_idx) = ( $info =~ /^(\d\|\d)/ );
        my  @allele_idx  = split /\|/, $allele_idx;
        # Hom phased, here is for diploid
        # how to update for triploid or more?
        my $is_Hom = ($allele_idx[0] == $allele_idx[1]);
        $V_Href->{VCFmutStat}->{i03_homPhaMut} ++ if( $is_Hom );
        # check
        if( $V_Href->{haploCount} != scalar(@allele_idx) ){
            warn_and_exit "<ERROR>\tencounter phased allele(s) ($chr:$pos) having count not $V_Href->{haploCount} (option '-hapct').\n";
        }
        # test chr-inner pos-sorted
        if(    $last_pos >  $pos
            && $last_chr eq $chr
        ){
            warn_and_exit "<ERROR>\tthe phased_vcf is not coordinates sorted.\n";
        }
        # update
        $last_pos = $pos;
        $last_chr = $chr;
        # alleles
        my @allele = ($refSeq, split /,/,$alt);
        # want InDel?
        my $is_InDel = ( sum( map {length($_)} @allele ) != scalar(@allele) );
        if( $is_InDel && !$V_Href->{use_InDel} ){
            $V_Href->{VCFmutStat}->{i04_skipINDEL} ++;
            next;
        }
        # allow hom-InDel keep to filter close-InDel
        if( !$is_InDel && $is_Hom ){
            next;
        }
        # accept this phased-alt, and start to record
        $V_Href->{VCFmutStat}->{i05_iniPhaMut} ++;
        $V_Href->{VCFmutStat}->{i05_hetPhaMut} ++ if( !$is_Hom );
        my $pMutNO = $V_Href->{VCFmutStat}->{i05_iniPhaMut};
        # create PhasedMut_OB
        my $PhasedMut_OB = BioFuse::BioInfo::Objects::Allele::PhasedMut_OB->new( chr=>$chr, pos=>$pos, phMutNO=>$pMutNO );
        # ref-info
        my $refBase = substr($refSeq,0,1);
        my $refTail = substr($refSeq,1); # might have
        # load different alleles
        my $phasedType = '';
        my $hasRefType = 0;
        for my $haplo_NO ( 1 .. $V_Href->{haploCount} ){
            my $allele_idx = $allele_idx[$haplo_NO-1];
            my $allele = $allele[$allele_idx];
            my ($alleleType, $alleleSeq) = ('null', $allele); # initialize
            # is ref
            if( $allele_idx eq 0 ){
                $alleleType = 'ref';
                $alleleSeq  = $refBase; # always store the ref-base at this position
                $hasRefType = 1;
            }
            # snv, keep the first base in $allele
            elsif( length($allele) == length($refSeq) ){ # length is might not one!
                $alleleType = 'snv';
                $alleleSeq = substr($alleleSeq,0,1);
            }
            # insertion behind this position
            # take the inserted-seq enclosed/prefixed in $allele
            elsif( length($allele) > length($refSeq) ){
                $alleleType = 'ins';
                $alleleSeq =~ s/^$refBase//; # first one is the ref-base, not in inserted-seq
                $alleleSeq =~ s/$refTail$// if( $refTail );
            }
            # deletion behind this position
            # take the deleted-seq enclosed/prefixed in $refBase
            elsif( length($allele) < length($refSeq) ){
                $alleleType = 'del';
                $alleleSeq = $refTail; # reset as $refTail
                my $alleleTail = substr($allele,1); # might have
                $alleleSeq =~ s/$alleleTail$// if( $alleleTail );
            }
            # to PhasedMut_OB
            $PhasedMut_OB->load_haplotype_allele( haplo_ID=>"h$haplo_NO", alleleType => $alleleType, alleleSeq => $alleleSeq );
            # phased-allele scenarios
            $phasedType .= $haplo_NO . uc(substr($alleleType,0,1));
        }
        # if all haplotypes are not ref, load ref
        $PhasedMut_OB->load_ref_allele( alleleSeq => $refBase ) unless( $hasRefType );
        # phased-allele scenarios
        $V_Href->{VCFmutStat}->{"i06_hetPha_$phasedType"} ++ if( !$is_Hom );
        # record
        my $PosIdx = Pos2Idx(pos => $pos, winSize => $V_Href->{phaMutPosWinSize});
        push @{$V_Href->{PhasedMut}->{$chr}->{$PosIdx}}, $PhasedMut_OB;
    }
    close VCF;
    # check
    if( $V_Href->{VCFmutStat}->{i05_hetPhaMut} == 0 ){
        warn_and_exit "<ERROR>\tno valid phased mutations found in VCF.\n".Dumper($V_Href->{VCFmutStat});
    }

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tload phased mutations from VCF OK.\n" unless $V_Href->{skipS01Report};

    # InDel things
    if( $V_Href->{use_InDel} ){
        # read reference context, determine read-edge distance to accept allele (has InDel)
        read_fasta_file(FaFile=>$V_Href->{GenomeRefFa}, subrtRef=>\&GetPhaseMutEdgeDist);
        # inform
        stout_and_sterr "[INFO]\t".`date`
                             ."\tdetermine read-edge distance for InDel alleles OK.\n" unless $V_Href->{skipS01Report};
        # discard close phased InDels
        &filterClosePhasedInDel;
    }

    # write report
    unless($V_Href->{skipS01Report}){
        delete $V_Href->{VCFmutStat}->{i05_iniPhaMut};
        $V_Href->{phasedMut_report} = GetPath(filekey => 'phasedMut_report');
        open (VCFRT, Try_GZ_Write($V_Href->{phasedMut_report})) || die "fail write VCF_load.phased_mut.report: $!\n";
        print VCFRT "$_:\t$V_Href->{VCFmutStat}->{$_}\n" for sort keys %{$V_Href->{VCFmutStat}};
        close VCFRT;
    }
}

#--- determine read-edge distance to accept allele (has InDel) based on reference context ---
## default is alleleSeqLen + delta, the delta is the InDelSeq's own-unit repTime (>=2), which should represent the local repeatability
## insertion: p3_Dist = default + refRepPartLen
## deletion:  p3_Dist = default + refRepPartLen - delLen
sub GetPhaseMutEdgeDist{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chr = $parm{segName};
    my $chrSeq_Sref = $parm{segSeq_Sref};

    return unless (exists $V_Href->{PhasedMut}->{$chr});

    for my $PosIdx (sort {$a<=>$b} keys %{$V_Href->{PhasedMut}->{$chr}}){
        for my $PhasedMut_OB (@{$V_Href->{PhasedMut}->{$chr}->{$PosIdx}}){
            next unless $PhasedMut_OB->has_indel; # any one is OK
            my $refContext = substr($$chrSeq_Sref, $PhasedMut_OB->get_pos, 500);
            my $maxEdgeDist = 0;
            for my $type ('ins', 'del'){
                next unless $PhasedMut_OB->has_alleleType(typeAref=>[$type]);
                my $hapAref = $PhasedMut_OB->get_alleleTypeToHapAref(type=>$type);
                for my $hapID (@$hapAref){
                    my $alleleSeq = $PhasedMut_OB->get_hapInfo(haplo_ID=>$hapID)->[1];
                    my ($repTime, $repUnit) = getStrRepUnit(str=>$alleleSeq, mode=>'w');
                    my $p3_Dist = length($alleleSeq) + max($repTime, 2);
                    if( $refContext =~ /^($repUnit)+/i ){
                        $p3_Dist += length($&);
                        $p3_Dist -= length($alleleSeq) if( $type eq 'del' );
                    }
                    $maxEdgeDist = max( $maxEdgeDist, $p3_Dist );
                }
            }
            # record available edgeDist once larger than user setting
            if($maxEdgeDist > $V_Href->{min_distToRedge}){
                $PhasedMut_OB->load_judgeReadEdgeDist(readEdgeDist=>$maxEdgeDist);
                # warn join(',', $PhasedMut_OB->get_chr, $PhasedMut_OB->get_pos, $maxEdgeDist)."\n$refContext\n";
                $V_Href->{VCFmutStat}->{i07_INDEL_rEdge} ++ unless $PhasedMut_OB->is_hom;
            }
        }
    }
}

#--- release phaseMut_OB space for saving memory ---
sub release_phaseMut_OB{
    for my $chr ( keys %{$V_Href->{PhasedMut}} ){
        for my $Aref (values %{$V_Href->{PhasedMut}->{$chr}}){
            $_->release_memory for @$Aref;
        }
    }
}

#--- discard close phased InDels ---
# hom-InDels are also discarded
sub filterClosePhasedInDel{

    # output only in 1st ('splitBam') step
    my $closeInDelReport = GetPath(filekey => 'closeInDelReport');
    open (CLIDRT, Try_GZ_Write($closeInDelReport)) || die "fail write VCF_load.phased_mut.closeInDel.list: $!\n" unless $V_Href->{skipS01Report};

    for my $chr ( sort keys %{$V_Href->{PhasedMut}} ){
        my %PhaMutFilter;
        my $PhaMutOB_l = undef; # globally on this chr
        my $PosIdx_l = undef;
        my $ArrayIdx_l = undef;
        # each PosIdx
        my $chrHref = $V_Href->{PhasedMut}->{$chr};
        for my $PosIdx_n ( sort {$a<=>$b} keys %$chrHref ){
            my $PhaMutOB_Aref = $chrHref->{$PosIdx_n};
            my $count = scalar @$PhaMutOB_Aref;
            # find alleles have InDel
            my @InDel_i = grep $PhaMutOB_Aref->[$_]->has_indel, ( 0 .. $count-1 );
            next unless( scalar @InDel_i ); # no InDel
            # find the first allele has InDel (globally on this chr)
            unless( defined $PhaMutOB_l ){
                $PhaMutOB_l = $PhaMutOB_Aref->[$InDel_i[0]];
                $PosIdx_l = $PosIdx_n;
                $ArrayIdx_l = shift @InDel_i;
            }
            # test distance on Haplotype
            for my $ArrayIdx_n ( @InDel_i ){
                my $PhaMutOB_n = $PhaMutOB_Aref->[$ArrayIdx_n];
                # have close distance
                if( ($PhaMutOB_n->get_pos - $PhaMutOB_l->get_pos) <= ($V_Href->{min_InDelDist} + $PhaMutOB_l->get_judgeReadEdgeDist) ){
                    my $filterBool = 1;
                    # whether has InDel on same Haplotype
                    if( $V_Href->{closeInDelSameHap} ){
                        my %InDelHap_l = map {($_,1)} ( @{$PhaMutOB_l->get_alleleTypeToHapAref(type=>'ins')}, @{$PhaMutOB_l->get_alleleTypeToHapAref(type=>'del')} );
                        my %InDelHap_i = map {($_,1)} ( @{$PhaMutOB_n->get_alleleTypeToHapAref(type=>'ins')}, @{$PhaMutOB_n->get_alleleTypeToHapAref(type=>'del')} );
                        my $ovHaplo = first { exists $InDelHap_i{$_} } keys %InDelHap_l;
                        $filterBool = ( defined $ovHaplo ? 1 : 0 ); # same haplotype (at least one)
                    }
                    if( $filterBool ){
                        $PhaMutFilter{$PosIdx_l}{$ArrayIdx_l} = 1 if( !exists $PhaMutFilter{$PosIdx_l}{$ArrayIdx_l} );
                        $PhaMutFilter{$PosIdx_n}{$ArrayIdx_n} = 1;
                    }
                }
                # update
                $PhaMutOB_l = $PhaMutOB_n;
                $PosIdx_l = $PosIdx_n;
                $ArrayIdx_l = $ArrayIdx_n;
            }
        }
        # discard close-InDels
        for my $PosIdx (sort {$a<=>$b} keys %PhaMutFilter){
            my @delPhaMutOB;
            for my $ArrayIdx (sort {$b<=>$a} keys %{$PhaMutFilter{$PosIdx}}){
                unshift @delPhaMutOB, splice( @{$chrHref->{$PosIdx}}, $ArrayIdx, 1 );
                $V_Href->{VCFmutStat}->{i08_INDEL_close} ++ unless $delPhaMutOB[0]->is_hom; # only count Het-InDel
            }
            # alert, output only in 1st ('splitBam') step
            for my $delPhaMutOB (@delPhaMutOB){
                print CLIDRT join("\t", @{$delPhaMutOB->get_infoSummary}) . "\n" unless $V_Href->{skipS01Report};
            }
            # become empty!
            if( scalar(@{$chrHref->{$PosIdx}}) == 0 ){
                delete $chrHref->{$PosIdx};
            }
        }
        # discard hom-InDel
        for my $PosIdx ( sort {$a<=>$b} keys %$chrHref ){
            my $PhaMutOB_Aref = $chrHref->{$PosIdx};
            my $count = scalar @$PhaMutOB_Aref;
            for my $ArrayIdx (reverse 0 .. $count-1 ){
                if( $PhaMutOB_Aref->[$ArrayIdx]->is_hom ){
                    splice( @$PhaMutOB_Aref, $ArrayIdx, 1 );
                }
            }
            # become empty!
            if( scalar(@{$chrHref->{$PosIdx}}) == 0 ){
                delete $chrHref->{$PosIdx};
            }
        }
    }
    # output only in 1st ('splitBam') step
    close CLIDRT unless $V_Href->{skipS01Report};

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tfilter close-InDel alleles with hom-InDel alleles OK.\n" unless $V_Href->{skipS01Report};
}

#--- find PhasedMut(s) by given pos/region ---
sub PosToPhasedMut{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chr = $parm{chr};
    my $pos = $parm{pos};
    my $ext = $parm{ext} || 0;

    if(exists $V_Href->{PhasedMut}->{$chr}){
        return FindOverlapIdxRegion(
                 IdxItvHref => $V_Href->{PhasedMut}->{$chr},
                 regionAref => $ext >= 0 ? [$pos, $pos+$ext] : [$pos+$ext, $pos],
                 winSize => $V_Href->{phaMutPosWinSize},
                 queryMode => 'ob_pos'
               );
    }
    else{
        return [];
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
