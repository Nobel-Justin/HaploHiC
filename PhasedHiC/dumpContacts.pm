package HaploHiC::PhasedHiC::dumpContacts;

use strict;
use warnings;
use Data::Dumper;
use File::Spec::Functions qw/ catfile /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Array qw/ binarySearch /;
use BioFuse::Util::GZfile qw/ Try_GZ_Write /;
use BioFuse::Util::Index qw/ Pos2Idx /;
use HaploHiC::LoadOn;
use HaploHiC::Extensions::FRAG2Gene qw/ load_enzyme_site_list /;
use HaploHiC::PhasedHiC::phasedPEtoContact qw/ load_phasedPE_contacts phasePE_contacts_to_count /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              dump_contacts
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::dumpContacts';
#----- version --------
$VERSION = "0.03";
$DATE = '2019-01-27';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        dump_contacts
                        get_contacts_idx
                        contacts_output
                     /;

#--- dump contacts ---
sub dump_contacts{
    # start from step after current step
    return if $V_Href->{stepToStart} > 5;

    # group merged bam files belonging to each tag (h[x]Intra and hInter)
    my %tagToMergeBam;
    for my $pairBamHref ( @{$V_Href->{PairBamFiles}} ){
        push @{$tagToMergeBam{$_}}, $pairBamHref->{mergeBam}->{$_} for keys %{$pairBamHref->{mergeBam}};
    }
    # use real number of bin_size
    $V_Href->{dumpBinSize} = $V_Href->{dump_allowBinSize}->{$V_Href->{dumpMode}}->{uc($V_Href->{dumpBinSize})};
    # load enzyme sites list in FRAG mode
    load_enzyme_site_list if($V_Href->{dumpMode} eq 'FRAG');

    # all OR specific tags, and check empty
    my @tag = sort keys %tagToMergeBam;
       @tag = grep /$V_Href->{dumpHapComb}/, @tag if defined $V_Href->{dumpHapComb};
    if(@tag == 0){
        warn_and_exit "<ERROR>\tNo Haplo-Comb selected to dump contacts by the option '-dphpcmb' ($V_Href->{dumpHapComb}).\n";
    }
    # deal with bam files of each tag
    for my $tag (@tag){
        $V_Href->{phasePEcontact} = {}; # reset contact container
        for my $mergeBam (@{$tagToMergeBam{$tag}}){
            # read phased bam
            my @subrtOpt = (subrtRef => \&load_phasedPE_contacts, subrtParmAref => [idxFunc => \&get_contacts_idx]);
            $mergeBam->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC', @subrtOpt);
            # inform
            my $mark = $mergeBam->get_tag;
            stout_and_sterr "[INFO]\t".`date`
                                 ."\tload contacts from $mark bam OK.\n";
        }
        # contacts to count
        ## here, do de-dup phased reads (optional)
        phasePE_contacts_to_count(tag => $tag);
        # output contacts
        &contacts_output(tag => $tag);
        # inform
        stout_and_sterr "[INFO]\t".`date`
                             ."\tdump contacts (m:$V_Href->{dumpMode};b:$V_Href->{dumpBinSize}) from $tag phased Hi-C PE-reads OK.\n";
    }

    # stop at current step
    exit(0) if $V_Href->{stepToStop} == 5;
}

#--- get contact idx in given mode ---
sub get_contacts_idx{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $get_remainder = $parm{get_remainder} || 0;

    if($V_Href->{dumpMode} eq 'BP'){
        return   $get_remainder
               ? $parm{pos} % $V_Href->{dumpBinSize}
               : Pos2Idx(pos => $parm{pos}, winSize => $V_Href->{dumpBinSize});
    }
    elsif($V_Href->{dumpMode} eq 'FRAG'){
        unless(exists $V_Href->{chr2enzymePos}->{$parm{chr}}){
            warn_and_exit "<ERROR>\tCannot find enzyme site of chromosome $parm{chr}:$parm{pos}\n";
        }
        my $FragIdx = binarySearch(query => $parm{pos}, array => $V_Href->{chr2enzymePos}->{$parm{chr}});
        return   $get_remainder
               ? $FragIdx % $V_Href->{dumpBinSize}
               : int($FragIdx / $V_Href->{dumpBinSize});
    }
}

#--- contacts output and bin log ---
sub contacts_output{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $tag = $parm{tag};

    my $dumpOutDir = catfile($V_Href->{outdir}, $V_Href->{dumpSubDir});
    (-d $dumpOutDir || `mkdir -p $dumpOutDir`);
    my $output_prefix = "dumpContacts.$V_Href->{dumpMode}.$V_Href->{dumpBinSize}";
    my $output_header = "##dumpMode: $V_Href->{dumpMode}, dumpBinSize: $V_Href->{dumpBinSize}\n"
                       .($V_Href->{dumpMode} eq 'FRAG' ? "##enzyme: $V_Href->{enzyme_type}\n" : '')
                       ."##haploCount: $V_Href->{haploCount}\n"
                       ."##".`date`;

    my @chr_Href = sort { $a->{turn} <=> $b->{turn} } values %{$V_Href->{ChrThings}};
    # get chr-idx interval
    ## chr-bin-range log
    if(!defined $V_Href->{dumpBinLog}){
        $V_Href->{dumpBinLog} = catfile($dumpOutDir, "$output_prefix.chrBinRange.txt");
        open (CHRBR, Try_GZ_Write($V_Href->{dumpBinLog})) || die "cannot write chr-bin-range file: $!\n";
        print CHRBR $output_header;
        print CHRBR join("\t", '#chr', 'chr_len', 'chrBinIdx_st', 'chrBinIdx_ed', 'wgBinIdx_st', "wgBinIdx_ed\n");
        my $preChrBinSum = 0;
        for my $chrHref (@chr_Href){
            my $chr = $chrHref->{chr};
            my $chrlen = $chrHref->{len};
            my $thisChrFirstIdx =  &get_contacts_idx(chr => $chr, pos => 1);
            my $thisChrLastIdx  =  &get_contacts_idx(chr => $chr, pos => $chrlen);
            my $thisChrBinCount =  &get_contacts_idx(chr => $chr, pos => $chrlen, get_remainder => 1)
                                 ? $thisChrLastIdx + 1
                                 : $thisChrLastIdx;
            print CHRBR join("\t", $chr,
                                   $chrlen,
                                   $thisChrFirstIdx,
                                   $thisChrLastIdx,
                                   $preChrBinSum+$thisChrFirstIdx,
                                   $preChrBinSum+$thisChrLastIdx) . "\n";
            # record preChrBinSum of this chr
            $chrHref->{preChrBinSum} = $preChrBinSum;
            # update
            $preChrBinSum += $thisChrBinCount;
        }
        close CHRBR;
    }

    # output contacts
    $V_Href->{dumpOutput} = catfile($dumpOutDir, "$output_prefix.$tag.contacts.txt.gz");
    open (CNTC, Try_GZ_Write($V_Href->{dumpOutput})) || die "cannot write contacts file: $!\n";
    print CNTC $output_header;
    print CNTC join("\t", '#hap_i:chr_i', 'chrBinIdx_i', "wgBinIdx_i", 'hap_j:chr_j', 'chrBinIdx_j', "wgBinIdx_j", "contacts\n");
    for my $chr_i_Href (@chr_Href){
        my $chr_i = $chr_i_Href->{chr};
        next unless exists $V_Href->{phasePEcontact}->{$chr_i};
        my $chr_i_preChrBinSum = $chr_i_Href->{preChrBinSum};
        for my $chr_j_Href (@chr_Href){
            my $chr_j = $chr_j_Href->{chr};
            next unless exists $V_Href->{phasePEcontact}->{$chr_i}->{$chr_j};
            my $chr_j_preChrBinSum = $chr_j_Href->{preChrBinSum};
            my $pIdxComb_Href = $V_Href->{phasePEcontact}->{$chr_i}->{$chr_j};
            for my $chr_pIdx_i (sort {$a<=>$b} keys %$pIdxComb_Href){
                my $wg_pIdx_i = $chr_i_preChrBinSum + $chr_pIdx_i;
                for my $chr_pIdx_j (sort {$a<=>$b} keys %{$pIdxComb_Href->{$chr_pIdx_i}}){
                    my $wg_pIdx_j = $chr_j_preChrBinSum + $chr_pIdx_j;
                    for my $hapComb (sort keys %{$pIdxComb_Href->{$chr_pIdx_i}->{$chr_pIdx_j}}){
                        my ($hap_i, $hap_j) = split /,/, $hapComb;
                        print CNTC join("\t", "$hap_i:$chr_i",
                                              $chr_pIdx_i,
                                              $wg_pIdx_i,
                                              "$hap_j:$chr_j",
                                              $chr_pIdx_j,
                                              $wg_pIdx_j,
                                              $pIdxComb_Href->{$chr_pIdx_i}->{$chr_pIdx_j}->{$hapComb})."\n";
                    }
                }
            }
        }
    }
    close CNTC;
}

#--- 
1; ## tell the perl script the successful access of this module.
