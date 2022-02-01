package HaploHiC::PhasedHiC::dumpContacts;

use strict;
use warnings;
use Data::Dumper;
use File::Spec::Functions qw/ catfile /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Array qw/ binarySearch /;
use BioFuse::Util::GZfile qw/ Try_GZ_Write Try_GZ_Read /;
use BioFuse::Util::Index qw/ Pos2Idx /;
use HaploHiC::LoadOn;
use HaploHiC::PhasedHiC::phasedPEtoContact qw/ load_phasedPE_contacts phasePE_contacts_to_count /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              dump_contacts
              load_enzyme_site_list
              get_contacts_idx
              get_dumpHeader
              write_dumpBinLog
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::dumpContacts';
#----- version --------
$VERSION = "0.10";
$DATE = '2019-03-27';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        dump_contacts
                        prepare
                        load_enzyme_site_list
                        get_contacts_idx
                        get_dumpHeader
                        write_dumpBinLog
                        contacts_output
                     /;

#--- dump contacts ---
sub dump_contacts{
    # start from step after current step
    return if $V_Href->{stepToStart} > 5;

    for my $dumpOptAf (@{$V_Href->{dump}}){
        # dump options
        # folder, prefix, header
        # real dumpBinSize
        # enzyme_site (if FRAG)
        &prepare(dumpOptAf => $dumpOptAf);
        # group merged bam files belonging to each tag (h[x]Intra and hInter)
        my %tagToMergeBam;
        for my $pairBamHref ( @{$V_Href->{PairBamFiles}} ){
            push @{$tagToMergeBam{$_}}, $pairBamHref->{mergeBam}->{$_} for keys %{$pairBamHref->{mergeBam}};
        }
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
                my $mark = $mergeBam->tag;
                # read phased bam
                my @lastChrPair = ('__NA__', '__NA__', $mark); # takes $mark by the way
                my @subrtOpt = (subrtRef => \&load_phasedPE_contacts, subrtParmAref => [idxFunc => \&get_contacts_idx, lastChrPairAf => \@lastChrPair]);
                $mergeBam->smartBam_PEread(samtools => $V_Href->{samtools}, readsType => 'HiC', deal_peOB_pool => 1, @subrtOpt);
                # contacts to count for last chr-pair, release memory
                ## here, do de-dup phased reads (in single run) (optional)
                phasePE_contacts_to_count(tag => "$mark " . join(',', @lastChrPair[0,1])) if $lastChrPair[0] ne '__NA__';
                # inform
                stout_and_sterr "[INFO]\t".`date`
                                     ."\tload contacts from $mark bam OK.\n";
            }
            # output contacts
            &contacts_output(tag => $tag);
            # inform
            stout_and_sterr "[INFO]\t".`date`
                                 ."\tdump contacts (m:$V_Href->{dumpMode};b:$V_Href->{dumpBinSize}) from $tag phased Hi-C PE-reads OK.\n";
        }
    }

    # stop at current step
    exit(0) if $V_Href->{stepToStop} == 5;
}

#--- prepare work ---
sub prepare{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $dumpOptAf = $parm{dumpOptAf};

    # dump options
    my ($dumpMode, $dumpBinSize, $dumpHapComb) = @$dumpOptAf;
    $V_Href->{dumpMode} = $dumpMode;
    $V_Href->{dumpBinSize} = $dumpBinSize;
    $V_Href->{dumpHapComb} = $dumpHapComb if $dumpHapComb ne 'all';
    # reset global variant
    $V_Href->{dumpBinLog} = undef;
    # folder
    my $dumpOutDir = catfile($V_Href->{outdir}, $V_Href->{dumpSubDir});
    (-d $dumpOutDir || `mkdir -p $dumpOutDir`);
    # file-name-prefix (full-path)
    $V_Href->{dumpFilePrefix} = catfile($dumpOutDir, "dumpContacts.$V_Href->{dumpMode}.$V_Href->{dumpBinSize}");
    # use real number of bin_size
    $V_Href->{dumpBinSize} = $V_Href->{dump_allowBinSize}->{$V_Href->{dumpMode}}->{uc($V_Href->{dumpBinSize})};
    # header of dump output file
    &get_dumpHeader;
    # load enzyme sites list in FRAG mode
    &load_enzyme_site_list if($V_Href->{dumpMode} eq 'FRAG');
}

#--- load enzyme site position list ---
sub load_enzyme_site_list{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $keep_all = $parm{keep_all} || 0;

    # read enzyme site file from juicer database
    open (EMLIST, Try_GZ_Read($V_Href->{enzyme_site})) || die "fail read enzyme site list: $!\n";
    while(<EMLIST>){
        next if(/^#/);
        my @info = split;
        my $chr = shift @info;
        # skip chr out of wanted list loaded if ever, if allowed
        if(    !$keep_all
            &&  scalar(keys %{$V_Href->{ChrThings}})
            && !exists($V_Href->{ChrThings}->{$chr})
        ){
            next;
        }
        # record array for bisect-serach
        $V_Href->{chr2enzymePos}->{$chr} = \@info;
    }
    close EMLIST;

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tread enzyme site DONE.\n" unless $V_Href->{skipS01Report};
}

#--- get contact idx in given mode ---
sub get_contacts_idx{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
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

#--- prepare header of dump file ---
sub get_dumpHeader{
    $V_Href->{dumpHeader}  = "##dumpMode: $V_Href->{dumpMode}, dumpBinSize: $V_Href->{dumpBinSize}\n";
    $V_Href->{dumpHeader} .= "##normMethod: $V_Href->{norm_method}\n";
    $V_Href->{dumpHeader} .= "##enzyme: $V_Href->{enzyme_type}\n"    if $V_Href->{dumpMode} eq 'FRAG';
    $V_Href->{dumpHeader} .= "##haploCount: $V_Href->{haploCount}\n" if $V_Href->{haploCount};
    $V_Href->{dumpHeader} .= "##dumpHapComb: $V_Href->{dumpHapComb}\n" if $V_Href->{dumpHapComb};
    $V_Href->{dumpHeader} .= "##".`date`;
}

#--- write chr-bin-range log of each refseg ---
sub write_dumpBinLog{
    open (CHRBR, Try_GZ_Write($V_Href->{dumpBinLog})) || die "cannot write chr-bin-range file: $!\n";
    print CHRBR $V_Href->{dumpHeader};
    print CHRBR join("\t", '#chr', 'chr_len', 'chrBinIdx_st', 'chrBinIdx_ed', 'wgBinIdx_st', "wgBinIdx_ed\n");
    my $preChrBinSum = 0;
    my @chr_Href = sort { $a->{turn} <=> $b->{turn} } values %{$V_Href->{ChrThings}};
    for my $chr_Href (@chr_Href){
        my $chr = $chr_Href->{chr};
        my $chrlen = $chr_Href->{len};
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
        $chr_Href->{preChrBinSum} = $preChrBinSum;
        # update
        $preChrBinSum += $thisChrBinCount;
    }
    close CHRBR;
}

#--- contacts output and bin log ---
sub contacts_output{
    # options
    # shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $tag = $parm{tag};

    # get chr-idx interval
    ## chr-bin-range log, only first time
    if(!defined $V_Href->{dumpBinLog}){
        $V_Href->{dumpBinLog} = "$V_Href->{dumpFilePrefix}.chrBinRange.txt";
        &write_dumpBinLog;
    }

    # output contacts
    $V_Href->{dumpOutput} = "$V_Href->{dumpFilePrefix}.$tag.contacts.txt.gz";
    open (CNTC, Try_GZ_Write($V_Href->{dumpOutput})) || die "cannot write contacts file: $!\n";
    print CNTC $V_Href->{dumpHeader};
    print CNTC join("\t", '#hap_i:chr_i', 'chrBinIdx_i', "wgBinIdx_i", 'hap_j:chr_j', 'chrBinIdx_j', "wgBinIdx_j", "contacts\n");
    my @chr_Href = sort { $a->{turn} <=> $b->{turn} } values %{$V_Href->{ChrThings}};
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
                    my $hapComb_Href = $pIdxComb_Href->{$chr_pIdx_i}->{$chr_pIdx_j};
                    # sum inter-haplotype same wg-idx contacts, e.g., Win2(h1,h2) + Win2(h2,h1)
                    my %skip_hapComb;
                    for my $hapComb (sort keys %$hapComb_Href){
                        my ($hap_i, $hap_j) = split /,/, $hapComb;
                        my $contact_count = $hapComb_Href->{$hapComb};
                        # add pair hapComb if it is inter-haplotype same wg-idx
                        if(    $wg_pIdx_i == $wg_pIdx_j
                            && $hap_i ne $hap_j
                        ){
                            next if exists $skip_hapComb{$hapComb};
                            my $pair_hapComb = "$hap_j,$hap_i";
                            $contact_count += ($hapComb_Href->{$pair_hapComb} || 0);
                            $skip_hapComb{$pair_hapComb} = 1;
                        }
                        # output
                        print CNTC join("\t", "$hap_i:$chr_i",
                                              $chr_pIdx_i,
                                              $wg_pIdx_i,
                                              "$hap_j:$chr_j",
                                              $chr_pIdx_j,
                                              $wg_pIdx_j,
                                              $contact_count)."\n";
                    }
                }
            }
        }
    }
    close CNTC;
}

#--- 
1; ## tell the perl script the successful access of this module.
