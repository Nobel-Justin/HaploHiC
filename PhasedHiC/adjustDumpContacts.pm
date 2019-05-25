package HaploHiC::PhasedHiC::adjustDumpContacts;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ max sum /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use HaploHiC::LoadOn;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              adjustDumpContacts
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::adjustDumpContacts';
#----- version --------
$VERSION = "0.01";
$DATE = '2019-02-02';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        adjustDumpContacts
                        check_header
                        prepare_fh
                        compare_chrPair
                        adjustHaploHiCdump
                        loadJuicerDPtoAdjust
                        dealJuicerDPinfo
                        dealHaploHiCRemainDP
                        dealJuicerRemainDP
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} adjustDump <[Options]>

     Options:

       # adjust HaploHiC dump contacts based on juicer dump contacts #

       # Inputs and Outputs #
        -db_hph  [s]  dump contacts file from step-5 of func 'haplo_div'. <required>
        -db_jci  [s]  dump contacts file from func 'juicerDump'. <required>

        -h|help       Display this help info.

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
    $V_Href->{ $_->[0] } = $_->[1] for
        map {
            if( !exists $V_Href->{$_->[0]} ){
                ( $_ );
            }
            else{
                warn_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
            }
        }
        (
            # input/output
            ## dumpContact->{HaploHiC}
            ## dumpContact->{Juicer}
            ## dumpContact->{HaploHiC_A} # adjust
            ## dumpContact->{HaploHiC_U} # uniq
            ## dumpContact->{Juicer_U} # uniq
            [ dumpContact => {} ],
            [ dumpContactJuicer => undef ],

            # intermediate variants
            [ dumpContactHeader => {} ],
            [ dumpContactFH => {} ],
            [ uniqPairChr => {} ],
            [ pairChrContact => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['dumpContactHaploHiC'],
                                  ['dumpContactJuicer']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-db_hph:s" => \$V_Href->{dumpContact}->{HaploHiC},
        "-db_jci:s" => \$V_Href->{dumpContact}->{Juicer},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{dumpContact}->{HaploHiC})
             || !file_exist(filePath=>$V_Href->{dumpContact}->{Juicer})
            );
}

#--- adjust HaploHiC dump contacts ---
sub adjustDumpContacts{

    # basic
    &check_header;
    &prepare_fh;

    # compare chrPair
    &compare_chrPair;

    # adjust HaploHiC dump contact
    &adjustHaploHiCdump;
}

#--- check header of dump contacts ---
sub check_header{
    my %dumpAttr;
    for my $tool (qw/ HaploHiC Juicer /){
        open (DP, Try_GZ_Read($V_Href->{dumpContact}->{$tool})) || die "fail read $tool dumpContact: $!\n";
        while(<DP>){
            last if $. >= 10;
            # record header
            $V_Href->{dumpContactHeader}->{$tool} .= $_ if /^#/;
            # attributes
            ## dumpMode: BP, dumpBinSize: 1000000
            if(/dumpMode:\s*(\S+),\s*dumpBinSize:\s*(\S+)/){
                $dumpAttr{dumpMode}{$tool} = $1;
                $dumpAttr{dumpBinSize}{$tool} = $2;
            }
            if(/normMethod:\s*(\S+)/){
                $dumpAttr{normMethod}{$tool} = $1;
            }
            if(/enzyme:\s*(\S+)/){
                $dumpAttr{enzyme}{$tool} = $1;
            }
        }
        close DP;
    }

    # compare dump-attributes
    if(    $dumpAttr{dumpMode}{HaploHiC}     ne $dumpAttr{dumpMode}{Juicer}
        || $dumpAttr{dumpBinSize}{HaploHiC}  ne $dumpAttr{dumpBinSize}{Juicer}
        || $dumpAttr{normMethod}{HaploHiC}   ne $dumpAttr{normMethod}{Juicer}
        || (   $dumpAttr{dumpMode}{HaploHiC} eq 'FRAG'
            && $dumpAttr{enzyme}{HaploHiC}   ne $dumpAttr{enzyme}{Juicer}
           )
    ){
        warn_and_exit "<ERROR>\tDifferent attributes between dump contacts files.\n";
    }
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tdump contacts file header checking DONE\n";
}

#--- prepare file-handle ---
sub prepare_fh{
    # input: HaploHiC dump
    open ($V_Href->{dumpContactFH}->{HaploHiC}, Try_GZ_Read($V_Href->{dumpContact}->{HaploHiC})) || die "fail read HaploHiC dumpContact: $!\n";
    # input: juicer dump
    open ($V_Href->{dumpContactFH}->{Juicer},   Try_GZ_Read($V_Href->{dumpContact}->{Juicer}))   || die "fail read Juicer dumpContact: $!\n";

    # output: HaploHiC dump adjusted
    $V_Href->{dumpContact}->{HaploHiC_A} = $V_Href->{dumpContact}->{HaploHiC} . ".adjustWithJuicer.gz";
    open ($V_Href->{dumpContactFH}->{HaploHiC_A}, Try_GZ_Write($V_Href->{dumpContact}->{HaploHiC_A})) || die "fail write adjustWithJuicer: $!\n";
    print {$V_Href->{dumpContactFH}->{HaploHiC_A}} $V_Href->{dumpContactHeader}->{HaploHiC};
    # output: HaploHiC dump unique
    $V_Href->{dumpContact}->{HaploHiC_U} = $V_Href->{dumpContact}->{HaploHiC} . ".adjustWithJuicer.HaploHiC_unique.gz";
    open ($V_Href->{dumpContactFH}->{HaploHiC_U}, Try_GZ_Write($V_Href->{dumpContact}->{HaploHiC_U})) || die "fail write adjustWithJuicer.HaploHiC_unique: $!\n";
    print {$V_Href->{dumpContactFH}->{HaploHiC_U}} $V_Href->{dumpContactHeader}->{HaploHiC};
    # output: juicer dump unique
    $V_Href->{dumpContact}->{Juicer_U}   = $V_Href->{dumpContact}->{HaploHiC} . ".adjustWithJuicer.Juicer_unique.gz";
    open ($V_Href->{dumpContactFH}->{Juicer_U}, Try_GZ_Write($V_Href->{dumpContact}->{Juicer_U})) || die "fail write adjustWithJuicer.Juicer_unique: $!\n";
    print {$V_Href->{dumpContactFH}->{Juicer_U}} $V_Href->{dumpContactHeader}->{Juicer};
}

#--- compare chr pair between two dump contact ---
sub compare_chrPair{
    my %chrPair;
    # get chrPair of HaploHiC dump
    open (DPHL, Try_GZ_Read($V_Href->{dumpContact}->{HaploHiC})) || die "fail read HaploHiC dumpContact: $!\n";
    while(<DPHL>){
        next if /^#/;
        my ($hapChr_i, $hapChr_j) = (split)[0,3];
        my ($chr_i) = ($hapChr_i =~ /h\d+:(.+)/);
        my ($chr_j) = ($hapChr_j =~ /h\d+:(.+)/);
        $chrPair{HaploHiC}{"$chr_i,$chr_j"} = 1;
    }
    close DPHL;
    # get chrPair of juicer dump
    open (DPJC, Try_GZ_Read($V_Href->{dumpContact}->{Juicer})) || die "fail read Juicer dumpContact: $!\n";
    while(<DPJC>){
        next if /^#/;
        my ($chr_i, $chr_j) = (split)[0,3];
        $chrPair{Juicer}{"$chr_i,$chr_j"} = 1;
    }
    close DPJC;
    # compare
    $V_Href->{uniqPairChr}->{HaploHiC}->{$_} = 1 for grep !exists $chrPair{Juicer}{$_},   keys %{$chrPair{HaploHiC}};
    $V_Href->{uniqPairChr}->{Juicer}->{$_}   = 1 for grep !exists $chrPair{HaploHiC}{$_}, keys %{$chrPair{Juicer}};
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tget unique and shared chrPair DONE\n";
}

#--- adjust dump contacts ---
sub adjustHaploHiCdump{
    my $lastChrPair = '_UNDEF_';
    my $DPHLFH = $V_Href->{dumpContactFH}->{HaploHiC};
    while(<$DPHLFH>){
        next if /^#/;
        my ($hapChr_i, $chrBinIdx_i, $wgBinIdx_i, $hapChr_j, $chrBinIdx_j, $wgBinIdx_j, $contacts) = (split);
        my ($hap_i, $chr_i) = ($hapChr_i =~ /(h\d+):(.+)/);
        my ($hap_j, $chr_j) = ($hapChr_j =~ /(h\d+):(.+)/);
        my $thisChrPair = "$chr_i,$chr_j";
        # HaploHiC uniq chrPair?
        if(exists $V_Href->{uniqPairChr}->{HaploHiC}->{$thisChrPair}){
            print {$V_Href->{dumpContactFH}->{HaploHiC_U}} $_;
            next;
        }
        # to adjust shared chrPair
        if(    $thisChrPair ne $lastChrPair
            && $lastChrPair ne '_UNDEF_'
        ){
            &loadJuicerDPtoAdjust(chrPair => $lastChrPair);
            &dealHaploHiCRemainDP(chrPair => $lastChrPair);
            # inform
            stout_and_sterr "[INFO]\t".`date`
                                 ."\tadjust $lastChrPair dump contacts DONE\n";
        }
        # record
        $V_Href->{pairChrContact}->{HaploHiC}->{$thisChrPair}->{$chrBinIdx_i}->{$chrBinIdx_j}->{"$hap_i,$hap_j"} = [$wgBinIdx_i, $wgBinIdx_j, $contacts];
        $lastChrPair = $thisChrPair;
    }
    # last chrPair
    &loadJuicerDPtoAdjust(chrPair => $lastChrPair);
    &dealHaploHiCRemainDP(chrPair => $lastChrPair);
    # deal with remained Juicer dump contacts
    &dealJuicerRemainDP;
    # close file-handles
    close $V_Href->{dumpContactFH}->{$_} for keys %{$V_Href->{dumpContactFH}};
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tadjust all dump contacts DONE\n";
}

#--- load juicer dump contacts to adjust ---
sub loadJuicerDPtoAdjust{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chrPair = $parm{chrPair};

    my $find_bool = 0;
    # last recorded
    if(exists $V_Href->{pairChrContact}->{Juicer}){
        &dealJuicerDPinfo(chrPair => $chrPair, JCdpText => $V_Href->{pairChrContact}->{Juicer}, find_bool_Sf => \$find_bool);
        delete $V_Href->{pairChrContact}->{Juicer};
    }
    # go on read juicer dump contacts
    my $DPJCFH = $V_Href->{dumpContactFH}->{Juicer};
    while(<$DPJCFH>){
        next if /^#/;
        last if &dealJuicerDPinfo(chrPair => $chrPair, JCdpText => $_, find_bool_Sf => \$find_bool);
    }
}

#--- use juicer dump contacts to adjust ---
sub dealJuicerDPinfo{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chrPair = $parm{chrPair};
    my $JCdpText = $parm{JCdpText};
    my $find_bool_Sf = $parm{find_bool_Sf};

    my ($chr_i, $chrBinIdx_i, $wgBinIdx_i, $chr_j, $chrBinIdx_j, $wgBinIdx_j, $contacts) = (split /\s+/, $JCdpText);
    my $thisChrPair = "$chr_i,$chr_j";
    # juicer uniq chrPair?
    if(exists $V_Href->{uniqPairChr}->{Juicer}->{$thisChrPair}){
        print {$V_Href->{dumpContactFH}->{Juicer_U}} $JCdpText;
        return 0;
    }
    # match?
    if($chrPair ne $thisChrPair){
        if($$find_bool_Sf){
            $V_Href->{pairChrContact}->{Juicer} = $JCdpText;
            return 1;
        }
        else{
            warn_and_exit "<ERROR>\tdifferent order of chrPair in dump contacts files.\n";
        }
    }
    else{
        $$find_bool_Sf = 1;
    }
    # shared chrPair to adjust
    if(exists $V_Href->{pairChrContact}->{HaploHiC}->{$thisChrPair}->{$chrBinIdx_i}->{$chrBinIdx_j}){
        my $chrBinPairHf = $V_Href->{pairChrContact}->{HaploHiC}->{$thisChrPair}->{$chrBinIdx_i}->{$chrBinIdx_j};
        my @hapComb = sort {$chrBinPairHf->{$a}->[2] <=> $chrBinPairHf->{$b}->[2]} keys %$chrBinPairHf;
        my $contacts_sum = sum( map {($chrBinPairHf->{$_}->[2])} @hapComb );
        my $contacts_adjSum = 0;
        for my $hapComb (@hapComb){
            my $contact_adj = max(1, int($contacts * $chrBinPairHf->{$hapComb}->[2] / $contacts_sum));
            push @{$chrBinPairHf->{$hapComb}}, $contact_adj;
            $contacts_adjSum += $contact_adj;
        }
        # not reach juicer contacts?
        if($contacts_adjSum < $contacts){
            $chrBinPairHf->{$hapComb[-1]}->[3] += ($contacts - $contacts_adjSum);
        }
        # output
        for my $hapComb (sort keys %$chrBinPairHf){
            my ($hap_i, $hap_j) = split /,/, $hapComb;
            my ($wgBinIdx_i, $wgBinIdx_j, $contacts_orig, $contacts_adj) = @{$chrBinPairHf->{$hapComb}};
            print {$V_Href->{dumpContactFH}->{HaploHiC_A}} join("\t", "$hap_i:$chr_i", $chrBinIdx_i, $wgBinIdx_i,
                                                                      "$hap_j:$chr_j", $chrBinIdx_j, $wgBinIdx_j,
                                                                      $contacts_adj)."\n";
        }
        # delete HaploHiC chrBinPairHf
        delete $V_Href->{pairChrContact}->{HaploHiC}->{$thisChrPair}->{$chrBinIdx_i}->{$chrBinIdx_j};
    }
    else{ # uniq bin-pair
        print {$V_Href->{dumpContactFH}->{Juicer_U}} $JCdpText;
    }
    # go on
    return 0;
}

#--- deal with remained HaploHiC dump contacts ---
sub dealHaploHiCRemainDP{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chrPair = $parm{chrPair};

    my ($chr_i, $chr_j) = (split /,/, $chrPair);
    my $chrPairHf = $V_Href->{pairChrContact}->{HaploHiC}->{$chrPair};
    for my $chrBinIdx_i (sort {$a<=>$b} keys %$chrPairHf){
        for my $chrBinIdx_j (sort {$a<=>$b} keys %{$chrPairHf->{$chrBinIdx_i}}){
            my $chrBinPairHf = $chrPairHf->{$chrBinIdx_i}->{$chrBinIdx_j};
            for my $hapComb (sort keys %$chrBinPairHf){
                my ($hap_i, $hap_j) = split /,/, $hapComb;
                my ($wgBinIdx_i, $wgBinIdx_j, $contacts_orig) = @{$chrBinPairHf->{$hapComb}};
                print {$V_Href->{dumpContactFH}->{HaploHiC_U}} join("\t", "$hap_i:$chr_i", $chrBinIdx_i, $wgBinIdx_i,
                                                                          "$hap_j:$chr_j", $chrBinIdx_j, $wgBinIdx_j,
                                                                          $contacts_orig)."\n";
            }
        }
    }
    # sweep
    delete $V_Href->{pairChrContact}->{HaploHiC}->{$chrPair};
}

#--- deal with remained Juicer dump contacts ---
sub dealJuicerRemainDP{
    # remained juicer dump contacts
    if(exists $V_Href->{pairChrContact}->{Juicer}){
        print {$V_Href->{dumpContactFH}->{Juicer_U}} $V_Href->{pairChrContact}->{Juicer};
        delete $V_Href->{pairChrContact}->{Juicer};
    }
    my $DPJCFH = $V_Href->{dumpContactFH}->{Juicer};
    while(<$DPJCFH>){
        next if /^#/;
        print {$V_Href->{dumpContactFH}->{Juicer_U}} $_;
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
