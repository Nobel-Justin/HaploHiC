package HaploHiC::PhasedHiC::plotPhasedLocReg;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Spec::Functions qw/ catfile /;
use List::Util qw/ max sum any /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Visual::Objects::Axis;
use BioFuse::Visual::Objects::BiAxis;
use BioFuse::Visual::Objects::Histogram;
use BioFuse::Visual::Objects::GradColor;
use BioFuse::Visual::SVG_Util::SVGWork qw/ initialize_SVG_obj output_SVG_file /;
use HaploHiC::LoadOn;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              plotPhasedLocReg
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::PhasedHiC::plotPhasedLocReg';
#----- version --------
$VERSION = "0.04";
$DATE = '2019-04-06';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        plotPhasedLocReg
                        prepare
                        loadLocRegInfo
                        filter_hapDivTag
                        filter_chrPair
                        filter_hapComb
                        plotLocRegInfo
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} plotLocReg <[Options]>

     Options:

       # plot phased local region information #

       # Inputs and Outputs #
        -dir     [s]  folder stores results of func 'haplo_div'. <required>
        -svg     [s]  svg figure output. <required>

       # Options #
        -max_lg  [s]  maximum size of local region. [2E7]
        -tmr_lg  [f]  ratio for unilateral trimming from largest local region size. [0.01]
        -tmr_pc  [f]  ratio for unilateral trimming from phased link counts. [0.01]
        -amt_fc  [s]  manually set the top of amount. [0]
        -rgb     [s]  rgb of heatmap color. ['255,0,0']
        -add_r   [f]  uniform addition ratio in phased local region. [0]<=0.1

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
            [ source_dir => undef ],
            [ svg_file => undef ],

            # options
            [ maxLocRegSize => 2E7 ],
            [ TrimRato => {UT=>0.01, LC=>0.01} ],
            [ color_rgb => '255,0,0' ],
            [ keep_hapDivTag_regex => [] ],
            [ keep_chr_regex => [] ],
            [ keep_hap_regex => [] ],
            [ keep_hInter => 0 ],
            [ forceAmountTop => 0 ],

            # intermediate variants
            [ statFiles => [] ],
            [ localRegionUnit => undef ],
            [ uniformAddRatio => 0 ],
            [ locRegInfo => {} ],
            [ maxLocRegUnitTime => 0 ],
            [ maxPhasedLinkC => 0 ],
            [ maxAmount => 0 ],
            [ axisLabel_x => undef ],
            [ axisLabel_y => undef ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['source_dir'],
                                  ['svg_file']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-dir:s" => \$V_Href->{source_dir},
        "-svg:s" => \$V_Href->{svg_file},
        # option
        "-max_lg:s" => \$V_Href->{maxLocRegSize},
        "-tmr_lg:f" => \$V_Href->{TrimRato}->{UT},
        "-tmr_pc:f" => \$V_Href->{TrimRato}->{LC},
        "-amt_fc:s" => \$V_Href->{forceAmountTop},
        "-rgb:s"    => \$V_Href->{color_rgb},
        "-trgx:s"   => \@{$V_Href->{keep_hapDivTag_regex}}, # hidden option
        "-crgx:s"   => \@{$V_Href->{keep_chr_regex}}, # hidden option
        "-hrgx:s"   => \@{$V_Href->{keep_hap_regex}}, # hidden option
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !defined $V_Href->{source_dir}   || !-d $V_Href->{source_dir}
             || !defined $V_Href->{svg_file}
             || $V_Href->{TrimRato}->{UT} >= 1
             || $V_Href->{TrimRato}->{LC} >= 1
             || $V_Href->{color_rgb} !~ /^(\d+),(\d+),(\d+)$/
             || $V_Href->{forceAmountTop} < 0
            );
}

#--- plot phased local region information ---
sub plotPhasedLocReg{

    # basic
    &prepare;

    # load local region information
    &loadLocRegInfo;

    # plot
    &plotLocRegInfo;
}

#--- prepare ---
sub prepare{
    # gather stat files
    @{$V_Href->{statFiles}} = glob catfile($V_Href->{source_dir}, '*.statOfPhasedLocReg.gz');
    if(@{$V_Href->{statFiles}} == 0){
        warn_and_exit "<ERROR>\tcannot find any statOfPhasedLocReg.gz file.\n";
    }
    # get basic info from header
    for my $key (qw/ localRegionUnit uniformAddRatio /){
        chomp(my $header = `zcat $V_Href->{statFiles}->[0] | grep -i '^##$key' | head -1`);
        if($header =~ /$key:\s+(\S+)/i){
            $V_Href->{$key} = $1;
        }
        else{
            warn_and_exit "<ERROR>\tcannot get $key from header of statOfPhasedLocReg.gz file.\n";
        }
    }
    # color
    if($V_Href->{color_rgb} =~ /^(\d+),(\d+),(\d+)$/){
        $V_Href->{color_rgb} = [$1, $2, $3];
    }
    # hap-div tag to keep
    if(scalar(@{$V_Href->{keep_hapDivTag_regex}}) == 0){
        delete $V_Href->{keep_hapDivTag_regex};
    }
    # chr-pair to keep
    if(scalar(@{$V_Href->{keep_chr_regex}}) == 0){
        delete $V_Href->{keep_chr_regex};
    }
    # hap-comb to keep
    if(scalar(@{$V_Href->{keep_hap_regex}}) != 0){
        if(any {$_ eq 'hInter'} @{$V_Href->{keep_hap_regex}}){
            $V_Href->{keep_hInter} = 1;
        }
    }
    else{
        delete $V_Href->{keep_hap_regex};
    }
    # label of axis
    my @labelTag_x = ("unit=$V_Href->{localRegionUnit}");
    for my $keepOpt (qw/ keep_hapDivTag_regex keep_chr_regex keep_hap_regex /){
        next unless exists $V_Href->{$keepOpt};
        push @labelTag_x, join('+',map {"'$_'"} @{$V_Href->{$keepOpt}});
    }
    $V_Href->{axisLabel_x} = "phased local region size (".join('; ',@labelTag_x).")";
    $V_Href->{axisLabel_y} = "phased contacts count (log10)";
}

#--- load local region information ---
sub loadLocRegInfo{
    my %CountForTrim;
    for my $statFile (@{$V_Href->{statFiles}}){
        open (STAT, Try_GZ_Read($statFile)) || die "cannot open statFile: $!\n";
        while(<STAT>){
            next if /^#/;
            my ($hapDivTag, $chrPair, $locRegSize, $phasedLinkC, $hapCombCount, $Amount) = (split)[0,1,2,3,4,5];
            # local region size
            next if $locRegSize > $V_Href->{maxLocRegSize};
            # hap-div tag
            next if exists $V_Href->{keep_hapDivTag_regex} && &filter_hapDivTag(hapDivTag => $hapDivTag);
            # chr-pair
            next if exists $V_Href->{keep_chr_regex} && &filter_chrPair(chrPair => $chrPair);
            # add_ratio reduce
            if($V_Href->{uniformAddRatio} != 0){
                my @hapCombC = split /;/, $hapCombCount;
                # set as before adding
                $phasedLinkC /= (1 + $V_Href->{uniformAddRatio} * scalar(@hapCombC));
                my $addC = $phasedLinkC * $V_Href->{uniformAddRatio};
                $hapCombCount = '';
                for my $hapCombC (@hapCombC){
                    if($hapCombC =~ /(h\d+,h\d+):([\d\.]+)/){
                        my $origC = $2 - $addC;
                        $hapCombCount .= "$1:$origC;" if $origC != 0;
                    }
                }
            }
            # hap-comb
            if(exists $V_Href->{keep_hap_regex}){
                my $matchHapBool = 0;
                for my $hapCombC (split /;/, $hapCombCount){
                    if($hapCombC =~ /(h\d+,h\d+):([\d\.]+)/){
                        $matchHapBool ||= &filter_hapComb(hapComb => $1);
                    }
                    last if $matchHapBool;
                }
                next unless $matchHapBool;
            }
            # record
            my $unitTime = sprintf "%.2f", $locRegSize / $V_Href->{localRegionUnit};
            $V_Href->{locRegInfo}->{$unitTime}->{$phasedLinkC} += $Amount;
            $CountForTrim{UT}{$unitTime} += $Amount;
            $CountForTrim{LC}{$phasedLinkC} += $Amount;
        }
        close STAT;
    }
    # check
    if(scalar(keys %CountForTrim) == 0){
        warn_and_exit "<ERROR>\tno data loaded.\n";
    }

    # which keys to trim
    my %KeyToTrim;
    for my $type (keys %CountForTrim){
        next if $V_Href->{TrimRato}->{$type} == 0;
        my $AllCount = sum( values %{$CountForTrim{$type}} );
        my $TrimCount = $AllCount * $V_Href->{TrimRato}->{$type};
        my $SumCount = 0;
        for my $key (sort {$b<=>$a} keys %{$CountForTrim{$type}}){
            if($SumCount + $CountForTrim{$type}{$key} <= $TrimCount){
                $SumCount += $CountForTrim{$type}{$key};
                $KeyToTrim{$type}{$key} = 1;
            }
            else{
                last;
            }
        }
    }
    # trim
    for my $unitTime (keys %{$V_Href->{locRegInfo}}){
        if(exists $KeyToTrim{UT}{$unitTime}){
            delete $V_Href->{locRegInfo}->{$unitTime};
            next;
        }
        for my $phasedLinkC (keys %{$V_Href->{locRegInfo}->{$unitTime}}){
            delete $V_Href->{locRegInfo}->{$unitTime}->{$phasedLinkC} if exists $KeyToTrim{LC}{$phasedLinkC};
        }
        if(scalar(keys %{$V_Href->{locRegInfo}->{$unitTime}}) == 0){
            delete $V_Href->{locRegInfo}->{$unitTime};
        }
    }

    # get maximum
    $V_Href->{maxLocRegUnitTime} = max(keys %{$V_Href->{locRegInfo}});
    for my $Hf (values %{$V_Href->{locRegInfo}}){
        for my $phasedLinkC (keys %$Hf){
            $V_Href->{maxPhasedLinkC} = max($phasedLinkC, $V_Href->{maxPhasedLinkC});
            $V_Href->{maxAmount} = max($Hf->{$phasedLinkC}, $V_Href->{maxAmount});
        }
    }

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tload local region data OK.\n";
}

#--- filter hap-div tag ---
sub filter_hapDivTag{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    if(any {$parm{hapDivTag}=~/\b$_\b/} @{$V_Href->{keep_hapDivTag_regex}}){
        return 0;
    }
    else{
        return 1;
    }
}

#--- filter chr-pair ---
sub filter_chrPair{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    if(any {$parm{chrPair}=~/\b$_\b/} @{$V_Href->{keep_chr_regex}}){
        return 0;
    }
    else{
        return 1;
    }
}

#--- filter hap-comb ---
sub filter_hapComb{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $hapComb = $parm{hapComb};

    # hap-Inter general
    if($V_Href->{keep_hInter}){
        if($hapComb =~ /(h\d+),(h\d+)/){
           return 0 if $1 ne $2;
        }
    }
    # specific hap-Comb
    if(any {$hapComb=~/\b$_\b/} @{$V_Href->{keep_hap_regex}}){
        return 0;
    }
    else{
        return 1;
    }
}

#--- plot ---
sub plotLocRegInfo{
    my $x_axisLen = 500;
    my $y_axisLen = 300;
    my $zero_x = 50;
    my $zero_y = $y_axisLen * 1.2;
    # X
    my $x_axis = BioFuse::Visual::Objects::Axis->new(origP_X=>$zero_x+2, origP_Y=>$zero_y+2, axisLen=>$x_axisLen, headAng=>90);
    $x_axis->set_tic(len=>3, width=>1, clockwise=>1);
    $x_axis->add_resol(minValue=>0, maxValue=>$V_Href->{maxLocRegUnitTime}+1);
    $x_axis->set_stub(fontsize=>10);
    $x_axis->set_label(text=>$V_Href->{axisLabel_x});
    # Y
    my $y_axis = BioFuse::Visual::Objects::Axis->new(origP_X=>$zero_x, origP_Y=>$zero_y, axisLen=>$y_axisLen, headAng=>0 );
    $y_axis->set_tic(len=>3, width=>1, clockwise=>0);
    $y_axis->add_resol(minValue=>0, maxValue=>(sprintf "%.2f", log($V_Href->{maxPhasedLinkC})));
    $y_axis->set_stub(verticalToBL=>1, fontsize=>10);
    $y_axis->set_label(text=>$V_Href->{axisLabel_y});
    # X-Y
    my $biAxis = BioFuse::Visual::Objects::BiAxis->new(axis_1=>$x_axis, axis_2=>$y_axis);
    # histogram
    my $x_resol = $x_axis->get_resol->[0]->{resol};
    my $histogram = BioFuse::Visual::Objects::Histogram->new(bi_axis => $biAxis, isStack=>1, pillarWid=>1/$x_resol);
    $histogram->set_attr(strokeCol=>['none'], strokeWid=>[0]);
    # GradColor
    my ($gc_x, $gc_y) = $x_axis->extendCoord(coordAf => [$x_axis->valueToSVGcoord(value=>$x_axis->get_stub->[-1])], relateAng=>90, distance=>20);
    my $gc_edValue = $V_Href->{forceAmountTop} || $V_Href->{maxAmount};
    my $amountGC = BioFuse::Visual::Objects::GradColor->new(stValue=>0, edValue=>$gc_edValue, anchorX=>$gc_x, anchorY=>$gc_y);
    $amountGC->set_axis(label=>'amout', clockwise=>1, stubText=>($V_Href->{forceAmountTop}?'+':''));
    # data
    for my $LocRegUnitTime (keys %{$V_Href->{locRegInfo}}){
        my @phasedLinkC = sort {$a<=>$b} keys %{$V_Href->{locRegInfo}->{$LocRegUnitTime}};
        my %fill;
        my @phasedLinkC_log;
        my $phasedLinkC_last = 0;
        my $i = 0;
        for my $phasedLinkC (@phasedLinkC){
            if($phasedLinkC_last && $phasedLinkC_last+1 != $phasedLinkC){
                push @phasedLinkC_log, log($phasedLinkC-1)-log($phasedLinkC_last+1);
                $fill{$i} = 'rgb('.join(',',@{$amountGC->get_stRGB}).')';
                $i++;
            }
            push @phasedLinkC_log, log($phasedLinkC) - log($phasedLinkC-1);
            my $Amount = $V_Href->{locRegInfo}->{$LocRegUnitTime}->{$phasedLinkC};
            $fill{$i} = 'rgb('.join(',',@{$amountGC->valueToRGB(value=>$Amount,allowOutOfRange=>1)}).')';
            $i++;
            $phasedLinkC_last = $phasedLinkC;
        }
        $histogram->load_data(vA=>$LocRegUnitTime, vB=>\@phasedLinkC_log, fill=>\%fill);
    }
    # SVG
    my $svg = initialize_SVG_obj;
    $histogram->draw(svg_obj=>$svg);
    $amountGC->draw_quick(svg_obj=>$svg);
    output_SVG_file(svg_obj=>$svg, svg_file=>$V_Href->{svg_file});

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tdraw local region stat OK.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
