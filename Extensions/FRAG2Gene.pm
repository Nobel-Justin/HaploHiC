package HaploHiC::Extensions::FRAG2Gene;

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use List::Util qw/ max min sum /;
use Data::Dumper;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist trible_run_for_success /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Array qw/ binarySearch /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen /;
use HaploHiC::LoadOn;
use HaploHiC::GetPath qw/ GetPath /;
use HaploHiC::Extensions::JuicerDump qw/ load_chr_Things /;
use HaploHiC::PhasedHiC::dumpContacts qw/ load_enzyme_site_list /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              DumpFRAGtoGeneLevel
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::Extensions::FRAG2Gene';
#----- version --------
$VERSION = "0.09";
$DATE = '2019-01-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        DumpFRAGtoGeneLevel
                        FRAG2gene_contact
                        convert_frag_to_gene_contacts
                        load_region_list
                        load_gene_info
                        load_gene_psl
                        load_gene_list
                        check_files
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} FRAG_gene <[Options]>
     
     Options:

       # obtain gene-level contacts based on results of juicerDump func, only for FRAG:1  #

       # Inputs and Outputs #
        -fragct [s]  FRAG mode contacts result. <required>
                      Note: 1) bin_size must be 1.
                            2) gene-level contacts will output in same folder.
                            3) merged.txt result from 'juicerDump' function with FRAG mode.
        -db_dir [s]  database folder made by 'juicer_db' function. <required>
        -ref_v  [s]  version of reference genome, see contents under 'db_dir'. <required>
        -output [s]  gene-level contacts output. [optional]
        -glist  [s]  gene list file which you want to deal with ONLY. [optional]
                      Note: default to find contacts between genes in this list.
        -regbed [s]  one-based BED list file of region to take into contact map. [optional]
                      Note: use this option to add concerned region(s), e.g., enhancer.

       # Options #
        -enzyme [s]  enzyme type. <required>
        -gtype  [s]  gene types (column NO.13) in gene-PSL file. ['protein_coding,lincRNA,miRNA']
        -gextd  [s]  to extend gene region bilaterally. [0]
                      Note: you may consider the promoter and terminator as the gene body, such as 1E3.
        -gextr  [s]  extending gene region to such ceiling length. [0]
                      Note: 1) default is 0, means disabled.
                            2) only works on extended region, not gene body (whatever length).
        -glsgsd      to find contacts of gene_list (-glist) genes with all genes. [disabled]
        -mindis [i]  minimum distance between two genes/regions to calculate contacts. [1000]
                      Note: gene/region pairs having overlap are filtered out.
        -hInter      denote the FRAG contacts are inter-haplotype, default is intra-haplotype.

        -h|help      Display this help info.

     Version:
        $VERSION on $DATE

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
            [ frag_contact => undef ],
            [ db_dir => undef ],
            [ g2g_contact_output => undef ],
            [ gene_list => undef ],
            [ user_region_bed => undef ],

            # options
            [ enzyme_type => undef ],
            [ gene_type => 'protein_coding,lincRNA,miRNA' ],
            [ gzip_output => 1 ],
            [ gene_extent => 0 ],
            [ gene_extent_to_reach => 0 ],
            [ gene_list_for_single_side => 0 ],
            [ minimum_cmp_distance => 1000 ],
            [ is_hInter => 0 ],
            [ isFrom_juicerDump => 0 ],
            [ isFrom_dumpContacts => 0 ],
            [ haploCount => undef ],
            [ hapComb => {} ],

            # intermediate variants
            [ frag_contact_dir => undef ],
            [ gene_psl => undef ],
            [ enzyme_site => undef ],
            [ GeneOrUserReg_FRAG_info_report => undef ],
            [ user_region => {} ],
            [ required_gene => {} ],
            [ need_chr => {} ],
            [ gene_info => {} ],
            [ chr2enzymePos => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['frag_contact'],
                                  ['db_dir'],
                                  ['g2g_contact_output'],
                                  ['gene_list'],
                                  ['juicer_dir']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-fragct:s" => \$V_Href->{frag_contact},
        "-db_dir:s" => \$V_Href->{db_dir},
        "-output:s" => \$V_Href->{g2g_contact_output},
        "-glist:s"  => \$V_Href->{gene_list},
        "-regbed:s" => \$V_Href->{user_region_bed},
        # options
        "-ref_v:s"  => \$V_Href->{ref_version},
        "-enzyme:s" => \$V_Href->{enzyme_type},
        "-gtype:s"  => \$V_Href->{gene_type},
        "-gextd:s"  => \$V_Href->{gene_extent},
        "-gextr:s"  => \$V_Href->{gene_extent_to_reach},
        "-glsgsd"   => \$V_Href->{gene_list_for_single_side},
        "-mindis:i" => \$V_Href->{minimum_cmp_distance},
        # "-gzip"     => \$V_Href->{gzip_output},
        "-hInter"   => \$V_Href->{is_hInter},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || ! file_exist(filePath=>$V_Href->{frag_contact})
             || !defined $V_Href->{db_dir} || !-d $V_Href->{db_dir}
             || !defined $V_Href->{ref_version}
             || !defined $V_Href->{enzyme_type}
             || $V_Href->{gene_extent_to_reach} < 0
             || $V_Href->{gene_extent} < 0
             || $V_Href->{minimum_cmp_distance} <= 0
            );
}

#--- calculate gene-level contacts from dump(FRAG) result ---
sub DumpFRAGtoGeneLevel{
    # basic
    &check_files;
    load_chr_Things;

    # load restriction enzyme site
    load_enzyme_site_list;

    # load concern region if user provided
    &load_region_list;

    # load gene information from psl
    &load_gene_info;

    # load FRAG contacts and convert to gene/region-level
    ## output gene/region-to-gene/region contacts
    &FRAG2gene_contact;

    # output gene/region's enzyme FRAG information
    &report_GeneOrUserReg_FRAGidx;
}

#--- write gene/region FRAG idx ---
sub report_GeneOrUserReg_FRAGidx{

    # variants
    my $userReg_Href = $V_Href->{user_region};
    my $GeneInfo_Href = $V_Href->{gene_info};

    my %type2Href = ( 'gene_list' => $GeneInfo_Href,
                      'user_reg' => $userReg_Href   );

    open (GURFIR, Try_GZ_Write($V_Href->{GeneOrUserReg_FRAG_info_report})) || die "cannot write gene/region FRAG info: $!\n";
    print GURFIR '#'.join("\t", 'name',
                                'chr',
                                'st_pos',
                                'ed_pos',
                                'st_posExt',
                                'ed_posExt',
                                'strand',
                                'gene_type',
                                'st_FRAG(idx,edge)',
                                'ed_FRAG(idx,edge)',
                                'cover_ratio'      )."\n";
    for my $type ( reverse sort keys %type2Href ){
        my $Href = $type2Href{$type};
        for my $chr (sort keys %$Href){
            for my $tag (sort keys %{$Href->{$chr}}){
                my $tagInfo_Href = $Href->{$chr}->{$tag};
                my $st_pos = $tagInfo_Href->{stpos};
                my $ed_pos = $tagInfo_Href->{edpos};
                my $st_posExt = $tagInfo_Href->{stpos_ext} || $st_pos;
                my $ed_posExt = $tagInfo_Href->{edpos_ext} || $ed_pos;
                my $gtype = $tagInfo_Href->{gtype} || 'user_region';
                my $strand = $tagInfo_Href->{strand} || 'NA';
                # $V_Href->{chr2enzymePos}->{$chr} = \@info;
                my $st_posExt_FRAG_idx = $tagInfo_Href->{frag_stI};
                my $ed_posExt_FRAG_idx = $tagInfo_Href->{frag_edI};
                my $st_posExt_FRAG_idx_Uedge = $V_Href->{chr2enzymePos}->{$chr}->[max(0,$st_posExt_FRAG_idx-1)];
                my $ed_posExt_FRAG_idx_Dedge = $V_Href->{chr2enzymePos}->{$chr}->[$ed_posExt_FRAG_idx];
                # sometimes, exp, chrM gene!
                if( !defined $ed_posExt_FRAG_idx_Dedge ){
                    $ed_posExt_FRAG_idx_Dedge = $V_Href->{chr2enzymePos}->{$chr}->[-1];
                    warn "<WARN>\tset ed_posExt_FRAG_idx_Dedge of gene/region ($tag) as endpos of its chr ($chr).\n";
                }
                my $FRAG_lenSum = $ed_posExt_FRAG_idx_Dedge - $st_posExt_FRAG_idx_Uedge + 1;
                my $cover_ratio = sprintf "%.3f", ($ed_posExt - $st_posExt) / $FRAG_lenSum;
                print GURFIR join("\t", $tag,
                                        $chr,
                                        $st_pos,
                                        $ed_pos,
                                        $st_posExt,
                                        $ed_posExt,
                                        $strand,
                                        $gtype,
                                        $st_posExt_FRAG_idx.',U:'.$st_posExt_FRAG_idx_Uedge,
                                        $ed_posExt_FRAG_idx.',D:'.$ed_posExt_FRAG_idx_Dedge,
                                        $cover_ratio )."\n";
            }
        }
    }
    close GURFIR;

    `gzip -f $V_Href->{GeneOrUserReg_FRAG_info_report}` if($V_Href->{gzip_output});

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\treport gene/region FRAG infomation DONE\n"
                         ."\t$V_Href->{GeneOrUserReg_FRAG_info_report}\n";
}

#--- load FRAG contacts and convert to gene-level ---
sub FRAG2gene_contact{

    open (GCONT, Try_GZ_Write($V_Href->{g2g_contact_output})) || die "cannot write gene-level contact: $!\n";

    my %frag_contact;
    my ($last_chr_a, $last_chr_b) = ('N/A', 'N/A');
    my @splitIdx = $V_Href->{isFrom_dumpContacts} ? (0,1,3,4,6) : (0..4);
    open (FRAG, Try_GZ_Read($V_Href->{frag_contact})) || die "cannot read FRAG contacts file: $!\n";
    while(<FRAG>){
        next if(/^#/);
        my ($chr_a, $idx_a, $chr_b, $idx_b, $contact) = (split)[@splitIdx];
        # extract hap
        my ($hap_a, $hap_b) = ('hx', 'hx');
        if($V_Href->{isFrom_dumpContacts}){
            ($hap_a, $chr_a) = ($chr_a =~ /^(h\d+):(\S+)$/);
            ($hap_b, $chr_b) = ($chr_b =~ /^(h\d+):(\S+)$/);
        }
        # filter chr
        if(    !exists $V_Href->{need_chr}->{$chr_a}
            || !exists $V_Href->{need_chr}->{$chr_b}
        ){
            next;
        }
        # testing chr turn (ascending sorted)
        if(    (    $chr_a ne $chr_b
                 && $V_Href->{ChrThings}->{$chr_a}->{turn} > $V_Href->{ChrThings}->{$chr_b}->{turn}
               )
            || (    $chr_a ne $last_chr_a
                 && $last_chr_a ne 'N/A'
                 && $V_Href->{ChrThings}->{$chr_a}->{turn} < $V_Href->{ChrThings}->{$last_chr_a}->{turn}
               )
            || (    $chr_b ne $last_chr_b && $chr_a eq $last_chr_a
                 && $last_chr_a ne 'N/A'
                 && $V_Href->{ChrThings}->{$chr_b}->{turn} < $V_Href->{ChrThings}->{$last_chr_b}->{turn}
               )
        ){
            warn_and_exit "<ERROR>\tchromosme turn in FRAG contact file is wrong.\n"
                                ."\tlast_chr_a=$last_chr_a, "."last_chr_b=$last_chr_b\n"
                                ."\tcurr_chr_a=$chr_a, "     ."curr_chr_b=$chr_b\n";
        }
        # testing idx turn (ascending sorted)
        if(    $chr_a eq $chr_b
            && $idx_a  > $idx_b
        ){
            warn_and_exit "<ERROR>\tindex turn in FRAG contact file is wrong.\n"
                                ."\tchr_a=$chr_a,idx_a=$idx_a,chr_b=$chr_b,idx_b=$idx_b\n";
        }
        # new chr-pair comes
        if(    $last_chr_a ne $chr_a
            || $last_chr_b ne $chr_b
        ){
            &convert_frag_to_gene_contacts(chr_a=>$last_chr_a, chr_b=>$last_chr_b, FRAG_Contact_Href=>\%frag_contact);
            %frag_contact = ();
            # update
            $last_chr_a = $chr_a;
            $last_chr_b = $chr_b;
        }
        # record
        my $hapComb = "$hap_a,$hap_b";
        $V_Href->{hapComb}->{$hapComb} ++;
        $frag_contact{$idx_a}{$idx_b}{$hapComb} = $contact;
    }
    close FRAG;
    # last chr-chr pair
    &convert_frag_to_gene_contacts(chr_a=>$last_chr_a, chr_b=>$last_chr_b, FRAG_Contact_Href=>\%frag_contact);
    %frag_contact = ();

    close GCONT;

    `gzip -f $V_Href->{g2g_contact_output}` if($V_Href->{gzip_output});

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tconvert frag contacts to gene-level contacts DONE\n"
                         ."\t$V_Href->{g2g_contact_output}\n";
}

#--- convert frag contacts to gene contacts ---
sub convert_frag_to_gene_contacts{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chr_a = $parm{chr_a};
    my $chr_b = $parm{chr_b};
    my $FRAG_Contact_Href = $parm{FRAG_Contact_Href};

    # variants
    my $userReg_Href = $V_Href->{user_region};
    my $GeneInfo_Href = $V_Href->{gene_info};

    my %type2Href = ( 'gene_list' => $GeneInfo_Href,
                      'user_reg' => $userReg_Href   );

    my $inform_ab_chr_bool = 0;
    my $is_same_chr_bool = ( $chr_a eq $chr_b );
    my %used_tag_a;
    # 'left' side
    for my $type_a ( sort keys %type2Href ){
        my $Href_a = $type2Href{$type_a};
        next if( !exists $Href_a->{$chr_a} );
        my $chrHref_a = $Href_a->{$chr_a};
        # 'right' side
        for my $type_b ( sort keys %type2Href ){
            my $Href_b = $type2Href{$type_b};
            next if( !exists $Href_b->{$chr_b} );
            my $chrHref_b = $Href_b->{$chr_b};
            # 'left' gene/region, ascending order frag_index sort
            for my $tag_a (sort {$chrHref_a->{$a}->{frag_stI} <=> $chrHref_a->{$b}->{frag_stI}} keys %$chrHref_a){
                # required genes?
                my $is_tag_a_required_gene = exists $V_Href->{required_gene}->{$tag_a};
                # filter required genes
                if(    $type_a eq 'gene_list'
                    && $type_b eq 'gene_list'
                    && defined $V_Href->{gene_list}
                    && !$V_Href->{gene_list_for_single_side}
                    && !$is_tag_a_required_gene
                ){
                    next;
                }
                # FRAG interval_a
                my $infoHref_a = $chrHref_a->{$tag_a};
                my $frag_stI_a = $infoHref_a->{frag_stI};
                my $frag_edI_a = $infoHref_a->{frag_edI};
                # 'right' gene/region, ascending order frag_index sort
                for my $tag_b (sort {$chrHref_b->{$a}->{frag_stI} <=> $chrHref_b->{$b}->{frag_stI}} keys %$chrHref_b){
                    # not same gene/region
                    next if( !$V_Href->{is_hInter} && $tag_a eq $tag_b );
                    # if it is ever used in 'left', never use it at 'right'
                    next if( exists $used_tag_a{$tag_b} );
                    # required genes?
                    my $is_tag_b_required_gene = exists $V_Href->{required_gene}->{$tag_b};
                    # filter required genes
                    if(    $type_a eq 'gene_list'
                        && $type_b eq 'gene_list'
                        && defined $V_Href->{gene_list}
                        && !( $V_Href->{gene_list_for_single_side} && ($is_tag_a_required_gene || $is_tag_b_required_gene) )
                        && !(!$V_Href->{gene_list_for_single_side} && ($is_tag_a_required_gene && $is_tag_b_required_gene) )
                    ){
                        next;
                    }
                    # chck overlap (distance)
                    ## ignore inter-haplotype neighbor genes
                    my $infoHref_b = $chrHref_b->{$tag_b};
                    if(    !$V_Href->{is_hInter} # allows inter-haplotype same tag
                        && $is_same_chr_bool
                        && Get_Two_Seg_Olen( $infoHref_a->{stpos_ext} || $infoHref_a->{stpos},
                                             $infoHref_a->{edpos_ext} || $infoHref_a->{edpos},
                                             $infoHref_b->{stpos_ext} || $infoHref_b->{stpos},
                                             $infoHref_b->{edpos_ext} || $infoHref_b->{edpos},
                                             $V_Href->{minimum_cmp_distance}
                                            )
                    ){
                        next;
                    }
                    # FRAG interval_b
                    my $frag_stI_b = $infoHref_b->{frag_stI};
                    my $frag_edI_b = $infoHref_b->{frag_edI};
                    # prepare compare region pairs
                    my $reverse = 0;
                    my $cmp_itv;
                    if(    $is_same_chr_bool
                        && $frag_stI_a > $frag_stI_b # must change to small -> large
                    ){
                        $cmp_itv = [ $frag_stI_b, $frag_edI_b, $frag_stI_a, $frag_edI_a ];
                        $reverse = 1;
                    }
                    else{
                        $cmp_itv = [ $frag_stI_a, $frag_edI_a, $frag_stI_b, $frag_edI_b ];
                    }
                    ## avoid shared frags in the second interval
                    if(    !$V_Href->{is_hInter} # allows inter-haplotype overlap tag
                        && $is_same_chr_bool
                        && $cmp_itv->[2] <= $cmp_itv->[1]
                    ){
                        warn "<WARN>\tcomparing $tag_a [$frag_stI_a, $frag_edI_a] vs $tag_b [$frag_stI_b, $frag_edI_b], "
                                    ."shrink [$cmp_itv->[2], $cmp_itv->[3]] due to overlap.\n";
                        $cmp_itv->[2]  = $cmp_itv->[1] + 1;
                    }
                    # accumulate frag contacts
                    my %hapCombToContact;
                    for my $I_1 ( $cmp_itv->[0] .. $cmp_itv->[1] ){
                        next if !exists $FRAG_Contact_Href->{$I_1};
                        for my $I_2 ( $cmp_itv->[2] .. $cmp_itv->[3] ){
                            next if !exists $FRAG_Contact_Href->{$I_1}->{$I_2};
                            for my $hapComb (sort keys %{$V_Href->{hapComb}}){
                                my $hapCombToRecord = $reverse ? join(',', (split /,/, $hapComb)[1,0]) : $hapComb;
                                $hapCombToContact{$hapCombToRecord} += ($FRAG_Contact_Href->{$I_1}->{$I_2}->{$hapComb} || 0);
                            }
                        }
                    }
                    # output
                    print GCONT join("\t", $tag_a,
                                           $tag_b,
                                           map {( "$_:" . ($hapCombToContact{$_} || 0) )} sort keys %{$V_Href->{hapComb}}
                                    )."\n";
                    # mark used as at 'left' side
                    # to avoid a<->b and b<->a once a and b is on the same chr
                    $used_tag_a{$tag_a} = 1;
                    # mark to inform
                    $inform_ab_chr_bool = 1;
                }
            }
        }
    }

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tconvert frag contacts to gene-level contacts ($chr_a vs $chr_b) DONE\n" if($inform_ab_chr_bool);
}

#--- load region bed list ---
sub load_region_list{

    # variants
    my $userReg_Href = $V_Href->{user_region};

    return if( !defined $V_Href->{user_region_bed} );

    open (REGBED, Try_GZ_Read($V_Href->{user_region_bed})) || die "fail read user_region_bed: $!\n";
    while(<REGBED>){
        next if(/^#/);
        my ($chr, $stpos, $edpos) = (split)[0,1,2];
        next if( !exists $V_Href->{chr2enzymePos}->{$chr} );
        my $region_tag = join('-', $chr, $stpos, $edpos);
        # map pos to FRAG index
        my $frag_stI = binarySearch(query => $stpos, array => $V_Href->{chr2enzymePos}->{$chr}); 
        my $frag_edI = binarySearch(query => $edpos, array => $V_Href->{chr2enzymePos}->{$chr});
        # record
        $V_Href->{need_chr}->{$chr} = 1;
        $userReg_Href->{$chr}->{$region_tag} = { chr=>$chr, stpos=>$stpos, edpos=>$edpos,
                                                 frag_stI=>$frag_stI, frag_edI=>$frag_edI };
    }
    close REGBED;

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tread user region bed list DONE\n";
}

#--- load gene's info ---
sub load_gene_info{

    # load gene_list if has
    &load_gene_list;

    # load gene info from gene-PSL file
    &load_gene_psl;
}

#--- load gene PSL ---
sub load_gene_psl{

    # variants
    my $GeneInfo_Href = $V_Href->{gene_info};

    # prepare required gene types
    for my $key (qw/ gene_type /){
        my %value = map {($_,1)} split /,+/, $V_Href->{$key};
        $V_Href->{$key} = \%value;
    }

    # read gene-psl file
    open (GPSL, Try_GZ_Read($V_Href->{gene_psl})) || die "fail read gene list: $!\n";
    while(<GPSL>){
        next if(/^#/);
        my ($gname, $strand, $gtype, $chr, $stpos, $edpos) = (split)[9,8,12,13,15,16];
        next if( !exists $V_Href->{chr2enzymePos}->{$chr} );
        next if( !exists $V_Href->{gene_type}->{$gtype} );
        if(    defined $V_Href->{gene_list}
            # && !defined $V_Href->{user_region_bed}
            && !$V_Href->{gene_list_for_single_side}
            && !exists $V_Href->{required_gene}->{$gname}
        ){
            next;
        }
        # pos process
        $stpos ++;
        my $geneLen = $edpos - $stpos + 1;
        my $extend = $V_Href->{gene_extent};
        # extent to reach
        if(    $V_Href->{gene_extent_to_reach} > 0
            && $geneLen+2*$extend > $V_Href->{gene_extent_to_reach}
        ){
            $extend = int( max( ($V_Href->{gene_extent_to_reach} - $geneLen) / 2, 0) );
        }
        my $stpos_ext = $stpos - $extend;
        my $edpos_ext = $edpos + $extend;
        # map gene pos to FRAG index
        my $frag_stI = binarySearch(query => $stpos_ext, array => $V_Href->{chr2enzymePos}->{$chr});
        my $frag_edI = binarySearch(query => $edpos_ext, array => $V_Href->{chr2enzymePos}->{$chr});
        # record
        $V_Href->{need_chr}->{$chr} = 1;
        $GeneInfo_Href->{$chr}->{$gname} = { chr=>$chr, stpos=>$stpos, edpos=>$edpos,
                                             stpos_ext=>$stpos_ext, edpos_ext=>$edpos_ext,
                                             gname=>$gname, gtype=>$gtype, strand=>$strand,
                                             frag_stI=>$frag_stI, frag_edI=>$frag_edI
                                            };
    }
    close GPSL;

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tread gene psl DONE\n";
}

#--- load gene list ---
sub load_gene_list{

    return if( !defined $V_Href->{gene_list} );

    open (GLIST, Try_GZ_Read($V_Href->{gene_list})) || die "fail read gene list: $!\n";
    while(<GLIST>){
        next if(/^#/);
        my ($gname) = (split)[0];
        $V_Href->{required_gene}->{$gname} = 1;
    }
    close GLIST;

    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tread gene list DONE\n";
}

#--- check dump(FRAG) result ---
sub check_files{
    # database file
    $V_Href->{chrLenFile} = GetPath(filekey => 'chrLenFile');
    $V_Href->{gene_psl} = GetPath(filekey => 'gene_psl');
    $V_Href->{enzyme_site} = GetPath(filekey => 'enzyme_site');
    if(    !-e $V_Href->{chrLenFile}
        || !-e $V_Href->{gene_psl}
        || !-e $V_Href->{enzyme_site}
    ){
        warn_and_exit "<ERROR>\tPlease make sure the existence of files below.\n"
                            ."\t$V_Href->{chrLenFile}\n"
                            ."\t$V_Href->{gene_psl}\n"
                            ."\t$V_Href->{enzyme_site}\n";
    }
    # match juicer's output folder structure
    $V_Href->{frag_contact_dir} = dirname($V_Href->{frag_contact});
    $V_Href->{isFrom_juicerDump} = $V_Href->{frag_contact_dir} =~ /aligned\b/ ? 1 : 0;
    # FRAG_binSize is 1 ?
    if(    $V_Href->{isFrom_juicerDump}
        && $V_Href->{frag_contact_dir} !~ /FRAG1$/
    ){
        warn_and_exit "<ERROR>\tPlease make sure this file is from 'juicerDump' function (FRAG mode) with bin_size as 1.\n"
                            ."\t$V_Href->{frag_contact}\n";
    }
    else{
        open (FRAG, Try_GZ_Read($V_Href->{frag_contact})) || die "cannot read FRAG contacts file: $!\n";
        while(<FRAG>){
            if(/^##dumpMode:\s(\S+),\sdumpBinSize:\s(\S+)/){
                $V_Href->{isFrom_dumpContacts} = 1;
                warn_and_exit "<ERROR>\tPlease make sure contacts file is from 'FRAG' mode with binSize as 1.\n"
                                    ."\t$V_Href->{frag_contact}\n" if($1 ne 'FRAG' || $2 != 1);
            }
            if(/^##enzyme:\s(\S+)/){
                warn_and_exit "<ERROR>\tEnzyme ($1) of contacts file is not input option ($V_Href->{enzyme_type}).\n"
                                    ."\t$V_Href->{frag_contact}\n" if($V_Href->{enzyme_type} ne $1);
            }
            if(/^##haploCount:\s(\S+)/){
                $V_Href->{haploCount} = $1;
                if($V_Href->{is_hInter}){
                    for my $i (1 .. $V_Href->{haploCount}){
                        $V_Href->{hapComb}->{"h$i,h$_"} = 0 for grep $_ != $i, (1 .. $V_Href->{haploCount});
                    }
                }
            }
            last if $.==10;
        }
        close FRAG;
    }
    # prepare outputs
    unless( defined $V_Href->{g2g_contact_output} ){
        $V_Href->{g2g_contact_output} = $V_Href->{frag_contact};
        $V_Href->{g2g_contact_output} =~ s/.gz$//;
        $V_Href->{g2g_contact_output} =~ s/.txt$//;
        $V_Href->{g2g_contact_output} = $V_Href->{g2g_contact_output} . ".gene-level.gextd$V_Href->{gene_extent}.gextr$V_Href->{gene_extent_to_reach}.contacts.txt";
    }
    else{
        $V_Href->{g2g_contact_output} =~ s/.gz$//;
    }
    $V_Href->{GeneOrUserReg_FRAG_info_report} = $V_Href->{g2g_contact_output}.'.GeneOrUserReg.FRAG_info';
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\tfiles checking DONE\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
