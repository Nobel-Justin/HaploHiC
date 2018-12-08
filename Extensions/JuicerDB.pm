package HaploHiC::Extensions::JuicerDB;

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw/ catfile /;
use Cwd qw/ abs_path /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use HaploHiC::LoadOn;
use HaploHiC::GetPath qw/ GetPath /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              JuicerMakeDBwork
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::Extensions::JuicerDB';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-14';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        JuicerMakeDBwork
                        ref_genome_part
                        gene_anno_part
                        enzyme_site_part
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} juicer_db <[Options]>
     
     Options:

       # run juicer_tools dump func in FRAG mode and merge all chromosomes to one final file #

       # Inputs and Outputs #
        -juicer [s]  path of juicer folder, which has sub-folder 'misc'. <required>
        -db_dir [s]  existing folder to store database files. <required>
        -ref_fa [s]  fasta file of reference genome. <required>
        -ref_v  [s]  version of reference genome, as sub-folder name in 'db_dir'. <required>
        -gpsl   [s]  gene psl file. <required>

       # Tools and Database #
        -bwa    [s]  BWA. <required>
        -samt   [s]  SAMtools. <required>

       # Options #
        -enzyme [s]  enzyme type. <required>
        -skip_rfg    skip reference genome work. [disabled]
        -skip_gan    skip gene annotation work. [disabled]
        -skip_ezy    skip enzyme site work. [disabled]
        -ref_idx     set if the reference genome has been indexed by BWA. [disabled]
                      Note: if set, just creates soft-links to index files.

        -h|help      Display this help info.

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
            [ juicer_dir => undef ],
            [ db_dir => undef ],
            [ GenomeRefFa => undef ],
            [ gene_psl => undef ],

            # software and database
            [ bwa => undef ],
            [ samtools => undef ],

            # options
            [ enzyme_type => undef ],
            [ refHasIdx => 0 ],
            [ skipRefGenome => 0 ],
            [ skipGeneAnno => 0 ],
            [ skipEnzymeSite => 0 ],

            # intermediate variants
            [ chrLenFile => undef ],
            [ HeaderForSam => undef ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['juicer_dir'],
                                  ['db_dir'],
                                  # ['GenomeRefFa'],
                                  ['gene_psl'],
                                  ['bwa'],
                                  ['samtools']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-juicer:s" => \$V_Href->{juicer_dir},
        "-db_dir:s" => \$V_Href->{db_dir},
        "-ref_fa:s" => \$V_Href->{GenomeRefFa},
        "-ref_v:s"  => \$V_Href->{ref_version},
        "-gpsl:s"   => \$V_Href->{gene_psl},
        # tools
        "-bwa:s"    => \$V_Href->{bwa},
        "-samt:s"   => \$V_Href->{samtools},
        # options
        "-enzyme:s" => \$V_Href->{enzyme_type},
        "-ref_idx"  => \$V_Href->{refHasIdx},
        "-skip_rfg" => \$V_Href->{skipRefGenome},
        "-skip_gan" => \$V_Href->{skipGeneAnno},
        "-skip_ezy" => \$V_Href->{skipEnzymeSite},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !defined $V_Href->{juicer_dir} || !-d $V_Href->{juicer_dir}
             || !defined $V_Href->{db_dir} || !-d $V_Href->{db_dir}
             || !defined $V_Href->{GenomeRefFa} || !-e $V_Href->{GenomeRefFa}
             || !defined $V_Href->{ref_version}
             || !defined $V_Href->{gene_psl} || !-e $V_Href->{gene_psl}
             || !defined $V_Href->{enzyme_type}
             || !defined $V_Href->{bwa} || !-e $V_Href->{bwa}
             || !defined $V_Href->{samtools} || !-e $V_Href->{samtools}
            );
}

#--- run juicer_tools dump func in FRAG mode and merge all chromosomes ---
sub JuicerMakeDBwork{
    &ref_genome_part;

    &gene_anno_part;

    &enzyme_site_part;
}

#--- database files of reference genome ---
sub ref_genome_part{
    return if $V_Href->{skipRefGenome};

    # update to full_path
    $V_Href->{GenomeRefFa} = catfile( abs_path(dirname $V_Href->{GenomeRefFa}), basename $V_Href->{GenomeRefFa} );
    # create folder and link fasta
    my $RefFaInDB = GetPath(filekey => 'GenomeRefFa');
    my $BWA_Idx_Folder = dirname $RefFaInDB;
    `mkdir -p $BWA_Idx_Folder`;
    `ln -sf $V_Href->{GenomeRefFa} $RefFaInDB`;
    # index files
    my $cmd;
    if($V_Href->{refHasIdx}){
        # link index files
        `ln -sf $V_Href->{GenomeRefFa}.$_ $RefFaInDB.$_` for qw/ amb ann bwt pac sa /;
    }
    else{
        # do index, current BWA 0.7.13 -a [auto]
        $cmd = "$V_Href->{bwa} index $RefFaInDB";
        trible_run_for_success($cmd, 'BAW index of reference genome');
    }
    # do faidx
    $cmd = "$V_Href->{samtools} faidx $RefFaInDB";
    trible_run_for_success($cmd, 'faidx of reference genome');
    # chrLenFile
    my $chrLenFile = GetPath(filekey => 'chrLenFile');
    `cut -f1,2 $RefFaInDB.fai > $chrLenFile`;
    # HeaderForSam
    my $HeaderForSam = GetPath(filekey => 'HeaderForSam');
    $cmd = 'echo -e "@HD\tVN:1.5\tGO:none\tSO:none" > '.$HeaderForSam;
    `$cmd`;
    $cmd = 'awk \'{print "@SQ\tSN:"$1"\tLN:"$2}\' '."$RefFaInDB.fai >> $HeaderForSam";
    `$cmd`;
    # inform
    stout_and_sterr "[INFO]\tPrepare database files of reference genome OK.\n"
                         ."\t$BWA_Idx_Folder\n"
                         ."\t$chrLenFile\n"
                         ."\t$HeaderForSam\n";
}

#--- gene annotation ---
sub gene_anno_part{
    return if $V_Href->{skipGeneAnno};

    my $genePSLinDB = GetPath(filekey => 'gene_psl');
    # create folder
    my $geneAnno_Folder = dirname $genePSLinDB;
    `mkdir -p $geneAnno_Folder`;
    `ln -sf $V_Href->{gene_psl} $genePSLinDB`;
    # inform
    stout_and_sterr "[INFO]\tPrepare database files of gene annotation OK.\n"
                         ."\t$genePSLinDB\n";
}

#--- enzyme site ---
sub enzyme_site_part{
    return if $V_Href->{skipEnzymeSite};

    my $RefFaInDB = GetPath(filekey => 'GenomeRefFa');
    my @enzyme_type = split /,/, $V_Href->{enzyme_type};
    for my $enzyme_type (@enzyme_type){
        $V_Href->{enzyme_type} = $enzyme_type; # update
        my $enzyme_site = GetPath(filekey => 'enzyme_site');
        # folder
        my $enzymeFolder = dirname $enzyme_site;
        `mkdir -p $enzymeFolder`;
        # create enzyme site
        my $cmd = "cd $enzymeFolder; python $V_Href->{juicer_dir}/misc/generate_site_positions.py $enzyme_type $V_Href->{ref_version} $RefFaInDB";
        trible_run_for_success($cmd, "get $enzyme_type site");
        # inform
        stout_and_sterr "[INFO]\tPrepare enzyme site list of $enzyme_type OK.\n"
                             ."\t$enzyme_site\n";
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
