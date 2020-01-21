package HaploHiC::Extensions::ShellForJuicer;

use strict;
use warnings;
use Getopt::Long;
use Cwd qw/ abs_path /;
use File::Basename;
use FindBin qw/ $RealBin /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use HaploHiC::LoadOn;
use HaploHiC::GetPath qw/ GetPath /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              shell_to_run_juicer
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::Extensions::ShellForJuicer';
#----- version --------
$VERSION = "0.07";
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
                        shell_to_run_juicer
                        check_files_parm
                        submit_to_flux
                        write_shell
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} run_juicer <[Options]>
     
     Options:

       # run juicer from FASTQ files, and could also do dump work according to user settings #

       # Inputs and Outputs #
        -juicer [s]  path of juicer folder, which has sub-folder, e.g., 'CPU'. <required>
        -db_dir [s]  database folder made by 'juicer_db' function. <required>
        -outdir [s]  output directory where topdir will be created. <required>
                       Note: Please make the output directory before running.
        -fqlist [s]  sample's fastq list, one line one file. <required>
        -sample [s]  sample ID. <required>

       # Options #

        # for juicer.sh
        -ref_v  [s]  version of reference genome, see contents under 'db_dir'. <required>
        -enzyme [s]  enzyme type. <required>

        # for juicer_tools dump
        # following parameters could use values concatenated by commas.
        -hic_fl [s]  prefix of '.hic' file: 'inter_30' or 'inter'. ['inter_30']
        -ctype  [s]  'observed' or 'oe'. ['observed']
        -norm   [s]  kind of normalization, check juicer_tools dump usage. ['NONE']
        -bin_BP [s]  bin size for BP mode. ['1MB,100KB']
                      Allowed: 2.5MB, 1MB, 500KB, 250KB, 100KB, 50KB, 25KB, 10KB, 5KB
        -bin_FG [s]  bin size for FRAG mode. ['500,1']
                      Allowed: 500, 200, 100, 50, 20, 5, 2, 1

        # for flux system PBS
        -no_flux     do not append flux job header. [disabled]
        -descp  [s]  description of experiment, enclosed in single quotes.
        -flux_a [s]  your Flux allocation to submit jobs. [none]
        -flux_q [s]  your Flux queue to submit jobs. [fluxod]
        -jname  [s]  name the job on flux system. ['juicer_{sampleID}']
        -email  [s]  email to inform when job's status changes. [disabled]
        -mact   [s]  action to email, (a) fails (b) starts (e) ends. [abe]
        -rtime  [s]  walltime to run your job, format: d:hh:mm:ss. [3:00:00:00]
                       Note: reckon the Walltime by { expect + 10-15% }.
        -submit      submit PBS to flux, or else, manually. [disabled]

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
            [ juicer_dir => undef ],
            [ db_dir => undef ],
            [ outdir => undef ],
            [ fqlist => undef ],
            [ sample => undef ],

            # settings shared
            [ skip_flux => 0 ],
            [ submit => 0 ],
            [ Flux_allocation => 'none' ], # indikar_fluxod
            [ Flux_queue => 'fluxod' ],
            [ job_name => undef ],
            [ email => undef ],
            [ action_to_email => 'abe' ],
            [ walltime_to_run => '3:00:00:00' ],

            # setting for 1st step
            [ enzyme_type => undef ],
            [ exp_descp => '' ],
            [ bwaMEMthreads => 4 ],

            # settings for 2nd step
            [ dump_hic_type => 'inter_30' ],
            [ dump_count_type => 'observed' ],
            [ dump_norm_method => 'NONE' ],
            [ dump_bin_size_BP => '1MB,100KB' ],
            [ dump_bin_size_FRAG => '500,1' ],

            # intermediate variants
            [ JuicerTopdir => undef ],
            [ CPU_dir => undef ],
            [ juicer_sh => undef ],
            [ HaploHiCperlBin => undef ],
            [ GenomeRefFa => undef ],
            [ chrLenFile => undef ],
            [ HeaderForSam => undef ],
            [ enzyme_site => undef ],
            [ JuicerFastqDir => undef ],
            [ JuicerShell => undef ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['outdir'],
                                  ['db_dir'],
                                  ['juicer_dir']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-juicer:s" => \$V_Href->{juicer_dir},
        "-db_dir:s" => \$V_Href->{db_dir},
        "-outdir:s" => \$V_Href->{outdir},
        "-fqlist:s" => \$V_Href->{fqlist},
        "-sample:s" => \$V_Href->{sample},
        # options
        ## juicer.sh
        "-ref_v:s"  => \$V_Href->{ref_version},
        "-enzyme:s" => \$V_Href->{enzyme_type},
        "-descp:s"  => \$V_Href->{exp_descp},
        ## juicer dump
        "-hic_fl:s" => \$V_Href->{hic_type},
        "-bin_BP:s" => \$V_Href->{dump_bin_size_BP},
        "-bin_FG:s" => \$V_Href->{dump_bin_size_FRAG},
        "-ctype:s"  => \$V_Href->{count_type},
        "-norm:s"   => \$V_Href->{norm_method},
        ## flux
        "-no_flux"  => \$V_Href->{skip_flux},
        "-flux_a:s" => \$V_Href->{Flux_allocation},
        "-flux_q:s" => \$V_Href->{Flux_queue},
        "-jname:s"  => \$V_Href->{job_name},
        "-email:s"  => \$V_Href->{email},
        "-mact:s"   => \$V_Href->{action_to_email},
        "-rtime:s"  => \$V_Href->{walltime_to_run},
        "-submit"   => \$V_Href->{submit},
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
             || !defined $V_Href->{ref_version}
             || !defined $V_Href->{outdir} || !-d $V_Href->{outdir}
             || !defined $V_Href->{fqlist} || !-e $V_Href->{fqlist}
             || !defined $V_Href->{sample}
             || !defined $V_Href->{enzyme_type}
            );
}

#--- run juicer from start ---
sub shell_to_run_juicer{
    &check_files_parm;

    &write_shell;

    &submit_to_flux;
}

#--- check files
sub check_files_parm{
    # requried juicer/database files to check
    $V_Href->{CPU_dir} = GetPath(filekey => 'CPU_dir');
    $V_Href->{juicer_sh} = GetPath(filekey => 'juicer_sh');
    $V_Href->{HaploHiCperlBin} = GetPath(filekey => 'HaploHiCperlBin');
    $V_Href->{chrLenFile}  = GetPath(filekey => 'chrLenFile');
    $V_Href->{HeaderForSam}= GetPath(filekey => 'HeaderForSam');
    $V_Href->{enzyme_site} = GetPath(filekey => 'enzyme_site');
    $V_Href->{GenomeRefFa} = GetPath(filekey => 'GenomeRefFa');

    if(    !-e $V_Href->{juicer_sh}
        || !-e $V_Href->{HaploHiCperlBin}
        || !-e $V_Href->{chrLenFile}
        || !-e $V_Href->{GenomeRefFa}
        || !-e $V_Href->{HeaderForSam}
        || !-e $V_Href->{enzyme_site}
    ){
        warn "Please make sure the existence of files/folder below.\n";
        warn "  $V_Href->{juicer_sh}\n";
        warn "  $V_Href->{HaploHiCperlBin}\n";
        warn "  $V_Href->{chrLenFile}\n";
        warn "  $V_Href->{GenomeRefFa}\n";
        warn "  $V_Href->{enzyme_site}\n";
        exit;
    }

    # topdir: workspace
    $V_Href->{JuicerTopdir} = GetPath(filekey => 'JuicerTopdir');

    # check qsub job if has, avoid mistakenly rerun
    my @shell_qsubInfo = glob( GetPath(filekey => 'JuicerShellGlob') );
    if( scalar(@shell_qsubInfo) != 0 ){
        my %jobid = map{ my ($jobid) = (`cat $_` =~ /^(\d+)/); ($jobid, $_) } @shell_qsubInfo;
        for my $jobid (sort {$a<=>$b} keys %jobid){
            my $qstat_returnWCL = `qstat $jobid 2>/dev/null | wc -l`;
            if( $qstat_returnWCL != 0 ){
                warn "The previous job ($jobid) is still on flux system.\n";
                warn "  check file $jobid{$jobid}\n";
                exit;
            }
        }
    }
    # set up fastq
    $V_Href->{JuicerFastqDir} = GetPath(filekey => 'JuicerFastqDir');
    # check and alert
    if( -e $V_Href->{JuicerFastqDir} ){
        warn "find existing 'fastq' folder in topdir:\n";
        warn "  $V_Href->{JuicerTopdir}\n";
        warn "please properly and CAREfully move or delete the 'fastq' folder!\n";
        exit;
    }
    # read fastq list
    my @fastq;
    open (FQLIST, $V_Href->{fqlist}) || die "fail read fqlist: $!\n";
    while(<FQLIST>){
        next if(/^#/);
        my ($fqFile) = (split)[0];
        if( !-e $fqFile ){
            warn "Please check the existence of fastq file in fqlist.\n";
            warn "  $fqFile\n";
            exit;
        }
        push @fastq, abs_path( $fqFile );
    }
    close FQLIST;
    # check fastq
    if( scalar(@fastq) == 0 ){
        warn "cannot find any available fastq file in fqlist:\n";
        warn "  $V_Href->{fqlist}\n";
        exit;
    }

    # parameters for juicer dump
    for my $key (qw/ hic_type dump_bin_size_FRAG dump_bin_size_BP count_type norm_method /){
        my %value = map {($_,1)} split /,+/, $V_Href->{$key};
        $V_Href->{$key} = \%value;
    }
    for my $BP_bin_size ( sort keys %{$V_Href->{dump_bin_size_BP}} ){
        if( !exists $V_Href->{dump_allowBinSize}->{BP}->{$BP_bin_size} ){
            warn "Only allows BP_bin_size: 2.5MB, 1MB, 500KB, 250KB, 100KB, 50KB, 25KB, 10KB, 5KB\n";
            exit;
        }
    }
    for my $FRAG_bin_size ( sort keys %{$V_Href->{dump_bin_size_FRAG}} ){
        if( !exists $V_Href->{dump_allowBinSize}->{FRAG}->{$FRAG_bin_size} ){
            warn "Only allows FRAG_bin_size: 500, 200, 100, 50, 20, 5, 2, 1\n";
            exit;
        }
    }

    # make workspace
    `mkdir -p $V_Href->{JuicerTopdir}`; # may already exist

    # soft link of fqfile
    `mkdir -p $V_Href->{JuicerFastqDir}`;
    `ln -s $_ $V_Href->{JuicerFastqDir}/` for @fastq;

    # sweep workspace
    `rm -rf $V_Href->{JuicerTopdir}/[^f]*`;

    # soft link of chrom.sizes file
    ## Comment on 2018-04-02
    ## finnaly, we make sure if we use 'relative path' of chrom.sizes file,
    ## the juicerbox desktop will recognize the genome version of .hic file.
    `ln -s $V_Href->{chrLenFile} $V_Href->{JuicerTopdir}/`;
    $V_Href->{chrLenFile} = basename($V_Href->{chrLenFile});

    # PBS shell
    $V_Href->{job_name} ||= 'juicer_'.$V_Href->{sample};
    $V_Href->{JuicerShell} = GetPath(filekey => 'JuicerShell');

    # inform
    print "files/folder checking DONE\n";
}

#--- write shell
sub write_shell{
    open (PBS, ">$V_Href->{JuicerShell}") || die "cannot write pbs shell file\n";

    unless( $V_Href->{skip_flux} ){
        # header for pbs
        print PBS "#   what account and queue to use.\n";
        print PBS "#PBS -A $V_Href->{Flux_allocation}\n";
        print PBS "#PBS -q $V_Href->{Flux_queue}\n";
        print PBS "#   job name - change to whatever you need.\n";
        print PBS "#PBS -N $V_Href->{job_name}\n";
        if( defined $V_Href->{email} ){
            print PBS "#   email and conditions: (a) fails (b) starts (e) ends\n";
            print PBS "#PBS -M $V_Href->{email}\n";
            print PBS "#PBS -m $V_Href->{action_to_email}\n";
        }
        print PBS "#   sends all environment variables on the login node.\n";
        print PBS "#PBS -V\n";
        print PBS "#   number of nodes and processors. max ppn = 12. Walltime: expect + 10-15% d:hh:mm:ss\n";
        print PBS "#PBS -l nodes=1:ppn=4,mem=15gb,qos=flux,walltime=$V_Href->{walltime_to_run}\n";
        print PBS "#   Output folder\n";
        print PBS "#PBS -j oe\n";
        print PBS "#PBS -o $V_Href->{job_name}.oe\n";
        print PBS "\n";
        print PBS "#   print nodes\n";
        print PBS "\n";
        print PBS 'if [ -e "$PBS_NODEFILE" ] ; then'."\n";
        print PBS '  echo "Running on"'."\n";
        print PBS '  cat $PBS_NODEFILE'."\n";
        print PBS "fi\n";
        print PBS 'if [ -d "$PBS_O_WORKDIR" ] ; then'."\n";
        print PBS '  cd $PBS_O_WORKDIR'."\n";
        print PBS "fi\n";
        print PBS "\n";
        print PBS "#  Put your job commands here:\n";
        print PBS "\n";
        print PBS 'echo "Running from $(pwd)";'."\n";
        print PBS "\n";
        # set error action
        print PBS "# exit immediately if any command exits with a non-zero status.\n";
        print PBS "set -e;\n\n";
        # load modules
        print PBS "# load modules\n";
        ## note that when bwa updates, make sure its previous indexed genome ref is still available,
        ## or else, use new version of bwa to make index again. --- by Wenlong Jia
        print PBS "module load bwa/0.7.15;\n"; # check on 2018-03-17
        print PBS "module load samtools/1.3.1;\n"; # check on 2018-03-17
        print PBS "\n";
    }

    # juicer
    print PBS "# run juicer\n";
    my $descp_parm = length($V_Href->{exp_descp}) == 0 ? '' : "-a '$V_Href->{exp_descp}'";
    print PBS "( bash $V_Href->{juicer_sh} -t $V_Href->{bwaMEMthreads} -z $V_Href->{GenomeRefFa} -p $V_Href->{chrLenFile} -y $V_Href->{enzyme_site} -d $V_Href->{JuicerTopdir} -D $V_Href->{CPU_dir} -s $V_Href->{enzyme_type} $descp_parm ) && ( date && echo juicer completes )\n";
    print PBS "\n";
    # sam to bam
    print PBS "# convert sam files to bam\n";
    print PBS "ls $V_Href->{JuicerTopdir}/*/*.sam | sed 's/.sam\$//' | while read samFilePref; do set -o pipefail; ( ( cat $V_Href->{HeaderForSam}; cat \$samFilePref.sam ) | samtools view -b -o \$samFilePref.bam ) && ( rm -rf \$samFilePref.sam ) && ( echo convert to bam OK: \$samFilePref.sam ); done\n";
    print PBS "# compress all large txt files (>100M)\n";
    print PBS "find $V_Href->{JuicerTopdir}/ -size +100M | grep '.txt\$' | while read largeTxt; do ( gzip \$largeTxt ) && ( echo compress largeTxt OK: \$largeTxt ); done\n";

    # dump from .hic file to contacts count txt
    for my $hic_type ( sort keys %{$V_Href->{hic_type}} ){
        for my $count_type ( sort keys %{$V_Href->{count_type}} ){
            for my $norm_method ( sort keys %{$V_Href->{norm_method}} ){
                # BP mode
                for my $BP_bin_size ( sort keys %{$V_Href->{dump_bin_size_BP}} ){
                    my $info_str = "juicer dump BP hic_file: $hic_type, BP_bin_size: $BP_bin_size, count_type: $count_type, norm_method: $norm_method";
                    print PBS "\n# $info_str\n";
                    print PBS "( perl $V_Href->{HaploHiCperlBin} juicerDump -dpmode BP -juicer $V_Href->{juicer_dir} -topdir $V_Href->{JuicerTopdir} -db_dir $V_Href->{db_dir} -ref_v $V_Href->{ref_version} -hic_fl $hic_type -bin $BP_bin_size -ctype $count_type -norm $norm_method) && ( echo $info_str OK )\n";
                }
                # FRAG mode
                for my $FRAG_bin_size ( sort keys %{$V_Href->{dump_bin_size_FRAG}} ){
                    my $info_str = "juicer dump FRAG hic_file: $hic_type, FRAG_bin_size: $FRAG_bin_size, count_type: $count_type, norm_method: $norm_method";
                    print PBS "\n# $info_str\n";
                    print PBS "( perl $V_Href->{HaploHiCperlBin} juicerDump -dpmode FRAG -enzyme $V_Href->{enzyme_type} -juicer $V_Href->{juicer_dir} -topdir $V_Href->{JuicerTopdir} -db_dir $V_Href->{db_dir} -ref_v $V_Href->{ref_version} -hic_fl $hic_type -bin $FRAG_bin_size -ctype $count_type -norm $norm_method) && ( echo $info_str OK )\n";
                }
            }
        }
    }

    close PBS;

    # inform
    print "write PBS file:\n  $V_Href->{JuicerShell}\n";
}

#--- submit to flux
sub submit_to_flux{
    return unless( $V_Href->{submit} );
    # inform
    print "submit PBS to flux system:\n";
    system("cd $V_Href->{JuicerTopdir}; qsub $V_Href->{JuicerShell} | tee -a $V_Href->{JuicerShell}.qsub_info");
}

#--- 
1; ## tell the perl script the successful access of this module.
