package HaploHiC::GetPath;

use strict;
use warnings;
use File::Spec::Functions qw/ catfile /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use HaploHiC::LoadOn;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
	          GetPath
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HaploHiC::GetPath';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--- get path ---
sub GetPath{
	# options
	shift if (@_ && $_[0] =~ /$MODULE_NAME/);
	my %parm = @_;
	my $filekey = $parm{filekey};

	# HaploHiC files
	if($filekey eq 'func_json'){
		return catfile($V_Href->{RealBin},'functions.json');
	}
	if($filekey eq 'HaploHiCperlBin'){
		return catfile($V_Href->{RealBin},'HaploHiC.pl');
	}
	# database file
	if($filekey eq 'chrLenFile'){
		return catfile($V_Href->{db_dir}, $V_Href->{ref_version}, $V_Href->{ref_version}.'.chrom.sizes');
	}
	if($filekey eq 'GenomeRefFa'){
		return catfile($V_Href->{db_dir}, $V_Href->{ref_version}, 'BWA_index', $V_Href->{ref_version}.'.fa');
	    ## Comment on 2018-03-23
	    ## Homo_sapiens_assembly19.fasta, it is required by juicer.
	    ## fixed name for online juicerbox: to recognize the genome version for further annotations.
	    ## (my $ref_version_number = $V_Href->{ref_version}) =~ s/\D//g;
	    ## $V_Href->{GenomeRefFa} = catfile( $V_Href->{db_dir}, $V_Href->{ref_version}, 'BWA_index', 'Homo_sapiens_assembly'.$ref_version_number.'.fasta' );
	}
	if($filekey eq 'enzyme_site'){
		return catfile($V_Href->{db_dir}, $V_Href->{ref_version}, 'enzyme_site', $V_Href->{enzyme_type}, $V_Href->{ref_version}.'_'.$V_Href->{enzyme_type}.'.txt');
	}
	if($filekey eq 'HeaderForSam'){
		return catfile($V_Href->{db_dir}, $V_Href->{ref_version}, $V_Href->{ref_version}.'.sam_header');
	}
	if($filekey eq 'juicer_tool_Jar'){
		return catfile($V_Href->{juicer_dir}, 'CPU', 'scripts', 'common', 'juicer_tools.jar');
	}
	if($filekey eq 'gene_psl'){
		return catfile($V_Href->{db_dir}, $V_Href->{ref_version}, 'Genes', 'gene.psl');
	}
	# juicer tools path
	if($filekey eq 'CPU_dir'){
		return catfile($V_Href->{juicer_dir}, 'CPU');
	}
	if($filekey eq 'juicer_sh'){
		return catfile($V_Href->{CPU_dir}, 'juicer.sh');
	}
	# juicer workspace
	if($filekey eq 'JuicerTopdir'){
		return catfile($V_Href->{outdir}, 'topdir_'.$V_Href->{sample});
	}
	if($filekey eq 'JuicerFastqDir'){
		return catfile($V_Href->{JuicerTopdir}, 'fastq');
	}
	if($filekey eq 'JuicerShellGlob'){
		return catfile($V_Href->{JuicerTopdir}, '*.qsub_info');
	}
	if($filekey eq 'JuicerShell'){
		return catfile($V_Href->{JuicerTopdir}, $V_Href->{job_name}.'.pbs');
	}
	if($filekey eq 'hic_file'){
		return catfile($V_Href->{JuicerTopdir}, 'aligned', $V_Href->{hic_type}.'.hic');
	}
	if($filekey eq 'JuicerAlignedDir'){
		return catfile($V_Href->{JuicerTopdir}, 'aligned');
	}
	# report on reading phased VCF
	if($filekey eq 'phasedMut_report'){
		return catfile($V_Href->{outdir}, 'VCF_load.phased_mut.report');
	}
	if($filekey eq 'closeInDelReport'){
		return catfile($V_Href->{outdir}, 'VCF_load.phased_mut.closeInDel.list');
	}
	# if($filekey eq ''){
	# 	return catfile();
	# }

	# reach here, not found
	warn_and_exit "<ERROR>\tunrecognized filekey $filekey in $MODULE_NAME.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
