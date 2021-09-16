# PLEASE use absolute path
# TOOLS
# HaploHiC
# https://github.com/deepomicslab/HaploHiC
HaploHiC=$PWD/../../HaploHiC.pl;

# enzyme type, e.g., MboI
enzyme=__PLEASE_SET__THE_ENZYME_TYPE__;

# juicer folder, which has sub-folder 'misc'
Juicer_folder=__PLEASE_SET__THE_JUICER_FOLDER_PATH__
# BWA, version: 0.7.17
BWA=__PLEASE_SET__THE_BWA_PATH__
# SAMtools, version: 1.9
SAMtools=__PLEASE_SET__THE_SAMTOOLS_PATH__;

# DATABASE
# reference genome fasta file
ref_fa=__PLEASE_SET__THE_REF_GENOME_FA_PATH__;
# e.g., 'hg19'
ref_id=__PLEASE_SET__THE_REF_GENOME_ID__;
# gene psl file
# see BioFuse 'get_gpsl' function, https://github.com/Nobel-Justin/BioFuse
# see SOAPfuse blog for PSL format
#     https://sourceforge.net/p/soapfuse/blog/2016/01/new-psl-file-format-applied-by-soapfuse/
gpsl=__PLEASE_SET__THE_GPSL_PATH__;
