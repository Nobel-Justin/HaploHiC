#verion=v1.0
#Wenlong JIA, 2021-09-15

# load variables
source ./setting.sh;

pm_dir=`dirname $HaploHiC`;
pm_dir=$pm_dir/..;
if [ ! -e $pm_dir ] || [ ! -e "$pm_dir/HaploHiC" ]; then
  echo "cannot find perl module folder: $pm_dir";
  exit 1;
fi

# add module
PERL5LIB=$pm_dir:$PERL5LIB; export PERL5LIB;

# path check
for file in HaploHiC db_dir SAMtools tabix; do
  if [ ! -e $(eval echo \$$file) ]; then
    echo "please use correct absolute path for '$file' in setting.sh";
    exit;
  fi;
done

enzyme=MboI;

bam_dir=$PWD/bam/;
result_dir=$PWD/result;
phsvcf=$PWD/database/NA12878/NA12878_diploid_2017_jan7.vcf.bgz;

mkdir -p $result_dir;

$HaploHiC haplo_div \
 -phsvcf $phsvcf \
 -bam $bam_dir/demo_R1.sortN.bam,$bam_dir/demo_R2.sortN.bam \
 -outdir $result_dir \
 -samt $SAMtools \
 -db_dir $db_dir \
 -ref_v hg19 \
 -enzyme $enzyme \
 -use_indel \
 -use_sp \
 -ucrfut 5E4 \
 -fork 10 \
 -dump FRAG:1
