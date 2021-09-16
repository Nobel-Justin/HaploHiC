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
for file in HaploHiC Juicer_folder SAMtools BWA ref_fa gpsl; do
  if [ ! -e $(eval echo \$$file) ]; then
    echo "please use correct absolute path for '$file' in setting.sh";
    exit;
  fi;
done

db_dir=$PWD/db;

mkdir -p $db_dir;

$HaploHiC juicer_db \
 -juicer $Juicer_folder \
 -db_dir $db_dir \
 -ref_fa $ref_fa \
 -ref_v $ref_id \
 -bwa $BWA \
 -samt $SAMtools \
 -enzyme $enzyme \
 -gpsl $gpsl
