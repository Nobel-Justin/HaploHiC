This demo case is from the simulation data in our submitted paper.
More cases could be found in the paper.

The genome reference is GRCh37.
The enzyme is MboI.

< Content from dropbox >

The alignment and database files are too large to store in GitHub.
Please download the 'bam' and 'database' folders from dropbox sharing:
  https://www.dropbox.com/sh/g2rk6sjud8scuqk/AADRAqUoRF66xJDTgN_IB2p5a?dl=0

Save these two folders in the current folder.

1. bam folder
   a) The alignemnt of sequencing reads of the simulated contacts.
   b) alignment reference is GRCh37.
   c) bam files are sorted by reads-id.

2. database folder of NA12878 sample
   a) VCF file of phased mutations.
   b) from the AlleleSeq database (version: Jan-7-2017).
      http://alleleseq.gersteinlab.org/downloads.html


< Operation >

1. set variables in the 'setting.sh' file
2. run the 'run.sh' file in this folder
    bash run.sh 1>olog 2>elog


Wenlong Jia
wenlongkxm@gmail.com
2021-09-15
