# FRAG_gene
`FRAG_gene` loads contact matrix of **FRAG** level, and outputs the **GENE** level contacts.

## Usage

```
Usage:   perl HaploHiC.pl FRAG_gene <[Options]>

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
   -glsgsd      use genes in gene_list (-glist) in single-side of gene-pair. [disabled]
   -mindis [i]  minimum distance between two genes/regions to calculate contacts. [1000]
                 Note: gene/region pairs having overlap are filtered out.
   -hInter      note the FRAG contacts are inter-haplotype, default is intra-haplotype.

   -h|help      Display this help info.

Version:
   0.09 on 2019-01-31

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```
