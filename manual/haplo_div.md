# haplo_div
`haplo_div` loads alignment file(s) of Hi-C data and haplotype information, and assigns haplotypes to Hi-C PE-reads of unknown parental origin according to contacts information (**local contacts ratio**) in local genomic regions with dynamic size.

## Usage

```
Usage:   perl HaploHiC.pl haplo_div <[Options]>

Options:

  # assign Hi-C reads to paternal and maternal haplotypes via Phased-Het-Mutation (PhaHetMut).  #
  # allocate haplotype-unknown Hi-C PE reads (HaploUKPE) based on Local Hap-Divide-Ratio (LHDR) #

  # Inputs and Outputs #
   -phsvcf  [s]  phased sorted VCF file. <required>
   -outdir  [s]  folder to store results and temporary files. <required>
   -bamlist [s]  file path list of paired SE bam or sortN PE bam. <required>
                  Note: 1) PE Hi-C reads are always aligned seperatedly.
                        2) list multiple runs, one row one run.
                        3) instance run_c: run_c_R1.bam,run_c_R2.bam
                        4) instance run_k: run_k.sortN.bam

  # Tools and Database #
   -samt    [s]  SAMtools. <required>
   -db_dir  [s]  database folder made by 'juicer_db' function. <required>
   -ref_v   [s]  version of reference genome, see contents under 'db_dir'. <required>

  # General Options #
   -enzyme  [s]  enzyme type. <required>
   -hapct   [i]  N-ploid of the sample. [2]
                  Note: default 2 is for homo sapiens.
   -use_indel    use phased InDels. [disabled]
   -min_idd [i]  minimum distance to accept adjacent indels. [5]
   -st_step [i]  start wrokflow from certain step. [1]
   -ed_step [i]  stop workflow at certain step. [5]
                  Note: step-NO. list:  1: splitBam;   2: sEndUKtoHap;  3: dEndUKtoHap;
                                        4: readsMerge; 5: dumpContacts

  # Options of step NO.1 #
   -flt_pe       discard PE-reads once either read has unqualitified alignment. [disabled]
                  Note: default to discard unqualitified alignment only.
   -use_sp       use PE-reads having supplementary alignments, suggested to enable. [disabled]
   -use_sd       use PE-reads having secondary alignments. [disabled]
   -use_mm       use multiple mapped reads. [disabled]
   -min_mq  [i]  the minimum mapping quality. [30]
   -min_re  [i]  the minimum distance to bilateral reads-edges to judge one het-allele. [5]
   -min_bq  [i]  the minimum base quality to accept mut-allele on one reads. [20]
   -min_ad  [i]  the minimum distance to nearby alteration on reads to accept one het-allele. [3]
   -qual    [i]  base quality offset. [33]
                  Note: '33' for Sanger and Illumina 1.8+; '64' for Solexa, Illumina 1.3+ and 1.5+.
   -max_rd  [s]  the maximum distance to determine close alignments. [1E3]
   -use_caln     accept PE-reads whose two ends are close alignments ('-max_rd'). [disabled]
                  Note: default to put them in 'invalid PE'.

  # Options of step NO.2-3 #
   -ucrfut  [s]  the unit size of unilateral local flanking region, minimum: 5E3. [1E4]
   -ucrfmx  [s]  the maximum size of unilateral local flanking region when extends. [1E7]
   -mpwr    [f]  ratio of '-ucrfut' to set as window size to store phased contacts. (0, [0.5]]
                  Note: the less '-mpwr' is, the more memory and cpu-time consumed.
   -min_ct  [i]  the minimum contact count from phased pairs to confirm phased local regions. >=[5]

  # Options of step NO.5 #
   -dump    [s]  how to dump. [BP:1MB]
                  Note: 1) format is DumpMode[:BinSize:HapComb], e.g., BP, FRAG:1, BP:500KB:h2Intra
                        2) can be used multiple times, e.g., -dump BP:500KB -dump FRAG:1 -dump BP
                        3) DumpMode: BP or FRAG.
                        4) BP   BinSize: 2.5MB, [1MB], 500KB, 250KB, 100KB, 50KB, 25KB, 10KB, 5KB
                        5) FRAG BinSize: 500, 200, 100, 50, 20, 5, 2, [1]
                        6) HapComb is to select any haplo combination matched this string. default all.
                           e.g., use 'Intra' to get intra-haplotype contacts of all chromosomes.
                                 use 'Inter' to get inter-haplotype contacts of all chromosomes.
                                 use 'h2Intra' to get only intra-h2 contacts of all chromosomes.

   -h|help       Display this help info.

Version:
   0.32 on 2019-04-05

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```

### Inputs

### Steps









