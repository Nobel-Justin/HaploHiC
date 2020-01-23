# haplo_div
`haplo_div` loads alignment file(s) of Hi-C data and haplotype information, and assigns haplotypes to Hi-C PE-reads of unknown parental origin according to contacts information (**local contacts ratio**) in local genomic regions with dynamic size. Currently, this functions ONLY works on **diploid**.

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

- sorted VCF file recording the phased mutations. (`-phsvcf` option, **required**)

HaploHiC utilizes phased mutations in this VCF file to distinguish parental origin of Hi-C reads. For instance, the VCF file recording phased germline mutations of NA12878 from the AlleleSeq database (version: Jan-7-2017): [link](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/).

HaploHiC initially loads ONLY phased SNVs from the VCF file. Use option `-use_indel` to enable phased InDels loading. Besides, considering alignment at the mutations site might suffer the affect of adjacent indels, we set a minimum distance (5bp) to avoid loading adjacent indels. User could use option `-min_idd` to reset this minimum distance.

- file path list of paired SE bam or sortN PE bam. (`-bamlist` option, **required**)

HaploHiC loads Hi-C PE-reads from alignment bam files listed in this file. If one sample has several paried bam files, all of their paths could be written in this bam list (one row one pair). HaploHiC accepts alignment of two main manners: 1) `paired SE` bam, and 2) `sortN PE` bam.\\
In the `paired SE` manner, paired Hi-C reads are aligned seperatedly and stored into own bam files.  For instance, if one sample has three runs (a,b,c) of paired-end Hi-C sequencing data in the `paired SE` manner, its bam list could be written as below:
```
run_a_R1.bam,run_a_R2.bam
run_b_R1.bam,run_b_R2.bam
run_c_R1.bam,run_c_R2.bam
```
In the `sortN PE` manner, paired Hi-C reads are aligned together and stored into one single bam files. Note that the bam files are all sorted by read-name (`samtools sort -n`). For instance, if one sample has four runs (m,n,i,k) of paired-end Hi-C sequencing data aligned in the `sortN PE` manner, its bam list could be written as below:
```
run_m.sortN.bam
run_n.sortN.bam
run_i.sortN.bam
run_k.sortN.bam
```
We suggest to use GATK to [realign](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md) the bam files with the phased VCF (`-phsvcf`). In the realignment process, user needs to firstly sort the bam file by coordinates (`samtools sort`), then do GATK realign, and finally sort the realigned bam file by read-name (`samtools sort -n`).

### Tools and Databas

- samtools path. (`-samt` option, **required**)

HaploHiC uses samtools to load input bam files and write results in bam format.

- database folder. (`-db_dir` option, **required**)

HaploHiC needs the database generated by the [`juicer_db`](./juicer_db.md) command, and automatically loads files under this database folder.

- reference version. (`-ref_v` option, **required**)

Files are kept in sub-folder under the database folder (`-db_dir`). The sub-folder is named as reference version, such as 'hg19', 'hg38', or 'mm9'. To use database file, user needs to state the reference version via this option.

### Steps
`haplo_div` has five steps: `splitBam`, `sEndUKtoHap`, `dEndUKtoHap`, `readsMerge`, and `dumpContacts`. User could set step interval to proceed via `-st_step` and `-ed_step` options.

- `splitBam` distributes Hi-C PE-reads into nine categories.

- `sEndUKtoHap` deals with only 'sEnd-h[x]' category, and assigns haplotype to 'UK' end in 'sEnd-h[x]' PE-reads based on local contacts ratio.

- `dEndUKtoHap` deals with only 'unknown' category, and assigns dual-side unknown PE to to certain haplotype-contacts based on local contacts ratio.

- `readsMerge`  merges 'h[x]Intra' and 'hInter' of previous steps, respectively.

- `dumpContacts` provides BP and FRAG mode summary of pairwise windows raw contacts.

### Outputs
