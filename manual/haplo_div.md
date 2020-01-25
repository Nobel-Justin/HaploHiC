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

#### sorted VCF file recording the phased mutations. (`-phsvcf` option, **required**)

HaploHiC utilizes phased mutations in this VCF file to distinguish parental origin of Hi-C reads. For instance, the VCF file recording phased germline mutations of NA12878 from the AlleleSeq database (version: Jan-7-2017): [link](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/).

HaploHiC initially loads ONLY phased SNVs from the VCF file. Use option `-use_indel` to enable phased InDels loading. Besides, considering alignment at the mutations site might suffer the affect of adjacent indels, we set a minimum distance (5bp) to avoid loading adjacent indels. User could use option `-min_idd` to reset this minimum distance.

#### file path list of paired SE bam or sortN PE bam. (`-bamlist` option, **required**)

HaploHiC loads Hi-C PE-reads from alignment bam files listed in this file. If one sample has several paried bam files, all of their paths could be written in this bam list (one row one pair). HaploHiC accepts alignment of two main manners: 1) `paired SE` bam, and 2) `sortN PE` bam.

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

#### output folder. (`-outdir` option, **required**)
All results will be stored under this folder (see [instance](#Outputs)). For each sequencing run, related summary files and subfolder will be created with specific prefix. The prefix is the bam file basename after removing the '.bam' postfix. The subfolders will store intermediate files for each sequencing run.

### Tools and Database

#### samtools path. (`-samt` option, **required**)

HaploHiC uses samtools to load input bam files and write results in bam format.

#### database folder. (`-db_dir` option, **required**)

HaploHiC needs the database generated by the [`juicer_db`](./juicer_db.md) command, and automatically loads files under this database folder.

#### reference version. (`-ref_v` option, **required**)

Files are kept in sub-folder under the database folder (`-db_dir`). The sub-folder is named as reference version, such as 'hg19', 'hg38', or 'mm9'. To use database file, user needs to state the reference version via this option.

#### enzyme type. (`-enzyme` option, **required**)

The restriction enzyme used in Hi-C fragment library construction.

### Outputs
Here is an instance, files under the output folder after all five steps (see [next](#Steps)) finished. Files and subfolder with the '**RunID**' prefix just represent one set of results of an sequencing run. If this sample has more data from other sequencing runs, there will be more sets of such results under the output folder. Note that the '**RunID**' is 
```
├── VCF_load.phased_mut.closeInDel.list   ──┬──> generated by 1st step 'splitBam'
├── VCF_load.phased_mut.report            ──┘
├── dumpContacts                  ──> folder stores results of 5th step 'dumpContacts'
├── RunID.merge.h1Intra.bam       ──┐
├── RunID.merge.h2Intra.bam       ──┤
├── RunID.merge.hInter.bam        ──┼──> generated by 4th step 'readsMerge'
├── RunID.PEtoHaplotype.report    ──┤
├── RunID.statOfPhasedLocReg.gz   ──┘
└── RunID-workspace
    ├── RunID.discarded.bam              ──┐
    ├── RunID.invalid.bam                ──┤
    ├── RunID.phMut-dEnd-h1.bam          ──┤
    ├── RunID.phMut-dEnd-h2.bam          ──┤
    ├── RunID.phMut-dEnd-hInter.bam      ──┼──> generated by 1st step 'splitBam'
    ├── RunID.phMut-sEnd-h1.bam          ──┤
    ├── RunID.phMut-sEnd-h2.bam          ──┤
    ├── RunID.phMut-sEnd-hInter.bam      ──┤
    ├── RunID.unknown.bam                ──┘
    ├── RunID.phMut-sEnd-h1.h1Intra.bam             ──┐
    ├── RunID.phMut-sEnd-h1.hInter.bam              ──┤
    ├── RunID.phMut-sEnd-h1.statOfPhasedLocReg.gz   ──┼──> generated by 2nd step 'sEndUKtoHap'
    ├── RunID.phMut-sEnd-h2.h2Intra.bam             ──┤
    ├── RunID.phMut-sEnd-h2.hInter.bam              ──┤
    ├── RunID.phMut-sEnd-h2.statOfPhasedLocReg.gz   ──┘
    ├── RunID.unknown.h1Intra.bam             ──┐
    ├── RunID.unknown.h2Intra.bam             ──┤
    ├── RunID.unknown.hInter.bam              ──┼──> generated by 3rd step 'dEndUKtoHap'
    └── RunID.unknown.statOfPhasedLocReg.gz   ──┘
```

### Steps
`haplo_div` has five steps: `splitBam`, `sEndUKtoHap`, `dEndUKtoHap`, `readsMerge`, and `dumpContacts`. User could set step interval to proceed via `-st_step` and `-ed_step` options.

- splitBam
  Firstly, HaploHiC loads phased mutations from the VCF input, and generates two files.
    - **VCF_load.phased_mut.report**
      basic statistics of phased mutations, including hom/het status, skipped mutations, allelic type of heterozygous mutation, close adjacent InDels, and tandem InDels.
    - **VCF_load.phased_mut.closeInDel.list**
      list of InDels avoided in following analysis due to too close distance (the `-min_idd` option).

  Then, HaploHiC loads bam file of each sequencing run, and distributes Hi-C PE-reads into different categories according to their mutations coverage.
    - *RunID*.discarded.bam
      Hi-C PE-reads filtered out due to `UM`/`MM`/`SD`/`SP`/`LM` depending on options: `-flt_pe`, `-use_sp`, `-use_sd`, `-use_mm`, `-min_mq`.
    - *RunID*.invalid.bam
      invalid Hi-C PE-reads: `Dangling`, `SelfCircle`, `DumpedPair-Forward/Reversed`, or `TooClose` alignment (see `-max_rd` and `-use_caln` options).
    - *RunID*.phMut-dEnd-h1.bam and *RunID*.phMut-dEnd-h2.bam
      Hi-C PE-reads whose two ends simutanously cover heterozygous mutations from same parental origin, which is the ONLY one this PE-reads support. The haplotype NO. ('h1' and 'h2') represent the index of the phased mutations in VCF file. For instance, from a phased site which has 'A|G' alleles, the 'A' allele is on the 'h1' origin, and 'G' allele is on the 'h2' origin.
    - *RunID*.phMut-dEnd-hInter.bam
      Hi-C PE-reads whose one or two ends cover heterozygous mutations from different parental origins.
    - *RunID*.phMut-sEnd-h1.bam and *RunID*.phMut-sEnd-h2.bam
      Hi-C PE-reads whose one end covers heterozygous mutations from ONLY one parental origin, while another end covers no heterozygous mutation (i.e., unknown parental origin, is given with the '*UK*' tag).
    - *RunID*.phMut-sEnd-hInter.bam
      Hi-C PE-reads whose one ends simutanously cover heterozygous mutations from different parental origins, while another end is '*UK*'.
    - *RunID*.unknown.bam
      Hi-C PE-reads whose two ends are '*UK*'.

- sEndUKtoHap
  In this step, HaploHiC deals with *RunID*.phMut-sEnd-h1.bam and *RunID*.phMut-sEnd-h2.bam. It assigns parental origin to the '*UK*' end of PE-reads based on calculation of the `local contacts ratio`.
    - PE-reads in *RunID*.phMut-sEnd-h1.bam are distributed into new files: *RunID*.phMut-sEnd-h1.h1Intra.bam and *RunID*.phMut-sEnd-h1.hInter.bam
    - PE-reads in *RunID*.phMut-sEnd-h2.bam are distributed into new files: *RunID*.phMut-sEnd-h2.h2Intra.bam and *RunID*.phMut-sEnd-h2.hInter.bam
    - The statistics files of `local contacts ratio` are generated: *RunID*.phMut-sEnd-h1.statOfPhasedLocReg.gz and *RunID*.phMut-sEnd-h2.statOfPhasedLocReg.gz

- dEndUKtoHap
  In this step, HaploHiC deals with *RunID*.unknown.bam. It assigns parental origin to the '*UK*' end of each PE-reads based on calculation of the `local contacts ratio`.
    - PE-reads in *RunID*.phMut-sEnd-h1.bam are distributed into three files: *RunID*.unknown.h1Intra.bam, *RunID*.unknown.h2Intra.bam, and *RunID*.unknown.hInter.bam
    - The statistics file of `local contacts ratio` is generated: *RunID*.unknown.statOfPhasedLocReg.gz.

- readsMerge
  In this step, HaploHiC merges 'h[x]Intra' and 'hInter' bams of previous steps, respectively.
    - *RunID*.merge.h1Intra.bam includes PE-reads from *RunID*.phMut-dEnd-h1.bam, *RunID*.phMut-sEnd-h1.h1Intra.bam, and *RunID*.unknown.h1Intra.bam
    - *RunID*.merge.h2Intra.bam includes PE-reads from *RunID*.phMut-dEnd-h2.bam, *RunID*.phMut-sEnd-h2.h2Intra.bam, and *RunID*.unknown.h2Intra.bam
    - *RunID*.merge.hInter.bam includes PE-reads from *RunID*.phMut-dEnd-hInter.bam, *RunID*.phMut-sEnd-h1.hInter.bam, *RunID*.phMut-sEnd-h2.hInter.bam, and *RunID*.unknown.hInter.bam
    - All statistics file of `local contacts ratio` are merged into *RunID*.statOfPhasedLocReg.gz

- dumpContacts
  In this step, HaploHiC provides BP and FRAG mode dump from merged results to summarize raw contacts of pairwise windows.

