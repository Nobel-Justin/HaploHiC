# run_juicer
`run_juicer` automatically sets up workspace and shell file for juicer operation.

## Usage

```bash
Usage:   perl HaploHiC.pl run_juicer <[Options]>

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
   0.07 on 2019-01-31

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```
