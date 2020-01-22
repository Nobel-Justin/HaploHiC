# juicerDump
`juicerDump` runs juicer `dump` on each pairwise inter/intra-chr chromosomes, and merge the results.

## Usage

```bash
Usage:   perl HaploHiC.pl juicerDump <[Options]>

Options:

  # run juicer_tools dump func and merge all chromosomes to one final file #

  # Inputs and Outputs #
   -juicer [s]  path of juicer folder, which has sub-folder, e.g., 'CPU'. <required>
   -topdir [s]  topdir stores juicer results. <required>
                   this func creates 'extract_matrix-<postfix>' folder under topdir.
   -db_dir [s]  database folder made by 'juicer_db' function. <required>
   -ref_v  [s]  version of reference genome, see contents under 'db_dir'. <required>

  # Options #
   -hic_fl [s]  prefix of '.hic' file: 'inter_30' or 'inter'. [inter_30]
   -dpmode [s]  mode of dump, 'BP' or 'FRAG'. [BP]
   -enzyme [s]  enzyme type. <required in FRAG mode>
   -bin    [s]  bin size. [1MB]
                  1)   BP mode: 2.5MB, 1MB, 500KB, 250KB, 100KB, 50KB, 25KB, 10KB, 5KB
                  2) FRAG mode: 500, 200, 100, 50, 20, 5, 2, 1
   -ctype  [s]  'observed' or 'oe'. [observed]
   -norm   [s]  kind of normalization, must be one of NONE/VC/VC_SQRT/KR. ['NONE']
                  Notes from juicer_tools dump usage:
                   1) VC is vanilla coverage,
                   2) VC_SQRT is square root of vanilla coverage;
                   3) KR is Knight-Ruiz or Balanced normalization.

   -h|help      Display this help info.

Version:
   0.11 on 2019-03-23

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```
