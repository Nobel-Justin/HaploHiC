# juicer_db
`juicer_db` automatically prepares database files. This database supports operation of both [juicer](https://github.com/aidenlab/juicer) and HaploHiC.

## Usage

```
Usage:   perl HaploHiC.pl juicer_db <[Options]>

Options:

  # run juicer_tools dump func in FRAG mode and merge all chromosomes to one final file #

  # Inputs and Outputs #
   -juicer [s]  path of juicer folder, which has sub-folder 'misc'. <required>
   -db_dir [s]  existing folder to store database files. <required>
   -ref_fa [s]  fasta file of reference genome. <required>
   -ref_v  [s]  version of reference genome, as sub-folder name in 'db_dir'. <required>
   -gpsl   [s]  gene psl file. <required>

  # Tools and Database #
   -bwa    [s]  BWA. <required>
   -samt   [s]  SAMtools. <required>

  # Options #
   -enzyme [s]  enzyme type. <required>
   -skip_rfg    skip reference genome work. [disabled]
   -skip_gan    skip gene annotation work. [disabled]
   -skip_ezy    skip enzyme site work. [disabled]
   -ref_idx     set if the reference genome has been indexed by BWA. [disabled]
                 Note: if set, just creates soft-links to index files.

   -h|help      Display this help info.

Version:
   0.01 on 2018-11-14

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```