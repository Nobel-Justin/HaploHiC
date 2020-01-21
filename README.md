# HaploHiC

Comprehensive haplotype division of Hi-C PE-reads based on local contacts ratio.

- Author: Wenlong Jia
- Email:  wenlongkxm@gmail.com

## Version
V0.32

## Installation

HaploHiC is written in PERL. All of its functions are packaged into a standalone PERL module. Besides, it also requires other additional modules.

To install HaploHiC, you need to download HaploHiC and add the current directory to the `PERL5LIB` path.
```bash
git clone https://github.com/Nobel-Justin/HaploHiC.git
PERL5LIB=$PERL5LIB:$PWD; export PERL5LIB
```
List of additional PERL modules required:
- [JSON](https://metacpan.org/pod/JSON)
- [SVG](https://metacpan.org/pod/SVG)
- [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager)
- [List::Util](https://metacpan.org/pod/List::Util)
- [BioFuse](https://github.com/Nobel-Justin/BioFuse)

If you encounter problems, please open an issue at the [project on Github](https://github.com/Nobel-Justin/HaploHiC/issues).

## Usage
HaploHiC provides series of commands for Hi-C data analysis.
```bash
Usage:   perl HaploHiC.pl <command>

Commands:

-- Haplotype work
   haplo_div      divide the diploid haplotype (P/M) Hi-C reads from alignment.
   plotLocReg     plot phased local region information.
   adjustDump     adjust HaploHiC dump contacts based on juicer dump contacts.

-- Extensions
   juicer_db      make database files for juicer and HaploHiC operations.
   run_juicer     set up workspace and shell file to run juicer.
   juicerDump     run juicer dump, and merge the inter-chr and intra-chr.
   FRAG_gene      convert FRAG contacts to gene level contacts.
```
