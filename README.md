# HaploHiC

Comprehensive haplotype division of Hi-C PE-reads based on local contacts ratio.

- Author: Wenlong Jia
- Email:  wenlongkxm@gmail.com

## Check.pm
### HaploHiC::Check
### VERSION = "0.02"
- check

## PhasedHiC/sEndSoloHapConfirm.pm
### HaploHiC::PhasedHiC::sEndSoloHapConfirm
### VERSION = "0.06"
- sEndhx_get_HapLink
- confirm_sEndSoloHapPE_HapLink
- prepareGetHapBamObj
- getTODOpairBamHrefArray
- startWriteGetHapBam
- assign_sEndUKend_haplotype
- findContactAnchor

## PhasedHiC/dumpContacts.pm
### HaploHiC::PhasedHiC::dumpContacts
### VERSION = "0.04"
- dump_contacts
- load_enzyme_site_list
- get_contacts_idx
- get_dumpHeader
- write_dumpBinLog
- contacts_output

## PhasedHiC/splitPairBam.pm
### HaploHiC::PhasedHiC::splitPairBam
### VERSION = "0.18"
- divide_pairBam
- prepareSplitBamObj
- sourceToSplitBam
- startWriteSplitBam
- loadrOBfromSourceBam
- rOB_to_peOB
- PEreads_to_splitBam
- judge_on_rEnd
- judge_reads_alignment
- judge_reads_by_phMut
- judge_on_PE
- selectOneFromMultiHap
- sortSplitBamByPEcontact
- write_peOB_to_chrPairBam
- chrPairBamToSortSplitBam
- write_sort_peOB_to_newSplitBam

## PhasedHiC/HaploDivideReads.pm
### HaploHiC::PhasedHiC::HaploDivideReads
### VERSION = "0.16"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- DivideHiCreadsToHaplotypes
- check_files
- prepare

## PhasedHiC/phasedPEtoContact.pm
### HaploHiC::PhasedHiC::phasedPEtoContact
### VERSION = "0.10"
- phasePE_to_contactCount
- load_phasedPE_contacts
- mPosToWinIdx
- phasePE_contacts_to_count
- get_rOBpair_HapLinkCount

## PhasedHiC/adjustDumpContacts.pm
### HaploHiC::PhasedHiC::adjustDumpContacts
### VERSION = "0.01"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- adjustDumpContacts
- check_header
- prepare_fh
- compare_chrPair
- adjustHaploHiCdump
- loadJuicerDPtoAdjust
- dealJuicerDPinfo
- dealHaploHiCRemainDP
- dealJuicerRemainDP

## PhasedHiC/mergeHaploReads.pm
### HaploHiC::PhasedHiC::mergeHaploReads
### VERSION = "0.03"
- merge_haplo_reads
- mergeReadsOfEachHapComb

## PhasedHiC/dEndUkHapConfirm.pm
### HaploHiC::PhasedHiC::dEndUkHapConfirm
### VERSION = "0.16"
- dEndUK_get_HapLink
- confirm_dEndUkHapPE_HapLink
- assign_dEndUKend_haplotype

## PhasedHiC/phasedMutWork.pm
### HaploHiC::PhasedHiC::phasedMutWork
### VERSION = "0.14"
- load_phased_VCF
- GetPhaseMutEdgeDist
- release_phaseMut_OB
- filterClosePhasedInDel
- PosToPhasedMut

## GetPath.pm
### HaploHiC::GetPath
### VERSION = "0.01"
- GetPath

## HaploHiC.pm
### HaploHiC::HaploHiC

## LoadOn.pm
### HaploHiC::LoadOn
### VERSION = "0.50"
- load_variants_dict

## Extensions/JuicerDump.pm
### HaploHiC::Extensions::JuicerDump
### VERSION = "0.10"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- JuicerDumpWrok
- check_juicer_files
- make_workspace
- load_chr_Things
- getPairChrDump
- mergePairChrDump

## Extensions/ShellForJuicer.pm
### HaploHiC::Extensions::ShellForJuicer
### VERSION = "0.07"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- shell_to_run_juicer
- check_files_parm
- write_shell
- submit_to_flux

## Extensions/JuicerDB.pm
### HaploHiC::Extensions::JuicerDB
### VERSION = "0.01"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- JuicerMakeDBwork
- ref_genome_part
- gene_anno_part
- enzyme_site_part

## Extensions/FRAG2Gene.pm
### HaploHiC::Extensions::FRAG2Gene
### VERSION = "0.09"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- DumpFRAGtoGeneLevel
- report_GeneOrUserReg_FRAGidx
- FRAG2gene_contact
- convert_frag_to_gene_contacts
- load_region_list
- load_gene_info
- load_gene_psl
- load_gene_list
- check_files

## RunFunc.pm
### HaploHiC::RunFunc
### VERSION = "0.53"
- options_alert_and_run
- func_run
- load_functions
- return_HELP_INFO
- HelpInfo_ExtractOpt
- para_alert

