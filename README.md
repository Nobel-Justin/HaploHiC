# HaploHiC

Comprehensive haplotype division of Hi-C PE-reads based on local contacts ratio.

- Author: Wenlong Jia
- Email:  wenlongkxm@gmail.com

## Check.pm
### HaploHiC::Check
### VERSION = "0.01"
- check

## PhasedHiC/sEndSoloHapConfirm.pm
### HaploHiC::PhasedHiC::sEndSoloHapConfirm
### VERSION = "0.05"
- sEndhx_get_HapLink
- confirm_sEndSoloHapPE_HapLink
- prepareGetHapBamObj
- getTODOpairBamHrefArray
- startWriteGetHapBam
- assign_sEndUKend_haplotype
- findContactAnchor

## PhasedHiC/dumpContacts.pm
### HaploHiC::PhasedHiC::dumpContacts
### VERSION = "0.02"
- dump_contacts
- get_contacts_idx
- contacts_output

## PhasedHiC/splitPairBam.pm
### HaploHiC::PhasedHiC::splitPairBam
### VERSION = "0.17"
- divide_pairBam
- prepareSplitBamObj
- SourceToSplitBam
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
### VERSION = "0.15"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- DivideHiCreadsToHaplotypes
- check_files
- prepare

## PhasedHiC/phasedPEtoContact.pm
### HaploHiC::PhasedHiC::phasedPEtoContact
### VERSION = "0.08"
- phasePE_to_contactCount
- load_phasedPE_contacts
- mPosToWinIdx
- phasePE_contacts_to_count
- get_rOBpair_HapLinkCount

## PhasedHiC/mergeHaploReads.pm
### HaploHiC::PhasedHiC::mergeHaploReads
### VERSION = "0.03"
- merge_haplo_reads
- mergeReadsOfEachHapComb

## PhasedHiC/dEndUkHapConfirm.pm
### HaploHiC::PhasedHiC::dEndUkHapConfirm
### VERSION = "0.15"
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
### VERSION = "0.51"
- load_variants_dict

## Extensions/ShellForJuicer.pm
### HaploHiC::Extensions::ShellForJuicer
### VERSION = "0.06"
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

## Extensions/JuicerDumpFRAG.pm
### HaploHiC::Extensions::JuicerDumpFRAG
### VERSION = "0.04"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- JuicerDumpFRAGWrok
- get_merged_dump_frag

## Extensions/JuicerDumpBP.pm
### HaploHiC::Extensions::JuicerDumpBP
### VERSION = "0.09"
- return_HELP_INFO
- Load_moduleVar_to_pubVarPool
- Get_Cmd_Options
- para_alert
- JuicerDumpBPWrok
- check_juicer_files
- make_workspace
- load_chr_Things
- get_interChr
- get_intraChr
- convert_intraChr_PosToBinIdx
- merge_inter_intra

## Extensions/FRAG2Gene.pm
### HaploHiC::Extensions::FRAG2Gene
### VERSION = "0.08"
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
- load_enzyme_site_list
- load_gene_psl
- load_gene_list
- check_files

## RunFunc.pm
### HaploHiC::RunFunc
### VERSION = "0.51"
- options_alert_and_run
- func_run
- load_functions
- return_HELP_INFO
- HelpInfo_ExtractOpt
- para_alert

