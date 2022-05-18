# TB_Pseudo_small_ncRNA
Pseudouridine RNA modifications on small non-coding RNAs in the human pathogen Trypanosoma brucei: Identification and functional implications



An UNIX shell and the R software are needed to use these scripts.

The "Unix_scripts" needs smalt, samtools, perl and Bedtools. You can use your own  alignment tools if you prefer, but be sure that the alignment output in SAM/BAM format. Bedtools is mandatory. These pipeline requires paired-end sequencing data.

**HydraPsiSeq pipeline**

HydraPsiSeq pipeline came from the repository - https://github.com/FlorianPichot/HydraPsiSeqPipeline with minor modifications

General pipeline:
 Alignment (Unix script) | --> | NormUCount (R script) | --> | Analysis with list (R script) |
 
**Ribomethseq pipeline**

Ribomethseq pipeline implements the RiboMethSeq scoring analysis as described in Birkedal et al. 2014
 
