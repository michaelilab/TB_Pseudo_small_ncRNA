#!/usr/bin/env tcsh

#Directory Structure:
#Main_Directory/RawData # raw fastq files
#Main_Directory/Analysis/RiboMethSeq_TB/rRNA #results files

set path1="Main_Directory/" #change to correct path
set path2="Scripts/"
set snrna="DB/TB_rRNA"
set file1="DB/rRNA.genome"

####################################################################################

foreach fastq_file ($path1/RawData/*_ME_R1.fastq.gz)
   echo $fastq_file
   set base=`basename $fastq_file|cut -f1,2 -d\_`
   echo $base
   set fastq_file2=`echo $fastq_file |sed -e 's/_R1/_R2/'` ;
   echo $fastq_file2

   smalt map $snrna $fastq_file $fastq_file2 > $base\_vs_rRNA.sam

   samtools view -f 0x02 -bS  $base\_vs_rRNA.sam > $base\_vs_rRNA_good_pairs.bam
   rm $base\_vs_rRNA.sam
end
###################################################################################
foreach bam_file ($path1/Analysis/RiboMethSeq_TB/rRNA/*.bam)
    echo $bam_file
    set prefix=`basename $bam_file|cut -f1,2 -d\_`
    ##BAM_TO_BED
    bedtools bamtobed -bedpe -i  $prefix\_vs_rRNA_good_pairs.bam  | awk '{print $1 "\t" $2  "\t" $6 "\t" $7 "\t" $8 "\t" $9}' | sort -k1,1 -k2,2n > $prefix\_vs_rRNA_good_pairs.sorted.bed
    ##GenomeCoverage
    bedtools genomecov -d -g $file1 -i $prefix\_vs_rRNA_good_pairs.sorted.bed  >  $prefix\_vs_rRNA_good_pairs.sorted.genomecov

    ##init
    cat $prefix\_vs_rRNA_good_pairs.sorted.bed  | perl $path2/CountInitiating.pl $file1   > $prefix\_vs_rRNA_good_pairs.sorted.init

    cat $prefix\_vs_rRNA_good_pairs.sorted.bed |perl $path2/Count3p.pl $file1 > $prefix\_vs_rRNA_good_pairs.sorted.3p
end

