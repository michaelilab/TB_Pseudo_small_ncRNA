#!/usr/bin/env tcsh

set path1="Sequencing_09Jan2022" # path to directory
set type="rRNA"
set path2="Scripts/" #path to Scripts
set snrna="DB/TB_rRNA"
set file1="DB/rRNA.genome"
set fasta="DB/TB_rRNA_chr2.fa"

foreach fastq_file ($path1/RawData/*ME*R1.fastq.gz)
   echo $fastq_file
   set base=`basename $fastq_file|cut -f1,2 -d\_`
   echo $base
   set fastq_file2=`echo $fastq_file |sed -e 's/_R1/_R2/'` ;
   echo $fastq_file2

   smalt map $snrna $fastq_file $fastq_file2 > $base\_vs_$type.sam

   samtools view -f 0x02 -bS  $base\_vs_$type.sam > $base\_vs_$type\_good_pairs.bam
   rm $base\_vs_$type.sam
end

foreach bam_file ($path1/Analysis/HydraPsiSeq_TB/$type/GoodPairs/*_good_pairs.bam)
    echo $bam_file
    set prefix=`basename $bam_file|cut -f1,2 -d\_`
    ##BAM_TO_BED
    bedtools bamtobed -bedpe -i  $prefix\_vs_$type\_good_pairs.bam  | awk '{print $1 "\t" $2  "\t" $6 "\t" $7 "\t" $8 "\t" $9}' | sort -k1,1 -k2,2n > $prefix\_vs_$type\_good_pairs.sorted.bed
    ##GenomeCoverage
    bedtools genomecov -d -g $file1 -i $prefix\_vs_$type\_good_pairs.sorted.bed  >  $prefix\_vs_$type\_good_pairs.sorted.genomecov

    ##init
    cat $prefix\_vs_$type\_good_pairs.sorted.bed  | perl $path2/CountInitiating.pl $file1   > $prefix\_vs_$type\_good_pairs.sorted.init

    cat $prefix\_vs_rRNA_good_pairs.sorted.bed |perl $path2/Count3p.pl $file1 > $prefix\_vs_rRNA_good_pairs.sorted.3p
end


foreach init_file ($path1/Analysis/HydraPsiSeq_TB/$type/GoodPairs/*.genomecov)
    set prefix=`basename $init_file|cut -f1,2 -d\_`
    perl $path2/AddBP2GenomeCov.pl $fasta $init_file > coverage_$prefix\_vs_$type\_good_pairs.txt
end
