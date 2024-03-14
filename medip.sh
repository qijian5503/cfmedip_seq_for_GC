ls raw_data/ |grep 1.fq.gz |sed s/"\_1.fq.gz"//g >name.txt
for i in $(cat name.txt);do
trim_galore -j 30 -q 20 --phred33 --stringency 1 --illumina  --length 20 -e 0.1 --paired raw_data/$i\_1.fq.gz raw_data/$i\_2.fq.gz --gzip -o bam_data/
bowtie2 -p 30 -x bowtie2_index/hg19_made/hg19 -1 bam_data/$i\_1_val_1.fq.gz -2 bam_data/$i\_2_val_2.fq.gz | samtools sort -O bam -@ 30 -o bam_data/$i\.bam
samtools sort -@ 30 -n -O bam -o bam_data/$i\_sorted.bam bam_data/$i\.bam
samtools fixmate -m -@ 30 -O bam bam_data/$i\_sorted.bam bam_data/$i\_sorted_fixmate.bam
samtools sort -@ 30 -O bam bam_data/$i\_sorted_fixmate.bam -o bam_data/$i\_mc_sorted.bam
samtools markdup -r -@ 30 bam_data/$i\_mc_sorted.bam bam_data/$i\_samtools_rmdup.bam
samtools sort -@ 30 bam_data/$i\_samtools_rmdup.bam -o bam_data/$i\_samtools_rmdup_sorted.bam
samtools index bam_data/$i\_samtools_rmdup_sorted.bam
rm bam_data/$i\_sorted.bam bam_data/$i\_samtools_rmdup.bam bam_data/$i\.bam bam_data/$i\_sorted_fixmate.bam bam_data/$i\_mc_sorted.bam

macs2 callpeak -t bam_data/$i\_samtools_rmdup_sorted.bam -m 10 30 -p 1e-5 -f BAM -g hs -n $i 2>$i\.masc2.log
grep -v \# $i\_peaks.xls  |sed '1,2d' >  $i\_peaks_xls.bed

done


ls *_peaks_xls.bed | xargs -i awk -v a={} 'BEGIN{OFS="\t";FS="\t"}{gsub("_peaks.xls.bed","",a);print $1,$2,$3,$5,$6,a}' {}  >> peaks_xls.txt

bedtools sort -i  peaks_xls.txt > all_peak_sorted.bed
bedtools merge -i all_peak_sorted.bed > all_peak_merged.bed
bedtools intersect -a all_peak_merged.bed -b $i\_peaks_xls.bed -wao > peaks_exp.bed


## new sample

bedtools sort -i  new_sample_peaks_xls.txt > new_sample_peak_sorted.bed
bedtools intersect -a all_peak_merged.bed -b new_sample_peak_sorted.bed -wao > new_sample_peaks_exp.bed



