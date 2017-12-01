#!/bin/bash -e 

# this script takes an output sam file from a lastz alignment to a target assembly, convert sam -> bam -> bed
# it creates a browser track from the bed, and spits a count of mapped read and how many were uniquely mapped
# it also cleans all the duplicated mappings from the bed file - for this, your thirds input should be some identifier 
# (preferably the sra ID)

if [ "$#" -ne 2 ]; then
	echo -e "\nLastZ mapping output processing\n  Usage:\n  $0 output_sam_file in_id"
	exit 1
fi

in_sam=$1
in_id=$2

base=`echo $in_sam | sed -e "s/.sam//g"`

echo -e "\n\t...converting: sam >> bam >> bed"
samtools view -bS ${base}.sam > ${base}.bam
bedtools bamtobed -bed12 -i ${base}.bam > ${base}.bed

echo -e "\n\t...Mapped read Summary:"
n_total_mapped=`cut -f4 ${base}.bed | sort -u | wc -l`
echo -e "\n\t\t...${n_total_mapped} total reads were mapped"
n_unique_mapped=`cut -f4 ${base}.bed| sort | uniq -c | awk '{if($1==1) print $0}' | wc -l`
echo -e "\n\t\t...${n_unique_mapped} reads were uniquely mapped"

echo -e "\n\t...Removing duplicated reads from ${base}.bed"
join -t$'\t' -1 1 -2 4 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10 2.11 2.12' \
	<(cut -f4 ${base}.bed | sort | uniq -c | awk '{if($1==1) print $2}' | sort -u) \
	<(sort -t$'\t' -k4,4 ${base}.bed) | bedSort stdin tmp.bed

cat <(echo -e "track name=${in_id}_lastz_mapped_reads_${base} description=\"${in_id} lastz mapped reads\" color=70,70,70") <(cat tmp.bed) > ${base}.bed

echo -e "\n\n\t...bed file (and browser track) ${base}.bed is ready\n\n"

rm -rf tmp.bed




