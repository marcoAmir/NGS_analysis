#!/bin/bash -e

# Amir Marcovitz
#
# A wrapper for indexing a genome+gtf for RNA-seq analysis with STAR
#
# Inputs:	genome assembly for which 2bit is availale (e.g., mm10)
#		relevant genome annotation file (gtf) from Ensembl
#
# for some assemblies (e.g., mm10, rn6, bosTau8, hg38) the script will download gtf
# from ensembl 
# before indexing, the script will make sure that chromosome names are consistent
# while removing chromosomes that have 'random' in their name
#
# requirements: 
#	kent source utilities	(https://github.com/ENCODE-DCC/kentUtils)
#	STAR v2.4 and up 	(https://github.com/alexdobin/STAR)

# defaults:
default_read_length=100								# modify if your RNA-seq library is different
default_avail_cores=8								# make sure it fites your system

path_to_genomes=/cluster/gbdb							# path to dir storing genome 2bits
genomeDir=/cluster/u/amirma/geneLoss/hg38/validations/rnaseq/indexed_genomes	# output dir for indexed genomes

if [ "$#" -eq 1 ]; then
	assembly=$1
	if [ ${assembly} = "mm10" ]; then
		gtf_file="ftp://ftp.ensembl.org/pub/release-86/gtf/mus_musculus/Mus_musculus.GRCm38.86.gtf.gz"
	elif [ ${assembly} = "rn6" ]; then
		gtf_file="ftp://ftp.ensembl.org/pub/release-86/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.86.gtf.gz"
	elif [ ${assembly} = "bosTau8" ]; then
		gtf_file="ftp://ftp.ensembl.org/pub/release-86/gtf/bos_taurus/Bos_taurus.UMD3.1.86.gtf.gz"
	elif [ ${assembly} = "hg38" ]; then
		gtf_file="ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz"
	else
		echo -e "\n\tError! genome-annotation file (gtf) or ftp path to gtf not found\n\n"
		exit 1
	fi
elif [ "$#" -eq 2 ]; then
	assembly=$1
	gtf_file=$2
else
	echo -e "\nA wrapper for indexing a genome+gtf for RNA-seq analysis with STAR\n  
  Usage:\n\t$0 assembly gtf_file\n  for some assemblies (e.g., mm10, rn6, bosTau8, hg38) gtf downloaded automatically\n"
	exit 1
fi

if [ ! -d "${genomeDir}/${assembly}" ]; then 
  mkdir ${genomeDir}/${assembly}
fi

# get the genome file (2bit->fa):
echo -e "\n\t...genome to fasta file"
twoBitToFa ${path_to_genomes}/${assembly}/${assembly}.2bit ${assembly}.tmp.fa

# filter out random chromosomes:
echo -e "\n\t...removing random chromosomes"
cat ${assembly}.tmp.fa | grep ">" | sed -e "s/>//g" | grep -v "random" > chroms.txt
faFilter -namePatList=chroms.txt ${assembly}.tmp.fa ${assembly}.fa
rm -rf ${assembly}.tmp.fa


# get gtf-file and process it:
if ls *.gtf*; then
	echo -e "\n\t...processing gtf file: ${gtf_file}"
else
	echo -e "\n\t...downloading gtf file: ${gtf_file}"
	wget ${gtf_file}
	gtf_file=`echo ${gtf_file} | awk -F'/' '{print $NF}'`
fi
echo -e "\n\t...making genome and and gtf chromosome names consistent"
zcat ${gtf_file} | egrep -v "^#" | awk -F'\t' \
'{if($1<=25 || $1=="X" || $1=="Y" || $1=="MT") {print "chr"$0} else {split($1,a,"."); \
	print a[1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}' | sed -e "s/chrMT/chrM/g" > tmp
comm -1 -2 <(cut -f1 tmp | sort -u) <(sed -e "s/chrUn_//g" chroms.txt | sort -u) | sort -u > chroms.gtf.txt
join -t$'\t' -1 1 -2 1 -o '1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9' <(sort -u tmp) <(sort -u chroms.gtf.txt) \
	| awk -F'\t' '{if($1~/GL/ || $1~/JH/) {print "chrUn_"$0} else {print $0}}' > ${assembly}.gtf
rm -rf tmp chroms.* ${gtf_file}
gtf_file=${assembly}.gtf
n_chroms=`cut -f1 ${gtf_file} | sort -u | wc -l`
echo -e "\t\t${n_chroms} references (chromosomes) in genome file ${assembly}.fa and annotation file ${gtf_file}\n\n"

# indexing genome+gtf with STAR
echo -e "\n\t...indexing genome+gtf with STAR\n"
if [ ${m_chroms} -gt 5000 ]; then
STAR --runThreadN ${default_avail_cores} --runMode genomeGenerate --genomeDir ${genomeDir}/${assembly} \
	--genomeFastaFiles ${assembly}.fa --sjdbGTFfile ${gtf_file} --sjdbOverhang ${default_read_length}
else 
STAR --runThreadN ${default_avail_cores} --runMode genomeGenerate --genomeDir ${genomeDir}/${assembly} \
	--genomeFastaFiles ${assembly}.fa --sjdbGTFfile ${gtf_file} --sjdbOverhang ${default_read_length} \
	--genomeChrBinNbits 18
fi
echo -e "\n\n"

rm -rf ${gtf_file} ${assembly}.fa









