#!/bin/bash

cd /home/lawrence
project="project"
mkdir -p ${project}
cd ${project}

ts=$(date +"%Y-%m-%d_%H-%M-%S")

project_dir=${PWD}
raw_reads_dir=${project_dir}/miRNA_RAW # miRNA fastq files need to be loaded into this file when the shell is running, this avoids issues with the root directory and permissions, do this as the genome is downloading
chr_dir="${project_dir}/genome/chromosomes"
genome_dir="${project_dir}/genome/whole_genome"
mirna_dir="${project_dir}/mirna"
genome_index_dir="${project_dir}/genome_index"
clip_reads_dir="${project_dir}/reads/clip"
raw_FastQC_dir="${project_dir}/fastqc/raw"
clip_fastq_dir="${project_dir}/fastqc/clip"
mirdeep_dir="${project_dir}/mirdeep"
mirdeep_oneW="${mirdeep_dir}/runs/oneW"
mirdeep_oneM="${mirdeep_dir}/runs/oneM"
mirdeep_threeM="${mirdeep_dir}/runs/threeM"
mirdeep_AD="${mirdeep_dir}/runs/AD"
mapper_dir="${mirdeep_dir}/mapper"

segregated_dir="${mirdeep_dir}/segregated"
one_week="${mirdeep_dir}/segregated/one_week"
one_month="${mirdeep_dir}/segregated/one_month"
three_month="${mirdeep_dir}/segregated/three_month"
AD="${mirdeep_dir}/segregated/AD"
standardise="${segregated_dir}/standardise"
AD_reads="${mapper_dir}/AD_reads"
oneW_reads="${mapper_dir}/oneW_reads"
oneM_reads="${mapper_dir}/oneM_reads"
threeM_reads="${mapper_dir}/threeM_reads"
AD_mapping="${mapper_dir}/AD_mapping"
oneW_mapping="${mapper_dir}/oneW_mapping"
oneM_mapping="${mapper_dir}/oneM_mapping"
threeM_mapping="${mapper_dir}/threeM_mapping"
mapper_output_1W="${mapper_dir}/mapper_output_1W"
mapper_output_1M="${mapper_dir}/mapper_output_1M"
mapper_output_3M="${mapper_dir}/mapper_output_3M"
mapper_output_AD="${mapper_dir}/mapper_output_AD"

#make appropriate directories
mkdir -p ${chr_dir}; mkdir -p ${genome_dir}; mkdir -p ${mirna_dir}; mkdir -p ${genome_index_dir}; mkdir -p ${mirdeep_dir}; mkdir -p ${clip_reads_dir}; mkdir -p ${raw_FastQC_dir}; 
mkdir -p ${clip_fastq_dir}; mkdir -p ${mapper_dir}; mkdir -p ${miRNA_RAW}; mkdir -p ${mirdeep_1W}; mkdir -p ${mirdeep_1M}; mkdir -p ${mirdeep_3M}; mkdir -p ${mirdeep_AD}; 
mkdir -p ${segregated_dir}; mkdir -p ${one_week}; mkdir -p ${one_month}; mkdir -p ${three_month}; mkdir -p ${AD}; mkdir -p ${standardise}; mkdir -p ${AD_reads}; mkdir -p ${oneW_reads};
mkdir -p ${oneM_reads}; mkdir -p ${threeM_reads}; mkdir -p ${AD_mapping}; mkdir -p ${oneW_mapping}; mkdir -p ${oneM_mapping}; mkdir -p ${threeM_mapping}; mkdir -p ${mapper_output_1W} ; 
mkdir -p ${mapper_output_1M}; mkdir -p ${mapper_output_3M}; mkdir -p ${mapper_output_AD}      

# Download sheep genome, unzip and edit header line for each chromosomal file, note the oar v1 is not the first version but is the fifth version, updated in Jan 2019, however it is the first version of the 'rambouillet' genome... groundbreaking stuff ladies and gentlemen, this is why we have so many issues with bioinformatics... ego and personal oddities
# Note the unplaced scaffolds have multiple FASTA headers so here we search and replace all header names (not just line 1) and number them sequentially so that each header is uniquely named.
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_Un/oar_ref_Oar_rambouillet_v1.0_chrUn.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chrUn.fa && awk '/^>/{$0=">chrUn_"(++i)}1' ${chr_dir}/oar_v1-0_chrUn.fa > ${chr_dir}/oar_v1-0_chrUn_temp.fa
mv ${chr_dir}/oar_v1-0_chrUn_temp.fa ${chr_dir}/oar_v1-0_chrUn.fa

wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_01/oar_ref_Oar_rambouillet_v1.0_chr1.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr1.fa && sed -i "1 s/^.*$/>chr01/" ${chr_dir}/oar_v1-0_chr1.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_02/oar_ref_Oar_rambouillet_v1.0_chr2.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr2.fa && sed -i "1 s/^.*$/>chr02/" ${chr_dir}/oar_v1-0_chr2.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_03/oar_ref_Oar_rambouillet_v1.0_chr3.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr3.fa && sed -i "1 s/^.*$/>chr03/" ${chr_dir}/oar_v1-0_chr3.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_04/oar_ref_Oar_rambouillet_v1.0_chr4.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr4.fa && sed -i "1 s/^.*$/>chr04/" ${chr_dir}/oar_v1-0_chr4.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_05/oar_ref_Oar_rambouillet_v1.0_chr5.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr5.fa && sed -i "1 s/^.*$/>chr05/" ${chr_dir}/oar_v1-0_chr5.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_06/oar_ref_Oar_rambouillet_v1.0_chr6.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr6.fa && sed -i "1 s/^.*$/>chr06/" ${chr_dir}/oar_v1-0_chr6.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_07/oar_ref_Oar_rambouillet_v1.0_chr7.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr7.fa && sed -i "1 s/^.*$/>chr07/" ${chr_dir}/oar_v1-0_chr7.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_08/oar_ref_Oar_rambouillet_v1.0_chr8.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr8.fa && sed -i "1 s/^.*$/>chr08/" ${chr_dir}/oar_v1-0_chr8.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_09/oar_ref_Oar_rambouillet_v1.0_chr9.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr9.fa && sed -i "1 s/^.*$/>chr09/" ${chr_dir}/oar_v1-0_chr9.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_10/oar_ref_Oar_rambouillet_v1.0_chr10.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr10.fa && sed -i "1 s/^.*$/>chr10/" ${chr_dir}/oar_v1-0_chr10.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_11/oar_ref_Oar_rambouillet_v1.0_chr11.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr11.fa && sed -i "1 s/^.*$/>chr11/" ${chr_dir}/oar_v1-0_chr11.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_12/oar_ref_Oar_rambouillet_v1.0_chr12.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr12.fa && sed -i "1 s/^.*$/>chr12/" ${chr_dir}/oar_v1-0_chr12.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_13/oar_ref_Oar_rambouillet_v1.0_chr13.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr13.fa && sed -i "1 s/^.*$/>chr13/" ${chr_dir}/oar_v1-0_chr13.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_14/oar_ref_Oar_rambouillet_v1.0_chr14.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr14.fa && sed -i "1 s/^.*$/>chr14/" ${chr_dir}/oar_v1-0_chr14.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_15/oar_ref_Oar_rambouillet_v1.0_chr15.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr15.fa && sed -i "1 s/^.*$/>chr15/" ${chr_dir}/oar_v1-0_chr15.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_16/oar_ref_Oar_rambouillet_v1.0_chr16.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr16.fa && sed -i "1 s/^.*$/>chr16/" ${chr_dir}/oar_v1-0_chr16.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_17/oar_ref_Oar_rambouillet_v1.0_chr17.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr17.fa && sed -i "1 s/^.*$/>chr17/" ${chr_dir}/oar_v1-0_chr17.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_18/oar_ref_Oar_rambouillet_v1.0_chr18.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr18.fa && sed -i "1 s/^.*$/>chr18/" ${chr_dir}/oar_v1-0_chr18.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_19/oar_ref_Oar_rambouillet_v1.0_chr19.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr19.fa && sed -i "1 s/^.*$/>chr19/" ${chr_dir}/oar_v1-0_chr19.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_20/oar_ref_Oar_rambouillet_v1.0_chr20.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr20.fa && sed -i "1 s/^.*$/>chr20/" ${chr_dir}/oar_v1-0_chr20.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_21/oar_ref_Oar_rambouillet_v1.0_chr21.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr21.fa && sed -i "1 s/^.*$/>chr21/" ${chr_dir}/oar_v1-0_chr21.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_22/oar_ref_Oar_rambouillet_v1.0_chr22.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr22.fa && sed -i "1 s/^.*$/>chr22/" ${chr_dir}/oar_v1-0_chr22.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_23/oar_ref_Oar_rambouillet_v1.0_chr23.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr23.fa && sed -i "1 s/^.*$/>chr23/" ${chr_dir}/oar_v1-0_chr23.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_24/oar_ref_Oar_rambouillet_v1.0_chr24.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr24.fa && sed -i "1 s/^.*$/>chr24/" ${chr_dir}/oar_v1-0_chr24.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_25/oar_ref_Oar_rambouillet_v1.0_chr25.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr25.fa && sed -i "1 s/^.*$/>chr25/" ${chr_dir}/oar_v1-0_chr25.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_26/oar_ref_Oar_rambouillet_v1.0_chr26.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chr26.fa && sed -i "1 s/^.*$/>chr26/" ${chr_dir}/oar_v1-0_chr26.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_MT/oar_ref_Oar_rambouillet_v1.0_chrMT.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chrMT.fa && sed -i "1 s/^.*$/>chrMT/" ${chr_dir}/oar_v1-0_chrMT.fa
wget -O - ftp://ftp.ncbi.nih.gov/genomes/Ovis_aries/CHR_X/oar_ref_Oar_rambouillet_v1.0_chrX.fa.gz | gunzip -c > ${chr_dir}/oar_v1-0_chrX.fa && sed -i "1 s/^.*$/>chrX/" ${chr_dir}/oar_v1-0_chrX.fa

	echo "[$ts] genome installed"

# The following steps are required to ensure the fasta file is in the correct miRDeep2 format
# Remove whitespace from fasta file
sed 's, ,_,g' -i ${chr_dir}/*.fa
# Replace U with T in the sequence (should not be neccesary)
sed -e '/^[^>]/s/U/T/g' -i ${chr_dir}/*.fa
# Replace non-stand sequence characters with N
sed -e '/^[^>]/s/[^ATGCatcg]/N/g' -i ${chr_dir}/*.fa
# Remove all blank lines 
sed '/^[[:space:]]*$/d' -i ${chr_dir}/*.fa 

# Combine all chromosomal assemblies into single genome fasta file
cat ${chr_dir}/*.fa > ${genome_dir}/oar_v1-0_genome.fa
 
# Rezip chromosomal files
gzip ${chr_dir}/*
 
# Download all mature miRNA sequences from miRBase 
wget -O - ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz | gunzip -c > ${mirna_dir}/mature_mirna.fa
# Linearise the miRNA file (one sequence per line after header) - this is so we can subset to select only oar as below.
awk -v ORS= '/^>/{$0=(NR==1?"":RS)$0RS} END {printf RS}1' ${mirna_dir}/mature_mirna.fa > ${mirna_dir}/mature_mirna_linear.fa
# Ensure FASTA format is applicable for miRDeep2
sed 's,[- ],_,g' -i ${mirna_dir}/mature_mirna_linear.fa
sed -e '/^[^>]/s/U/T/g' -i ${mirna_dir}/mature_mirna_linear.fa
sed -e '/^[^>]/s/[^ATGCatcg]/N/g' -i ${mirna_dir}/mature_mirna_linear.fa
sed '/^[[:space:]]*$/d' -i ${mirna_dir}/mature_mirna_linear.fa 
# Subset FASTA file to isolatee 'oar', 'bta', 'hsa' and 'mmu' (related species) - only works if sequences on one line after header id.
awk '/^>oar.*/{print;getline;print;}' ${mirna_dir}/mature_mirna_linear.fa > ${mirna_dir}/mature_mirna_oar.fa
awk '/^>bta.*/{print;getline;print;}' ${mirna_dir}/mature_mirna_linear.fa > ${mirna_dir}/mature_mirna_bta.fa
awk '/^>hsa.*/{print;getline;print;}' ${mirna_dir}/mature_mirna_linear.fa > ${mirna_dir}/mirna_hsa.fa
awk '/^>mmu.*/{print;getline;print;}' ${mirna_dir}/mature_mirna_linear.fa > ${mirna_dir}/mature_mirna_mmu.fa
 
# Download all hairpin miRNA sequences from miRBase and perform cleaning as above
wget -O - ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz | gunzip -c > ${mirna_dir}/hairpin_mirna.fa
awk -v ORS= '/^>/{$0=(NR==1?"":RS)$0RS} END {printf RS}1' ${mirna_dir}/hairpin_mirna.fa > ${mirna_dir}/hairpin_mirna_linear.fa
sed 's,[- ],_,g' -i ${mirna_dir}/hairpin_mirna_linear.fa
sed -e '/^[^>]/s/U/T/g' -i ${mirna_dir}/hairpin_mirna_linear.fa
sed -e '/^[^>]/s/[^ATGCatcg]/N/g' -i ${mirna_dir}/hairpin_mirna_linear.fa
sed '/^[[:space:]]*$/d' -i ${mirna_dir}/hairpin_mirna_linear.fa
awk '/^>oar.*/{print;getline;print;}' ${mirna_dir}/hairpin_mirna_linear.fa > ${mirna_dir}/hairpin_mirna_oar.fa 

	echo "[$ts] all files standardised"

# build a bowtie index with the newly downloaded sheep genome
bowtie-build ${genome_dir}/oar_v1-0_genome.fa ${genome_index_dir}/oar_v1-0
	echo "[$ts]-- Bowtie done"
 
# Run fastqc on all raw data and output to fastqc raw results directory.
fastqc -o ${raw_FastQC_dir} ${raw_reads_dir}/*

	echo "[$ts]-- Fastqc done"
 
# Forward adapter used for sequencing, if different adapters were used change the sequence here
adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
min_length=18
max_length=28
error_rate=0.3

# Extract sheep number between the second and third underscore '_' these numbers are not unique to one sheep but unique to each experimental group 
rna_seq_r1="${raw_reads_dir}/*_1.*"

cd $raw_reads_dir

for raw_fastq in ${rna_seq_r1[@]}; do

	echo "[$ts]-- Processing Raw fastq file $raw_fastq"

	no_ext=${raw_fastq%%.*}
	file_name=${no_ext##*/}

	interim=${file_name%_*}
	sheep_number=${interim##*_}
	buffer_underscore=$"_"

	#create unique sheep ID
	sheepID=${interim#*_}
	     
	# Assign designated experimental grouping (1W = 1 week 1M = 1 month 3M = 3 months AD = adult)
	if (( ${sheep_number} > 248 )) && (( ${sheep_number} < 256 )); then
		uniq_sheepID="1M${buffer_underscore}${sheepID}"
	elif (( ${sheep_number} > 262 )) && (( ${sheep_number} <271 )); then
		uniq_sheepID="3M${buffer_underscore}${sheepID}"
	elif (( ${sheep_number} > 270 )) && (( ${sheep_number} <279 )); then
		uniq_sheepID="AD${buffer_underscore}${sheepID}"
	else 
		uniq_sheepID="1W${buffer_underscore}${sheepID}"
	fi

	echo "$uniq_sheepID"
	clip_fastq_path=${clip_reads_dir}/${uniq_sheepID}_clip

	echo "[$ts]-- Running cutadapt."

	#run cutadapt and send the output file to the clipped fastq file directory
	cutadapt -a $adapter --trimmed-only -m $min_length -M $max_length -e $error_rate --max-n=0 -o $clip_fastq_path $raw_fastq -j 6
	 
	# If neccesary, gunzip file and remove .gz file extension
	if [[ ${clip_fastq_path} == *gz ]]; then
	    gunzip ${clip_fastq_path}
	    clip_fastq_path=${clip_fastq_path%.*}
	fi
done
	echo "[$ts]-- cutadapt complete, takes a while eh? good thing we use several cores"
#add fastq attachement to the file names in the directories, then move them into the appropiate directories for the different expermiental groups
cd ${clip_reads_dir}
for file in ${clip_reads_dir}; do
	#find -type f -name '*.fastq' | while read f; do mv "$f" "${f%.fastq}"; done
	find . -type f -exec mv '{}' '{}'.fastq \;
done

cd ${mapper_dir}
fastq_files="${clip_reads_dir}/*.fastq*"

for fastq in ${fastq_files}; do
	no_ex=${fastq%.*}
	x=${no_ex%%_*}
	file_n=${x##*/}

if [ "${file_n}" == "1W" ]; then mv ${fastq} ${one_week}
	elif [ "${file_n}" == "1M" ]; then mv ${fastq} ${one_month} 
	elif [ "${file_n}" == "3M" ]; then mv ${fastq} ${three_month} 
	else mv ${fastq} ${AD}
	fi 
echo "$file_n"
done

#remove all the files in the mapper directories if they are not already empty
rm ${oneW_reads}/*.*; rm ${oneM_reads}/*.*; rm ${threeM_reads}/*.* ; rm ${AD_reads}/*.*
rm ${oneW_mapping}/*.* ; rm ${oneM_mapping}/*.* ; rm ${threeW_mapping}/*.* ; rm ${AD_mapping}/*.*
rm -rf ${mapper_output_1W}/* ; rm -rf ${mapper_output_1M}/*; rm -rf ${mapper_output_3M}/*; rm -rf ${mapper_output_AD}/*  

##IMPORTANT, all files will now be concatonated in their respective experimental groups... if you don't want this to happen, change the code here. 
	echo "[$ts] -- files undergoing concatonation into respective experimental groups"
cd ${one_week}; cat ${one_week}/*.fastq > one_week.fastq
	echo "[$ts] -- 1 done"
cd ${one_month}; cat ${one_month}/*.fastq > one_month.fastq
	echo "[$ts] -- 2 done"
cd ${three_month}; cat ${three_month}/*.fastq > three_month.fastq
	echo "[$ts] -- 3 done"
cd ${AD}; cat ${AD}/*.fastq > AD_only.fastq
	echo "[$ts] -- all files concatonated"

#perform mapper.pl on the fastq files in their appropriate directories, the files used to be done one by one in the experimental groups so there is a bit of an artifact with the naming protocol, this should not be an issue as they all end up in their respective directories.  
#ONE WEEK
cd ${mapper_output_1W}
fastq_oneW="${one_week}/*one_week.fastq*"
n=1
oneW_r="${oneW_reads}/${n}_reads.fa"
oneW_mappings="${oneW_mapping}/${n}_genome_vs_reads.arf"

for fastq1W in ${fastq_oneW}; do
	oneW_r="${oneW_reads}/${n}_reads.fa"
	oneW_mappings="${oneW_mapping}/${n}_genome_vs_reads.arf"

	mapper.pl ${fastq1W} -e -v -o 6 -q -h -l 18 -i -j -m -g 01W -p ${genome_index_dir}/oar_v1-0 -r 5 -s ${oneW_r} -t ${oneW_mappings}
	let "n += 1"
done

#ONE MONTH
cd ${mapper_output_1M}
fastq_oneM="${one_month}/*one_month.fastq*"
n=1
oneM_r="${oneM_reads}/${n}_reads.fa"
oneM_mappings="${oneM_mapping}/${n}_genome_vs_reads.arf"

for fastq1M in ${fastq_oneM}; do 
	oneM_r="${oneM_reads}/${n}_reads.fa"
	oneM_mappings="${oneM_mapping}/${n}_genome_vs_reads.arf"
	mapper.pl ${fastq1M} -e -v -o 6 -q -l 18 -h -i -j -m -g 01M -p ${genome_index_dir}/oar_v1-0 -r 5 -s ${oneM_r} -t ${oneM_mappings}
	let "n += 1"
done

#THREE MONTHS
cd ${mapper_output_3M}
fastq_threeM="${three_month}/*three_month.fastq*"
n=1
threeM_r="${threeM_reads}/${n}_reads.fa"
threeM_mappings="${threeM_mapping}/${n}_genome_vs_reads.arf"

for fastq3M in ${fastq_threeM}; do 
	threeM_r="${threeM_reads}/${n}_reads.fa"
	threeM_mappings="${threeM_mapping}/${n}_genome_vs_reads.arf"
	mapper.pl ${fastq3M} -e -v -o 6 -q -l 18 -h -i -j -m -g 03M -p ${genome_index_dir}/oar_v1-0 -r 5 -s ${threeM_r} -t ${threeM_mappings}
	let "n += 1"
done

#ADULT
cd ${mapper_output_AD}
fastq_AD="${AD}/*AD_only.fastq*"
n=1
AD_r="${AD_reads}/${n}_reads.fa"
AD_mappings="${AD_mapping}/${n}_genome_vs_reads.arf"

for fastqAD in ${fastq_AD}; do 
	AD_r="${AD_reads}/${n}_reads.fa"
	AD_mappings="${AD_mapping}/${n}_genome_vs_reads.arf"
	mapper.pl ${fastqAD} -e -v -o 6 -q -l 18 -h -i -j -m -g 0AD -p ${genome_index_dir}/oar_v1-0 -r 5 -s ${AD_r} -t ${AD_mappings}
	let "n += 1"
done

###########################

#remove all whitespaces from the miRNA, .fa and .arf files		 
sed 's/ /_/g' -i ${mirna_dir}/*.*
sed '/^[[:space:]]*$/d' -i ${mirna_dir}/*.* 
	echo "[$ts] -- 1 done"

sed 's, ,_,g' -i ${AD_reads}/*.*
sed '/^[[:space:]]*$/d' -i ${AD_reads}/*.*
	echo "[$ts] -- 2 done "

sed 's, ,_,g' -i ${oneW_reads}/*.*
sed '/^[[:space:]]*$/d' -i ${oneW_reads}/*.*
	echo "[$ts] -- 3 done"

sed 's, ,_,g' -i ${oneM_reads}/*.*
sed '/^[[:space:]]*$/d' -i ${oneM_reads}/*.* 
	echo "[$ts] -- 4 done"

sed 's, ,_,g' -i ${threeM_reads}/*.*
sed '/^[[:space:]]*$/d' -i ${threeM_reads}/*.*
	echo "[$ts] -- 5 done"

sed 's, ,_,g' -i ${genome_dir}/*.*
sed '/^[[:space:]]*$/d' -i ${genome_dir}/*.* 

#move files into appropriate mirDeep directories
cd ${oneW_reads}; mv 1_reads.fa ${mirdeep_oneW}; cd ${oneW_mapping}; mv 1_genome_vs_reads.arf ${mirdeep_oneW}
cd ${oneM_reads}; mv 1_reads.fa ${mirdeep_oneM}; cd ${oneM_mapping}; mv 1_genome_vs_reads.arf ${mirdeep_oneM}
cd ${threeM_reads}; mv 1_reads.fa ${mirdeep_threeM}; cd ${threeM_mapping}; mv 1_genome_vs_reads.arf ${mirdeep_threeM}
cd ${AD_reads}; mv 1_reads.fa ${mirdeep_AD}; cd ${AD_mapping}; mv 1_genome_vs_reads.arf ${mirdeep_AD}

cd ${mirdeep_oneW}
miRDeep2.pl 1_reads.fa ${genome_dir}/oar_v1-0_genome.fa 1_genome_vs_reads.arf ${mirna_dir}/mature_mirna_oar.fa ${mirna_dir}/mirna_hsa.fa ${mirna_dir}/hairpin_mirna_oar.fa 2>report.log -v -P

cd ${mirdeep_oneM}
miRDeep2.pl 1_reads.fa ${genome_dir}/oar_v1-0_genome.fa 1_genome_vs_reads.arf ${mirna_dir}/mature_mirna_oar.fa ${mirna_dir}/mirna_hsa.fa ${mirna_dir}/hairpin_mirna_oar.fa 2>report.log -v -P

cd ${mirdeep_threeM}
miRDeep2.pl 1_reads.fa ${genome_dir}/oar_v1-0_genome.fa 1_genome_vs_reads.arf ${mirna_dir}/mature_mirna_oar.fa ${mirna_dir}/mirna_hsa.fa ${mirna_dir}/hairpin_mirna_oar.fa 2>report.log -v -P

cd ${mirdeep_AD}
miRDeep2.pl 1_reads.fa ${genome_dir}/oar_v1-0_genome.fa 1_genome_vs_reads.arf ${mirna_dir}/mature_mirna_oar.fa ${mirna_dir}/mirna_hsa.fa ${mirna_dir}/hairpin_mirna_oar.fa 2>report.log -v -P

echo "[$ts] MirDeep2 performed succesfully, script finished, please find results in the 'runs' subdirectory of the mirDeep directory"

exit
