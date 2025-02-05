#!/bin/bash

conda activate dnaseq
#List of compatible packages in dnaseq channel: 
#sratoolkit
#fastqc
#fastp
#cutadapt
#trimmomatic
#bwa
#samtools
#picard
#gatk4
#entrez-direct

#Spostare i file da linux a windows
for file in './4fastqc_output/*.html'; do
	cp $file /mnt/c/Users/Davide/Downloads/
done 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 0 download files 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Bioproject link: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1166965

#Experiment link case: https://www.ncbi.nlm.nih.gov/sra/SRR30834327
#Experiment link control: https://www.ncbi.nlm.nih.gov/sra/SRR30834326

prefetch SRR30834327
prefetch SRR30834326

fasterq-dump SRR30834327
fasterq-dump SRR30834326

# download hg19 chromosome fasta files
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 1 Fastqc: reads quality control 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Set up directories 
mkdir 0input_files 
mkdir 1fastqc_output/

for file in './*.fastq'; do
	mv $file 0input_files/
done 

#Run FastQC
fastqc 0input_files/* -o 1fastqc_output/

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 2 Clean fastq files

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Identify adapters 
for forward in 0input_files/*_1.fastq; do
	reverse=${forward/_1/_2}
	echo "Data analysis with fastp for file $forward and $reverse"
	fastp -i $forward -I $reverse --detect_adapter_for_pe
	#rm fastp.json fastp.html
done

#The following adapters were identified: 
# >Illumina TruSeq Adapter Read 1
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# >Illumina TruSeq Adapter Read 2
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# >Illumina TruSeq Adapter Read 1
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# >Illumina TruSeq Adapter Read 2
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

mkdir 2cut
for forward in 0input_files/*_1.fastq; do 
	reverse=${forward/_1/_2}
	basename_f=$(basename $forward)
	basename_r=$(basename $reverse)
	out_for=2cut/$basename_f
	out_rev=2cut/$basename_r

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o $out_for -p $out_rev  \
	--cores 10 \
	$forward $reverse
done

#Trim low quality reads
mkdir 3trimmed
out_dir="3trimmed"
for file_forward in `ls 2cut/*_1.fastq`; do
	file_reverse=${file_forward/_1./_2.}
	basename_f=$(basename $file_forward)
	basename_r=$(basename $file_reverse)	
	for_pai=$out_dir/${basename_f/_1/_1_PAI}
	rev_pai=$out_dir/${basename_r/_2/_2_PAI}
	for_unp=$out_dir/${basename_f/_1/_1_UNP}
	rev_unp=$out_dir/${basename_r/_2/_2_UNP}
	echo $for_pai
	trimmomatic PE -phred33 -threads 10 $file_forward $file_reverse $for_pai $for_unp $rev_pai $rev_unp \
	TRAILING:30 LEADING:30 MINLEN:75 AVGQUAL:30
done

#Check read quality after adapter removal and trimming
mkdir 4fastqc_output/
fastqc 3trimmed/*  -o 4fastqc_output/    

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 3 Alignment 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#ref tutorial: https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day1/Sequence%20Alignment_July2015_ShamithSamarajiwa.pdf

mkdir 5bwa
mkdir 5bwa/hg19 
#after unzipping in . 
mv chr*.fa 5/bwa/hg19

#concatenate chromosomes in single file 
cat *.fa > hg19.fa
rm chr*.fa

#create reference index 
bwa index -p hg19bwaidx -a bwtsw hg19.fa

for forward in 3trimmed/*_1_PAI.fastq; do 
	reverse=${forward/_1/_2}
	indexed_forward=5bwa/$(basename $forward _PAI.fastq).sai
	indexed_reverse=5bwa/$(basename $reverse _PAI.fastq).sai
	bwa aln -t 4 5bwa/hg19/hg19bwaidx $forward > $indexed_forward
	bwa aln -t 4 5bwa/hg19/hg19bwaidx $reverse > $indexed_reverse
	bwa sampe 5bwa/hg19/hg19bwaidx $indexed_forward $indexed_reverse $forward $reverse > ${indexed_forward/.sai/2_pe.sam}
done

mkdir 6samtools/
for sequence in 5bwa/*.sam; do 
	basename=$(basename $sequence .sam)
	#samtools view -bT hg19.fa sequence1.sam > sequence1.bam # when no header
	samtools view -bS $sequence > 6samtools/$basename.bam # when SAM header present
	#Could also sort and index alignment here but in this particular pipeline would be performed in STEP 4 
	#samtools sort -O bam -o 6samtools/$basename.sorted.bam -T temp 6samtools/$basename.bam # sort by coordinate to streamline data processing
	#samtools index 6samtools/$basename.sorted.bam # a position-sorted BAM file can also beindexed
done

#check alignment quality 
samtools stats file.bam 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 4 Processing of alignment

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

mkdir 7picard/

#To identify and remove PCR duplicates Read groups (RG) are necessary in the BAM files

#Here we use entrez-direct to retrieve and add RG info to alignments 
for alignment in 6samtools/*.bam; do
    run_basename=$(basename "$alignment" _12_pe.bam)  # Extract run name from BAM filename

    echo "Processing: $alignment (Run: $run_basename)"

    # Retrieve metadata for the current run numbers correspond to column name of data of interest e.g. $1 -> Run 
    metadata=$(esearch -db sra -query "$run_basename" | efetch -format runinfo | awk -F ',' 'NR>1 {print $1","$12","$19","$20","$26}')
    
    if [[ -z "$metadata" ]]; then
        echo "WARNING: No metadata found for $run_basename, skipping..."
        continue
    fi

    # Parse metadata into separate variables
    IFS=',' read -r Run LibraryName Platform Model Sample <<< "$metadata"

    # Check if the BAM file exists
    if [[ -f "$alignment" ]]; then
        echo "Adding read group to $alignment"
        picard AddOrReplaceReadGroups \
            -I "$alignment" \
            -O "7picard/${run_basename}_RG.bam" \
            -RGID "$Run" \
            -RGLB "$LibraryName" \
            -RGPL "$Platform" \
            -RGPU "$Model" \
            -RGSM "$Sample" \
			-VALIDATION_STRINGENCY LENIENT #This parameter is necessary to skip unmapped reads that are not crucial for the subsequent analysis
    else
        echo "ERROR: BAM file $alignment not found!"
    fi
done

for alignment in 6samtools/*.bam; do
	echo $alignment
	# picard MarkDuplicates \
	# 	-I $alignment \
	# 	-O 7picard/$(basename $alignment .bam)_marked_duplicates.bam \
	# 	-M 7picard/$(basename $alignment .bam)_marked_dup_metrics.txt
	# 	REMOVE_DUPLICATES=true
done

for case in 6samtools/*7*.bam; do 
	normal=${case/7/6}
	echo $case
	echo $normal 
	# picard MarkDuplicates \
	# 	-I sorted.bam \
	# 	-O marked_duplicates.bam \
	# 	-M marked_dup_metrics.txt
done

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 5 Copy number calling pipeline

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

conda deactivate 
conda activate cnvkit #necessary because python version in (dnaseq) is not the compatible with (cnvkit) 
#List of compatible packages in cnvkit:
#cnvkit 
#samtools


#check read length to download correct mappability file 
for sorted in 6samtools/*.sorted.bam; do
	samtools view $sorted | awk '{print length($10)}' | sort | uniq -c
done

#being majority of reads > 100 we'll use 100mer.bigwig 
#download: https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/

mkdir 7cnvkitin/
mkdir 8cnvkitout/
#Downloaded exome capture kit file, mappability file are moved in 7cnvkit/ 

#For better accuracy calculating mappability with cnvkit 
cnvkit.py access 5bwa/hg19/hg19.fa -o 7cnvkit/access-hg19.bed

#Cnvkit pipeline 
num_cores=$(nproc)
echo "Number of available cores: $num_cores"

for case in 6samtools/*7*.sorted.bam #SRR30834327 is tumor (case) outdir would be called as case 
do 
	sample_name=$(basename "$case" .sorted.bam)
	normal=${case/7/6}
	#Check that normal exists 
    if [ ! -f "$normal" ]; then
        echo "Error: normal file not found. Skipping $case."
		#could be flat reference with --output-reference "$sample_name.cnn"
        continue
    fi

	if ! [ -d "8cnvkitout/$sample_name" ]; then
		echo "Creating directory $sample_name"
		mkdir "8cnvkitout/$sample_name"
		echo "Launching batch command for $(basename $case)"
		cnvkit.py batch $case -n $normal --targets 7cnvkitin/S03723314_Covered.bed --fasta 5bwa/hg19/hg19.fa --access 7cnvkitin/access-hg19.bed --output-dir "8cnvkitout/$sample_name" --scatter -p $num_cores --drop-low-coverage 
	fi 
done


