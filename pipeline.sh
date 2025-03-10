#!/bin/bash

conda activate dnaseq
#List of compatible packages in dnaseq channel: (TS2 are on server) 
#sratoolkit
#fastqc
#fastp
#cutadapt
#trimmomatic
#bwa
#samtools (TS2)
#picard (TS2)
#gatk4 (TS2)
#bcftools (TS2)
#htslib (TS2 only)
#vcftools (TS2 only)
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

#Sort alignment, mark and remove duplicates, index output files 
for alignment in 7picard/*.bam; do
	echo "Processing file: $alignment"
	base_name=${alignment%.bam}

	echo "Sorting BAM file before marking duplicates..."
    samtools sort "$alignment" -o "${base_name}_sorted.bam"

	echo "Running picard MarkDuplicates..."
	picard MarkDuplicates \
		-I "${base_name}_sorted.bam" \
		-O ${base_name}_marked_duplicates.bam \
		-M ${base_name}_marked_dup_metrics.txt \
		-VALIDATION_STRINGENCY LENIENT \
		--REMOVE_DUPLICATES true

	echo "Indexing BAM file..."
    samtools index "${base_name}_marked_duplicates.bam"
done

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 5 Download and process db references

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#annovar download wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
#added annovar to ~/.bashrc 

#download 
#files:
# dbSNP, Mills, 1000Genomes https://data.broadinstitute.org/snowman/hg19/variant_calling/vqsr_resources/Exome/v2/
# gnomAD https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz from https://gnomad.broadinstitute.org/downloads#v2-core-dataset
# EXAC https://hgdownload.soe.ucsc.edu/gbdb/hg19/ExAC/ExAC.r0.3.sites.vep.hg19.vcf.gz

# gnomAD with annovar: annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome humandb/ => vcfconverter.sh => *.vcf

#ESP6500SI https://hgdownload.soe.ucsc.edu/gbdb/hg19/evs/
#download all files 
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/evs/ESP6500SI-V2-SSA137.updatedProteinHgvs.chr{1..22}.snps_indels.vcf.gz
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/evs/ESP6500SI-V2-SSA137.updatedProteinHgvs.chrX.snps_indels.vcf.gz
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/evs/ESP6500SI-V2-SSA137.updatedProteinHgvs.chrY.snps_indels.vcf.gz

#concatenate all chromosomes of ESP65000 in one file 
for file in ESP6500SI-V2-SSA137.updatedProteinHgvs.chr*.snps_indels.vcf.gz; do
    bcftools index "$file"
done

ls ESP6500SI-V2-SSA137.updatedProteinHgvs.chr*.snps_indels.vcf.gz | sort > vcf_list.txt

bcftools concat -a -f vcf_list.txt -O z -o ESP6500SI-V2-SSA137_merged.vcf.gz

#remove single files 
rm ESP6500SI-V2-SSA137.updatedProteinHgvs.chr*.snps_indels.vcf.gz
rm ESP6500SI-V2-SSA137.updatedProteinHgvs.chr*.snps_indels.vcf.gz.csi

#rename gnomeAD file
mv gnomad.exomes.r2.1.1.sites.vcf.bgz gnomad.exomes.r2.1.1.sites.vcf.gz

#create dict file for hg19 
gatk CreateSequenceDictionary \
    -R hg19/hg19.fa \
    -O hg19/hg19.dict

#convert b37 files in hg19 
#chain file downloaded from https://github.com/broadgsa/gatk/blob/master/public/chainFiles/b37tohg19.chain

###############
#Convert reference files from b37 to hg19 using CrossMap
###############

#new env for CrossMap (incompatible)
conda create --name crossmap
conda activate crossmap 

# Install crossmap 
conda install -c conda-forge compilers
conda install -c conda-forge cython
conda install -c conda-forge libgcc
conda install -c bioconda htslib
conda install pip 
pip install CrossMap

# Path configuration
CHAIN_FILE="hg19/b37tohg19.chain"          # Path to the chain file for b37 to hg19 conversion
REF_FASTA="hg19/hg19.fa"                   # hg19 reference FASTA file (must have .fai and .dict)
INPUT_DIR="refdb"                          # Directory containing original b37 files
REJECT_DIR="refdb/rejected"                # Directory for unmapped variants

# Check if required files exist
if [ ! -f "$CHAIN_FILE" ]; then
    echo "Error: Chain file $CHAIN_FILE not found"
    exit 1
fi

if [ ! -f "$REF_FASTA" ]; then
    echo "Error: Reference FASTA file $REF_FASTA not found"
    exit 1
fi

# Check if bgzip is installed
if ! command -v bgzip &> /dev/null; then
    echo "Error: bgzip not found. Please install htslib with: conda install -c bioconda htslib"
    exit 1
fi

# Create output directories
mkdir -p $REJECT_DIR

# Loop through all b37 VCF files
for file in $INPUT_DIR/*.b37.vcf.gz; do
  # Check if matching files exist
  if [ ! -f "$file" ]; then
    echo "No .b37.vcf.gz files found in $INPUT_DIR"
    exit 1
  fi
  
  # Extract base name (e.g., dbsnp_138.b37.vcf.gz → dbsnp_138)
  base_name=$(basename "$file" .b37.vcf.gz)
  
  # Define temporary uncompressed output file
  temp_output="$INPUT_DIR/${base_name}.hg19.vcf"
  # Define final bgzip-compressed output file
  #final_output="$INPUT_DIR/${base_name}.hg19.vcf.gz"
  
  echo "Converting $file to hg19 format..."
  
  # Run CrossMap conversion 
  CrossMap vcf --chromid a \
    $CHAIN_FILE \
    $file \
    $REF_FASTA \
    $temp_output \
    2>&1 | tee "$REJECT_DIR/${base_name}.log"
  
  # Check if conversion was successful
  if [ -f "$temp_output" ]; then
    echo "Compressing with bgzip..."
    # Compress with bgzip (block gzip) for proper genomic data handling
    bgzip -f "$temp_output"
    echo "Output file compressed with bgzip: ${temp_output}.gz"
    
    # Create tabix index for the compressed VCF
    # tabix -p vcf "${temp_output}.gz"
    # echo "Created tabix index: ${temp_output}.gz.tbi"
  else
    echo "Error: CrossMap did not create output file $temp_output"
    continue
  fi
  
  # Handle unmapped variants if they exist
  if [ -f "${temp_output}.unmap" ]; then
    mv "${temp_output}.unmap" "$REJECT_DIR/${base_name}.rejected.vcf"
    echo "Unmapped variants saved to $REJECT_DIR/${base_name}.rejected.vcf"
  fi
done

#move all b37 files in other directory 
for file in refdb/*b37.vcf.gz; do 
	mv $file notuseddb 
done

conda deactivate 
conda activate dnaseq
###############

#Contigs in Mills and dbsnp still are formatted as "ID=1" instead of "ID=chr1"
#Following script is to convert formatting to prevent gatk errors. 
#Create chr mapping
echo -e "1\tchr1\n2\tchr2\n3\tchr3\n4\tchr4\n5\tchr5\n6\tchr6\n7\tchr7\n8\tchr8\n9\tchr9\n10\tchr10\n11\tchr11\n12\tchr12\n13\tchr13\n14\tchr14\n15\tchr15\n16\tchr16\n17\tchr17\n18\tchr18\n19\tchr19\n20\tchr20\n21\tchr21\n22\tchr22\nX\tchrX\nY\tchrY\nMT\tchrM" > chr_map.txt

rename_and_compress_vcf() {
    local chr_map="$1"     
    local input_vcf="$2"   
    local output_vcf="$3"  

	#Check input files
    if [[ ! -f "$chr_map" ]]; then
        echo "Error: Chromosome mapping file '$chr_map' not found!" >&2
        return 1
    fi

    if [[ ! -f "$input_vcf" ]]; then
        echo "Error: Input VCF file '$input_vcf' not found!" >&2
        return 1
    fi

    #Change formatting 
    bcftools annotate --rename-chrs "$chr_map" "$input_vcf" -Oz -o "$output_vcf"

    #Check output file 
    if [[ -f "$output_vcf" ]]; then
        echo "Success: VCF file processed as '$output_vcf'"
    else
        echo "Error: Something went wrong, '$output_vcf' was not created!" >&2
        return 1
    fi
}

rename_and_compress_vcf chr_map.txt refdb/dbsnp_138.hg19.vcf.gz refdb/dbsnp_138.hg19.fixed.vcf.gz
rename_and_compress_vcf chr_map.txt refdb/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz refdb/Mills_and_1000G_gold_standard.indels.hg19.fixed.vcf.gz

#Index all files 
for file in refdb/*.vcf.gz; do 
	if [ ! -f "${file}.tbi" ]; then 
		tabix -p vcf "$file" 
	fi
done 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 6 Base recalibration 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

mkdir 8gatk/ 

#paths configuration 
REF_GENOME="hg19/hg19.fa"  # .fai also present in hg19/
OUTPUT_DIR="8gatk"
DBSNP="refdb/dbsnp_138.hg19.vcf.gz"
MILLS="refdb/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz"
G1000="refdb/1000G_phase1.snps.high_confidence.hg19.vcf.gz"

#Base Recalibration
for file in input_files/*marked_duplicates.bam; do 
sample_name=$(basename "$file" _marked_duplicates.bam)
	gatk BaseRecalibrator \
	-I $file \
	-R $REF_GENOME \
	--known-sites $DBSNP \
	--known-sites $MILLS \
	--known-sites $G1000 \
	-O ${OUTPUT_DIR}/${sample_name}_recal_data.table

	gatk ApplyBQSR \
	-I $file \
	-R $REF_GENOME \
	--bqsr-recal-file ${OUTPUT_DIR}/${sample_name}_recal_data.table \
	-O ${OUTPUT_DIR}/${sample_name}_tumor_recal.bam

	#samtools index ${OUTPUT_DIR}/${sample_name}_tumor_recal.bam
done

#Annotation of known variants with 

#MuTect2 
# Useful info: Mutect2 does not require a germline resource nor a panel of normals (PoN) to run, although both are recommended. 
# The tool prefilters sites for the matched normal and the PoN. 
# If a variant is absent from a given germline resource, then the value for --af-of-alleles-not-in-resource is used as an imputed allele frequency. 
# Below is an excerpt of a known variants resource with population allele frequencies
# https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2

for case in 7picard/*7*.bam; do #SRR30834327 is tumor (case) 
	normal=${case/7/6} 
	echo $case
	echo $normal 
	# gatk Mutect2 \
	# 	-R 5bwa/hg19/hg19.fa \
	# 	-I $case \
	# 	-I $normal \
	# 	-tumor tumor_sample_name \
	# 	-normal normal_sample_name \
	# 	-O 8gatk/somatic_variants.vcf
done

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP X Copy number calling pipeline

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
cnvkit.py access 5bwa/hg19/hg19.fa -o 7cnvkitin/access-hg19.bed

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


