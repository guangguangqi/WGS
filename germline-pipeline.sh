# Create a new conda environment for WGS analysis
# This isolates our tools from other projects and prevents conflicts
conda create -n wgs_analysis python=3.9
 
# Activate the environment - all subsequent installations will go here
conda activate wgs_analysis
 
# Add bioinformatics channels to conda
# bioconda: specialized channel for bioinformatics software
# conda-forge: community-maintained packages with high quality standards
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict  # Prevents package conflicts
 
# Install GATK and essential tools for the entire workflow
conda install -y \
    gatk4 \          # Main variant calling toolkit
    bwa \            # Read alignment tool
    samtools \       # SAM/BAM file manipulation
    picard \         # Java-based tools for BAM processing
    trim-galore \    # Adapter trimming and quality control
    fastqc \         # Sequencing quality assessment
    bcftools \       # VCF file manipulation and statistics
    snpeff \         # Variant annotation and effect prediction
    bedtools \       # Genomic interval operations
    vcftools \       # VCF filtering and format conversion
    tabix  \         # VCF indexing and querying
    tree             # Show directory structure
 
# Verify GATK installation and check version
gatk --version

# Create organized directory structure for the project
mkdir -p ~/wgs_analysis/{reference,data,results}
cd ~/wgs_analysis/reference
 
# Download human reference genome (GRCh38/hg38 - the current standard)
# This is the "map" we'll align our reads to
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
 
# Create BWA index - this allows BWA to quickly find where reads should align
# This step takes ~2-3 hours but only needs to be done once
bwa index Homo_sapiens_assembly38.fasta

# Download GATK resource bundle - these are "truth sets" of known variants
# Used for base quality score recalibration and variant filtering
 
# HapMap: High-quality SNPs used for training variant filters
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
 
# 1000 Genomes Omni: Another high-quality variant set for training
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
 
# 1000 Genomes high-confidence SNPs: Large collection of validated variants
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
 
# Mills and 1000G indels: High-quality insertions and deletions
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
 
# dbSNP: Database of known variants - helps identify novel vs. known variants
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi


#step 1
# Set up variables and create directory structure
REFERENCE="~/wgs_analysis/reference/Homo_sapiens_assembly38.fasta"
SAMPLE="tumor1"
THREADS=20
 
# Create organized directory structure for results
mkdir -p ~/wgs_analysis/results/{qc,trimmed,aligned,recal,var}/${SAMPLE}
 
# Run FastQC to assess the quality of raw sequencing data
# This generates HTML reports showing:
# - Base quality scores across read positions
# - GC content distribution
# - Adapter contamination levels
# - Overrepresented sequences
fastqc \
    ~/wgs_analysis/data/${SAMPLE}_R1.fastq.gz \
    ~/wgs_analysis/data/${SAMPLE}_R2.fastq.gz \
    -o ~/wgs_analysis/results/qc/${SAMPLE}


#step2 
# Trim adapters and low-quality bases using Trim Galore
# This removes:
# - Illumina sequencing adapters that can interfere with alignment
# - Low-quality bases (quality score < 20) from read ends
# - Very short reads (< 50bp) that may align poorly
trim_galore \
    --paired \           # Indicates we have paired-end reads (R1 and R2)
    --quality 20 \       # Remove bases with quality score < 20
    --length 50 \        # Discard reads shorter than 50bp after trimming
    --fastqc \           # Run FastQC on trimmed reads for comparison
    --output_dir ~/wgs_analysis/results/trimmed/${SAMPLE} \
    ~/wgs_analysis/data/${SAMPLE}_R1.fastq.gz \
    ~/wgs_analysis/data/${SAMPLE}_R2.fastq.gz


#step 3
# Align trimmed reads to the reference genome using BWA-MEM
# Align reads to reference genome using BWA-MEM algorithm
# BWA-MEM is specifically designed for 70bp-1Mbp Illumina reads
bwa mem \
    -t ${THREADS} \     # Use multiple CPU threads for faster processing
    -M \                # Mark shorter split hits as secondary (required for Picard compatibility)
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}_lib" \  # Read group info (required for GATK)
    ${REFERENCE} \      # Reference genome file
    ~/wgs_analysis/results/trimmed/${SAMPLE}/${SAMPLE}_R1_val_1.fq.gz \  # Trimmed forward reads
    ~/wgs_analysis/results/trimmed/${SAMPLE}/${SAMPLE}_R2_val_2.fq.gz \  # Trimmed reverse reads
    | samtools sort -@ ${THREADS} -o ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}.bam  # Sort reads by position
 
# Create index for the BAM file - this allows random access to specific regions
samtools index ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}.bam


#step 4 
# Mark duplicates using Picard to avoid PCR bias in variant calling
# Mark PCR and optical duplicates using GATK MarkDuplicates
# Duplicates arise from:
# - PCR amplification during library preparation
# - Optical duplicates from sequencing (reads from same DNA cluster)
# These can bias variant calling by making variants appear more frequent
gatk MarkDuplicates \
    -I ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}.bam \
    -O ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}_marked_duplicates.bam \
    -M ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}_duplicate_metrics.txt \  # File containing duplication statistics
    --CREATE_INDEX true  # Automatically create BAM index


#step 5
# Base Quality Score Recalibration (BQSR) to correct systematic errors
# Step 5-1: Generate recalibration data table
# This analyzes patterns of base quality score errors using known variant sites
gatk BaseRecalibrator \
    -I ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}_marked_duplicates.bam \
    -R ${REFERENCE} \
    --known-sites ~/wgs_analysis/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz \        # Known SNPs from dbSNP
    --known-sites ~/wgs_analysis/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \  # Known indels
    -O ~/wgs_analysis/results/recal/${SAMPLE}/${SAMPLE}_recal_data.table
 
# Step 5-2: Apply base quality score recalibration
# This updates the quality scores in your BAM file based on the patterns found above
gatk ApplyBQSR \
    -I ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}_marked_duplicates.bam \
    -R ${REFERENCE} \
    --bqsr-recal-file ~/wgs_analysis/results/recal/${SAMPLE}/${SAMPLE}_recal_data.table \
    -O ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bam

#step6
#alignment quantility assenment 
# Collect alignment summary metrics
# This provides statistics about how well your reads aligned to the reference
gatk CollectAlignmentSummaryMetrics \
    -R ${REFERENCE} \
    -I ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bam \
    -O ~/wgs_analysis/results/qc/${SAMPLE}/${SAMPLE}_alignment_summary.txt
 
# Collect insert size metrics (for paired-end data)
# Insert size is the distance between paired reads - important for structural variant detection
gatk CollectInsertSizeMetrics \
    -I ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bam \
    -O ~/wgs_analysis/results/qc/${SAMPLE}/${SAMPLE}_insert_size_metrics.txt \
    -H ~/wgs_analysis/results/qc/${SAMPLE}/${SAMPLE}_insert_size_histogram.pdf


#step7
# Variant calling using GATK HaplotypeCaller
# Call variants in GVCF (Genomic Variant Call Format) mode
# GVCF mode records information about ALL sites (variant and non-variant)
# This allows for accurate joint genotyping when analyzing multiple samples
gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}_recalibrated.bam \
    -O ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}.g.vcf.gz \
    -ERC GVCF \          # Emit Reference Confidence mode - creates GVCF format
    --dbsnp ~/wgs_analysis/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz  # Annotate with known variants


# Convert GVCF to regular VCF for single-sample analysis
# This step is REQUIRED before variant filtering
gatk GenotypeGVCFs \
    -R ${REFERENCE} \
    -V ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}.g.vcf.gz \
    -O ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_raw_variants.vcf.gz


#step8
#Variant Quality Score Recalibration (VQSR)
#SNP
# Build SNP recalibration model using multiple high-quality training sets
# Each training set has different confidence levels and purposes
gatk VariantRecalibrator \
    -R ${REFERENCE} \
    -V ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_raw_variants.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
        ~/wgs_analysis/reference/hapmap_3.3.hg38.vcf.gz \     # Highest quality training set
    --resource:omni,known=false,training=true,truth=false,prior=12.0 \
        ~/wgs_analysis/reference/1000G_omni2.5.hg38.vcf.gz \  # High quality, slightly lower confidence
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
        ~/wgs_analysis/reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz \  # Large training set
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
        ~/wgs_analysis/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz \  # Known variants (not for training)
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \  # Variant annotations to use
    -mode SNP \             # Process SNPs only
    -O ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_snps.recal \
    --tranches-file ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_snps.tranches
 
# Apply SNP recalibration filter
# This assigns PASS/FAIL based on the machine learning model
gatk ApplyVQSR \
    -R ${REFERENCE} \
    -V ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_raw_variants.vcf.gz \
    --recal-file ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_snps.recal \
    --tranches-file ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_snps.tranches \
    -mode SNP \
    --truth-sensitivity-filter-level 99.0 \  # Keep 99% of true variants (high sensitivity)
    -O ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_snps_recalibrated.vcf.gz

###Indel
# Build indel recalibration model
# Indels are harder to call accurately, so we use fewer annotations
gatk VariantRecalibrator \
    -R ${REFERENCE} \
    -V ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_snps_recalibrated.vcf.gz \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 \
        ~/wgs_analysis/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \  # High-quality indels
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
        ~/wgs_analysis/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -an QD -an ReadPosRankSum -an FS -an SOR \  # Fewer annotations for indels
    -mode INDEL \
    --max-gaussians 4 \     # Limit model complexity (indels have less training data)
    -O ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_indels.recal \
    --tranches-file ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_indels.tranches
 
# Apply indel recalibration
# Using 95% sensitivity (slightly more stringent than SNPs)
gatk ApplyVQSR \
    -R ${REFERENCE} \
    -V ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_snps_recalibrated.vcf.gz \
    --recal-file ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_indels.recal \
    --tranches-file ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_indels.tranches \
    -mode INDEL \
    --truth-sensitivity-filter-level 95.0 \
    -O ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz

#step 9
# Functional Annotation with SnpEff
# Annotate variants with functional consequences
# This predicts the effect of each variant on genes and proteins
snpEff ann \
    -Xmx32g \    # Allocate 32GB memory for large genomes
    -stats ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_annotation_stats.html \  # Generate stats report
    GRCh38.105 \  # Latest Ensembl-based database version
    ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz \
    > ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_annotated.vcf

# Convert VCF to a readable table format with key annotation fields
# This extracts the most important information into columns
gatk VariantsToTable \
    -V ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_annotated.vcf \
    -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN \  # Select specific fields to extract
    -O ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_variants_table.tsv


#step 10 Variant Statistics and Quality Assessment
# Generate comprehensive variant statistics using bcftools
# This provides detailed metrics about your variant calls
bcftools stats ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz > ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_variant_stats.txt
 
# Count different types of variants
# SNPs (single nucleotide polymorphisms)
bcftools view -v snps ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz | bcftools query -f '.\n' | wc -l > ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_snp_count.txt
 
# Indels (insertions and deletions)
bcftools view -v indels ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_filtered.vcf.gz | bcftools query -f '.\n' | wc -l > ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_indel_count.txt

#!/bin/bash
# Quick QC Summary Script for WGS Analysis
# This script extracts key quality metrics from all analysis steps
 
SAMPLE=tumor1
 
echo "========================================="
echo "WGS Quality Control Summary for ${SAMPLE}"
echo "========================================="
echo
 
echo "=== ALIGNMENT QUALITY ==="
echo -n "Mapping Rate: "
grep -A1 "FIRST_OF_PAIR" ~/wgs_analysis/results/qc/${SAMPLE}/${SAMPLE}_alignment_summary.txt | tail -1 | cut -f7
echo "  (Benchmark: >95%)"
echo
 
echo -n "Duplicate Rate: "
grep -A1 "LIBRARY" ~/wgs_analysis/results/aligned/${SAMPLE}/${SAMPLE}_duplicate_metrics.txt | tail -1 | cut -f9
echo "  (Benchmark: <30%)"
echo
 
echo -n "Mean Insert Size: "
grep -A1 "MEDIAN_INSERT_SIZE" ~/wgs_analysis/results/qc/${SAMPLE}/${SAMPLE}_insert_size_metrics.txt | tail -1 | cut -f1
echo "bp  (Benchmark: 300-500bp)"
echo
 
echo "=== VARIANT CALLING QUALITY ==="
echo -n "Total SNPs: "
cat ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_snp_count.txt
echo "  (Benchmark: 4-5 million)"
echo
 
echo -n "Total Indels: "
cat ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_indel_count.txt
echo "  (Benchmark: 0.5-0.8 million)"
echo
 
echo -n "Ti/Tv Ratio: "
# Extract Ti/Tv ratio from bcftools stats output
grep -v "^#" ~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_variant_stats.txt | grep "TSTV" | cut -f5
echo "  (Benchmark: 2.0-2.1)"
echo
 
echo "For detailed variant annotation statistics, open:"
echo "~/wgs_analysis/results/var/${SAMPLE}/${SAMPLE}_annotation_stats.html"
echo "========================================="

