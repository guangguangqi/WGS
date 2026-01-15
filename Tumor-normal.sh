#-----------------------------------------------
# STEP 1: Activate existing GATK environment
#-----------------------------------------------
 
# Activate the WGS data analysis environment from Part 1
# If you haven't completed Part 1, please follow that tutorial first
conda activate wgs_analysis

#-----------------------------------------------
# STEP 2a: Download somatic-specific reference files (hg38)
#-----------------------------------------------
 
# Create directory for somatic analysis references
mkdir -p ~/references/somatic_resources
cd ~/references/somatic_resources
 
echo "Downloading somatic analysis reference files..."
 
# Panel of Normals (PON) - contains common technical artifacts
echo "Downloading Panel of Normals..."
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
 
# Germline resource - population allele frequencies from gnomAD
echo "Downloading germline resource..."
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
 
# Common variants for contamination estimation
echo "Downloading contamination resource..."
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi
 
# Link to reference genome from Part 1
echo "Linking reference genome from Part 1..."
ln -s ~/wgs_analysis/reference/Homo_sapiens_assembly38.fasta ./
ln -s ~/wgs_analysis/reference/Homo_sapiens_assembly38.fasta.fai ./
ln -s ~/wgs_analysis/reference/Homo_sapiens_assembly38.dict ./
 
echo "✓ Reference files downloaded successfully!"

#-----------------------------------------------
# STEP 2b: Download somatic-specific reference files (hg19/b37)
#-----------------------------------------------
 
# Create directory for somatic analysis references
mkdir -p ~/references/somatic_resources_hg19
cd ~/references/somatic_resources_hg19
 
echo "Downloading somatic analysis reference files for hg19..."
 
# Panel of Normals (PON) - contains common technical artifacts
echo "Downloading Panel of Normals..."
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf.idx
 
# Germline resource - population allele frequencies from gnomAD
echo "Downloading germline resource..."
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx
 
# Common variants for contamination estimation
echo "Downloading contamination resource..."
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/small_exac_common_3.vcf
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/small_exac_common_3.vcf.idx
 
# Download and prepare hg19 reference genome
echo "Downloading hg19 reference genome..."
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict
 
echo "✓ Reference files downloaded successfully!"

#-----------------------------------------------
# STEP 3: Create organized project directory structure
#-----------------------------------------------
 
# Create main project directory
mkdir -p ~/somatic_analysis_matched
cd ~/somatic_analysis_matched
 
# Create subdirectories for different analysis stages
mkdir -p {input_data,raw_calls,filtered_calls,contamination_analysis}
mkdir -p {converted_tables,maf_files,qc_reports}
 


#-----------------------------------------------
# STEP 4: Link to processed BAM files from Part 1
#-----------------------------------------------
 
# Set paths based on your Part 1 analysis location
PART1_DIR="~/wgs_analysis/results/aligned"
SOMATIC_DIR="~/somatic_analysis_matched"
 
echo "Linking processed BAM files from Part 1..."
 
cd ${SOMATIC_DIR}/input_data
 
# Link tumor1 files (BAM and index)
echo "Linking tumor1 files..."
ln -s ${PART1_DIR}/tumor1/tumor1_recalibrated.bam tumor1_recalibrated.bam
ln -s ${PART1_DIR}/tumor1/tumor1_recalibrated.bai tumor1_recalibrated.bai
 
# Link normal1 files (BAM and index)
echo "Linking normal1 files..."
ln -s ${PART1_DIR}/normal1/normal1_recalibrated.bam normal1_recalibrated.bam
ln -s ${PART1_DIR}/normal1/normal1_recalibrated.bai normal1_recalibrated.bai
 
cd ${SOMATIC_DIR}
echo "✓ Input data preparation complete!"

#-----------------------------------------------
# STEP 5: Run Mutect2 for somatic variant detection
#-----------------------------------------------
 
# Set up variables for clarity and reusability
REFERENCE="~/references/somatic_resources/Homo_sapiens_assembly38.fasta"
GERMLINE_RESOURCE="~/references/somatic_resources/af-only-gnomad.hg38.vcf.gz"
PON="~/references/somatic_resources/1000g_pon.hg38.vcf.gz"
INPUT_DIR="${SOMATIC_DIR}/input_data"
OUTPUT_DIR="${SOMATIC_DIR}/raw_calls"
 
echo "Running Mutect2 to identify potential somatic mutations..."
 
# Run Mutect2 to identify potential somatic mutations
gatk Mutect2 \
    -R $REFERENCE \
    -I ${INPUT_DIR}/tumor1_recalibrated.bam \
    -I ${INPUT_DIR}/normal1_recalibrated.bam \
    -tumor tumor1 \
    -normal normal1 \
    --germline-resource $GERMLINE_RESOURCE \
    --panel-of-normals $PON \
    --f1r2-tar-gz ${OUTPUT_DIR}/tumor1_f1r2.tar.gz \
    -O ${OUTPUT_DIR}/tumor1_raw.vcf.gz \
    --native-pair-hmm-threads 8 \
    --max-reads-per-alignment-start 50
 
echo "✓ Mutect2 variant calling complete!"
 
# Generate basic statistics about the raw calls
echo "Generating call statistics..."
bcftools stats ${OUTPUT_DIR}/tumor1_raw.vcf.gz > ${OUTPUT_DIR}/tumor1_raw_stats.txt


#-----------------------------------------------
# STEP 6: Generate pileup summaries for contamination analysis
#-----------------------------------------------
 
COMMON_VARIANTS="~/references/somatic_resources/small_exac_common_3.hg38.vcf.gz"
CONTAM_DIR="${SOMATIC_DIR}/contamination_analysis"
 
echo "Generating pileup summaries for contamination assessment..."
 
# Generate pileup summary for tumor sample
echo "Processing tumor1 sample..."
gatk GetPileupSummaries \
    -I ${INPUT_DIR}/tumor1_recalibrated.bam \
    -V $COMMON_VARIANTS \
    -L $COMMON_VARIANTS \
    -O ${CONTAM_DIR}/tumor1_pileups.table
 
# Generate pileup summary for normal sample
echo "Processing normal1 sample..."
gatk GetPileupSummaries \
    -I ${INPUT_DIR}/normal1_recalibrated.bam \
    -V $COMMON_VARIANTS \
    -L $COMMON_VARIANTS \
    -O ${CONTAM_DIR}/normal1_pileups.table
 
echo "✓ Pileup summaries generated successfully!"

#-----------------------------------------------
# STEP 7: Calculate contamination estimates
#-----------------------------------------------
 
echo "Calculating contamination estimates..."
 
# Calculate contamination by comparing tumor vs normal allele frequencies
gatk CalculateContamination \
    -I ${CONTAM_DIR}/tumor1_pileups.table \
    -matched ${CONTAM_DIR}/normal1_pileups.table \
    -O ${CONTAM_DIR}/tumor1_contamination.table \
    --tumor-segmentation ${CONTAM_DIR}/tumor1_segments.table
 
echo "✓ Contamination analysis complete!"

#-----------------------------------------------
# STEP 8: Learn read orientation artifacts for filtering
#-----------------------------------------------
 
echo "Learning read orientation artifacts..."
 
# Analyze the read orientation data collected during Mutect2 calling
gatk LearnReadOrientationModel \
    -I ${OUTPUT_DIR}/tumor1_f1r2.tar.gz \
    -O ${OUTPUT_DIR}/tumor1_orientation_model.tar.gz
 
echo "✓ Read orientation model training complete!"

#-----------------------------------------------
# STEP 9: Apply comprehensive Mutect2 filtering
#-----------------------------------------------
 
echo "Applying FilterMutectCalls to remove artifacts..."
 
# Apply GATK's comprehensive filtering to remove false positive calls
gatk FilterMutectCalls \
    -R $REFERENCE \
    -V ${OUTPUT_DIR}/tumor1_raw.vcf.gz \
    --contamination-table ${CONTAM_DIR}/tumor1_contamination.table \
    --tumor-segmentation ${CONTAM_DIR}/tumor1_segments.table \
    --ob-priors ${OUTPUT_DIR}/tumor1_orientation_model.tar.gz \
    -O ${SOMATIC_DIR}/filtered_calls/tumor1_filtered.vcf.gz
 
echo "✓ FilterMutectCalls complete!"
 
# Generate filtering statistics
bcftools view -H ${SOMATIC_DIR}/filtered_calls/tumor1_filtered.vcf.gz | cut -f7 | sort | uniq -c | sort -nr


#-----------------------------------------------
# STEP 10: Apply additional quality filters for high-confidence calls
#-----------------------------------------------
 
echo "Applying additional quality filters..."
 
# Extract only variants that passed FilterMutectCalls
bcftools view -f PASS \
    ${SOMATIC_DIR}/filtered_calls/tumor1_filtered.vcf.gz \
    -O z \
    -o ${SOMATIC_DIR}/filtered_calls/tumor1_pass.vcf.gz
 
# Apply stringent quality filters for high-confidence calls
bcftools filter \
    -i 'FORMAT/AF[0:0] >= 0.05 && FORMAT/DP[0:0] >= 10 && INFO/TLOD >= 6.3 && (FORMAT/AF[0:1] <= 0.03 || FORMAT/AF[0:1] == ".")' \
    ${SOMATIC_DIR}/filtered_calls/tumor1_pass.vcf.gz \
    -O z \
    -o ${SOMATIC_DIR}/filtered_calls/tumor1_high_confidence.vcf.gz
 
# Index the final high-confidence VCF file
bcftools index -t ${SOMATIC_DIR}/filtered_calls/tumor1_high_confidence.vcf.gz
 
# Generate final statistics
raw_count=$(bcftools view -H ${OUTPUT_DIR}/tumor1_raw.vcf.gz | wc -l)
pass_count=$(bcftools view -H ${SOMATIC_DIR}/filtered_calls/tumor1_pass.vcf.gz | wc -l)
hc_count=$(bcftools view -H ${SOMATIC_DIR}/filtered_calls/tumor1_high_confidence.vcf.gz | wc -l)
 
echo ""
echo "Final filtering cascade results:"
echo "  Raw Mutect2 calls: $raw_count"
echo "  PASS calls: $pass_count"
echo "  High-confidence calls: $hc_count"
echo "  Final success rate: $(echo "scale=2; $hc_count * 100 / $raw_count" | bc)%"
 
echo "✓ Quality filtering complete"
 
# Statistcs for the tumor sample used in this tutorial
# Raw Mutect2 calls: 107878
# PASS calls: 1686
# High-confidence calls: 1453
# Final success rate: 1.34%


#-----------------------------------------------
# STEP 11: Convert VCF to human-readable tables
#-----------------------------------------------
 
echo "Converting VCF files to human-readable tables..."
 
# Use GATK's VariantsToTable for comprehensive data extraction
gatk VariantsToTable \
    -V ${SOMATIC_DIR}/filtered_calls/tumor1_high_confidence.vcf.gz \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
    -F TLOD -F NLOD -F ECNT \
    -GF GT -GF AD -GF AF -GF DP \
    -O ${SOMATIC_DIR}/converted_tables/tumor1_mutations.tsv
 
echo "✓ Human-readable table created!"


