#-----------------------------------------------
# STEP 6: Run tumor-only analysis with custom PON
#-----------------------------------------------
 
echo "Running tumor-only analysis with custom Panel of Normals..."
 
# Run Mutect2 in tumor-only mode using our custom PON
gatk Mutect2 \
    -R $REFERENCE \
    -I ${INPUT_DIR}/tumor1_recalibrated.bam \       # Tumor sample only
    -tumor tumor1 \                          # Tumor sample name
    --germline-resource $GERMLINE_RESOURCE \ # Population frequencies
    --panel-of-normals pon_creation/custom_pon.vcf.gz \  # Our custom PON
    --f1r2-tar-gz raw_calls/tumor1_custom_pon_f1r2.tar.gz \
    -O raw_calls/tumor1_custom_pon_raw.vcf.gz \
    --native-pair-hmm-threads 8
 
echo "Tumor-only calling with custom PON complete!"
 
# Generate statistics
custom_pon_variants=$(bcftools view -H raw_calls/tumor1_custom_pon_raw.vcf.gz | wc -l)
echo "Raw variants with custom PON: $custom_pon_variants"

#-----------------------------------------------
# STEP 9: Apply aggressive filtering for tumor-only analysis
#-----------------------------------------------
 
echo "Applying aggressive filtering for tumor-only analysis..."
 
# Contamination analysis (tumor-only - less reliable)
gatk GetPileupSummaries \
    -I ${INPUT_DIR}/tumor1_recalibrated.bam \
    -V $COMMON_VARIANTS \
    -L $COMMON_VARIANTS \
    -O contamination/tumor1_only_pileups.table
 
gatk CalculateContamination \
    -I contamination/tumor1_only_pileups.table \
    -O contamination/tumor1_only_contamination.table
 
# Learn read orientation model
gatk LearnReadOrientationModel \
    -I raw_calls/tumor1_only_f1r2.tar.gz \
    -O raw_calls/tumor1_only_orientation_model.tar.gz
 
# Apply FilterMutectCalls
gatk FilterMutectCalls \
    -R $REFERENCE \
    -V raw_calls/tumor1_only_raw.vcf.gz \
    --contamination-table contamination/tumor1_only_contamination.table \
    --ob-priors raw_calls/tumor1_only_orientation_model.tar.gz \
    -O filtered_calls/tumor1_only_filtered.vcf.gz
 
# Extract PASS variants
bcftools view -f PASS \
    filtered_calls/tumor1_only_filtered.vcf.gz \
    -O z \
    -o filtered_calls/tumor1_only_pass.vcf.gz
 
# Apply very stringent quality filters for tumor-only analysis
# Higher thresholds compensate for lack of normal sample
bcftools filter \
    -i 'FORMAT/AF[0:0] >= 0.10 && FORMAT/DP[0:0] >= 20 && INFO/TLOD >= 10.0 && INFO/POPAF < 0.001' \
    filtered_calls/tumor1_only_pass.vcf.gz \
    -O z \
    -o filtered_calls/tumor1_only_high_confidence.vcf.gz
 
# Index final VCF
bcftools index -t filtered_calls/tumor1_only_high_confidence.vcf.gz
 
echo "Aggressive filter criteria for tumor-only analysis:"
echo "  Tumor AF ≥ 10% (high frequency threshold)"
echo "  Tumor depth ≥ 20 reads (high confidence requirement)"
echo "  TLOD ≥ 10.0 (very strong statistical evidence)"
echo "  Population AF < 0.1% (aggressive germline filtering)"
 
# Final statistics
tumor_only_hc=$(bcftools view -H filtered_calls/tumor1_only_high_confidence.vcf.gz | wc -l)
echo "High-confidence tumor-only variants: $tumor_only_hc"
 
# Return to main analysis directory
cd ~/somatic_analysis_matched
 
echo "Tumor-only analysis complete!"