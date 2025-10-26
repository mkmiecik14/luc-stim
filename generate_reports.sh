#!/bin/bash

# Ensure we're in the project root directory
cd "$(dirname "$0")"

# Set output directory
OUTPUT_DIR="output/r-reports"

# Generate all 4 bandwise analysis report variants using quarto render
echo "Generating bandwise analysis reports..."

# Create output directory
mkdir -p "$OUTPUT_DIR"

# dB FDR-corrected
echo "Generating dB FDR-corrected report..."
quarto render reports/bandwise-analyses.qmd \
  --output bandwise-analyses-db-fdr.pdf \
  -P results_file:"output/r-analysis/bandwise-analyses-db.rds" \
  -P use_fdr:true \
  -P unit_type:"dB" \
  --to pdf
mv -f bandwise-analyses-db-fdr.pdf "$OUTPUT_DIR"/

# dB uncorrected
echo "Generating dB uncorrected report..."
quarto render reports/bandwise-analyses.qmd \
  --output bandwise-analyses-db-uncorrected.pdf \
  -P results_file:"output/r-analysis/bandwise-analyses-db.rds" \
  -P use_fdr:false \
  -P unit_type:"dB" \
  --to pdf
mv -f bandwise-analyses-db-uncorrected.pdf "$OUTPUT_DIR"/

# uV FDR-corrected
echo "Generating uV FDR-corrected report..."
quarto render reports/bandwise-analyses.qmd \
  --output bandwise-analyses-uv-fdr.pdf \
  -P results_file:"output/r-analysis/bandwise-analyses-uv.rds" \
  -P use_fdr:true \
  -P unit_type:"uV" \
  --to pdf
mv -f bandwise-analyses-uv-fdr.pdf "$OUTPUT_DIR"/

# uV uncorrected
echo "Generating uV uncorrected report..."
quarto render reports/bandwise-analyses.qmd \
  --output bandwise-analyses-uv-uncorrected.pdf \
  -P results_file:"output/r-analysis/bandwise-analyses-uv.rds" \
  -P use_fdr:false \
  -P unit_type:"uV" \
  --to pdf
mv -f bandwise-analyses-uv-uncorrected.pdf "$OUTPUT_DIR"/

echo "All reports generated successfully!"
echo "Reports saved to: $OUTPUT_DIR/"
ls -la "$OUTPUT_DIR"/bandwise-analyses-*.pdf