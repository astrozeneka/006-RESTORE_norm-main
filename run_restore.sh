#!/bin/bash

# RESTORE Normalization Runner Script
# This script runs the complete RESTORE normalization pipeline

set -e  # Exit on any error

echo "=== RESTORE Normalization Pipeline ==="
echo "Starting at $(date)"

# Configuration
CSV_FILE="sample_cell_data.csv"
MARKER_FILE="sample_markers.csv"
SAVE_DIR="./results"
FLOOR=50

# Check if input files exist
if [ ! -f "$CSV_FILE" ]; then
    echo "Error: $CSV_FILE not found!"
    echo "Please make sure you have the sample data file in the current directory."
    exit 1
fi

if [ ! -f "$MARKER_FILE" ]; then
    echo "Error: $MARKER_FILE not found!"
    echo "Please make sure you have the marker file in the current directory."
    exit 1
fi

# Check if Python scripts exist
if [ ! -f "calculate_threshold.py" ]; then
    echo "Error: calculate_threshold.py not found!"
    echo "Please make sure the RESTORE Python scripts are in the current directory."
    exit 1
fi

if [ ! -f "normalize.py" ]; then
    echo "Error: normalize.py not found!"
    echo "Please make sure the RESTORE Python scripts are in the current directory."
    exit 1
fi

# Create output directory structure
echo "Creating output directories..."
mkdir -p "$SAVE_DIR"
mkdir -p "$SAVE_DIR/img/threshs"
mkdir -p "$SAVE_DIR/thresh_dicts"

# Validate data (if validation script exists)
if [ -f "validate_sample_data.py" ]; then
    echo "Validating input data..."
    python validate_sample_data.py
    echo ""
fi

# Step 1: Calculate thresholds for all markers
echo "Step 1: Calculating thresholds for all markers..."
echo "This will run 6 parallel processes (one for each marker)"

# Array of marker indices (0-5 for CD68, CD31, CD45, CK19, CD3, CD20)
MARKER_INDICES=(0 1 2 3 4 5)
MARKER_NAMES=("CD68" "CD31" "CD45" "CK19" "CD3" "CD20")

# Start parallel threshold calculations
for i in "${MARKER_INDICES[@]}"; do
    echo "  Starting threshold calculation for marker index $i (${MARKER_NAMES[$i]})..."
    python calculate_threshold.py \
        --CSV "$CSV_FILE" \
        --mfname "$MARKER_FILE" \
        --pos_idx "$i" \
        --save_dir "$SAVE_DIR" \
        --floor "$FLOOR" &
done

# Wait for all background processes to complete
echo "  Waiting for all threshold calculations to complete..."
wait

echo "✓ All threshold calculations completed!"

# Check if threshold files were created
echo "Checking threshold calculation results..."
THRESHOLD_FILES_FOUND=0
for i in "${MARKER_INDICES[@]}"; do
    MARKER_NAME="${MARKER_NAMES[$i]}"
    THRESHOLD_FILE="$SAVE_DIR/thresh_dicts/sample_cell_data/sample_cell_data_${MARKER_NAME}_thresh_dict.pkl"
    if [ -f "$THRESHOLD_FILE" ]; then
        echo "  ✓ Found threshold file for $MARKER_NAME"
        ((THRESHOLD_FILES_FOUND++))
    else
        echo "  ✗ Missing threshold file for $MARKER_NAME: $THRESHOLD_FILE"
    fi
done

if [ $THRESHOLD_FILES_FOUND -eq 6 ]; then
    echo "✓ All threshold files found!"
else
    echo "✗ Only $THRESHOLD_FILES_FOUND/6 threshold files found. Check error logs."
    echo "Error logs are in: $SAVE_DIR/*_thresh_log.log"
    exit 1
fi

# Step 2: Normalize the data
echo ""
echo "Step 2: Normalizing the data..."
python normalize.py \
    --CSV "$CSV_FILE" \
    --mfname "$MARKER_FILE" \
    --save_dir "$SAVE_DIR"

# Check if normalized file was created
NORMALIZED_FILE="$SAVE_DIR/sample_cell_data_RESTORE.csv"
if [ -f "$NORMALIZED_FILE" ]; then
    echo "✓ Normalization completed successfully!"
    echo "✓ Normalized data saved to: $NORMALIZED_FILE"
else
    echo "✗ Normalization failed! Output file not found: $NORMALIZED_FILE"
    echo "Check error logs in: $SAVE_DIR/*_norm_log.log"
    exit 1
fi

# Display summary
echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "Completed at $(date)"
echo ""
echo "Output files:"
echo "  - Normalized data: $NORMALIZED_FILE"
echo "  - Threshold plots: $SAVE_DIR/img/threshs/"
echo "  - Threshold data: $SAVE_DIR/thresh_dicts/"
echo "  - Log files: $SAVE_DIR/*_log.log"
echo ""
echo "Next steps:"
echo "  1. Examine the normalized data: head $NORMALIZED_FILE"
echo "  2. Check threshold plots in: $SAVE_DIR/img/threshs/"
echo "  3. Compare before/after statistics"
echo ""
echo "Quick comparison:"
echo "Original data shape: $(head -1 $CSV_FILE | tr ',' '\n' | wc -l) columns, $(tail -n +2 $CSV_FILE | wc -l) rows"
echo "Normalized data shape: $(head -1 $NORMALIZED_FILE | tr ',' '\n' | wc -l) columns, $(tail -n +2 $NORMALIZED_FILE | wc -l) rows"

# Optional: Show first few lines of results
echo ""
echo "First 5 lines of normalized data:"
head -6 "$NORMALIZED_FILE"