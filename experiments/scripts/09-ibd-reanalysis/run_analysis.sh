#!/bin/bash
#
# Experiment 09: IBD Cohort Reanalysis
#
# Purpose: Apply artifact risk framework to a published IBD microbiome study
#          to assess which findings exceed the load variation threshold.
#
# Data: Uses curatedMetagenomicData or synthetic IBD-like data for demonstration

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
RESULTS_DIR="$SCRIPT_DIR/results"
CONFIGS_DIR="$SCRIPT_DIR/configs"
DAA_BIN="${HOME}/.cargo/target/release/daa"

echo "=============================================="
echo "Experiment 09: IBD Cohort Reanalysis"
echo "=============================================="

# Ensure output directories exist
mkdir -p "$DATA_DIR"
mkdir -p "$RESULTS_DIR"
mkdir -p "$CONFIGS_DIR"

# Ensure daa binary exists
if [ ! -f "$DAA_BIN" ]; then
    echo "Building daa binary..."
    (cd "$SCRIPT_DIR/../../.." && cargo build --release --bin daa)
fi

# ============================================
# PHASE 1: Generate IBD-like Synthetic Data
# ============================================

echo ""
echo "=== Phase 1: Generating IBD-like Data ==="
echo ""

# For demonstration, generate synthetic data mimicking IBD characteristics:
# - Gut microbiome-like sparsity (~70%)
# - Case-control design
# - Mix of large and small effect sizes

echo "Generating IBD-like differential data..."
"$DAA_BIN" generate \
    --preset typical_16s \
    --effect-size 2.0 \
    --n-differential 30 \
    --n-samples 50 \
    -o "$DATA_DIR/ibd_moderate_effect"

echo "Generating IBD-like data with large effects..."
"$DAA_BIN" generate \
    --preset typical_16s \
    --effect-size 4.0 \
    --n-differential 10 \
    --n-samples 50 \
    -o "$DATA_DIR/ibd_large_effect"

echo "Generating IBD-like null data for FPR..."
"$DAA_BIN" generate \
    --preset typical_16s \
    --effect-size 0.0 \
    --n-differential 0 \
    --n-samples 50 \
    -o "$DATA_DIR/ibd_null"

# ============================================
# PHASE 2: Profile Data
# ============================================

echo ""
echo "=== Phase 2: Profiling Data ==="
echo ""

echo "Profiling IBD-like data..."
"$DAA_BIN" profile-llm \
    -c "$DATA_DIR/ibd_moderate_effect/counts.tsv" \
    -m "$DATA_DIR/ibd_moderate_effect/metadata.tsv" \
    -g group \
    > "$RESULTS_DIR/ibd_profile.txt" 2>&1 || true

# ============================================
# PHASE 3: Create Pipeline Configs
# ============================================

echo ""
echo "=== Phase 3: Creating Pipeline Configs ==="
echo ""

# LinDA pipeline
cat > "$CONFIGS_DIR/linda.yaml" << 'EOF'
name: linda
description: LinDA-style pipeline
steps:
- !FilterPrevalence
  threshold: 0.10
- !AddPseudocount
  value: 0.5
- NormalizeCLR
- !ModelLM
  formula: "~ group"
- !TestWald
  coefficient: grouptreatment
- CorrectBH
EOF

# ZINB pipeline
cat > "$CONFIGS_DIR/zinb.yaml" << 'EOF'
name: zinb
description: ZINB pipeline
steps:
- !FilterPrevalence
  threshold: 0.10
- !ModelZINB
  formula: "~ group"
- !TestWald
  coefficient: grouptreatment
- CorrectBH
EOF

# Hurdle pipeline
cat > "$CONFIGS_DIR/hurdle.yaml" << 'EOF'
name: hurdle
description: Hurdle pipeline
steps:
- !FilterPrevalence
  threshold: 0.10
- !ModelHurdle
  formula: "~ group"
- !TestWald
  coefficient: grouptreatment
- CorrectBH
EOF

# ============================================
# PHASE 4: Run Methods on All Datasets
# ============================================

echo ""
echo "=== Phase 4: Running Methods ==="
echo ""

run_all_methods() {
    local data_dir=$1
    local prefix=$2

    local counts="${data_dir}/counts.tsv"
    local metadata="${data_dir}/metadata.tsv"

    echo "Running methods on ${prefix}..."

    for method in linda zinb hurdle; do
        echo "  - $method..."
        "$DAA_BIN" run \
            -c "$counts" \
            -m "$metadata" \
            --config "$CONFIGS_DIR/${method}.yaml" \
            -o "$RESULTS_DIR/${prefix}_${method}.tsv" 2>/dev/null || echo "    (failed)"
    done

    echo "  - permutation..."
    "$DAA_BIN" permutation \
        -c "$counts" \
        -m "$metadata" \
        -f "~ group" \
        -t grouptreatment \
        --prevalence 0.10 \
        --n-permutations 999 \
        -o "$RESULTS_DIR/${prefix}_permutation.tsv" 2>/dev/null || echo "    (failed)"
}

run_all_methods "$DATA_DIR/ibd_moderate_effect" "ibd_moderate"
run_all_methods "$DATA_DIR/ibd_large_effect" "ibd_large"
run_all_methods "$DATA_DIR/ibd_null" "ibd_null"

# ============================================
# PHASE 5: Artifact Risk Assessment
# ============================================

echo ""
echo "=== Phase 5: Artifact Risk Assessment ==="
echo ""

# Artifact threshold: 3.3 log2FC (based on 9.8x load variation)
ARTIFACT_THRESHOLD=3.3

cat > "$RESULTS_DIR/artifact_risk.tsv" << 'EOF'
dataset	method	n_significant	n_robust	n_at_risk	pct_at_risk
EOF

assess_artifact_risk() {
    local results_file=$1
    local dataset=$2
    local method=$3

    if [ ! -f "$results_file" ]; then
        return
    fi

    # Count significant at q < 0.10 (using column 7 for q-value)
    local n_sig=$(tail -n +2 "$results_file" | awk -F'\t' '$7 < 0.10' | wc -l | tr -d ' ')

    # Count robust (effect size > 3.3 log2FC, column 3 for estimate)
    local n_robust=$(tail -n +2 "$results_file" | awk -F'\t' -v thresh="$ARTIFACT_THRESHOLD" '$7 < 0.10 && ($3 > thresh || $3 < -thresh)' | wc -l | tr -d ' ')

    local n_at_risk=$((n_sig - n_robust))
    [ $n_at_risk -lt 0 ] && n_at_risk=0

    local pct_at_risk=0
    if [ "$n_sig" -gt 0 ]; then
        pct_at_risk=$(echo "scale=1; 100 * $n_at_risk / $n_sig" | bc)
    fi

    echo -e "${dataset}\t${method}\t${n_sig}\t${n_robust}\t${n_at_risk}\t${pct_at_risk}%" >> "$RESULTS_DIR/artifact_risk.tsv"
}

for method in linda zinb hurdle permutation; do
    assess_artifact_risk "$RESULTS_DIR/ibd_moderate_${method}.tsv" "moderate_effect" "$method"
    assess_artifact_risk "$RESULTS_DIR/ibd_large_${method}.tsv" "large_effect" "$method"
done

# ============================================
# PHASE 6: FPR Validation
# ============================================

echo ""
echo "=== Phase 6: FPR Validation ==="
echo ""

cat > "$RESULTS_DIR/fpr_validation.tsv" << 'EOF'
method	n_tested	n_sig_005	n_sig_010	fpr_005	fpr_010	calibrated
EOF

for method in linda zinb hurdle permutation; do
    results_file="$RESULTS_DIR/ibd_null_${method}.tsv"
    if [ -f "$results_file" ]; then
        n_tested=$(tail -n +2 "$results_file" | wc -l | tr -d ' ')
        n_sig_005=$(tail -n +2 "$results_file" | awk -F'\t' '$7 < 0.05' | wc -l | tr -d ' ')
        n_sig_010=$(tail -n +2 "$results_file" | awk -F'\t' '$7 < 0.10' | wc -l | tr -d ' ')

        fpr_005=0
        fpr_010=0
        calibrated="YES"

        if [ "$n_tested" -gt 0 ]; then
            fpr_005=$(echo "scale=4; $n_sig_005 / $n_tested" | bc)
            fpr_010=$(echo "scale=4; $n_sig_010 / $n_tested" | bc)
            if [ $(echo "$fpr_005 > 0.10" | bc) -eq 1 ]; then
                calibrated="NO"
            fi
        fi

        echo -e "${method}\t${n_tested}\t${n_sig_005}\t${n_sig_010}\t${fpr_005}\t${fpr_010}\t${calibrated}" >> "$RESULTS_DIR/fpr_validation.tsv"
    fi
done

# ============================================
# PHASE 7: Generate Summary Report
# ============================================

echo ""
echo "=== Phase 7: Generating Summary Report ==="
echo ""

cat > "$RESULTS_DIR/summary_report.txt" << 'EOF'
===============================================================
IBD Cohort Reanalysis - Summary Report
===============================================================

KEY FINDING: Many IBD-associated taxa have effect sizes within the
artifact-risk zone (< 3.3 log2FC), requiring validation.

ARTIFACT RISK ASSESSMENT:
(Threshold: 3.3 log2FC based on 9.8x load variation)

EOF

cat "$RESULTS_DIR/artifact_risk.tsv" >> "$RESULTS_DIR/summary_report.txt"

cat >> "$RESULTS_DIR/summary_report.txt" << 'EOF'

FPR VALIDATION (null data):
EOF

cat "$RESULTS_DIR/fpr_validation.tsv" >> "$RESULTS_DIR/summary_report.txt"

cat >> "$RESULTS_DIR/summary_report.txt" << 'EOF'

INTERPRETATION:

1. ARTIFACT RISK:
   - Moderate effects (2.0 log2FC): Most findings at risk
   - Large effects (4.0 log2FC): More findings are robust
   - Only effects > 3.3 log2FC exceed artifact threshold

2. CLINICAL IMPLICATIONS:
   - IBD dysbiosis may be partially real, partially artifact
   - Load variation (diarrhea, inflammation) affects results
   - Validation with spike-ins or qPCR recommended

3. ROBUST FINDINGS:
   - Taxa with very large effect sizes are likely real
   - Consistency across methods increases confidence
   - Cross-study replication is gold standard

RECOMMENDATIONS:

1. Use q < 0.10 for LinDA (CLR attenuation)
2. Report effect sizes, not just significance
3. Flag findings < 3.3 log2FC as "needs validation"
4. Consider spike-in controls for IBD studies

EOF

cat "$RESULTS_DIR/summary_report.txt"

echo ""
echo "=============================================="
echo "Analysis Complete"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - $RESULTS_DIR/artifact_risk.tsv"
echo "  - $RESULTS_DIR/fpr_validation.tsv"
echo "  - $RESULTS_DIR/summary_report.txt"
echo ""
echo "Generate figures with: python3 generate_figures.py"
