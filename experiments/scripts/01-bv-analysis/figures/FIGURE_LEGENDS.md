# Figure Legends - BV Compositional Analysis

## Figure 1: Compositional Closure Demonstration

**File**: `fig1_compositional_closure.png/pdf`

**Panel A**: Distribution of effect sizes (log2 fold change, healthy vs BV) for all 55 taxa tested. Blue histogram shows Hurdle model (count-based) estimates; purple shows LinDA (CLR-based) estimates. Vertical dashed line indicates zero (no change).

**Panel B**: Sum of all effect sizes across taxa. Hurdle model sum = -30.45, indicating net lower abundance in healthy samples. LinDA (CLR) sum = 0.0000 exactly, demonstrating compositional closure—a mathematical property of the CLR transformation that forces increases and decreases to balance. The red horizontal line and annotation highlight that compositional closure is not a biological finding but a mathematical constraint.

---

## Figure 2: Sensitivity Analysis Heatmap

**File**: `fig2_sensitivity_heatmap.png/pdf`

Effect of total bacterial load assumptions on biological interpretation. Rows show key taxa; columns show assumed BV load relative to healthy (1× = equal load, 20× = BV has 20-fold higher total bacteria). Cell values are load-corrected log2 fold changes calculated as: Observed_CLR + log2(Load_BV/Load_Healthy).

Color scale: Blue indicates higher in healthy; red indicates higher in BV. The yellow box highlights taxa whose interpretation changes substantially at 10× load (plausible for BV biofilm). At 10× load, classic "BV pathogens" like Prevotella (-0.5) and Megasphaera (-0.3) show near-zero absolute change, suggesting their apparent increase is a compositional artifact of Lactobacillus dilution rather than genuine pathogen bloom.

---

## Figure 3: Method Comparison

**File**: `fig3_method_comparison.png/pdf`

**Panel A**: Scatter plot comparing effect size estimates between Hurdle (count-based, x-axis) and LinDA (CLR-based, y-axis) models. Points colored by significance status: green = both significant at q<0.05; blue = Hurdle only; purple = LinDA only; gray = neither. Diagonal dashed line indicates 1:1 correspondence. Key taxa (L. crispatus, Prevotella, Gardnerella, L. iners) are labeled. The moderate correlation but substantial scatter demonstrates that method choice substantially affects individual taxon estimates.

**Panel B**: Significance agreement between methods. Bar chart showing number of taxa in each category. LinDA identifies more significant taxa (47) than Hurdle (28), partially due to compositional coupling inflating apparent effects when one dominant taxon changes.

---

## Figure 4: Forest Plot of Top Differential Taxa

**File**: `fig4_forest_plot.png/pdf`

Forest plot showing the 20 taxa with largest absolute effect sizes from LinDA/CLR analysis. Horizontal bars indicate log2 fold change (healthy vs BV) with 95% confidence intervals. Green points indicate taxa higher in healthy samples; red points indicate taxa higher in BV.

The largest effects are Lactobacillus species (positive, higher in healthy) and known BV-associated bacteria like Prevotella, Megasphaera, Dialister, and Sneathia (negative, higher in BV). Gray shading indicates |log2FC| < 1 zone where effects are modest.

---

## Figure 5: Conceptual Load Scenarios

**File**: `fig5_load_scenarios.png/pdf`

Schematic illustrating two scenarios that produce identical relative abundance profiles but different absolute abundance patterns.

**Left panel (Healthy)**: Reference state with ~80% Lactobacillus dominance at 100 units total load.

**Center panel (Scenario A)**: BV with constant total load. Lactobacillus dies and pathogens replace it. True Lactobacillus depletion with pathogen increase.

**Right panel (Scenario B)**: BV with 10× increased total load (biofilm). Lactobacillus absolute abundance unchanged (hatched baseline shows original level); pathogens bloom on top (cross-hatched bars show additional growth). Lactobacillus appears depleted only because it's diluted by the pathogen bloom.

**Key message** (bottom text): Both scenarios produce IDENTICAL relative abundance profiles. Standard 16S analysis cannot distinguish between them without absolute quantification.

---

## Figure 6: Load Correction Waterfall

**File**: `fig6_load_correction_waterfall.png/pdf`

Paired bar chart comparing observed CLR effect sizes (purple) with load-corrected estimates assuming 10× higher load in BV (orange). Gray shaded region indicates the "no meaningful effect" zone (|log2FC| < 1). Red asterisks mark taxa where load correction changes the biological interpretation (effect crosses from outside to inside the neutral zone, or vice versa).

Red curved arrows connect paired bars, visualizing the direction and magnitude of interpretation change. The key insight is that Prevotella, Megasphaera, Sneathia, and Atopobium—classic "BV pathogens"—all shift from clearly elevated in BV to near-neutral when load is considered, while Gardnerella remains elevated even after correction.

---

## Summary Figure

**File**: `fig_summary.png/pdf`

Three-panel summary combining key findings:

**Panel A**: Compositional closure demonstration (simplified from Fig. 1B).

**Panel B**: Key differential taxa showing the five taxa with largest effects.

**Panel C**: Load-corrected interpretation heatmap (same as Fig. 2 but with abbreviated taxon names).

This figure is suitable for graphical abstract or presentations.

---

## Technical Notes

- All figures generated at 300 DPI for publication quality
- PDF versions available for vector graphics (scalable)
- Color scheme: Green (#2ecc71) = healthy/protective; Red (#e74c3c) = BV/pathogenic; Blue (#3498db) = Hurdle model; Purple (#9b59b6) = LinDA model; Gray (#95a5a6) = neutral/non-significant
- Generated by `generate_figures.py` using matplotlib/seaborn
