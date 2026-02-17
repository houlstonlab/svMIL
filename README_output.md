**GENERATED USING AI - I HAVE CHECKED THAT IT SEEMS CORRECT, BUT DOUBLE CHECK ANYTHING YOU USE**
# svMIL Output Files Documentation

This directory contains the output files from the svMIL (Structural Variant Multiple Instance Learning) pipeline run on the Nunes colorectal cancer dataset. The pipeline identifies pathogenic structural variants (SVs) that disrupt topologically associating domains (TADs) and affect gene expression.

## Directory Structure

```
output_colorectal/
├── tadDisruptionsZScores/          # TAD disruption analysis and z-score calculations
├── linkedSVGenePairs/              # SV-gene pairs with computed features
├── patientGeneMutationPairs/       # Patient-specific mutation annotations
├── multipleInstanceLearning/       # Machine learning classification results
│   ├── similarityMatrices/         # Feature matrices and similarity calculations
│   └── leaveOnePatientOutCV/       # Cross-validation results per SV type
└── rocCurves/                      # ROC curve data for model performance
```

---

## 1. tadDisruptionsZScores/

This directory contains the statistical analysis of TAD disruptions and their impact on gene expression.

### zScores.txt

Contains z-scores for each patient-gene pair where the gene is in a TAD disrupted by a structural variant. Z-scores quantify differential expression between patients with and without TAD disruptions.

**Columns:**
1. **Patient_Gene**: Patient ID and gene symbol (format: `U####_GENE` or `UM###_GENE`)
2. **P-value**: Raw p-value from two-tailed test (assuming normal distribution: `stats.norm.sf(abs(z))*2`)
3. **Adjusted P-value**: Bonferroni corrected p-value for multiple testing correction
4. **Significant**: Boolean indicating if Bonferroni-adjusted p-value < 0.05 (True/False)
5. **P-value (duplicate)**: Same as column 2
6. **Z-score**: Standardized effect size calculated as `(patient_expression - mean_control) / std_control` (positive = upregulated, negative = downregulated)

**Notes:**
- `nan` values indicate genes without sufficient data for statistical analysis (e.g., no expression, insufficient samples)
- Pairs with `|z-score| > 1.5` are considered potentially pathogenic (positive training examples)
- Pairs excluded if they have coding mutations (SNVs/CNVs) to isolate non-coding effects

**Example:**
```
U0001_MAP2K7    1.0167849953292189e-05  0.0009252743457495891   True    1.0167849953292189e-05  -4.413572811267576
```
This shows patient U0001 has TAD disruption near MAP2K7 with significant downregulation (z-score = -4.41, p_adj < 0.001).

---

### tadPositiveAndNegativeSet.txt

Lists TADs and the patients affected (positive set) or unaffected (negative set) by SV disruptions.

**Columns:**
1. **TAD_coordinates**: TAD genomic coordinates (format: `chr#_start_end`)
2. **Positive_samples**: List of patient IDs with SVs disrupting this TAD
3. **Negative_samples**: List of patient IDs without disruptions in this TAD
4. **SV_types**: List of SV types causing disruptions (DEL=deletion, DUP=duplication, INV=inversion, ITX=translocation)

**Notes:**
- TADs are considered disrupted if an SV breakpoint crosses the TAD boundary
- SVs fully contained within a TAD are NOT considered disruptive
- Used to compute z-scores by comparing expression between positive and negative sets
- **All TADs from the reference file are included**, even those with no disruptions or with empty sample sets
- TADs with empty positive AND negative sets (`[] [] []`) occur when:
  - **The TAD contains no genes**
  - OR all genes in the TAD are filtered out due to coding mutations in ALL patients (both disrupted and non-disrupted)
  - OR genes in the TAD have no expression data in the cohort
- Patients are only added to negative set if the TAD contains at least one gene without mutations

**Example:**
```
chr1_0_960000   ['U0800', 'U3070']   ['U0008', 'U0012', ...]   ['INV', 'INV', 'ITX']
```
This TAD is disrupted by inversions/translocations in patients U0800 and U3070.

---

## 2. linkedSVGenePairs/

Contains feature matrices for each SV-gene pair used in machine learning classification.

### nonCoding_geneSVPairs.txt_

Master list of all non-coding SV-gene pairs identified in the cohort with their initial feature set.

**Format:** Tab-delimited text with 48 columns

**Columns:**
1. **SV_metadata**: Unique identifier with gene, SV coordinates, patient, and SV type
2-48. **Initial features**: Including eQTL gains/losses, regulatory element overlaps, and basic genomic features

**Notes:**
- Contains header row starting with `# SV_metadata`
- This is the intermediate output before z-score labeling and full feature computation
- Serves as input for generating pathogenic/non-pathogenic feature matrices

---

### nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt

Feature matrix for SV-gene pairs classified as potentially pathogenic based on z-score threshold (`|z| > 1.5`).

**Format:** 81 columns (no header row)

**Columns:**
1. **Pair_ID**: Unique identifier for the SV-gene pair (format: `GENE_chr#_pos1_pos1_chr#_pos2_pos2_PATIENT_...`)
2-81. **Features**: 80 genomic and regulatory features per SV-gene pair:

   **Element Loss/Gain Features (columns 2-47):**
   - For 23 regulatory element types, binary indicators of loss and gain (46 features total)
   - Element types: eQTL, enhancer, promoter, cpg, tf, h3k4me3, h3k27ac, h3k27me3, h3k4me1, dnaseI, rnaPol, CTCF, CTCF+Enhancer, CTCF+Promoter, Enhancer, Heterochromatin, Poised_Promoter, Promoter, Repeat, Repressed, Transcribed, superEnhancer, ctcf
   - Each element has: `{element}_loss` and `{element}_gain`

   **Strength Features (columns 48-61):**
   - For 7 element types, quantitative strength scores of lost/gained elements (14 features total)
   - Element types: enhancer, ctcf, rnaPol, h3k4me3, h3k27ac, h3k27me3, h3k4me1
   - Each element has: `{element}_strength_loss` and `{element}_strength_gain`
   
   **Additional Features (columns 62-81):**
   - Distance metrics, conservation scores, TAD disruption indicators, and other genomic context features (20 features)

**Notes:**
- These pairs have significant differential expression (z-score threshold met)
- Used as positive training examples in the classifier
- Values are normalized to [0-1000] range or binary (0/1)

**Example:**
```
CSF3R_chr1_36626418_36626418_chr1_38247969_38247969_U0800_2_2_1_1_INV_1000.0_1000.0_1000.0_538.0  0.0  0.0  0.0  0.0  1.0  ...
```

---

### nonCoding_geneSVPairs.txt_nonPathogenicPairsFeatures.txt

Feature matrix for SV-gene pairs classified as non-pathogenic (z-score below threshold or not significant).

**Column structure:** Same as pathogenicPairsFeatures.txt (see above)

**Notes:**
- These pairs show no significant differential expression
- Used as negative training examples in the classifier
- Larger in size due to more negative examples in the dataset

---

### normalizedBags.pkl

Python pickle file containing normalized feature bags for Multiple Instance Learning. Each "bag" is a collection of SV-gene pair instances for one patient-gene combination.

**Format:** Binary pickle file (use `pickle.load()` to read)

**Structure:**
- Dictionary with keys = patient-gene pair identifiers
- Values = arrays of normalized feature vectors

---

### bags.pkl

Python pickle file containing raw (un-normalized) feature bags before normalization.

**Format:** Binary pickle file

---

## 3. patientGeneMutationPairs/

Contains dictionaries mapping patients to genes affected by different mutation types. Used to filter out coding mutations from the analysis.

### cnvPatientsAmp.npy

NumPy dictionary file mapping genes to patients with copy number amplifications (copy number > 2.3).

**Format:** Binary NumPy file (use `np.load(..., allow_pickle=True).item()`)

**Structure:**
```python
{
    'GENE1': ['U0001', 'U0002', ...],
    'GENE2': ['U0005', 'U0010', ...],
    ...
}
```

---

### cnvPatientsDel.npy

NumPy dictionary file mapping genes to patients with copy number deletions (copy number < 1.7).

**Format:** Same as cnvPatientsAmp.npy

---

### snvPatients.npy

NumPy dictionary file mapping genes to patients with single nucleotide variants (SNVs) affecting coding regions.

**Format:** Same structure as CNV files

**Notes:**
- Includes all SNV effects annotated by snpEff (missense, nonsense, splice site, etc.)
- Used to exclude genes from z-score analysis if they have coding mutations

---

### svPatientsDel.npy / svPatientsDup.npy / svPatientsInv.npy / svPatientsItx.npy

NumPy dictionary files mapping genes to patients with different structural variant types:
- **svPatientsDel**: Deletions
- **svPatientsDup**: Duplications (tandem duplications)
- **svPatientsInv**: Inversions
- **svPatientsItx**: Translocations (inter-chromosomal)

**Format:** Same structure as CNV/SNV files

---

## 4. multipleInstanceLearning/

Contains machine learning classification results and model performance metrics.

### leaveOnePatientOutCV/

Cross-validation results using leave-one-patient-out (LOPO) strategy, with separate models for each SV type.

#### leaveOnePatientOutCV_{SVTYPE}.txt

Detailed cross-validation results for each SV type (DEL, DUP, INV, ITX).

**Columns:**
1. **Pair_ID**: SV-gene pair identifier (same format as feature files)
2. **True_label**: Ground truth label (0=non-pathogenic, 1=pathogenic)
3. **Predicted_label**: Model prediction after applying optimal threshold

**Notes:**
- Each row represents one test instance from one CV fold
- Model trained on all other patients, tested on held-out patient
- Predictions made using SV type-specific classifier

**Example:**
```
PSMB8_chr6_32670373_32670373_chr6_32806557_32806557_U0136_3_7_7_8_DEL_1000.0_1000.0_1000.0_1000.0   0   1
```
This deletion pair was predicted as pathogenic (1) but is actually non-pathogenic (0) - a false positive.

---

#### leaveOnePatientOutCV_{SVTYPE}_threshold.txt

Contains the optimal classification threshold for each SV type model.

**Format:** Single floating point value per line

**Notes:**
- Threshold optimized to maximize Youden's J statistic (sensitivity + specificity - 1)
- Used to convert continuous prediction scores to binary labels
- Different threshold per SV type due to different class distributions

**Example:**
```
0.8475508351488744
```

---

#### leaveOnePatientOutCV_{SVTYPE}_FINAL_AUC.txt

Performance metrics for the model across all cross-validation folds.

**Format:** 3 rows × 2 columns (tab-delimited)

**Row 1 - ROC AUC Statistics:**
1. **Mean AUC**: Average area under ROC curve across all CV folds
2. **Std Dev AUC**: Standard deviation of AUC across folds

**Row 2 - Precision-Recall Metrics:**
1. **AUPRC**: Area under precision-recall curve
2. **AP**: Average precision score

**Row 3 - Performance at Selected Threshold:**
1. **Precision**: Precision at the optimal threshold (from threshold.txt file)
2. **Recall**: Recall (sensitivity) at the optimal threshold

**Notes:**
- AUC ranges from 0 to 1 (0.5 = random, 1.0 = perfect classifier)
- AUPRC is particularly useful for imbalanced datasets (common in this analysis)
- Average precision summarizes the precision-recall curve
- Precision/Recall values show performance at the threshold optimized for Youden's J statistic

**Example:**
```
0.8475508351488744      0.3100805508133753
0.8580905649906915      0.8585398314552293
0.563302752293578       0.9967532467532467
```
- Mean AUC = 0.848 (±0.310 std dev) indicates good discriminative ability
- AUPRC = 0.858, AP = 0.859 show strong precision-recall performance
- At optimal threshold: Precision = 0.563, Recall = 0.997 (high sensitivity)

---

### similarityMatrices/

Contains pre-computed similarity matrices and filtered features.

#### lowVarianceIdx_{SVTYPE}.txt

Indices of features with zero or near-zero variance across all instances, which were removed before model training.

**Format:** One feature index per line (0-indexed)

**Notes:**
- Features with no variation provide no discriminative information
- Removed to reduce dimensionality and improve model efficiency
- Specific to each SV type due to different feature distributions

**Example:**
```
0.0
1.0
10.0
```
Features at indices 0, 1, and 10 were removed due to low variance.

---

## 5. rocCurves/

Contains data for plotting Receiver Operating Characteristic (ROC) curves to visualize model performance.

**Files:** Typically include true positive rates (TPR), false positive rates (FPR), and thresholds for each SV type model.

---

## Summary Statistics

Based on the current run:

- **Total samples analyzed**: 821 non-hypermutated colorectal cancer patients
- **TADs in reference**: 14,893 total
  - 12,701 with SV disruptions detected (before filtering)
  - 6,058 with disruptions and analyzable genes (in output file)
  - 993 without disruptions but with analyzable genes
  - 7,842 with no analyzable genes (empty in output file)
- **Unique TAD-disrupting SVs**: 374 (112 DEL, 43 DUP, 127 INV, 92 ITX)
- **Patient-gene pairs with z-scores**: 91 (with 11 having insufficient data = NaN)
- **Features per SV-gene pair**: 80 genomic/regulatory features

---

## Usage Notes

1. **Python files (.npy, .pkl)**: Load with NumPy or pickle
   ```python
   import numpy as np
   import pickle
   
   # Load NumPy dictionary
   cnv_amp = np.load('cnvPatientsAmp.npy', allow_pickle=True).item()
   
   # Load pickle file
   with open('normalizedBags.pkl', 'rb') as f:
       bags = pickle.load(f)
   ```

2. **Text files**: Load with pandas or NumPy
   ```python
   import pandas as pd
   import numpy as np
   
   # Load with pandas
   zscores = pd.read_csv('zScores.txt', sep='\t', header=None)
   
   # Load with NumPy
   cv_results = np.loadtxt('leaveOnePatientOutCV_DEL.txt', dtype='object')
   ```

3. **NaN handling**: NaN values in zScores.txt represent genes without valid statistics (no expression data or insufficient samples). These are automatically treated as non-significant by the pipeline.

4. **Model interpretation**: Higher z-scores (>1.5 or <-1.5) indicate genes significantly affected by TAD disruption. The machine learning model learns which genomic features best predict these significant effects.

---

## Questions or Issues

For questions about the pipeline or output files, refer to the main svMIL repository README or contact the pipeline maintainers.

**Pipeline version**: svMIL adapted for Nunes colorectal cancer data (HMF → Nunes)
**Run date**: January 30, 2026
**Dataset**: Nunes CRC (821 non-hypermutated tumor samples)