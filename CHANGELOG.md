# Changelog

All notable changes to the svMIL project are documented in this file.

## [Current] - 2026-02-17

### Commit: 36f8dc9 - "Updated for use with the Nunes cohort"

This major update adds complete support for the Nunes colorectal cancer cohort dataset and upgrades the pipeline to support hg38 genome reference for colorectal cancer.

---

## Added

### Data Support

#### New Genome Reference (hg38)
- **Primary genome support upgraded to hg38** throughout the pipeline
- Added [data/chromosomes/hg38Coordinates.txt](data/chromosomes/hg38Coordinates.txt) - chromosome coordinate mappings for hg38
- Added [data/genes/allGenesAndIdsHg38.txt](data/genes/allGenesAndIdsHg38.txt) - complete gene annotation for hg38 (78,655 genes)
- Added [data/genes/ensembl_protein_coding_genes.bed](data/genes/ensembl_protein_coding_genes.bed) - protein-coding gene coordinates (20,084 genes)
- Added [data/genes/CCGC_hg38.tsv](data/genes/CCGC_hg38.tsv) and [data/genes/Cosmic_CancerGeneCensus_v103_GRCh38.tsv](data/genes/Cosmic_CancerGeneCensus_v103_GRCh38.tsv) - cancer gene catalogs for hg38

#### Epigenomic Data (hg38-lifted)
- **CTCF binding sites**: [data/ctcf/coad/ENCFF230HWI_hg38.bed](data/ctcf/coad/ENCFF230HWI_hg38.bed) - 16,610 sites for colorectal tissue
- **DNase I hypersensitivity**: [data/dnase/coad/ENCFF882CST_hg38.bed](data/dnase/coad/ENCFF882CST_hg38.bed) - 50,682 open chromatin regions
- **H3K27ac (active enhancers/promoters)**: [data/h3k27ac/coad/ENCFF389IUW_h3k27ac_hg38.bed](data/h3k27ac/coad/ENCFF389IUW_h3k27ac_hg38.bed) - 110,953 peaks (85,817 clustered)
- **H3K27me3 (repressive marks)**: [data/h3k27me3/coad/ENCFF781AEI_h3k27me3_hg38.bed](data/h3k27me3/coad/ENCFF781AEI_h3k27me3_hg38.bed) - 109,957 peaks (79,617 clustered)
- **H3K4me1 (enhancer marks)**: [data/h3k4me1/coad/ENCFF005BSX_h3k4me1_hg38.bed](data/h3k4me1/coad/ENCFF005BSX_h3k4me1_hg38.bed) - 200,411 peaks (170,037 clustered)
- **H3K4me3 (promoter marks)**: [data/h3k4me3/coad/ENCFF893JJR_h3k4me3_hg38.bed](data/h3k4me3/coad/ENCFF893JJR_h3k4me3_hg38.bed) - 38,922 peaks (36,135 clustered)
- **RNA Pol II binding**: [data/rnapol/coad/ENCFF207KQM_hg38.bed](data/rnapol/coad/ENCFF207KQM_hg38.bed) - 4,685 sites
- **Transcription factors**: [data/tf/tf_experimentallyValidated_hg38.bed](data/tf/tf_experimentallyValidated_hg38.bed) - 1,429,312 TF binding sites (504,009 clustered)

#### Regulatory Elements
- **Enhancers**: [data/enhancers/coad/encoderoadmap_elasticnet.104.hg38_lifted.csv](data/enhancers/coad/encoderoadmap_elasticnet.104.hg38_lifted.csv) - 11,757 predicted enhancers
- **Promoters**: [data/promoters/human_epdnew_xsvji.bed](data/promoters/human_epdnew_xsvji.bed) - 21,071 promoter regions
- **Super enhancers**: [data/superEnhancers/coad/sigmoid_colon_hg38.bed](data/superEnhancers/coad/sigmoid_colon_hg38.bed) - 1,022 super enhancer regions
- **CpG islands**: [data/cpg/cpgIslandExt_hg38.txt](data/cpg/cpgIslandExt_hg38.txt) - 32,038 CpG island annotations
- **ChromHMM states**: [data/chromhmm/Sigmoid_Colon/E106_18_core_K27ac_hg38lift_dense.bed](data/chromhmm/Sigmoid_Colon/E106_18_core_K27ac_hg38lift_dense.bed) and mnemonics versions - 513,769+ chromatin state annotations

#### Metadata
- Added [metadata/metadata_nunes.tsv](metadata/metadata_nunes.tsv) - sample metadata for 821 Nunes cohort MSS samples
- Added [src/test/metadata/metadata_nunes_test.tsv](src/test/metadata/metadata_nunes_test.tsv) - test dataset metadata (11 samples)

#### Configuration
- Added [settings/settings_Colorectal/settings.py](settings/settings_Colorectal/settings.py) - complete settings configuration for colorectal cancer analysis with Nunes cohort

### Code Features

#### New SV Parsing Functions
- **`getSVs_nunes()`** in [src/linkSVsGenes/inputParser.py](src/linkSVsGenes/inputParser.py) - Parse structural variants from Nunes cohort VCF files
  - Reads metadata to identify samples by cancer type
  - Validates expression data availability before processing
  - Handles glob patterns: `*_sv_*{sampleId}-T*.vcf.gz`
  - Filters out `.SV.vcf.gz` files to avoid duplicates

- **`getSVsFromFile_nunes_single()`** in [src/linkSVsGenes/inputParser.py](src/linkSVsGenes/inputParser.py) - Comprehensive VCF breakend notation parser
  - Extracts SVCLASS from INFO field (deletion, tandem-duplication, inversion, translocation)
  - Maps SVCLASS to standard SV types: DEL, DUP, INV, ITX
  - Intelligent chromosome ordering for inter-chromosomal events
    - Handles numeric chromosomes (1-22) by numeric value
    - Special handling for sex chromosomes (X, Y) and mitochondrial (MT)
    - Maintains consistent coordinate ordering
  - Adds 'chr' prefix if missing from chromosome names
  - Prevents duplicate SV entries from paired breakends
  - Returns list of SV objects with normalized coordinates

#### CNV Analysis for Nunes Data
- **`getPatientsWithCNVGeneBased_nunes()`** in [src/tadDisruptionsZScores/determinePatientGeneMutationPairs.py](src/tadDisruptionsZScores/determinePatientGeneMutationPairs.py)
  - Reads gene-level copy number data from `*{sampleId}-T*_gene_copynumber.tsv` files
  - Classifies amplifications (CN > 2.3) and deletions (CN < 1.7)
  - Filters neutral copy number events (1.7 ≤ CN ≤ 2.3)
  - Returns per-patient gene amplification and deletion dictionaries

#### SNV Parsing for Nunes Data
- **`getPatientsWithSVs_nunes()`** in [src/tadDisruptionsZScores/determinePatientGeneMutationPairs.py](src/tadDisruptionsZScores/determinePatientGeneMutationPairs.py)
  - Parses SV data to identify genes directly overlapped by structural variants
  - Categorizes overlaps by SV type (DEL, DUP, INV, ITX)
  - For intrachromosomal SVs: identifies all genes between breakpoints
  - For interchromosomal SVs: identifies genes containing breakpoints
  - Validates expression data availability
  - Returns per-patient dictionaries for each SV type

#### Pipeline Integration
- **Nunes data source support** in [src/linkSVsGenes/main.py](src/linkSVsGenes/main.py)
  - Added conditional branch: `if settings.general['source'] == 'Nunes'`
  - Calls `InputParser().getSVs_nunes()` for Nunes cohort processing
  - Integrated with existing PCAWG and HMF pipelines

- **TAD disruption analysis** in [src/tadDisruptionsZScores/computeZScoresDisruptedTads.py](src/tadDisruptionsZScores/computeZScoresDisruptedTads.py)
  - Added Nunes data source handling
  - Calls appropriate SV parsing function based on source

- **Patient-gene mutation pairs** in [src/tadDisruptionsZScores/determinePatientGeneMutationPairs.py](src/tadDisruptionsZScores/determinePatientGeneMutationPairs.py)
  - Full Nunes cohort integration for CNV, SNV, and SV analysis
  - Cross-references mutation data with sample metadata

### Documentation
- **Added [README_output.md](README_output.md)** - Comprehensive documentation of output files
  - Detailed description of all output directory structures
  - Column-by-column explanations for each output file type
  - Example data rows with interpretations
  - Notes on data processing and filtering criteria
  - Covers:
    - TAD disruption z-scores and statistical significance
    - Positive/negative sample sets for disrupted TADs
    - Feature matrices for machine learning
    - Cross-validation results
    - ROC curve data
    - Z-score distribution plots

---

## Changed

### Data Processing

#### Input Parsing Improvements
- **Line-by-line whitespace handling** in [src/linkSVsGenes/inputParser.py](src/linkSVsGenes/inputParser.py)
  - Added `.strip()` to gene file parsing (line 505) to handle trailing whitespace
  - Handles missing or 'NA' values in start/end positions (line 531)
  - More robust parsing with empty string checks

#### Code Cleanup and Optimization
- **Removed redundant header skipping** across multiple parsers in [src/linkSVsGenes/inputParser.py](src/linkSVsGenes/inputParser.py):
  - `readNonCausalGeneFile()` - removed line counter, direct parsing
  - `readTADFile()` - removed 2-line header skip and counter
  - `getCTCFSites()` - removed header skip logic
  - `geteQTLs()` - removed header skip logic
  - `getEnhancers()` - removed header skip logic
  - `getPromoters()` - removed header skip logic
  - `getCpGIslands()` - removed header skip logic
  - `getTFs()` - removed header skip logic
  - `getHistones()` - removed header skip logic
  - `getDnaseISites()` - removed header skip logic
  - `getRnaPolSites()` - removed header skip logic
  
  **Impact**: Cleaner code with ~22 fewer lines across parsers, assumes BED files have no headers (standard format)

- **Removed unused variables**:
  - `nonCausalGeneNameDict` removed from `readNonCausalGeneFile()`
  - Multiple `lineCount` variables removed (11 instances)

### Script Updates

#### Preprocessing Pipeline
- **Modified [src/preprocess.sh](src/preprocess.sh)** - Updated shell script for new data paths and hg38 processing

#### Settings File Generation
- **Updated [src/createSettingsFiles.py](src/createSettingsFiles.py)** - Added Nunes cohort configuration templates

### Machine Learning Pipeline

#### Bug Fix in Bag Filtering
- **Fixed empty bag handling** in [src/multipleInstanceLearning/generateSimilarityMatrices.py](src/multipleInstanceLearning/generateSimilarityMatrices.py)
  - **Issue**: After variance filtering, some bags were left empty, causing numpy array creation errors
  - **Solution**: Added explicit filtering to remove empty bags before array conversion
  - Filter applied to both positive and negative bags (lines 302-317)
  - Also filters corresponding bag pair names to maintain alignment
  - Prevents downstream errors in similarity matrix computation

#### Gene Ranking
- **Updated [src/linkSVsGenes/geneRanking.py](src/linkSVsGenes/geneRanking.py)** - Enhanced ranking algorithm with 99 line changes to accommodate new features

#### Generic Clustering
- **Modified [src/DataProcessing/genericClustering.py](src/DataProcessing/genericClustering.py)** - 76 line changes for improved clustering performance

#### Derivative TAD Analysis
- **Added [src/linkSVsGenes/derivativeTADMaker.py](src/linkSVsGenes/derivativeTADMaker.py)** - New module (73 lines) for analyzing derivative TAD structures formed by SVs

---

## Summary Statistics

**Files Changed**: 60 files  
**Insertions**: +5,037,918 lines  
**Deletions**: -1,000,346 lines  
**Net Change**: +4,037,572 lines

**Key Additions by Category**:
- **Epigenomic data files**: ~3.5M lines (clustered histone marks, TF sites, ChromHMM states)
- **Gene annotations**: ~120,000 lines (hg38 gene lists, cancer gene catalogs)
- **Code**: ~800 lines (new parsers, bug fixes, optimizations)
- **Documentation**: ~400 lines (README_output.md)
- **Metadata**: ~833 lines (sample metadata, test data)

---

## Technical Notes

### Coordinate System
- All new data uses **1-based coordinates** (BED format standard)
- Chromosome naming: 'chr' prefix enforced throughout pipeline
- Inter-chromosomal events: chromosomes ordered numerically, with X/Y/MT handled specially

### Data Quality
- Unmapped regions tracked in `.unmapped` files during liftover
- Clustered files reduce redundancy in overlapping peaks
- Expression data validation prevents processing incomplete samples

### Compatibility
- Maintains backward compatibility with PCAWG and HMF data sources
- Source selection via `settings.general['source']` parameter
- Supported values: 'PCAWG', 'HMF', 'HMF_simple', 'Nunes'

---

## Previous Version

### Commit: b17f9b3 - "readme clarification"
- Minor documentation updates
- No functional code changes
