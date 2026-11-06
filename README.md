# Computational Biology Research Portfolio

This repository presents five reproducible computational biology projects spanning single-cell transcriptomics, spatial biology, cancer genomics, metagenomics and machine learning for cancer pharmacogenomics.

The analyses were developed by **Sriram Gummadi**, a bioinformatics researcher at La Trobe University whose PhD thesis has been submitted and is currently under examination. The work reflects a broader interest in combining biological interpretation with robust statistical modelling, reproducible computation and scientific visualisation.

## Research projects

### Single-Cell RNA-seq of Human Lung Disease

Donor-aware analysis of the GSE136831 human lung atlas, examining cellular composition and transcriptional remodelling in idiopathic pulmonary fibrosis and chronic obstructive pulmonary disease. The workflow includes pseudobulk differential expression, pathway enrichment, macrophage subgroup analysis, ATII heterogeneity and robustness diagnostics.

**Directory:** [`Single-Cell-RNAseq-Lung-Disease`](Single-Cell-RNAseq-Lung-Disease/)

### Spatial Transcriptomics of Breast Cancer

Spatial analysis of breast tumour tissue integrating histopathology, spatial gene expression, clustering, tumour microenvironment programmes, spatially variable genes and pathway enrichment. The code distinguishes descriptive spatial structure from patient-level biological inference.

**Directory:** [`Spatial-Transcriptomics-Breast-Cancer`](Spatial-Transcriptomics-Breast-Cancer/)

### Pan-Cancer Transcriptomic Survival Analysis

Cancer-specific and pan-cancer analysis of TCGA transcriptomic and survival data. The workflow includes cohort construction, expression-integrity checks, Cox proportional-hazards modelling, effect heterogeneity, clinical adjustment, proportional-hazards diagnostics, internal validation and pathway analysis.

**Directory:** [`Pan-Cancer-Transcriptomic-Survival`](Pan-Cancer-Transcriptomic-Survival/)

### Shotgun Metagenomics of Colorectal Cancer

Multi-cohort analysis of harmonised gut shotgun-metagenomic profiles. The workflow evaluates microbial diversity, compositional associations, cross-project synthesis, predicted functional potential, healthy-adenoma-CRC contrasts and independent external validation.

**Directory:** [`Shotgun-Metagenomics-Colorectal-Cancer`](Shotgun-Metagenomics-Colorectal-Cancer/)

### Explainable AI for Cancer Drug Response

Multi-omics modelling of cancer drug sensitivity using RNA expression, proteomics, copy-number and mutation features. The workflow includes leakage-conscious prediction, nested validation, SHAP-based model interpretation, cross-platform assessment and targeted MAPK-response analyses.

**Directory:** [`Explainable-AI-Cancer-Drug-Response`](Explainable-AI-Cancer-Drug-Response/)

## Reproducibility and interpretation

- Clinical samples or participants, rather than individual cells or measurements, are used as the biological replication unit where appropriate.
- Feature filtering, preprocessing and model selection are kept within training data for predictive analyses where supported by the original design.
- Multiple-testing correction, model diagnostics and sensitivity analyses are retained throughout the workflows.
- Predictive associations and model explanations are not presented as evidence of biological causality.
- Large source datasets, intermediate results and generated figures are intentionally excluded from this repository.

The scripts pass static syntax validation. Complete numerical reproduction requires the corresponding source datasets, documented software dependencies and an end-to-end execution of each workflow.

## Installation, execution and tests

Python and R dependencies are documented in [`DEPENDENCIES.md`](DEPENDENCIES.md), with package lists in [`requirements-python.txt`](requirements-python.txt) and [`requirements-r.txt`](requirements-r.txt).

Each project now has an explicit pipeline runner. These runners stop on the first failed script so downstream analyses cannot continue from incomplete prerequisites. Project 5 intentionally runs `03_rna_models.py` before `02_response_landscape.py` because the latter can consume the held-out prediction table produced by the RNA-model script; existing filenames are preserved for compatibility.

Use dry-run mode to inspect Python execution order without running analyses:

```text
python Spatial-Transcriptomics-Breast-Cancer/run_pipeline.py --dry-run
python Explainable-AI-Cancer-Drug-Response/run_pipeline.py --dry-run
```

Run lightweight automated checks with:

```text
python -m unittest discover -s tests -v
```

The tests validate interfaces and small synthetic fixtures only. They do not rerun scientific models or alter archived results.

## Languages and methods

**Languages:** R and Python

**Methods:** single-cell RNA sequencing, spatial transcriptomics, pseudobulk differential expression, survival analysis, pathway enrichment, microbial community analysis, multi-omics integration, machine learning, nested cross-validation and explainable AI.

## Research profiles

- [Google Scholar](https://scholar.google.com/citations?user=HEPzv0gAAAAJ&hl=en)
- [ORCID](https://orcid.org/0009-0002-5609-5371)
- [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Sriram+Gummadi%5BAuthor%5D)
- [ResearchGate](https://www.researchgate.net/profile/Sriram-Gummadi)
- [LinkedIn](https://au.linkedin.com/in/sriramgummadi)
- [Email](mailto:S.Gummadi@latrobe.edu.au)

## License and reuse

The analysis code is provided for research and portfolio review. Source datasets remain subject to their original licenses, access conditions and citation requirements.
