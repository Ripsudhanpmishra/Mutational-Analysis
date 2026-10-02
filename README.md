# GPCR Signaling Pathway Mutational Analysis

Computational pipeline for genomic variant analysis across a diabetes/metabolic-focused **GPCR signaling gene panel** — spanning G-protein-coupled receptors (GPCRs), their G-protein / downstream signaling / effector / transcription-factor partners, and GPCR accessory proteins (GRKs, arrestins, RAMPs, MRAPs).

Developed during the **RISE-UP Summer Internship, IIT Jammu** (May–July 2026), under the guidance of **Dr. Mithu Baidya**.

> 📄 A full, script-level methodology write-up is available in [`docs/Methodology.docx`](docs/Methodology.docx) (or wherever you place it in this repo) — this README is a navigable summary of that document.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Pipeline](#pipeline)
  - [1. GPCRs](#1-gpcrs)
  - [2. G-Proteins / Downstream Signaling / Effectors / Transcription Factors](#2-g-proteins--downstream-signaling--effectors--transcription-factors)
  - [3. Accessory Proteins](#3-accessory-proteins)
- [Data Sources](#data-sources)
- [Requirements](#requirements)
- [Usage](#usage)
- [Scripts Reference](#scripts-reference)
- [Manual / Non-Scripted Steps](#manual--non-scripted-steps)
- [Known Limitations](#known-limitations)
- [Acknowledgments](#acknowledgments)
- [Author](#author)
- [License](#license)

---

## Overview

This project screens a 30+ gene panel for coding variants relevant to GPCR-mediated metabolic signaling, using population-scale variant data (UK Biobank / T2DKP), GRCh37 coordinates, and Ensembl VEP annotation. Variants are structurally contextualized using **GPCRdb** — Ballesteros–Weinstein (BW) numbering, structural domain/region assignment, and G-protein coupling — then filtered down to functionally important positions and visualized.

The panel is organized into three tiers, each processed through a shared extraction → annotation → merge backbone, then tier-specific structural mapping and visualization:

| Tier | Genes | Structural context |
|---|---|---|
| **GPCRs** | Receptor panel (e.g. GLP1R, GLP2R, MC4R, SCTR, FFAR1–3, LPAR1/3/5/6, S1PR1/2/5, GPBAR, GPR142, GPR183, GPR75, ADRA1A, ADRB2/3, AGTR1, CASR, CNR1, GHSR, MTNR1B, NPY1R, PTGER3, SSTR2, SUCNR1, …) | BW numbering, GPCRdb domains/regions, G-protein coupling, experimental structure coverage |
| **G-proteins / downstream / effectors / TFs** | Gα subunits (GNA11, GNAI1, GNAO1, GNAQ, GNAS, GNA13, GNA12) and related signaling/effector/TF genes | Mutation type and burden per gene |
| **Accessory proteins** | GRKs, β-arrestins (ARRBs), RAMPs, MRAPs | Total + type-wise mutation counts, family-specific domain-wise analysis |


---

## Pipeline

### 1. GPCRs

1. **Per-gene extraction** — Bash/AWK scripts pull each receptor's variant records into a VCF and a parallel CSV from the source dataset (GRCh37, PLINK chromosome convention handled).
2. **Annotation** — Ensembl VEP (`--database` mode) annotates consequence, transcript/protein position, amino-acid change, canonical-transcript flag.
3. **Merge** — a generalized Python merge script combines annotated per-gene VCFs, with a ±3 bp indel-tolerance fallback for coordinate mismatches between callers.
4. **Visualization (position-level)** — lollipop and circular/radial bar plots of mutation distribution per receptor.
5. **Snake plot annotation** — GPCRdb snake plots per receptor, manually annotated with mutation positions.
6. **Domain-wise mutation CSV** — mutations mapped to GPCRdb structural domains (TM helices, loops, termini, H8).
7. **BW numbering** — `bw_classifier_v3.py` assigns Ballesteros–Weinstein position, region, and class to every residue, across both legacy and current GPCRdb export formats.
8. **G-protein coupling** — the same script attaches a curated, pharma-database-cross-verified G-protein coupling call per receptor.
9. **Compilation + functional-site filtering** — all per-receptor layers are merged into one master file, then manually filtered in Excel for mutations at functionally critical positions (BW motifs, interface contacts, binding-pocket residues).
10. **Structure coverage** — experimentally determined GPCR/G-protein structures identified by manually cross-referencing GPCRdb/PDB listings, then visualized.

### 2. G-Proteins / Downstream Signaling / Effectors / Transcription Factors

1. Per-gene VCF/CSV extraction (same Bash/AWK backbone).
2. VEP annotation + Python merge (Gα subunits use the generalized merge pipeline with indel tolerance).
3. Visualization of mutation type and burden per gene (consequence-type breakdowns, per-gene comparisons).

### 3. Accessory Proteins

1. Per-gene VCF/CSV extraction for GRKs, arrestins, RAMPs, MRAPs.
2. VEP annotation + Python merge.
3. Total and type-wise mutation compilation per protein, visualized comparatively across the panel.
4. **Family-wise domain analysis** — a separate script per family (GRK / arrestin / RAMP / MRAP) maps mutations onto family-specific functional domains.

---

## Data Sources

- **UK Biobank** and **T2DKP** — population-scale variant/association data (GRCh37)
- **Ensembl VEP** — variant consequence annotation
- **GPCRdb** — BW numbering tables, structural domain annotations, snake plots, structure listings
- **PDB** — experimental structure cross-reference
- A pharmacological interaction database (cite the specific source/version used) — cross-verification of G-protein coupling assignments

---

## Requirements

```
# Core
pandas
numpy
matplotlib
plotly
openpyxl

# Annotation
Ensembl VEP (with --database connectivity, or a local cache — see Known Limitations)
```

Bash/AWK extraction scripts require a standard POSIX shell environment (tested on WSL/Linux).

---

## Usage

```bash
# 1. Extract a gene's variants (edit CHR / coordinate window inside the script, or parameterize it)
bash scripts/extraction/extract_gene_vcf.sh
bash scripts/extraction/extract_gene_csv.sh

# 2. Annotate with VEP
vep -i data/vcf/<gene>_grch37_targets.vcf \
    -o data/annotated/<gene>_annotated.txt \
    --database --assembly GRCh37 --protein --symbol --canonical --tab

# 3. Merge annotation with association/summary data, filter to coding variants
python scripts/merge/merge_annotated_vcfs.py

# 4. Structural annotation (BW numbering + G-protein coupling)
python scripts/structural/bw_classifier_v3.py <variants_folder> <residue_table.xlsx> <output_folder>

# 5. Visualize
python scripts/visualization/lollipop_plot.py
python scripts/visualization/circular_radial_plot.py
python scripts/visualization/domain_mapper.py
```

---

## Scripts Reference

| Script | Stage | Tier |
|---|---|---|
| `extract_gene_vcf.sh` | Extraction (Bash/AWK) | All |
| `extract_gene_csv.sh` | Extraction (Bash/AWK) | All |
| VEP command (`--database` mode) | Annotation | All |
| `merge_annotated_vcfs.py` | Merge + coding-variant filter (±3 bp indel tolerance) | All |
| `bw_classifier_v3.py` | BW numbering + G-protein coupling annotation | GPCRs |
| `domain_mapper.py` | Domain-wise mutation distribution visualization | GPCRs |
| `lollipop_plot.py` | Lollipop plot (allele frequency × position, domain track) | GPCRs |
| `circular_radial_plot.py` | Circular bar plot (grouped by consequence type) | GPCRs |
| `mutation_spectrum_viz.py` | Mutation type/count summary charts | G-proteins/downstream/effectors/TFs |
| `accessory_mutation_summary.py` | Total + type-wise mutation compilation | Accessory proteins |
| `domain_analysis_grk.py` / `_arrestin.py` / `_ramp.py` / `_mrap.py` | Family-wise domain analysis | Accessory proteins |

---

## Manual / Non-Scripted Steps

A few steps in this pipeline were performed manually rather than by script — documented here for reproducibility:

- **Snake plot annotation** — GPCRdb snake plot images annotated by hand with mutation positions.
- **Functional-site filtering** — the compiled master GPCR file was filtered for mutations at functionally important positions using Excel filters/formulas, not a script.
- **Structure lookup** — experimentally determined structure status was determined by manually cross-referencing GPCRdb/PDB listings, not by script.

If you want these fully reproducible from the repo, consider adding the raw Excel filter criteria or a short written protocol alongside the corresponding data files.

---

## Known Limitations

- VEP annotation was run in `--database` mode (live connection to Ensembl's remote DB), not offline cache mode — note this if re-running, since it requires network access and a working Perl `DBD::mysql` setup.
- The annotation/merge step matches variants to VEP's `Location` field by exact string; indel positions reported as a range by VEP can fail this match and be silently dropped — spot-check indel counts pre/post-merge.
- `DOMAIN_PRESETS` in `lollipop_plot.py` only covers a subset of panel receptors; any gene without a preset renders with a plain backbone and no domain track.
- Color-coding for `splice_region_variant` differs between `circular_radial_plot.py` and `lollipop_plot.py` — reconcile if this consequence type appears in your data.

---

## Acknowledgments

- Dr. Mithu Baidya — research guide, IIT Jammu
- GPCRdb, for structural numbering, domain, and coupling reference data
- Mahajan et al. 2022 (*Nature Genetics*) and Amisten et al. GPCR atlas papers — background literature informing panel design

---

## Author

**Ripsudhan**
B.Tech Biotechnology, NIT Durgapur

---

## License

_Add a license (e.g., MIT) if you intend this repository to be reused — none specified yet._
