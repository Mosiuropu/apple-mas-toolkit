# Apple MAS Toolkit

**Marker-Assisted Selection Toolkit for Apple (*Malus domestica* Borkh.) Molecular Breeding**

A comprehensive, peer-reviewed-structured Python toolkit for chromosome-aware marker selection, genotyping primer design, and downstream statistical analysis in apple breeding programmes.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Repository Architecture](#repository-architecture)
3. [Installation](#installation)
4. [Formulating a Breeding Hypothesis](#formulating-a-breeding-hypothesis)
5. [Step-by-Step Research Workflow](#step-by-step-research-workflow)
6. [Marker Database Reference](#marker-database-reference)
7. [Analytics Engine Reference](#analytics-engine-reference)
8. [Apple Breeding Resources](#apple-breeding-resources)
9. [Contributing](#contributing)
10. [Citation](#citation)
11. [License](#license)

---

## Introduction

The Apple MAS Toolkit provides students and researchers with a structured framework for:

- **Mapping agronomic traits** to validated molecular markers (SSR, SCAR, InDel, SNP) across all 17 *Malus* chromosomes.
- **Querying markers** by chromosome, trait, or keyword to design chromosome-specific assay panels.
- **Generating genotyping primer manifests** with forward/reverse primer sequences, expected allele sizes, and PCR conditions.
- **Performing statistical analyses** including segregation distortion tests (chi-square), Polymorphism Information Content (PIC) calculation, and genotype-phenotype association evaluation.
- **Visualising marker data** with publication-quality box-and-whisker plots comparing trait variance across marker genotype classes.

The toolkit aligns with nomenclature and coordinate conventions from the **Genome Database for Rosaceae (GDR)** and references the GDDH13 v1.1 reference genome assembly (Daccord et al., 2014).

---

## Repository Architecture

```
apple-mas-toolkit/
├── README.md                          # This handbook
├── pyproject.toml                     # Package metadata and dependencies
├── LICENSE                            # MIT License
├── CONTRIBUTING.md                    # Contribution guidelines
├── data/
│   └── apple_marker_map.json          # Curated marker database (17 chromosomes)
├── src/
│   ├── marker_selector.py             # Interactive marker query tool
│   └── breeding_analytics.py          # Statistical analysis suite
├── tests/
│   └── test_all.py                    # Test suite
├── docs/
│   ├── getting_started.md             # Installation and first steps
│   └── marker_reference.md            # Detailed marker documentation
└── notebooks/
    └── 01_complete_analysis.ipynb     # Walkthrough notebook
```

---

## Installation

```bash
git clone https://github.com/Mosiuropu/apple-mas-toolkit-hidden.git
cd apple-mas-toolkit
pip install -e .
```

### Requirements

- Python ≥ 3.9
- pandas ≥ 1.5
- numpy ≥ 1.23
- matplotlib ≥ 3.6
- seaborn ≥ 0.12
- scipy ≥ 1.9

---

## Formulating a Breeding Hypothesis

A well-structured breeding hypothesis connects a phenotypic objective to specific molecular markers. Use the following framework:

1. **Identify the target trait.** What phenotype are you selecting for? (e.g., firm fruit with scab resistance)
2. **Select candidate markers.** Which validated markers control or are linked to this trait?
3. **Determine the population.** What type of segregating population will you genotype? (F1, BC1, BC2, F2, RIL)
4. **Predict segregation.** What genotypic ratios do you expect under Mendelian inheritance?
5. **Design the assay.** Which primer pairs, dye channels, and fragment sizes do you need?

**Example hypothesis:**

> "In a BC1 population derived from a cross between a firm, scab-susceptible cultivar and a soft, scab-resistant donor, seedlings carrying the ACS1-2 allele (341 bp, Md-ACS1, Chr 15) and the Vf resistance allele (200 bp, AL07, Chr 1) will exhibit significantly higher fruit firmness and field resistance to *Venturia inaequalis* than seedlings lacking one or both alleles."

---

## Step-by-Step Research Workflow

### Step 1: Explore the Marker Database

```python
from src.marker_selector import MarkerSelector

selector = MarkerSelector()

# View all 17 chromosomes
chromosomes = selector.list_chromosomes()
for c in chromosomes:
    print(f"Chr {c['chromosome']:2s} ({c['linkage_group']}): {c['description'][:60]}")

# List all available trait categories
for cat in selector.list_trait_categories():
    print(f"  {cat['category']:20s} — {cat['description'][:50]}")
```

### Step 2: Retrieve Markers by Chromosome

Query all validated markers on a specific linkage group to map structural traits:

```python
# All markers on Chromosome 10 (ACO1, PG1, Co)
chr10_markers = selector.markers_by_chromosome("10")
for m in chr10_markers:
    print(f"  {m['marker_name']:12s} | {m['target_trait']:30s} | {m['position_cM']} cM")
```

### Step 3: Select Markers by Trait

Choose target traits and retrieve associated markers:

```python
# Markers for fruit quality
quality = selector.markers_by_trait("fruit_quality")
for m in quality:
    print(f"  {m['marker_name']:12s} | {m['target_trait']}")

# Markers for scab resistance
scab = selector.markers_by_trait("scab_resistance")
for m in scab:
    print(f"  {m['marker_name']:12s} | Chr {m['chromosome']} | {m['target_trait']}")
```

### Step 4: Generate a Genotyping Primer Manifest

Produce a laboratory-ready Markdown document listing all primers, expected allele sizes, and PCR conditions:

```python
manifest = selector.genotyping_manifest(
    traits=["fruit_quality", "disease_resistance"],
    population_type="BC1",
    population_size=200,
    title="Firmness + Scab Resistance Genotyping Plan",
)
print(manifest)

# Save to file
selector.export_manifest(
    "output/priming_plan.md",
    traits=["fruit_quality", "disease_resistance"],
    population_type="BC1",
)
```

### Step 5: Import and Parse Genotype Data

Prepare your experimental data in a standardised CSV format:

```python
import pandas as pd

df = pd.read_csv("your_genotype_data.csv")
# Expected columns: sample_id, marker, genotype, trait, no_amplification
```

**Input format (long CSV):**

```csv
sample_id,marker,genotype,trait,no_amplification
S001,AL07,161:161,Apple Scab Resistance,False
S001,Md_ACS1,341:341,Fruit Firmness,False
S001,Ma1_INDEL,409:456,Fruit Acidity,False
S002,AL07,142:142,Apple Scab Resistance,False
S002,Md_ACS1,200:200,Fruit Firmness,False
S002,Ma1_INDEL,409:409,Fruit Acidity,False
```

### Step 6: Run Segregation Distortion Tests

Verify whether observed genotype ratios conform to expected Mendelian proportions:

```python
from src.breeding_analytics import (
    segregation_distortion_test,
    segregation_distortion_from_df,
    batch_segregation_test,
)

# Single marker test
result = segregation_distortion_test(
    observed_genotypes={"AA": 48, "AB": 52, "BB": 0},
    population_type="BC1",
)
print(result["interpretation"])

# From DataFrame
seg = segregation_distortion_from_df(df, marker="AL07", population_type="BC1")
print(seg["interpretation"])

# Batch test across all markers
seg_batch = batch_segregation_test(df, population_type="BC1")
print(seg_batch[["marker", "chi2", "p_value", "significant"]])
```

### Step 7: Genotype-Phenotype Mapping

Visualise trait variance across marker genotype classes:

```python
from src.breeding_analytics import genotype_phenotype_boxplot

genotype_phenotype_boxplot(
    genotypes=["AA", "AB", "BB", "AA", "BB", "AB"],
    phenotypic_values=[12.5, 14.2, 11.8, 13.1, 10.9, 14.0],
    marker_name="Md_ACS1",
    trait_name="Fruit Firmness (kg/cm²)",
    save_path="output/acs1_firmness_boxplot.png",
)
```

### Step 8: Calculate Polymorphism Metrics

Evaluate marker informativeness:

```python
from src.breeding_analytics import calculate_pic, calculate_pic_from_genotypes, batch_pic

# Single marker PIC
pic = calculate_pic({"196": 0.35, "168": 0.25, "210": 0.15})
print(f"PIC: {pic:.4f}")

# From genotype list
pic_direct = calculate_pic_from_genotypes(["390:390", "390:490", "490:490"])
print(f"PIC: {pic_direct:.4f}")

# Batch PIC for all markers
pic_batch = batch_pic(df)
print(pic_batch)
```

### Step 9: Generate a Comprehensive Report

```python
from src.breeding_analytics import generate_report

report_path = generate_report(df, population_type="BC1", output_dir="output/report")
print(f"Report saved to: {report_path}")
```

---

## Marker Database Reference

The curated database (`data/apple_marker_map.json`) contains validated markers across 17 chromosomes:

| Marker | Type | Chr | Trait | Key Alleles |
|--------|------|-----|-------|-------------|
| AL07 | SSR | 1 | Apple Scab Resistance | 161 (R), 142 (S) |
| Vf_SCAR | SCAR | 1 | Apple Scab Resistance | 200 (R), NB (S) |
| AM19 | SCAR | 1 | Apple Scab Resistance | 200 (R), NB (S) |
| Rvi2_SSR | SSR | 2 | Apple Scab Resistance | 135 (R) |
| Fb_Mr5 | SSR | 3 | Fire Blight Resistance | 180 (R) |
| Rvi4_SSR | SSR | 5 | Apple Scab Resistance | 125 (R) |
| CRISP_SSR | SSR | 5 | Fruit Crispness | 196 (crisp), 168 (mealy) |
| Ma3_SSR | SSR | 8 | Fruit Acidity (modifier) | 170 |
| BP16 | SSR | 8 | Bitter Pit Susceptibility | 203 (R), 368 (S) |
| BP13 | SSR | 8 | Bitter Pit Susceptibility | 224 (R), 232 (S) |
| MYB10 | SSR | 9 | Flesh Color | 390 (white), 490 (red) |
| RED_TE | SCAR | 9 | Skin Color | 750 (red), NB (non-red) |
| Md_ACO1 | SSR | 10 | Fruit Firmness / Ethylene | 237 (firm), 300 (standard) |
| Md_PG1 | SSR | 10 | Fruit Texture / Softening | 289/292 (firm), 298 (soft) |
| Co_INDEL | SCAR | 10 | Columnar Growth | 180 (columnar), 220 (standard) |
| Pl2_SSR | SSR | 11 | Powdery Mildew Resistance | 178 (R), 150 (S) |
| Rvi15_SSR | SSR | 11 | Apple Scab Resistance | 145 (R) |
| Pl5_SSR | SSR | 12 | Powdery Mildew Resistance | 132 (R) |
| Md_ACS1 | SSR | 15 | Fruit Firmness / Ethylene | 341/247 (firm), 200 (standard) |
| Ma1_SNP | SNP | 16 | Fruit Acidity | Ma1 (tart), ma1 (sweet) |
| Ma1_INDEL | InDel | 16 | Fruit Acidity | 409 (tart), 456 (sweet) |

**Trait categories:** `fruit_quality`, `disease_resistance`, `tree_architecture`, `acidity`, `texture`, `color`, `scab_resistance`, `mildew_resistance`

---

## Analytics Engine Reference

### Segregation Distortion Test

Tests whether observed genotype ratios deviate from expected Mendelian proportions using chi-square goodness-of-fit analysis.

```python
from src.breeding_analytics import segregation_distortion_test

result = segregation_distortion_test(
    observed_genotypes={"AA": 48, "AB": 52, "BB": 0},
    population_type="BC1",
    significance=0.05,
)
print(result["interpretation"])
```

### Genotype-Phenotype Boxplot

Generates box-and-whisker plots comparing trait variance across marker genotype classes (e.g., Homozygous Resistant vs. Heterozygous vs. Homozygous Susceptible). Includes Kruskal-Wallis H-test for non-parametric significance testing.

```python
from src.breeding_analytics import genotype_phenotype_boxplot

fig = genotype_phenotype_boxplot(
    genotypes=["R/R", "R/S", "S/S", "R/R", "S/S", "R/S"],
    phenotypic_values=[12.5, 14.2, 11.8, 13.1, 10.9, 14.0],
    marker_name="AL07",
    trait_name="Apple Scab Resistance Score",
    save_path="output/scab_boxplot.png",
)
```

### Polymorphism Information Content (PIC)

PIC quantifies the informativeness of a marker based on allele frequencies:

    PIC = 1 − Σ(pᵢ²) − Σᵢ<ⱼ(2 · pᵢ² · pⱼ²)

```python
from src.breeding_analytics import calculate_pic

pic = calculate_pic({"196": 0.35, "168": 0.25, "210": 0.15})
print(f"PIC: {pic:.4f}")
```

**Interpretation:** PIC > 0.5 = highly informative; PIC 0.25–0.5 = reasonably informative; PIC < 0.25 = slightly informative.

---

## Apple Breeding Resources

### Genomic Databases

| Resource | URL | Description |
|----------|-----|-------------|
| Genome Database for Rosaceae (GDR) | [rosaceae.org](https://www.rosaceae.org) | Genomic, genetic, and breeding data for Rosaceae. Marker nomenclature used in this toolkit follows GDR conventions. |
| GDDH13 Reference Genome | GDR | Golden Delicious doubled-haploid assembly used as the coordinate system for all markers. |
| Ensembl Plants (*Malus domestica*) | [plants.ensembl.org](https://plants.ensembl.org/Malus_domestica/) | Gene annotation, synteny, and variation data. |
| USDA National Plant Germplasm System | [npgsweb.ars-grin.gov](https://npgsweb.ars-grin.gov/) | Public repository of apple cultivars and wild species. |

### QTL Mapping and Genomic Selection Tools

| Tool | Language | Use | Repository |
|------|----------|-----|------------|
| qtl / qtl2 | R | QTL mapping in segregating populations | [github.com/rqtl/qtl](https://github.com/rqtl/qtl) |
| sommer | R | Genomic Selection via mixed models (GEBV) | [github.com/Covarrubias-Pazaran/sommer](https://github.com/Covarrubias-Pazaran/sommer) |
| rrBLUP | R | Genomic prediction using ridge regression | [CRAN](https://CRAN.R-project.org/package=rrBLUP) |
| GAPIT | R | GWAS and genomic prediction | [github.com/jiabocwang/GAPIT](https://github.com/jiabocwang/GAPIT) |

### Apple-Specific Research Repositories

| Repository | Focus | URL |
|------------|-------|-----|
| MylesLab/apple-aroma | Genetic architecture of apple aroma | [github.com/MylesLab/apple-aroma](https://github.com/MylesLab/apple-aroma) |
| CooperstoneLab/apple | Pedigree-based genotyping and quality | [github.com/CooperstoneLab/apple](https://github.com/CooperstoneLab/apple) |

### Key Publications

- Daccord et al., 2014, *Nature Genetics* — Golden Delicious genome
- Velasco et al., 2010, *Nature Genetics* — Apple genome resequencing
- Gianfranceschi et al., 1999, *Theor Appl Genet* — Vf scab resistance marker
- Bai et al., 2012, *Mol Genet Genomics* — Ma locus acidity
- Longhi et al., 2013, *BMC Plant Biology* — Fruit texture markers
- Espley et al., 2007, *Plant Journal* — MYB10 red flesh color
- Botstein et al., 1980, *Am J Hum Genet* — PIC definition
- Sunako et al., 1999, *Plant Physiology* — Md-ACS1 ethylene gene
- Costa et al., 2005, *Euphytica* — Md-ACO1 ethylene gene

---

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on code style, testing, and pull request procedures.

---

## Citation

If you use this toolkit in your research, please cite:

```bibtex
@software{apu2024_apple_mas,
  author    = {Md Mosiur Rahman Bhuyin Apu},
  title     = {Apple MAS Toolkit: Marker-Assisted Selection for Apple Breeding},
  year      = {2024},
  url       = {https://github.com/Mosiuropu/apple-mas-toolkit-hidden}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.
