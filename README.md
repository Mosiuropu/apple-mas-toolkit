# Apple MAS Toolkit

**Marker-Assisted Selection (MAS) Toolkit for Apple (*Malus × domestica*) Molecular Breeding**

A comprehensive Python toolkit for mapping agronomic traits to validated molecular markers, executing marker-assisted selection workflows, and performing statistical analyses on genotype data within apple breeding programs.

---

## Table of Contents

1. [Overview](#overview)
2. [Supported Markers](#supported-markers)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Workflow Guide](#workflow-guide)
6. [Marker Reference](#marker-reference)
7. [Analytics Reference](#analytics-reference)
8. [Command-Line Interface](#command-line-interface)
9. [Python API Reference](#python-api-reference)
10. [Input Data Formats](#input-data-formats)
11. [Apple Breeding Resources](#apple-breeding-resources)
12. [Contributing](#contributing)
13. [Citation](#citation)
14. [License](#license)

---

## Overview

The Apple MAS Toolkit provides researchers and students with a structured framework for:

- **Mapping agronomic characters** to validated molecular markers (SSR, SCAR, InDel, SNP).
- **Filtering and selecting markers** based on trait targets through an interactive selection interface.
- **Generating genotyping plans** with primer sequences, linkage groups, expected allele sizes, and references.
- **Performing statistical analyses** including segregation distortion tests, Polymorphism Information Content (PIC) calculation, and genotype-phenotype association evaluation.
- **Visualizing marker data** with publication-quality plots (allele frequencies, genotype distributions, heatmaps, boxplots).
- **Scoring and ranking** accessions for breeding suitability using customizable trait profiles.

The toolkit aligns its nomenclature and data structures with standards established by the **Genome Database for Rosaceae (GDR)** and references the Golden Delicious v1.0 reference genome (Daccord et al., 2014).

---

## Supported Markers

The database includes validated markers across six trait categories:

| Marker | Trait | Chr | LG | Type | Key Alleles |
|--------|-------|-----|----|------|-------------|
| **Fruit Quality** | | | | | |
| MA_INDEL | Fruit Acidity (Ma locus) | 16 | LG16 | SCAR | 409 (high acid), 456 (low acid) |
| ACS | Ethylene / Firmness (Md-ACS1) | 10 | LG10 | SCAR | 200 (standard), 341 (firm) |
| ACO | Ethylene / Firmness (Md-ACO1) | 10 | LG10 | SCAR | 237 (firm), 300 (standard) |
| MD_PG1 | Firmness (Polygalacturonase) | 10 | LG10 | SSR | 289, 292 (firm), 298 (soft) |
| CRISP_SSR | Fruit Crispness | 5 | LG05 | SSR | 196 (crisp), 168 (mealy) |
| **Color** | | | | | |
| MYB10 | Flesh Color | 9 | LG09 | SSR | 390 (white), 490 (red) |
| RED_TE | Skin Color | 9 | LG09 | SCAR | 750 (red), NB (non-red) |
| **Disease Resistance** | | | | | |
| Vf_SCAR | Apple Scab (Vf/Rvi6) | 1 | LG01 | SCAR | 200 (resistant), NB (susceptible) |
| AL07 | Apple Scab (Vf-linked SSR) | 1 | LG01 | SSR | 161 (Vf-linked), 142 (susceptible) |
| PL2_SSR | Powdery Mildew (Pl2) | 11 | LG11 | SSR | 178 (resistant), 150 (susceptible) |
| BP16 | Bitter Pit Susceptibility | 8 | LG08 | SSR | 203 (resistant), 368 (susceptible) |
| BP13 | Bitter Pit Susceptibility | 8 | LG08 | SSR | 224 (low risk), 232 (high risk) |
| **Growth Habit** | | | | | |
| CO_INDEL | Columnar Growth (Co locus) | 10 | LG10 | SCAR | 180 (columnar), 220 (standard) |

---

## Installation

```bash
# Clone the repository
git clone https://github.com/Mosiuropu/apple-mas-toolkit.git
cd apple-mas-toolkit

# Install in development mode
pip install -e .

# Or install with development dependencies
pip install -e ".[dev]"
```

### Requirements

- Python >= 3.9
- pandas >= 1.5
- numpy >= 1.23
- matplotlib >= 3.6
- seaborn >= 0.12
- scipy >= 1.9

---

## Quick Start

```python
from apple_mas import MarkerSelector, MarkerDatabase
from apple_mas.analytics import segregation_distortion_from_df, batch_pic

# 1. Explore available markers and traits
selector = MarkerSelector()
traits = selector.list_available_traits()
print(traits)

# 2. Build a genotyping plan for a specific breeding objective
plan = selector.build_genotyping_plan(
    target_traits=["fruit_quality", "disease_resistance"],
    population_type="BC1",
    population_size=200,
)
print(selector.plan_summary(plan))

# 3. Export the plan
selector.export_plan(plan, "my_genotyping_plan.json")
selector.export_plan(plan, "marker_table.csv", fmt="csv")
```

---

## Workflow Guide

### Step 1: Define Breeding Objectives

Before selecting markers, clearly define the target traits for your breeding program. The toolkit supports six trait categories:

| Category | Markers | Description |
|----------|---------|-------------|
| `fruit_quality` | MA_INDEL, ACS, ACO, MD_PG1, CRISP_SSR | Acidity, firmness, crispness |
| `disease_resistance` | BP16, BP13, Vf_SCAR, AL07, PL2_SSR | Apple scab, powdery mildew, bitter pit |
| `color` | MYB10, RED_TE | Skin and flesh pigmentation |
| `texture` | ACS, ACO, MD_PG1, CRISP_SSR | Firmness, crispness, ethylene |
| `acidity` | MA_INDEL | Malic acid content |
| `growth_habit` | CO_INDEL | Columnar tree architecture |

```python
from apple_mas import MarkerSelector

selector = MarkerSelector()
traits = selector.list_available_traits()
print(traits)
```

### Step 2: Select Markers

Use the `MarkerSelector` to query markers by trait, keyword, or chromosome:

```python
# List all markers
all_markers = selector.list_available_markers()

# Search by keyword
scab_markers = selector.search_markers("scab")

# Get markers for a specific trait
disease_markers = selector.get_markers_for_trait("disease_resistance")

# Get primer information for a specific marker
primers = selector.get_marker_primer_info("Vf_SCAR")
print(primers)
```

### Step 3: Generate a Genotyping Plan

Construct a complete genotyping plan specifying your target traits and population type:

```python
plan = selector.build_genotyping_plan(
    target_traits=["fruit_quality", "disease_resistance"],
    population_type="BC1",
    population_size=200,
    notes="Backcross population for fruit quality improvement",
)

# View summary
print(selector.plan_summary(plan))

# Export to JSON or CSV
selector.export_plan(plan, "plan.json")
selector.export_plan(plan, "markers.csv", fmt="csv")
```

### Step 4: Import and Parse Genotype Data

The toolkit supports multiple input formats:

```python
from apple_mas import DataParser

parser = DataParser()

# From CSV (long format)
df = parser.from_csv("genotypes.csv", sample_col="sample_id",
                     marker_col="marker", genotype_col="genotype")

# From wide-format CSV
df = parser.from_wide_csv("genotypes_wide.csv", sample_col="sample_id")

# From Excel
df = parser.from_excel("genotypes.xlsx")

# From Python dictionary
data = {"S1": {"MYB10": "390:390", "Vf_SCAR": "200"}, ...}
df = parser.from_dict(data)

# Generate sample data for testing
df, raw = parser.create_sample_dataset()
```

### Step 5: Analyze Marker Data

Perform statistical analyses on your genotype data:

```python
from apple_mas import AlleleAnalyzer
from apple_mas.analytics import segregation_distortion_from_df, batch_pic

analyzer = AlleleAnalyzer()

# Allele frequencies
freq = analyzer.allele_frequencies(df, marker="Vf_SCAR")

# PIC calculation
pic = batch_pic(df)
print(pic)

# Segregation distortion test (F1 population)
seg_result = segregation_distortion_from_df(df, marker="ACS", population_type="F1")
print(seg_result["interpretation"])

# Hardy-Weinberg equilibrium test
hwe = analyzer.hardy_weinberg_test(df, marker="ACS")

# Cross-marker diversity summary
summary = analyzer.cross_marker_summary(df)
```

### Step 6: Visualize Results

Generate publication-quality figures:

```python
from apple_mas import MarkerVisualizer

viz = MarkerVisualizer()

# Allele frequency bar plot
viz.allele_frequency_plot(df, marker="Vf_SCAR", save_path="vf_freq.png")

# Genotype distribution pie chart
viz.genotype_distribution(df, marker="ACS", save_path="acs_geno.png")

# Multi-marker comparison
viz.multi_marker_frequency(df, save_path="all_freq.png")

# Genotype heatmap
viz.genotype_heatmap(df, save_path="heatmap.png")

# Genotype-phenotype boxplot (requires external phenotype data)
from apple_mas.analytics import genotype_phenotype_boxplot
genotype_phenotype_boxplot(
    genotypes=["AA", "AB", "BB", "AA", "BB", ...],
    phenotypic_values=[12.5, 14.2, 11.8, 13.1, 10.9, ...],
    marker_name="ACS",
    trait_name="Fruit Firmness",
    save_path="acs_firmness_boxplot.png",
)
```

### Step 7: Score and Rank Parents

Rank accessions using predefined or custom breeding profiles:

```python
from apple_mas import ParentScorer

scorer = ParentScorer()

# List available profiles
profiles = scorer.list_profiles()
print(profiles)

# Rank by premium table profile
ranked = scorer.rank_accessions(df, profile="premium_table", top_n=20)
print(ranked)

# Find complementary parent pairs
pairs = scorer.find_complementary_parents(df, n_pairs=5, profile="balanced")
print(pairs)

# Visualize rankings
viz.ranking_barplot(ranked, save_path="ranking.png")
```

---

## Marker Reference

### Apple Scab Resistance (Vf/Rvi6)

The Vf locus on chromosome 1 (LG01) was introgressed from *Malus micromalus* and confers broad-spectrum resistance to *Venturia inaequalis*. Two markers are available:

- **Vf_SCAR**: Dominant SCAR marker detecting presence/absence of the resistance allele.
- **AL07**: Co-dominant SSR marker linked to the Vf locus, enabling zygosity determination.

For durable resistance, pyramiding with additional Rvi loci (Rvi2, Rvi4, Rvi5, Rvi15) is recommended.

### Powdery Mildew Resistance (Pl2)

The Pl2 locus on chromosome 11 (LG11) was identified from *Malus zumi* and confers moderate to high resistance to *Podosphaera leucotricha*. The PL2_SSR marker enables selection at the seedling stage.

### Columnar Growth (Co)

The Co locus on chromosome 10 (LG10) controls columnar growth habit, originating from the 'Wijcik McIntosh' mutation. The CO_INDEL marker distinguishes homozygous columnar, heterozygous semi-columnar, and standard growth habit genotypes.

### Fruit Acidity (Ma locus)

The Ma locus on chromosome 16 (LG16) encodes aluminum-activated malate transporters (ALMT) controlling malic acid accumulation. The MA_INDEL marker identifies high-acidity (409 bp) and low-acidity (456 bp) alleles.

### Fruit Texture and Firmness

Three markers on chromosome 10 (LG10) provide comprehensive texture profiling:

- **ACS** (Md-ACS1): Ethylene biosynthesis; 341 bp allele associated with reduced ethylene and extended shelf life.
- **ACO** (Md-ACO1): Ethylene oxidation; 237 bp allele associated with enhanced firmness.
- **MD_PG1** (Polygalacturonase): Cell wall degradation; 289/292 bp alleles associated with firm texture.

### Fruit Crispness

The CRISP_SSR marker on chromosome 5 (LG05) is associated with crispness retention during cold storage. The 196 bp allele is linked to delayed mealiness onset.

---

## Analytics Reference

### Segregation Distortion Test

Tests whether observed genotype ratios deviate from expected Mendelian proportions using chi-square goodness-of-fit analysis.

```python
from apple_mas.analytics import segregation_distortion_test

result = segregation_distortion_test(
    observed_genotypes={"AA": 48, "AB": 52, "BB": 0},
    population_type="F1",
    significance=0.05,
)
print(result["interpretation"])
```

### Genotype-Phenotype Boxplot

Evaluates phenotypic variance grouped by marker genotype class (e.g., AA vs. AB vs. BB). Includes Kruskal-Wallis H-test for non-parametric significance testing.

```python
from apple_mas.analytics import genotype_phenotype_boxplot

fig = genotype_phenotype_boxplot(
    genotypes=["AA", "AB", "BB", ...],
    phenotypic_values=[12.5, 14.2, 11.8, ...],
    marker_name="ACS",
    trait_name="Fruit Firmness (kg/cm2)",
    save_path="boxplot.png",
)
```

### Polymorphism Information Content (PIC)

PIC quantifies the informativeness of a marker based on allele frequencies. For multi-allelic markers:

    PIC = 1 - sum(p_i^2) - sum_{i<j} (2 * p_i^2 * p_j^2)

```python
from apple_mas.analytics import calculate_pic

pic = calculate_pic({"196": 0.35, "168": 0.25, "210": 0.15})
print(f"PIC: {pic:.4f}")
```

---

## Command-Line Interface

```bash
# Show marker database summary
apple-mas info

# Show details for a specific marker
apple-mas info --marker Vf_SCAR

# List all markers
apple-mas info --list

# Search markers by keyword
apple-mas info --search "scab"

# Select markers and build a genotyping plan
apple-mas select --traits fruit_quality disease_resistance --plan-out plan.json

# Select markers for columnar habit
apple-mas select --traits growth_habit --population BC1 --size 200

# Analyze genotype data
apple-mas analyze --input data.csv

# Run segregation distortion tests
apple-mas analytics --input data.csv --test segregation

# Calculate PIC for all markers
apple-mas analytics --input data.csv --test pic

# Rank accessions
apple-mas rank --input data.csv --profile premium_table --top 20

# Generate visualizations
apple-mas visualize --input data.csv --output plots/ --format png

# Generate comprehensive report
apple-mas report --input data.csv --output report/

# Create sample dataset
apple-mas sample-data --output sample_data/
```

---

## Python API Reference

### Core Classes

| Class | Module | Description |
|-------|--------|-------------|
| `MarkerDatabase` | `apple_mas.marker_db` | Query the marker database |
| `MarkerSelector` | `apple_mas.selector` | Interactive trait-to-marker selection workflow |
| `DataParser` | `apple_mas.data_parser` | Import and standardize genotype data |
| `AlleleAnalyzer` | `apple_mas.allele_analysis` | Allele frequency and diversity analysis |
| `ParentScorer` | `apple_mas.parent_scorer` | Score and rank accessions for breeding |
| `MarkerVisualizer` | `apple_mas.visualization` | Generate publication-quality plots |

### Analytics Functions

| Function | Module | Description |
|----------|--------|-------------|
| `segregation_distortion_test()` | `apple_mas.analytics` | Chi-square test for segregation distortion |
| `genotype_phenotype_boxplot()` | `apple_mas.analytics` | Boxplot of phenotype by marker genotype |
| `calculate_pic()` | `apple_mas.analytics` | PIC calculation from allele frequencies |
| `calculate_pic_from_df()` | `apple_mas.analytics` | PIC calculation from genotype DataFrame |
| `batch_segregation_test()` | `apple_mas.analytics` | Segregation tests for all markers |
| `batch_pic()` | `apple_mas.analytics` | PIC values for all markers |

---

## Input Data Formats

### Long Format (CSV)

```csv
sample_id,marker,genotype
SAMPLE_001,MYB10,390:390
SAMPLE_001,Vf_SCAR,200
SAMPLE_001,ACS,200:341
SAMPLE_002,MYB10,390:490
SAMPLE_002,Vf_SCAR,NB
SAMPLE_002,ACS,341:341
```

### Wide Format (CSV)

```csv
sample_id,MYB10,Vf_SCAR,ACS,CO_INDEL
SAMPLE_001,390:390,200,200:341,180:220
SAMPLE_002,390:490,NB,341:341,220:220
```

### Python Dictionary

```python
data = {
    "SAMPLE_001": {"MYB10": "390:390", "Vf_SCAR": "200", "ACS": "200:341"},
    "SAMPLE_002": {"MYB10": "390:490", "Vf_SCAR": "NB", "ACS": "341:341"},
}
df = parser.from_dict(data)
```

---

## Apple Breeding Resources

This toolkit is designed to complement existing tools and databases used in apple genomics and breeding research. Below are key resources that every apple breeder and student should know.

### Genomic Databases

| Resource | URL | Description |
|----------|-----|-------------|
| **Genome Database for Rosaceae (GDR)** | [rosaceae.org](https://www.rosaceae.org) | Comprehensive genomic, genetic, and breeding data for Rosaceae species including apple. Marker nomenclature and chromosome conventions used in this toolkit follow GDR standards. |
| **Apple Genome (Daccord et al., 2014)** | *Nature Genetics* | Golden Delicious v1.0 reference genome assembly used as the coordinate system for all markers in this database. |
| **Apple REFPOP** | [maxapress.com](https://www.maxapress.com/data/article/frures/preview/pdf/FruRes-2023-0027.pdf) | Standard reference population for genomics-assisted apple breeding. |
| **USDA National Plant Germplasm System** | [npgsweb.ars-grin.gov](https://npgsweb.ars-grin.gov/) | Public repository of apple cultivars and wild species for breeding programs. |

### QTL Mapping and Genomic Selection Tools

| Tool | Language | Primary Use | Repository |
|------|----------|-------------|------------|
| **qtl / qtl2** | R | QTL mapping in segregating populations | [github.com/rqtl/qtl](https://github.com/rqtl/qtl) |
| **sommer** | R | Genomic Selection via mixed models (GEBV calculation) | [github.com/Covarrubias-Pazaran/sommer](https://github.com/Covarrubias-Pazaran/sommer) |
| **rrBLUP** | R | Genomic prediction using ridge regression | [CRAN](https://CRAN.R-project.org/package=rrBLUP) |
| **GAPIT** | R | GWAS and genomic prediction | [github.com/jiabocwang/GAPIT](https://github.com/jiabocwang/GAPIT) |

### Apple-Specific Research Repositories

| Repository | Focus | URL |
|------------|-------|-----|
| **MylesLab/apple-aroma** | Genetic architecture of apple aroma compounds | [github.com/MylesLab/apple-aroma](https://github.com/MylesLab/apple-aroma) |
| **CooperstoneLab/apple** | Pedigree-based genotyping and apple quality traits | [github.com/CooperstoneLab/apple](https://github.com/CooperstoneLab/apple) |

### Key Publications

- **Daccord et al., 2014** — *Nature Genetics* — Golden Delicious genome sequence
- **Velasco et al., 2010** — *Nature Genetics* — Apple genome resequencing
- **Gianfranceschi et al., 1999** — *Theor Appl Genet* — Vf scab resistance marker
- **Bai et al., 2012** — *Mol Genet Genomics* — Ma locus acidity
- **Longhi et al., 2013** — *BMC Plant Biology* — Fruit texture markers
- **Espley et al., 2007** — *Plant Journal* — MYB10 red flesh color
- **Chagné et al., 2013** — *Plant Physiology* — Skin color markers

### Recommended Workflow for Students

1. **Define your breeding objective** — Choose target traits (e.g., firmness + disease resistance)
2. **Select markers** — Use `apple-mas select` to generate a genotyping plan
3. **Design crosses** — Use the marker database to plan crosses that combine desirable alleles
4. **Genotype seedlings** — Run PCR with selected markers and score alleles
5. **Analyze your data** — Use `apple-mas analyze` and `apple-mas analytics` for allele frequencies, PIC, and segregation tests
6. **Visualize results** — Generate publication-quality plots with `apple-mas visualize`
7. **Rank candidates** — Use `apple-mas rank` to identify the best individuals for your breeding objective
8. **Integrate with QTL mapping** — Export your data for use with qtl/qtl2, sommer, or GAPIT for genome-wide analyses

---

## External Data Standards

This toolkit aligns with standards from the following resources:

- **GDR (Genome Database for Rosaceae)**: Marker naming conventions and chromosome linkage group numbering follow GDR standards based on the Golden Delicious v1.0 reference genome.
- **Community marker aggregation conventions**: Data import schemas are modeled after established community tools that aggregate disparate marker names and accession thesauri.
- **MAS validation frameworks**: Genotypic/phenotypic data linkage blocks follow structures used in validated apple marker-assisted seedling selection models.

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
  url       = {https://github.com/Mosiuropu/apple-mas-toolkit}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Acknowledgments

- Genome Database for Rosaceae (GDR) for nomenclature standards and reference genome resources.
- The broader apple genomics community for establishing marker validation frameworks and sharing open-source tools.
- University of Western Australia, School of Agriculture and Environment.
