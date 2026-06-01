# Getting Started

This guide walks you through setting up and using the Apple MAS Toolkit.

## Prerequisites

- Python 3.9 or higher
- pip or conda package manager
- Git (for cloning)

## Installation

### Option 1: Install from source (recommended)

```bash
git clone https://github.com/mosiurrahmanapu/apple-mas-toolkit.git
cd apple-mas-toolkit
pip install -e .
```

### Option 2: Install with development dependencies

```bash
pip install -e ".[dev,notebooks]"
```

## Quick Start

### 1. Create a Sample Dataset

```bash
apple-mas sample-data --output sample_data/
```

This creates a realistic synthetic sample dataset with representative apple marker genotypes.

### 2. Explore the Marker Database

```bash
# Show database summary
apple-mas info

# List all markers
apple-mas info --list

# Get details for a specific marker
apple-mas info --marker MYB10

# Search markers by keyword
apple-mas info --search "texture"

# View PCR methodology
apple-mas info --pcr
```

### 3. Analyze Your Data

```bash
# Analyze all markers
apple-mas analyze --input sample_data/sample_genotypes.csv

# Analyze a specific marker
apple-mas analyze --input sample_data/sample_genotypes.csv --marker ACS
```

### 4. Rank Accessions for Breeding

```bash
# Rank with default balanced profile
apple-mas rank --input sample_data/sample_genotypes.csv

# Rank with specific profile
apple-mas rank --input sample_data/sample_genotypes.csv --profile premium_table

# Find complementary parent pairs
apple-mas rank --input sample_data/sample_genotypes.csv --pairs 10

# Save results
apple-mas rank --input sample_data/sample_genotypes.csv --output rankings.csv
```

### 5. Generate Visualizations

```bash
# Generate all plots
apple-mas visualize --input sample_data/sample_genotypes.csv --output plots/

# Generate specific plot types
apple-mas visualize --input sample_data/sample_genotypes.csv --type frequency
apple-mas visualize --input sample_data/sample_genotypes.csv --type genotype

# Generate plots for a specific marker
apple-mas visualize --input sample_data/sample_genotypes.csv --marker MYB10

# Change output format
apple-mas visualize --input sample_data/sample_genotypes.csv --format pdf
```

### 6. Generate a Complete Report

```bash
apple-mas report --input sample_data/sample_genotypes.csv --profile balanced --output report/
```

## Using the Python API

### Basic Workflow

```python
from apple_mas import DataParser, AlleleAnalyzer, ParentScorer, MarkerVisualizer

# Initialize
parser = DataParser()
analyzer = AlleleAnalyzer()
scorer = ParentScorer()
viz = MarkerVisualizer()

# Load data (from CSV, Excel, or create sample)
df, _ = parser.create_sample_dataset()

# Analyze alleles
freq = analyzer.allele_frequencies(df, marker="MYB10")
print(freq)

# Get diversity metrics
cross = analyzer.cross_marker_summary(df)
print(cross)

# Rank accessions
ranked = scorer.rank_accessions(df, profile="premium_table")
print(ranked.head(10))

# Find parent pairs
pairs = scorer.find_complementary_parents(df, n_pairs=5)
print(pairs)

# Generate plots
viz.allele_frequency_plot(df, marker="MYB10", save_path="freq.png")
viz.ranking_barplot(ranked, save_path="ranking.png")
```

### Loading Your Own Data

```python
# From CSV (long format)
df = parser.from_csv("my_data.csv", sample_col="Accession", marker_col="Locus", allele_col="Genotype")

# From Excel
df = parser.from_excel("my_data.xlsx", sheet_name="Sheet1")

# From wide-format CSV (one column per marker)
df = parser.from_wide_csv("my_data_wide.csv", sample_col="Accession")

# From a dictionary
data = {
    "APPLE_001": {"MYB10": "390:390", "RED_TE": "750", "ACS": "200:341"},
    "APPLE_002": {"MYB10": "390:490", "RED_TE": "NB", "ACS": "341:341"},
}
df = parser.from_dict(data)
```

### Customizing Breeding Profiles

```python
# Use a predefined profile
ranked = scorer.rank_accessions(df, profile="juice_cider")

# List available profiles
print(scorer.list_profiles())

# Access profile descriptions
from apple_mas.parent_scorer import BREEDING_PROFILES
for name, info in BREEDING_PROFILES.items():
    print(f"{name}: {info['description']}")
```

### Querying the Marker Database

```python
from apple_mas import MarkerDatabase

db = MarkerDatabase()

# List all markers
print(db.list_markers())

# Get marker details
info = db.get_marker("MYB10")
print(info["trait"])
print(info["allele_definitions"])

# Find markers for a trait
color_markers = db.get_markers_by_trait("color")

# Get chromosome locations
chr9_markers = db.get_markers_by_chromosome("9")

# Search markers
results = db.search_markers("bitter pit")

# Get breeding recommendations
notes = db.get_breeding_notes("MYB10")
```

## Understanding the Output

### Allele Frequency Output

| allele | count | frequency |
|--------|-------|-----------|
| 390    | 95    | 0.950     |
| 490    | 5     | 0.050     |

### Genotype Distribution Output

| genotype | count | frequency | percentage |
|----------|-------|-----------|------------|
| 390:390  | 47    | 0.959     | 95.9%      |
| 390:490  | 2     | 0.041     | 4.1%       |

### Diversity Metrics

- **PIC (Polymorphism Information Content)**: 0 = monomorphic, 1 = maximally informative
- **Observed Heterozygosity**: Proportion of heterozygous individuals
- **Expected Heterozygosity**: Expected proportion under Hardy-Weinberg equilibrium

### Parent Rankings

Accessions are scored based on:
1. **Marker scores**: Each genotype receives a score based on its desirability
2. **Trait weighting**: Traits are weighted according to the selected breeding profile
3. **No-amplification penalty**: Accessions with missing data are penalized

## Troubleshooting

### Common Issues

**"Column not found" error**
- Check that your CSV column names match the expected format
- Use `sample_col`, `marker_col`, and `genotype_col` parameters to specify column names

**"Marker not found" error**
- Ensure marker names match exactly: MYB10, RED_TE, MA_INDEL, BP16, BP13, ACS, ACO, MD_PG1
- Check for typos or different naming conventions

**Plots not displaying**
- In Jupyter notebooks, use `%matplotlib inline` before plotting
- For non-interactive environments, plots are saved to files

**Import errors**
- Ensure you've installed the package: `pip install -e .`
- Check that you're using Python 3.9+
