# Apple MAS Toolkit

**Marker-Assisted Selection (MAS) Toolkit for Apple Breeding Programs**

A comprehensive Python toolkit for analyzing molecular markers (SSR and SCAR) in apple (*Malus domestica*) germplasm to support parent selection for fruit quality traits.

---

## Overview

This toolkit provides researchers and breeders with tools to:

- **Query marker information** for 8 key markers associated with fruit quality traits
- **Parse genotype data** from GeneMapper, CSV, Excel, or structured formats
- **Analyze allele frequencies** and genotype distributions with statistical metrics
- **Score and rank accessions** based on customizable breeding profiles
- **Generate publication-quality visualizations** for analysis results
- **Find complementary parent pairs** for crossing programs

## Supported Markers

| Marker | Trait | Chromosome | Type | Key Alleles |
|--------|-------|------------|------|-------------|
| MYB10 | Flesh Color | 9 | SSR | 390 (white), 490 (red) |
| RED_TE | Skin Color | 9 | SCAR | 750 (red), NB (non-red) |
| MA_INDEL | Fruit Acidity | 16 | SCAR | 409 (high), 456 (low) |
| BP16 | Bitter Pit | 8 | SSR | 203 (resistant), 368 (susceptible) |
| BP13 | Bitter Pit | 8 | SSR | 224-260 (variable) |
| ACS | Fruit Texture | 10 | SCAR | 200, 341 (firmness) |
| ACO | Fruit Texture | 10 | SCAR | 237, 300 (ethylene) |
| MD_PG1 | Fruit Texture | 10 | SSR | 289, 292, 298 (firmness) |

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/mosiurrahmanapu/apple-mas-toolkit.git
cd apple-mas-toolkit

# Install in development mode
pip install -e .

# Or install with development dependencies
pip install -e ".[dev]"
```

### Python API

```python
from apple_mas import DataParser, AlleleAnalyzer, ParentScorer, MarkerVisualizer, MarkerDatabase

# Initialize components
db = MarkerDatabase()
parser = DataParser()
analyzer = AlleleAnalyzer()
scorer = ParentScorer()
viz = MarkerVisualizer()

# Load sample data
df, raw = parser.create_sample_dataset()

# Analyze a specific marker
freq = analyzer.allele_frequencies(df, marker="MYB10")
print(freq)

# Get diversity metrics
pic = analyzer.polymorphism_information_content(df, marker="MYB10")
print(f"PIC: {pic:.4f}")

# Rank accessions by breeding profile
ranked = scorer.rank_accessions(df, profile="premium_table")
print(ranked.head(10))

# Find complementary parent pairs
pairs = scorer.find_complementary_parents(df, n_pairs=5)
print(pairs)

# Generate visualizations
viz.allele_frequency_plot(df, marker="MYB10", save_path="myb10_freq.png")
viz.genotype_distribution(df, marker="MYB10", save_path="myb10_geno.png")
viz.multi_marker_frequency(df, save_path="all_markers.png")
viz.ranking_barplot(ranked, save_path="ranking.png")
```

### Command Line Interface

```bash
# Show marker database summary
apple-mas info

# Show details for a specific marker
apple-mas info --marker MYB10

# List all markers
apple-mas info --list

# Show PCR methodology
apple-mas info --pcr

# Search markers by keyword
apple-mas info --search "texture"

# Analyze genotype data from CSV
apple-mas analyze --input my_data.csv

# Analyze a specific marker
apple-mas analyze --input my_data.csv --marker ACS

# Rank accessions with a breeding profile
apple-mas rank --input my_data.csv --profile premium_table --top 20

# Find complementary parent pairs
apple-mas rank --input my_data.csv --pairs 10

# Generate visualization plots
apple-mas visualize --input my_data.csv --output plots/ --format png

# Generate a comprehensive report
apple-mas report --input my_data.csv --profile balanced --output report/

# Create sample dataset for testing
apple-mas sample-data --output sample_data/
```

## Breeding Profiles

The toolkit includes predefined breeding profiles with optimized trait weightings:

| Profile | Description | Optimized Traits |
|---------|-------------|------------------|
| `premium_table` | Premium table apple | Texture (2x), Bitter Pit (1.5x), Acidity (1.5x), Skin (1x) |
| `juice_cider` | Juice/cider production | Acidity (2.5x), Others (0.5x) |
| `red_flesh_novelty` | Novelty red-fleshed apple | Flesh Color (3x), Texture (1.5x), Others (1x) |
| `storage_shelf_life` | Long storage/shelf life | Texture (3x), Bitter Pit (2x), Others (1x) |
| `balanced` | Equal weighting | All traits (1x) |

## Input Data Formats

### Long Format (CSV/Excel)

```csv
sample_id,marker,genotype
ANABP_001,MYB10,390:390
ANABP_001,RED_TE,750
ANABP_001,ACS,200:341
ANABP_002,MYB10,390:490
ANABP_002,RED_TE,NB
ANABP_002,ACS,341:341
```

### Wide Format (CSV)

```csv
sample_id,MYB10,RED_TE,ACS,ACO,MD_PG1,BP16,BP13,MA_INDEL
ANABP_001,390:390,750,200:341,300:300,289:292,203:368,232:250,409:456
ANABP_002,390:490,NB,341:341,237:300,298:298,203:203,224:232,456:456
```

### Python Dictionary

```python
data = {
    "ANABP_001": {"MYB10": "390:390", "RED_TE": "750", "ACS": "200:341"},
    "ANABP_002": {"MYB10": "390:490", "RED_TE": "NB", "ACS": "341:341"},
}
df = parser.from_dict(data)
```

## Analysis Metrics

### Diversity Metrics
- **Polymorphism Information Content (PIC)**: Measures marker informativeness (0-1)
- **Observed Heterozygosity**: Proportion of heterozygous genotypes
- **Expected Heterozygosity**: 1 - Σ(pi²) based on allele frequencies
- **Heterozygosity Deficit**: Difference between expected and observed

### Statistical Tests
- **Hardy-Weinberg Equilibrium**: Chi-squared test for biallelic markers
- **Allele Frequency Analysis**: Count and frequency calculations
- **Genotype Distribution**: Proportion analysis across accessions

## Visualization Types

| Plot Type | Method | Description |
|-----------|--------|-------------|
| Allele Frequency | `allele_frequency_plot()` | Bar plot of allele frequencies |
| Genotype Distribution | `genotype_distribution()` | Pie chart of genotype proportions |
| Multi-Marker Comparison | `multi_marker_frequency()` | Grid of allele frequency plots |
| Genotype Heatmap | `genotype_heatmap()` | Heatmap of genotype frequencies |
| Sample Composition | `sample_composition_chart()` | Stacked bar chart |
| Parent Ranking | `ranking_barplot()` | Horizontal bar chart of ranked accessions |
| Trait Radar | `trait_score_radar()` | Spider/radar chart comparing accessions |

## Documentation

- [Getting Started](docs/getting_started.md)
- [Marker Reference](docs/marker_reference.md)
- [Contributing Guide](CONTRIBUTING.md)

## Data Sources

This toolkit incorporates marker information and allele definitions from:

- Australian National Apple Breeding Program (ANABP) germplasm analysis
- Published literature on apple molecular markers (see individual marker references)

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](docs/contributing.md) for guidelines.

## Citation

If you use this toolkit in your research, please cite:

```bibtex
@software{apu2024_apple_mas,
  author = {Md Mosiur Rahman Bhuyin Apu},
  title = {Apple MAS Toolkit: Marker-Assisted Selection for Apple Breeding},
  year = {2024},
  url = {https://github.com/mosiurrahmanapu/apple-mas-toolkit}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- Australian National Apple Breeding Program (ANABP) for germplasm access
- Assoc. Prof. Michael Considine and Dr. Sultan Mia for supervision
- University of Western Australia, School of Agriculture and Environment
