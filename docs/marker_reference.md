# Marker Reference

Complete reference for all molecular markers included in the Apple MAS Toolkit.

---

## Flesh Color: MYB10

**Full Name:** MdMYB10 (R2R3-MYB Transcription Factor)
**Chromosome:** 9
**Marker Type:** SSR
**Dye:** HEX

### Primer Sequences
- Forward: `/5HEX/GGAGGGGAATGAAGAAGAGG`
- Reverse: `TCCACAGAAGCAAACACTGAC`

### Allele Definitions

| Allele (bp) | Phenotype | Frequency (ANABP) | Desirability |
|-------------|-----------|-------------------|--------------|
| 390 | White/green flesh | 97.92% | Common |
| 490 | Red/pink flesh | 2.08% | Desirable |

### Genotype Interpretation

| Genotype | Phenotype | Notes |
|----------|-----------|-------|
| 390:390 | White/green flesh | Most common; low anthocyanin |
| 390:490 | Red/pink flesh | Heterozygous; moderate anthocyanin |
| 490:490 | Deep red flesh | Homozygous; high anthocyanin |

### Biological Context

MYB10 is a member of the R2R3-MYB transcription factor family, essential for regulating anthocyanin biosynthesis in apple flesh. The 490 bp allele activates genes involved in red pigmentation. Red-fleshed apples are rare in commercial germplasm, making accessions carrying this allele valuable breeding resources.

### Breeding Notes

- The 490 bp allele is crucial for introducing red flesh color
- Environmental conditions (light exposure, temperature) influence MYB10 expression
- The co-dominant nature aids in selecting offspring with desired traits
- Red-fleshed apples appeal to consumers for appearance and higher antioxidant content

### Key References

1. Espley et al., 2007, Plant Journal
2. Chagné et al., 2013, Plant Physiology
3. Bars-Cortina et al., 2017, J Agric Food Chem

---

## Skin Color: RED TE

**Full Name:** Red Transposable Element
**Chromosome:** 9
**Marker Type:** SCAR
**Probe Type:** Non-fluorescent

### Primer Sequences
- Forward: `/5HEX/GGTCACCCAACCCACACTGGGCCTTG`
- Reverse: `CGGCCGCAATCGCAAGACGCAGA`

### Allele Definitions

| Allele | Phenotype | Frequency (ANABP) | Desirability |
|--------|-----------|-------------------|--------------|
| 750 | Red skin | 64.06% | Desirable |
| NB | Non-red skin | 33.33% | Varies |

### Biological Context

RED TE targets a transposable element insertion within the MdMYB10 promoter region. The insertion boosts MdMYB10 expression, increasing anthocyanin production for red skin pigmentation. This is a dominant marker producing binary (present/absent) results.

### Limitations

- Cannot distinguish homozygous (750/750) from heterozygous (750/NB) genotypes
- Does not predict intensity of red pigmentation
- Environmental factors (light, temperature) affect final skin color expression

### Breeding Notes

- Simple, cost-effective, reliable for large-scale screening
- For zygosity information, use dual-band markers (Red_TE_Pos_F1 / Red_TE_Null_F2)
- High temperatures can reduce anthocyanin accumulation even with 750 bp allele

### Key References

1. Chagné et al., 2013
2. Zhang et al., 2019, Nature Communications
3. Moriya et al., 2017, Euphytica

---

## Fruit Acidity: MA_INDEL

**Full Name:** Ma Locus Insertion-Deletion
**Chromosome:** 16
**Marker Type:** SCAR
**Dye:** HEX

### Primer Sequences
- Forward: `/5HEX/AAATGAACAGGACCCAGACG`
- Reverse: `ACGACTCCAATCCAACATCC`

### Allele Definitions

| Allele (bp) | Phenotype | Acidity Effect | Frequency (ANABP) | Desirability |
|-------------|-----------|----------------|-------------------|--------------|
| 409 | High acidity | High malic acid | 8.0% | Desirable (tart) |
| 456 | Low acidity | Low malic acid | 10.0% | Desirable (sweet) |
| 452 | Variable | Intermediate | 5.21% | Needs characterization |
| 387 | Rare | Variable | 6.25% | Needs characterization |
| 402 | Rare | Variable | 6.25% | Needs characterization |

### Genotype Interpretation

| Genotype | Phenotype | Flavor Profile |
|----------|-----------|----------------|
| 409:409 | High acidity | Tart |
| 409:456 | Intermediate | Balanced |
| 456:456 | Low acidity | Sweet |
| 409:452 | Intermediate | Variable |
| 387:402 | Variable (rare) | Unknown |

### Biological Context

The Ma locus is defined in a 150 kb region containing 44 predicted genes on chromosome 16. Encodes aluminum-activated malate transporters (ALMT) that control malic acid accumulation. A natural truncation mutation in one of these genes is associated with low fruit acidity.

### Breeding Notes

- 409 bp allele: Higher malic acid (tart flavor, good for cider)
- 456 bp allele: Lower acidity (sweet)
- Rare alleles (387, 402, 452) need phenotypic validation
- Environmental factors (temperature, light, soil) can modify expression

### Key References

1. Bai et al., 2012, Mol Genet Genomics
2. Xu et al., 2012, Mol Breeding
3. Bianco et al., 2016, Plant Journal
4. Zheng et al., 2023, Plant Physiology

---

## Bitter Pit: BP16

**Full Name:** BP16-indel
**Chromosome:** 8
**Marker Type:** SSR
**Dye:** FAM

### Primer Sequences
- Forward: `/5HEX/TGCGATGTTAGAACGAGAGC`
- Reverse: `TTGGAAACTGGATTTTATTTGC`

### Allele Definitions

| Allele (bp) | Phenotype | Frequency (ANABP) | Desirability |
|-------------|-----------|-------------------|--------------|
| 203 | Lower susceptibility (resistance) | 22.0% | Desirable |
| 206 | Intermediate susceptibility | 6.0% | Acceptable |
| 368 | Higher susceptibility | 45.0% | Undesirable |

### Genotype Interpretation

| Genotype | Phenotype | Risk Level |
|----------|-----------|------------|
| 203:203 | Lower susceptibility | Low |
| 203:206 | Lower susceptibility | Low |
| 203:368 | Moderate susceptibility | Moderate |
| 206:206 | Moderate susceptibility | Moderate |
| 206:368 | Moderate-high susceptibility | Moderate-high |
| 368:368 | High susceptibility | High |

### Breeding Notes

- 368:368 genotype requires enhanced calcium management
- 203:203 genotype provides the best genetic resistance
- Combine with BP13 analysis for comprehensive assessment
- Field validation recommended

### Key References

1. Buti et al., 2015, Mol Breeding
2. Thapa et al., 2021, Plant Genome

---

## Bitter Pit: BP13

**Full Name:** BP13-SSR
**Chromosome:** 8
**Marker Type:** SSR
**Dye:** FAM

### Primer Sequences
- Forward: `/5HEX/CCTGTTCGCAAACAAGAAGG`
- Reverse: `GCGTAGTCAATCAAAACATTCG`

### Allele Definitions

| Allele (bp) | Phenotype | Frequency (ANABP) | Desirability |
|-------------|-----------|-------------------|--------------|
| 224 | Reduced risk | 10.42% | Desirable |
| 232 | Higher susceptibility | 30.21% | Undesirable |
| 250 | Variable | 15.0% | Intermediate |
| 252 | Variable | 5.0% | Intermediate |

### Breeding Notes

- Use in combination with BP16 for comprehensive assessment
- Non-232 bp combinations suggest reduced risk
- Wide allele range (222-260 bp) indicates substantial genetic variation

### Key References

1. Buti et al., 2015, Mol Breeding
2. Thapa et al., 2021, Plant Genome

---

## Fruit Texture: ACS

**Full Name:** Md-ACS1 (1-aminocyclopropane-1-carboxylic acid synthase)
**Chromosome:** 10
**Marker Type:** SCAR

### Primer Sequences
- Forward: `/5HEX/CGAGGTTGACTCAAATCAAAAC`
- Reverse: `GCTGATGAATGAGTCGTTGC`

### Allele Definitions

| Allele (bp) | Phenotype | Effect | Frequency (ANABP) | Desirability |
|-------------|-----------|--------|-------------------|--------------|
| 200 | Md-ACS1-1 | Standard ethylene | 38.0% | Standard |
| 341 | Md-ACS1-2 | Reduced ethylene | 38.0% | Desirable |

### Genotype Interpretation

| Genotype | Phenotype | Shelf Life |
|----------|-----------|------------|
| 200:200 | Standard ethylene | Standard |
| 200:341 | Intermediate | Slightly extended |
| 341:341 | Reduced ethylene | Extended |

### Breeding Notes

- 341 bp allele is favorable for extended shelf life
- Combine with ACO and Md-PG1 for comprehensive texture profiling

### Key References

1. Sunako et al., 1999, Plant Physiology
2. Costa et al., 2005, Euphytica

---

## Fruit Texture: ACO

**Full Name:** Md-ACO1 (1-aminocyclopropane-1-carboxylic acid oxidase)
**Chromosome:** 10
**Marker Type:** SCAR

### Primer Sequences
- Forward: `/5HEX/TGTTACCAATACATTTCAATTCTCG`
- Reverse: `GCTGATGAATGAGTCGTTGC`

### Allele Definitions

| Allele (bp) | Phenotype | Effect | Frequency (ANABP) | Desirability |
|-------------|-----------|--------|-------------------|--------------|
| 237 | ACO1-1 | Reduced ethylene oxidation | 23.0% | Desirable |
| 300 | ACO1-2 | Standard ethylene oxidation | 67.0% | Standard |

### Genotype Interpretation

| Genotype | Phenotype | Firmness |
|----------|-----------|----------|
| 237:237 | Reduced ethylene | Enhanced |
| 237:300 | Intermediate | Intermediate |
| 300:300 | Standard ethylene | Softer |

### Breeding Notes

- 237 bp allele is favorable for firm apples
- Best used in combination with ACS

### Key References

1. Costa et al., 2005, Euphytica
2. Zhu & Barritt, 2008, Tree Genetics & Genomes

---

## Fruit Texture: MD_PG1

**Full Name:** Md-PG1-10kDa SSR (Polygalacturonase)
**Chromosome:** 10
**Marker Type:** SSR
**Dye:** FAM

### Primer Sequences
- Forward: `/56-FAM/TTTCTTCCTTGGGTTTTTGG`
- Reverse: `CGAAGCCACTCTCCTTCTCC`

### Allele Definitions

| Allele (bp) | Phenotype | Frequency (ANABP) | Desirability |
|-------------|-----------|-------------------|--------------|
| 289 | Favorable for firmness | 20.0% | Desirable |
| 292 | Favorable for firmness | 20.0% | Desirable |
| 298 | Softer texture | 25.0% | Less desirable |

### Genotype Interpretation

| Genotype | Phenotype | Cell Wall Status |
|----------|-----------|-----------------|
| 289:292 | Firm fruit | Intact |
| 289:298 | Moderate firmness | Variable |
| 292:298 | Moderate firmness | Variable |
| 292:292 | Firm fruit | Intact |
| 298:298 | Softer texture | Degraded |

### Breeding Notes

- Alleles 289 and 292 are favorable for firmness
- 298:298 genotype should be avoided for firm-textured cultivars
- Best combined with ACS and ACO

### Key References

1. Longhi et al., 2013, BMC Plant Biology
2. Costa et al., 2010, J Exp Botany

---

## Marker-Trait Summary

| Marker | Trait | Chr | Type | Primary Alleles | Key Reference |
|--------|-------|-----|------|-----------------|---------------|
| MYB10 | Flesh Color | 9 | SSR | 390, 490 | Espley et al., 2007 |
| RED_TE | Skin Color | 9 | SCAR | 750, NB | Chagné et al., 2013 |
| MA_INDEL | Acidity | 16 | SCAR | 409, 456 | Bai et al., 2012 |
| BP16 | Bitter Pit | 8 | SSR | 203, 368 | Buti et al., 2015 |
| BP13 | Bitter Pit | 8 | SSR | 224-260 | Buti et al., 2015 |
| ACS | Texture | 10 | SCAR | 200, 341 | Sunako et al., 1999 |
| ACO | Texture | 10 | SCAR | 237, 300 | Costa et al., 2005 |
| MD_PG1 | Texture | 10 | SSR | 289, 292, 298 | Longhi et al., 2013 |
