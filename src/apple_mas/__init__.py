"""
Apple MAS Toolkit — Marker-Assisted Selection for Apple Breeding.

A comprehensive toolkit for analysing molecular markers (SSR, SNP, SCAR,
and InDel) in apple (*Malus domestica* Borkh.) to support marker-assisted
selection for fruit quality, disease resistance, and growth habit traits.

Modules
-------
- ``apple_mas.marker_db`` — Query the legacy marker database.
- ``apple_mas.selector`` — Interactive trait-to-marker selection workflow.
- ``apple_mas.data_parser`` — Import and standardise genotype data.
- ``apple_mas.allele_analysis`` — Allele frequency and diversity analysis.
- ``apple_mas.parent_scorer`` — Score and rank accessions for breeding.
- ``apple_mas.visualization`` — Generate publication-quality plots.
- ``apple_mas.analytics`` — Segregation distortion, PIC, and boxplots.

Standalone modules (imported from ``src/``):
- ``marker_selector`` — Chromosome-aware marker query and primer manifest.
- ``breeding_analytics`` — Statistical analysis and visualisation suite.
"""

__version__ = "3.0.0"
__author__ = "Md Mosiur Rahman Bhuyin Apu"

from apple_mas.marker_db import MarkerDatabase
from apple_mas.allele_analysis import AlleleAnalyzer
from apple_mas.parent_scorer import ParentScorer
from apple_mas.visualization import MarkerVisualizer
from apple_mas.data_parser import DataParser
from apple_mas.selector import MarkerSelector
from apple_mas.analytics import (
    segregation_distortion_test,
    segregation_distortion_from_df,
    genotype_phenotype_boxplot,
    calculate_pic,
    calculate_pic_from_df,
    batch_segregation_test,
    batch_pic,
)

__all__ = [
    "MarkerDatabase",
    "AlleleAnalyzer",
    "ParentScorer",
    "MarkerVisualizer",
    "DataParser",
    "MarkerSelector",
    "segregation_distortion_test",
    "segregation_distortion_from_df",
    "genotype_phenotype_boxplot",
    "calculate_pic",
    "calculate_pic_from_df",
    "batch_segregation_test",
    "batch_pic",
]
