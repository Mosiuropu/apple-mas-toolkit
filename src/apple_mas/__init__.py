"""
Apple MAS Toolkit - Marker-Assisted Selection for Apple Breeding Programs

A comprehensive toolkit for analyzing molecular markers (SSR, SNP, SCAR,
and InDel) in apple (Malus x domestica) germplasm to support marker-assisted
selection for fruit quality, disease resistance, and growth habit traits.
"""

__version__ = "2.0.0"
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
