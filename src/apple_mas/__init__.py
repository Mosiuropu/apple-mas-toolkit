"""
Apple MAS Toolkit - Marker-Assisted Selection for Apple Breeding Programs

A comprehensive toolkit for analyzing molecular markers (SSR and SCAR) in apple
germplasm to support parent selection for fruit quality traits.
"""

__version__ = "1.0.0"
__author__ = "Md Mosiur Rahman Bhuyin Apu"

from apple_mas.marker_db import MarkerDatabase
from apple_mas.allele_analysis import AlleleAnalyzer
from apple_mas.parent_scorer import ParentScorer
from apple_mas.visualization import MarkerVisualizer
from apple_mas.data_parser import DataParser

__all__ = [
    "MarkerDatabase",
    "AlleleAnalyzer",
    "ParentScorer",
    "MarkerVisualizer",
    "DataParser",
]
