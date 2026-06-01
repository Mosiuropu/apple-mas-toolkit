"""
Selector - Interactive Marker-Assisted Selection Workflow.

Provides a structured workflow interface for students and researchers to:
- Map agronomic traits to validated molecular markers
- Query the marker database for primer sequences, linkage groups, and allele sizes
- Generate tailored genotyping plan configuration files
"""

import json
import os
from collections import OrderedDict
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from apple_mas.marker_db import MarkerDatabase


class MarkerSelector:
    """Interactive workflow interface for selecting markers by target trait.

    Guides the user through trait selection, marker retrieval, and genotyping
    plan generation. Designed for students and breeders beginning a
    marker-assisted selection project.

    Parameters
    ----------
    db_path : str, optional
        Path to a custom marker database JSON. Uses the built-in database
        if None.

    Examples
    --------
    >>> selector = MarkerSelector()
    >>> plan = selector.build_genotyping_plan(
    ...     target_traits=["fruit_acidity", "disease_resistance"],
    ...     population_type="F1",
    ... )
    >>> selector.export_plan(plan, "my_genotyping_plan.json")
    """

    def __init__(self, db_path: Optional[str] = None):
        self._db = MarkerDatabase(db_path)
        self._trait_map = self._build_trait_index()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def list_available_traits(self) -> pd.DataFrame:
        """Return a summary table of all available trait categories.

        Returns
        -------
        pd.DataFrame
            Columns: ['category', 'marker_count', 'markers', 'description'].
        """
        rows = []
        for cat_name, cat_info in self._db._db.get("trait_groups", {}).items():
            markers = cat_info.get("markers", [])
            rows.append({
                "category": cat_name,
                "marker_count": len(markers),
                "markers": ", ".join(markers),
                "description": cat_info.get("description", ""),
            })
        return pd.DataFrame(rows)

    def list_available_markers(self) -> pd.DataFrame:
        """Return a summary table of all markers in the database.

        Returns
        -------
        pd.DataFrame
            Columns: ['marker', 'full_name', 'trait', 'chromosome', 'lg',
            'marker_type', 'n_alleles'].
        """
        rows = []
        for marker_name in self._db.list_markers():
            info = self._db.get_marker(marker_name)
            rows.append({
                "marker": marker_name,
                "full_name": info.get("full_name", ""),
                "trait": info.get("trait", ""),
                "chromosome": info.get("chromosome", ""),
                "lg": info.get("lg", ""),
                "marker_type": info.get("marker_type", ""),
                "n_alleles": len(info.get("allele_sizes_bp", [])),
            })
        return pd.DataFrame(rows)

    def search_markers(self, query: str) -> List[str]:
        """Search markers by keyword.

        Parameters
        ----------
        query : str
            Search term (case-insensitive). Matches marker names, traits,
            notes, and breeding notes.

        Returns
        -------
        list of str
            Matching marker identifiers.
        """
        return self._db.search_markers(query)

    def get_markers_for_trait(self, trait_category: str) -> Dict[str, Any]:
        """Retrieve all markers associated with a trait category.

        Parameters
        ----------
        trait_category : str
            One of: 'color', 'texture', 'acidity', 'disease_resistance',
            'growth_habit', 'fruit_quality'.

        Returns
        -------
        dict
            Marker details keyed by marker name.

        Raises
        ------
        KeyError
            If the trait category is not found.
        """
        marker_names = self._db.get_markers_by_trait(trait_category)
        return {name: self._db.get_marker(name) for name in marker_names}

    def get_marker_primer_info(self, marker_name: str) -> Dict[str, Any]:
        """Retrieve primer sequences and PCR-relevant details for a marker.

        Parameters
        ----------
        marker_name : str
            Marker identifier.

        Returns
        -------
        dict
            Keys: 'marker', 'primer_sequences', 'dye', 'allele_sizes_bp',
            'marker_type'.
        """
        info = self._db.get_marker(marker_name)
        return {
            "marker": marker_name,
            "primer_sequences": info.get("primer_sequences", {}),
            "dye": info.get("dye", "N/A"),
            "allele_sizes_bp": info.get("allele_sizes_bp", []),
            "marker_type": info.get("marker_type", ""),
        }

    def get_expected_alleles(self, marker_name: str) -> Dict[str, Any]:
        """Retrieve allele definitions with phenotype mappings.

        Parameters
        ----------
        marker_name : str
            Marker identifier.

        Returns
        -------
        dict
            Allele definitions keyed by allele size.
        """
        info = self._db.get_marker(marker_name)
        return info.get("allele_definitions", {})

    def get_marker_plan_table(self, marker_names: List[str]) -> pd.DataFrame:
        """Build a formatted table summarizing markers for a genotyping plan.

        Parameters
        ----------
        marker_names : list of str
            Markers to include in the plan.

        Returns
        -------
        pd.DataFrame
            Columns: ['marker', 'trait', 'chromosome', 'lg', 'marker_type',
            'primer_forward', 'primer_reverse', 'dye', 'allele_range_bp',
            'key_alleles'].
        """
        rows = []
        for name in marker_names:
            try:
                info = self._db.get_marker(name)
                primers = info.get("primer_sequences", {})
                alleles = info.get("allele_sizes_bp", [])
                allele_defs = info.get("allele_definitions", {})
                key_alleles = ", ".join(
                    f"{k}: {v.get('phenotype', '')}"
                    for k, v in allele_defs.items()
                    if v.get("desirability") in ("desirable", "desirable_for_tart", "desirable_for_sweet", "desirable_for_high_density")
                )
                rows.append({
                    "marker": name,
                    "trait": info.get("trait", ""),
                    "chromosome": info.get("chromosome", ""),
                    "lg": info.get("lg", ""),
                    "marker_type": info.get("marker_type", ""),
                    "primer_forward": primers.get("forward", ""),
                    "primer_reverse": primers.get("reverse", ""),
                    "dye": info.get("dye", "N/A"),
                    "allele_range_bp": f"{min(alleles)}-{max(alleles)}" if alleles else "N/A",
                    "key_alleles": key_alleles,
                })
            except KeyError:
                continue
        return pd.DataFrame(rows)

    def build_genotyping_plan(
        self,
        target_traits: List[str],
        population_type: str = "F1",
        population_size: Optional[int] = None,
        additional_markers: Optional[List[str]] = None,
        notes: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Construct a genotyping plan for a breeding project.

        Parameters
        ----------
        target_traits : list of str
            Trait categories to include. Examples: 'fruit_quality',
            'disease_resistance', 'growth_habit', 'acidity', 'color',
            'texture'.
        population_type : str
            Type of segregating population. One of: 'F1', 'BC1', 'BC2',
            'F2', 'RIL', 'open_pollinated'.
        population_size : int, optional
            Expected number of individuals to genotype.
        additional_markers : list of str, optional
            Extra marker names to include beyond those mapped to traits.
        notes : str, optional
            Free-text notes for the plan.

        Returns
        -------
        dict
            Complete genotyping plan with marker details, primer info,
            and methodology recommendations.

        Raises
        ------
        ValueError
            If no markers are found for the specified traits.
        """
        marker_names = set()
        for trait in target_traits:
            try:
                names = self._db.get_markers_by_trait(trait)
                marker_names.update(names)
            except KeyError:
                continue

        if additional_markers:
            marker_names.update(additional_markers)

        if not marker_names:
            raise ValueError(
                f"No markers found for traits: {target_traits}. "
                f"Available categories: {list(self._db._db.get('trait_groups', {}).keys())}"
            )

        sorted_markers = sorted(marker_names)

        marker_details = {}
        for name in sorted_markers:
            info = self._db.get_marker(name)
            primers = info.get("primer_sequences", {})
            marker_details[name] = {
                "full_name": info.get("full_name", ""),
                "trait": info.get("trait", ""),
                "trait_category": info.get("trait_category", ""),
                "chromosome": info.get("chromosome", ""),
                "lg": info.get("lg", ""),
                "marker_type": info.get("marker_type", ""),
                "primer_forward": primers.get("forward", ""),
                "primer_reverse": primers.get("reverse", ""),
                "dye": info.get("dye", "N/A"),
                "allele_sizes_bp": info.get("allele_sizes_bp", []),
                "allele_definitions": info.get("allele_definitions", {}),
                "references": info.get("references", []),
                "breeding_notes": info.get("breeding_notes", ""),
            }

        methodology = self._db.get_pcr_methodology()

        plan = {
            "plan_title": "Genotyping Plan - Apple MAS Toolkit",
            "version": "2.0",
            "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "organism": "Malus domestica Borkh.",
            "population": {
                "type": population_type,
                "size": population_size,
                "segregation_notes": self._segregation_notes(population_type),
            },
            "target_traits": target_traits,
            "markers": marker_details,
            "marker_count": len(marker_details),
            "chromosomes_covered": sorted(
                set(str(m.get("chromosome", "")) for m in marker_details.values())
            ),
            "methodology": methodology,
            "notes": notes or "",
        }

        return plan

    def export_plan(
        self,
        plan: Dict[str, Any],
        output_path: str,
        fmt: str = "json",
    ) -> str:
        """Export a genotyping plan to file.

        Parameters
        ----------
        plan : dict
            Genotyping plan from ``build_genotyping_plan()``.
        output_path : str
            Destination file path.
        fmt : str
            Output format: 'json' or 'csv' (CSV exports the marker table).

        Returns
        -------
        str
            Absolute path to the exported file.
        """
        os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)

        if fmt == "csv":
            marker_names = list(plan.get("markers", {}).keys())
            table = self.get_marker_plan_table(marker_names)
            table.to_csv(output_path, index=False)
        else:
            with open(output_path, "w") as fh:
                json.dump(plan, fh, indent=2, default=str)

        return os.path.abspath(output_path)

    def plan_summary(self, plan: Dict[str, Any]) -> str:
        """Return a human-readable summary of a genotyping plan.

        Parameters
        ----------
        plan : dict
            Genotyping plan from ``build_genotyping_plan()``.

        Returns
        -------
        str
            Formatted summary string.
        """
        lines = [
            "=" * 60,
            "  GENOTYPING PLAN SUMMARY",
            "=" * 60,
            f"  Generated: {plan.get('generated', 'N/A')}",
            f"  Population: {plan.get('population', {}).get('type', 'N/A')}",
            f"  Population size: {plan.get('population', {}).get('size', 'Not specified')}",
            f"  Target traits: {', '.join(plan.get('target_traits', []))}",
            f"  Markers selected: {plan.get('marker_count', 0)}",
            f"  Chromosomes covered: {', '.join(plan.get('chromosomes_covered', []))}",
            "",
        ]

        for name, detail in plan.get("markers", {}).items():
            lines.append(f"  {name}")
            lines.append(f"    Trait: {detail.get('trait', 'N/A')}")
            lines.append(f"    Chr: {detail.get('chromosome', '?')}  LG: {detail.get('lg', '?')}")
            lines.append(f"    Type: {detail.get('marker_type', 'N/A')}  Dye: {detail.get('dye', 'N/A')}")
            fwd = detail.get("primer_forward", "")
            rev = detail.get("primer_reverse", "")
            if fwd and rev:
                lines.append(f"    Forward: {fwd}")
                lines.append(f"    Reverse: {rev}")
            alleles = detail.get("allele_sizes_bp", [])
            lines.append(f"    Allele range: {min(alleles)}-{max(alleles)} bp" if alleles else "    Allele range: N/A")
            lines.append("")

        lines.append("=" * 60)
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_trait_index(self) -> Dict[str, List[str]]:
        """Build an internal index mapping trait keywords to marker names."""
        index: Dict[str, List[str]] = {}
        for group_name, group_info in self._db._db.get("trait_groups", {}).items():
            index[group_name] = group_info.get("markers", [])
        return index

    @staticmethod
    def _segregation_notes(population_type: str) -> str:
        """Return expected segregation ratios for a population type."""
        notes = {
            "F1": "All loci should be heterozygous if parents are homozygous for different alleles. Useful for confirming parentage.",
            "BC1": "Expected 1:1 segregation ratio for each locus (homozygous vs. heterozygous). Suitable for backcross breeding.",
            "BC2": "Expected 3:1 ratio (homozygous recurrent parent : heterozygous). Marker-assisted backcross selection.",
            "F2": "Expected 1:2:1 genotypic ratio for co-dominant markers. Suitable for mapping and selection experiments.",
            "RIL": "Recombinant inbred lines; expected near-homozygous. Fixed genotypes suitable for replicated phenotyping.",
            "open_pollinated": "Segregation ratios depend on parental allele frequencies. Population structure analysis recommended.",
        }
        return notes.get(population_type, "Standard Mendelian segregation expected.")
