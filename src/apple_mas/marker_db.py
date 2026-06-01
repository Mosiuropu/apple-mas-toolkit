"""
Marker Database - Access and query apple molecular marker information.

Loads the comprehensive marker database and provides methods to query
markers by trait, chromosome, or name. Includes allele definitions,
breeding recommendations, and PCR methodology.
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Any


class MarkerDatabase:
    """Access and query the apple molecular marker database.

    Provides methods to retrieve marker information, filter by trait category,
    look up allele definitions, and access breeding recommendations.

    Examples
    --------
    >>> db = MarkerDatabase()
    >>> db.list_markers()
    ['MYB10', 'RED_TE', 'MA_INDEL', 'BP16', 'BP13', 'ACS', 'ACO', 'MD_PG1']
    >>> info = db.get_marker("MYB10")
    >>> print(info["trait"])
    Flesh Color
    """

    def __init__(self, db_path: Optional[str] = None):
        """Initialize the marker database.

        Parameters
        ----------
        db_path : str, optional
            Path to a custom marker database JSON file. If None, uses the
            built-in database shipped with the package.
        """
        if db_path is None:
            db_path = str(
                Path(__file__).parent.parent.parent
                / "data"
                / "marker_database"
                / "apple_markers.json"
            )
        self._db_path = db_path
        self._db = self._load_database()

    def _load_database(self) -> Dict[str, Any]:
        """Load the marker database from JSON."""
        with open(self._db_path, "r") as f:
            return json.load(f)

    @property
    def version(self) -> str:
        """Database version string."""
        return self._db.get("version", "unknown")

    def list_markers(self) -> List[str]:
        """Return all marker names in the database.

        Returns
        -------
        list of str
            Marker identifiers (e.g., ['MYB10', 'RED_TE', ...]).
        """
        return list(self._db["markers"].keys())

    def get_marker(self, marker_name: str) -> Dict[str, Any]:
        """Get full information for a single marker.

        Parameters
        ----------
        marker_name : str
            Marker identifier (e.g., 'MYB10', 'ACS').

        Returns
        -------
        dict
            Complete marker information including alleles, references, and notes.

        Raises
        ------
        KeyError
            If the marker name is not found in the database.
        """
        markers = self._db["markers"]
        if marker_name not in markers:
            available = ", ".join(markers.keys())
            raise KeyError(f"Marker '{marker_name}' not found. Available: {available}")
        return markers[marker_name]

    def get_markers_by_trait(self, trait_category: str) -> List[str]:
        """Get all markers associated with a specific trait category.

        Parameters
        ----------
        trait_category : str
            One of: 'color', 'texture', 'acidity', 'disease_resistance'.

        Returns
        -------
        list of str
            Marker names associated with the given trait category.
        """
        groups = self._db.get("trait_groups", {})
        if trait_category not in groups:
            available = ", ".join(groups.keys())
            raise KeyError(
                f"Trait category '{trait_category}' not found. Available: {available}"
            )
        return groups[trait_category]["markers"]

    def get_markers_by_chromosome(self, chromosome: str) -> List[str]:
        """Get all markers located on a specific chromosome.

        Parameters
        ----------
        chromosome : str
            Chromosome number as string (e.g., '8', '9', '10', '16').

        Returns
        -------
        list of str
            Marker names on the given chromosome.
        """
        chrom_map = self._db.get("chromosomal_locations", {})
        return chrom_map.get(str(chromosome), [])

    def get_trait_categories(self) -> List[str]:
        """List all available trait categories.

        Returns
        -------
        list of str
            Trait category names.
        """
        return list(self._db.get("trait_groups", {}).keys())

    def get_trait_description(self, trait_category: str) -> str:
        """Get description for a trait category.

        Parameters
        ----------
        trait_category : str
            Trait category name.

        Returns
        -------
        str
            Description of the trait category.
        """
        groups = self._db.get("trait_groups", {})
        if trait_category not in groups:
            raise KeyError(f"Trait category '{trait_category}' not found.")
        return groups[trait_category]["description"]

    def get_allele_info(self, marker_name: str, allele_size: Any) -> Dict[str, Any]:
        """Get information for a specific allele of a marker.

        Parameters
        ----------
        marker_name : str
            Marker identifier.
        allele_size : str or int
            Allele size (e.g., 390, '750', 'NB').

        Returns
        -------
        dict
            Allele definition including phenotype and frequency.

        Raises
        ------
        KeyError
            If marker or allele is not found.
        """
        marker = self.get_marker(marker_name)
        alleles = marker.get("allele_definitions", {})
        allele_key = str(allele_size)
        if allele_key not in alleles:
            available = ", ".join(str(k) for k in alleles.keys())
            raise KeyError(
                f"Allele '{allele_size}' not found for marker '{marker_name}'. "
                f"Available alleles: {available}"
            )
        return alleles[allele_key]

    def get_genotype_phenotype(self, marker_name: str, genotype: str) -> Dict[str, Any]:
        """Get phenotype interpretation for a specific genotype.

        Parameters
        ----------
        marker_name : str
            Marker identifier.
        genotype : str
            Genotype string (e.g., '390:390', '750', '200:341').

        Returns
        -------
        dict
            Phenotype interpretation for the genotype.
        """
        marker = self.get_marker(marker_name)
        interpretations = marker.get("genotype_interpretation", {})
        if genotype not in interpretations:
            available = ", ".join(interpretations.keys())
            raise KeyError(
                f"Genotype '{genotype}' not found for marker '{marker_name}'. "
                f"Available: {available}"
            )
        return interpretations[genotype]

    def get_breeding_notes(self, marker_name: str) -> str:
        """Get breeding recommendations for a marker.

        Parameters
        ----------
        marker_name : str
            Marker identifier.

        Returns
        -------
        str
            Breeding notes and recommendations.
        """
        marker = self.get_marker(marker_name)
        return marker.get("breeding_notes", "No breeding notes available.")

    def get_pcr_methodology(self) -> Dict[str, Any]:
        """Get the standard PCR methodology used for marker analysis.

        Returns
        -------
        dict
            PCR conditions, thermal cycling parameters, and fragment analysis details.
        """
        return self._db.get("methodology", {})

    def search_markers(self, query: str) -> List[str]:
        """Search markers by keyword in name, trait, or notes.

        Parameters
        ----------
        query : str
            Search term (case-insensitive).

        Returns
        -------
        list of str
            Matching marker names.
        """
        query_lower = query.lower()
        results = []
        for name, info in self._db["markers"].items():
            searchable = " ".join([
                name,
                info.get("full_name", ""),
                info.get("trait", ""),
                info.get("notes", ""),
                info.get("breeding_notes", ""),
            ]).lower()
            if query_lower in searchable:
                results.append(name)
        return results

    def summary(self) -> str:
        """Print a summary of the marker database.

        Returns
        -------
        str
            Formatted summary string.
        """
        lines = [
            "=" * 60,
            "  Apple Molecular Marker Database",
            f"  Version: {self.version}",
            f"  Organism: {self._db.get('organism', 'Unknown')}",
            "=" * 60,
            "",
            f"  Total markers: {len(self.list_markers())}",
            "",
            "  Markers by trait category:",
        ]
        for category, info in self._db.get("trait_groups", {}).items():
            markers_str = ", ".join(info["markers"])
            lines.append(f"    {category:20s} -> {markers_str}")
            lines.append(f"    {'':20s}    {info['description'][:60]}...")
        lines.append("")
        lines.append("  Chromosomal locations:")
        for chrom, markers in self._db.get("chromosomal_locations", {}).items():
            lines.append(f"    Chr {chrom:2s}: {', '.join(markers)}")
        lines.append("")
        lines.append("=" * 60)
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"MarkerDatabase(version='{self.version}', "
            f"markers={len(self.list_markers())})"
        )
