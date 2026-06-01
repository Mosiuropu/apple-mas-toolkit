"""
Marker Selector - Interactive student query tool for apple molecular markers.

Enables students and researchers to filter markers by chromosome, trait,
or keyword, and generate formatted genotyping primer manifests for
experimental assay design.

References:
    - GDDH13 v1.1 reference genome (Daccord et al., 2014)
    - Genome Database for Rosaceae (GDR) nomenclature conventions
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional
from datetime import datetime


class MarkerSelector:
    """Query the apple marker database and generate genotyping plans.

    Provides methods to filter markers by chromosome, trait, or keyword,
    retrieve primer sequences, and produce formatted primer manifests
    suitable for laboratory assay design.

    Parameters
    ----------
    db_path : str or None
        Path to a custom marker database JSON.  If ``None`` the built-in
        ``data/apple_marker_map.json`` shipped with the package is used.

    Examples
    --------
    >>> selector = MarkerSelector()
    >>> selector.list_chromosomes()
    >>> markers = selector.markers_by_chromosome("10")
    >>> manifest = selector.genotyping_manifest(
    ...     traits=["fruit_quality", "disease_resistance"],
    ...     population_type="BC1",
    ... )
    """

    def __init__(self, db_path: Optional[str] = None) -> None:
        if db_path is None:
            db_path = str(
                Path(__file__).resolve().parent.parent
                / "data"
                / "apple_marker_map.json"
            )
        self._db_path = db_path
        self._db = self._load_database()

    # ------------------------------------------------------------------
    # Database loading
    # ------------------------------------------------------------------

    def _load_database(self) -> Dict[str, Any]:
        with open(self._db_path, "r") as fh:
            return json.load(fh)

    @property
    def version(self) -> str:
        return self._db.get("version", "unknown")

    @property
    def organism(self) -> str:
        return self._db.get("organism", "Malus domestica")

    @property
    def reference_genome(self) -> str:
        return self._db.get("reference_genome", "GDDH13 v1.1")

    # ------------------------------------------------------------------
    # Chromosome queries
    # ------------------------------------------------------------------

    def list_chromosomes(self) -> List[Dict[str, Any]]:
        """Return metadata for all 17 apple chromosomes.

        Returns
        -------
        list of dict
            Each entry contains ``chromosome``, ``linkage_group``,
            ``approximate_length_Mb``, ``major_loci``, and ``description``.
        """
        results = []
        for chrom_num, info in sorted(
            self._db.get("chromosomes", {}).items(),
            key=lambda x: int(x[0]),
        ):
            results.append({
                "chromosome": chrom_num,
                "linkage_group": info.get("linkage_group", ""),
                "approximate_length_Mb": info.get("approximate_length_Mb", 0),
                "major_loci": info.get("major_loci", []),
                "description": info.get("description", ""),
            })
        return results

    def get_chromosome_info(self, chromosome: str) -> Dict[str, Any]:
        """Return detailed information for a single chromosome.

        Parameters
        ----------
        chromosome : str
            Chromosome number as a string, e.g. ``"1"`` or ``"16"``.

        Returns
        -------
        dict
            Chromosome metadata including linkage group, length, and loci.

        Raises
        ------
        KeyError
            If the chromosome number is not found.
        """
        chroms = self._db.get("chromosomes", {})
        if chromosome not in chroms:
            available = ", ".join(sorted(chroms.keys(), key=int))
            raise KeyError(
                f"Chromosome '{chromosome}' not found. Available: {available}"
            )
        info = chroms[chromosome]
        info["chromosome"] = chromosome
        return info

    # ------------------------------------------------------------------
    # Marker queries
    # ------------------------------------------------------------------

    def list_markers(self) -> List[str]:
        """Return all marker names in the database."""
        return sorted(self._db.get("markers", {}).keys())

    def get_marker(self, marker_name: str) -> Dict[str, Any]:
        """Retrieve full information for a single marker.

        Parameters
        ----------
        marker_name : str
            Marker identifier, e.g. ``"AL07"`` or ``"Ma1_SNP"``.

        Returns
        -------
        dict
            Complete marker record.

        Raises
        ------
        KeyError
            If the marker is not in the database.
        """
        markers = self._db.get("markers", {})
        if marker_name not in markers:
            available = ", ".join(sorted(markers.keys()))
            raise KeyError(
                f"Marker '{marker_name}' not found. Available: {available}"
            )
        return markers[marker_name]

    def markers_by_chromosome(self, chromosome: str) -> List[Dict[str, Any]]:
        """Return all validated markers located on a given chromosome.

        Parameters
        ----------
        chromosome : str
            Chromosome number (``"1"``–``"17"``).

        Returns
        -------
        list of dict
            Marker records sorted by genetic position (cM).
        """
        results = []
        for name, info in self._db.get("markers", {}).items():
            if str(info.get("chromosome", "")) == str(chromosome):
                entry = dict(info)
                entry["marker_name"] = name
                results.append(entry)
        results.sort(key=lambda x: x.get("position_cM", 0))
        return results

    def markers_by_trait(self, trait_category: str) -> List[Dict[str, Any]]:
        """Return all markers associated with a trait category.

        Parameters
        ----------
        trait_category : str
            One of the categories defined in ``trait_groups``:
            ``fruit_quality``, ``disease_resistance``, ``tree_architecture``,
            ``acidity``, ``texture``, ``color``, ``scab_resistance``,
            ``mildew_resistance``.

        Returns
        -------
        list of dict
            Marker records for the requested trait category.

        Raises
        ------
        KeyError
            If the trait category is not defined.
        """
        groups = self._db.get("trait_groups", {})
        if trait_category not in groups:
            available = ", ".join(sorted(groups.keys()))
            raise KeyError(
                f"Trait category '{trait_category}' not found. "
                f"Available: {available}"
            )
        marker_names = groups[trait_category].get("markers", [])
        return [
            {**self.get_marker(m), "marker_name": m}
            for m in marker_names
            if m in self._db.get("markers", {})
        ]

    def search_markers(self, query: str) -> List[Dict[str, Any]]:
        """Search markers by keyword across names, traits, and notes.

        Parameters
        ----------
        query : str
            Case-insensitive search term.

        Returns
        -------
        list of dict
            Matching marker records.
        """
        query_lower = query.lower()
        results = []
        for name, info in self._db.get("markers", {}).items():
            searchable = " ".join([
                name,
                info.get("full_name", ""),
                info.get("target_trait", ""),
                info.get("locus", ""),
                info.get("notes", ""),
                info.get("trait_category", ""),
            ]).lower()
            if query_lower in searchable:
                entry = dict(info)
                entry["marker_name"] = name
                results.append(entry)
        return results

    def list_trait_categories(self) -> List[Dict[str, str]]:
        """Return all defined trait categories with descriptions.

        Returns
        -------
        list of dict
            Each entry has ``category``, ``description``, and ``markers``.
        """
        results = []
        for cat_name, cat_info in self._db.get("trait_groups", {}).items():
            results.append({
                "category": cat_name,
                "description": cat_info.get("description", ""),
                "markers": cat_info.get("markers", []),
            })
        return results

    # ------------------------------------------------------------------
    # Genotyping manifest generation
    # ------------------------------------------------------------------

    def build_marker_table(
        self,
        marker_names: List[str],
    ) -> List[Dict[str, Any]]:
        """Build a formatted table summarizing markers for a genotyping plan.

        Parameters
        ----------
        marker_names : list of str
            Markers to include.

        Returns
        -------
        list of dict
            Rows with ``Marker_Name``, ``Type``, ``Chromosome``,
            ``Linkage_Group``, ``Position_cM``, ``Forward_Primer``,
            ``Reverse_Primer``, ``Dye``, ``Target_Trait``, and
            ``Allele_Range``.
        """
        rows = []
        for name in marker_names:
            try:
                info = self.get_marker(name)
            except KeyError:
                continue
            alleles = info.get("allele_sizes_bp", [])
            allele_range = (
                f"{min(alleles)}–{max(alleles)} bp"
                if alleles
                else "N/A (SNP/KASP)"
            )
            rows.append({
                "Marker_Name": name,
                "Type": info.get("type", ""),
                "Chromosome": info.get("chromosome", ""),
                "Linkage_Group": info.get("linkage_group", ""),
                "Position_cM": info.get("position_cM", ""),
                "Forward_Primer": info.get("forward_primer", ""),
                "Reverse_Primer": info.get("reverse_primer", ""),
                "Dye": info.get("dye", ""),
                "Target_Trait": info.get("target_trait", ""),
                "Allele_Range": allele_range,
            })
        return rows

    def genotyping_manifest(
        self,
        traits: Optional[List[str]] = None,
        chromosomes: Optional[List[str]] = None,
        marker_names: Optional[List[str]] = None,
        population_type: str = "F1",
        population_size: Optional[int] = None,
        title: str = "Genotyping Primer Manifest",
        notes: str = "",
    ) -> str:
        """Generate a formatted Markdown primer manifest.

        This method produces a complete, laboratory-ready document listing
        all primers, expected allele sizes, and PCR conditions for the
        selected markers.

        Parameters
        ----------
        traits : list of str, optional
            Trait categories to include (e.g. ``["fruit_quality"]``).
        chromosomes : list of str, optional
            Chromosome numbers to include (e.g. ``["1", "10"]``).
        marker_names : list of str, optional
            Explicit marker names to include.
        population_type : str
            Population type: ``F1``, ``BC1``, ``BC2``, ``F2``, ``RIL``,
            or ``open_pollinated``.
        population_size : int, optional
            Expected number of individuals.
        title : str
            Title for the manifest document.
        notes : str
            Free-text notes to append.

        Returns
        -------
        str
            Markdown-formatted primer manifest.
        """
        # Collect unique marker names from all sources
        selected = set()

        if traits:
            for trait in traits:
                try:
                    for m in self.markers_by_trait(trait):
                        selected.add(m["marker_name"])
                except KeyError:
                    continue

        if chromosomes:
            for chrom in chromosomes:
                for m in self.markers_by_chromosome(chrom):
                    selected.add(m["marker_name"])

        if marker_names:
            selected.update(marker_names)

        if not selected:
            # Default: all markers
            selected = set(self.list_markers())

        marker_table = self.build_marker_table(sorted(selected))

        # Build manifest
        seg_notes = self._segregation_notes(population_type)
        lines = [
            f"# {title}",
            "",
            f"**Date:** {datetime.now().strftime('%Y-%m-%d')}",
            f"**Organism:** {self.organism}",
            f"**Reference Genome:** {self.reference_genome}",
            f"**Population Type:** {population_type}",
        ]
        if population_size:
            lines.append(f"**Population Size:** {population_size}")
        lines.append(f"**Markers Selected:** {len(marker_table)}")
        lines.append("")
        lines.append("---")
        lines.append("")

        # Marker table
        lines.append("## Primer Summary")
        lines.append("")
        lines.append(
            "| Marker | Type | Chr | LG | cM | Forward Primer | Reverse Primer | "
            "Dye | Trait | Allele Range |"
        )
        lines.append(
            "|--------|------|-----|----|----|----------------|----------------|"
            "-----|-------|-------------|"
        )
        for row in marker_table:
            lines.append(
                f"| {row['Marker_Name']} | {row['Type']} | {row['Chromosome']} "
                f"| {row['Linkage_Group']} | {row['Position_cM']} "
                f"| `{row['Forward_Primer']}` | `{row['Reverse_Primer']}` "
                f"| {row['Dye']} | {row['Target_Trait']} | {row['Allele_Range']} |"
            )
        lines.append("")

        # Allele definitions
        lines.append("## Allele Definitions")
        lines.append("")
        for row in marker_table:
            try:
                info = self.get_marker(row["Marker_Name"])
            except KeyError:
                continue
            alleles = info.get("associated_alleles", {})
            if not alleles:
                continue
            lines.append(f"### {row['Marker_Name']} ({row['Target_Trait']})")
            lines.append("")
            lines.append("| Allele | Phenotype | Desirability |")
            lines.append("|--------|-----------|--------------|")
            for allele, defn in alleles.items():
                pheno = defn.get("phenotype", "")
                desir = defn.get("desirability", "")
                lines.append(f"| {allele} | {pheno} | {desir} |")
            lines.append("")

        # PCR conditions
        lines.append("## PCR Conditions")
        lines.append("")
        pcr = self._db.get("methodology", {}).get("PCR_conditions", {})
        if pcr:
            lines.append(f"- **Reaction Volume:** {pcr.get('total_volume_ul', 25)} µL")
            comps = pcr.get("components", {})
            for k, v in comps.items():
                label = k.replace("_", " ").replace("uM", " µM").replace("mM", " mM")
                lines.append(f"  - {label}: {v}")
            tc = pcr.get("thermal_cycling", {})
            if tc:
                lines.append(f"- **Cycles:** {tc.get('cycles', 30)}")
                init = tc.get("initial_denaturation", {})
                lines.append(
                    f"  - Initial denaturation: {init.get('temp_C', 95)}°C / "
                    f"{init.get('time_min', 5)} min"
                )
                denat = tc.get("denaturation", {})
                lines.append(
                    f"  - Denaturation: {denat.get('temp_C', 94)}°C / "
                    f"{denat.get('time_sec', 30)} sec"
                )
                ann = tc.get("annealing", {})
                lines.append(
                    f"  - Annealing: {ann.get('temp_C', 55)}°C / "
                    f"{ann.get('time_sec', 30)} sec"
                )
                ext = tc.get("extension", {})
                lines.append(
                    f"  - Extension: {ext.get('temp_C', 72)}°C / "
                    f"{ext.get('time_sec', 60)} sec"
                )
                fin = tc.get("final_extension", {})
                lines.append(
                    f"  - Final extension: {fin.get('temp_C', 72)}°C / "
                    f"{fin.get('time_min', 10)} min"
                )
        lines.append("")

        # Segregation notes
        lines.append("## Segregation Expectations")
        lines.append("")
        lines.append(f"**Population Type:** {population_type}")
        lines.append(f"{seg_notes}")
        lines.append("")

        if notes:
            lines.append("## Additional Notes")
            lines.append("")
            lines.append(notes)
            lines.append("")

        lines.append("---")
        lines.append("")
        lines.append(
            "*Generated by the Apple MAS Toolkit — "
            "see https://github.com/Mosiuropu/apple-mas-toolkit-hidden*"
        )

        return "\n".join(lines)

    def export_manifest(
        self,
        output_path: str,
        **kwargs: Any,
    ) -> str:
        """Write a genotyping manifest to a Markdown file.

        Parameters
        ----------
        output_path : str
            Destination file path (``.md`` recommended).
        **kwargs
            Forwarded to :meth:`genotyping_manifest`.

        Returns
        -------
        str
            Absolute path to the written file.
        """
        manifest = self.genotyping_manifest(**kwargs)
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        Path(output_path).write_text(manifest, encoding="utf-8")
        return str(Path(output_path).resolve())

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _segregation_notes(population_type: str) -> str:
        notes = {
            "F1": (
                "All loci should be heterozygous if parents are homozygous "
                "for different alleles. Useful for confirming parentage."
            ),
            "BC1": (
                "Expected 1:1 segregation ratio per locus "
                "(homozygous vs. heterozygous)."
            ),
            "BC2": (
                "Expected 3:1 ratio "
                "(homozygous recurrent parent : heterozygous)."
            ),
            "F2": (
                "Expected 1:2:1 genotypic ratio for co-dominant markers. "
                "Suitable for mapping and selection experiments."
            ),
            "RIL": (
                "Recombinant inbred lines; expected near-homozygous. "
                "Fixed genotypes suitable for replicated phenotyping."
            ),
            "open_pollinated": (
                "Segregation ratios depend on parental allele frequencies. "
                "Population structure analysis recommended."
            ),
        }
        return notes.get(population_type, "Standard Mendelian segregation expected.")

    def __repr__(self) -> str:
        n_markers = len(self.list_markers())
        n_chroms = len(self._db.get("chromosomes", {}))
        return (
            f"MarkerSelector(version='{self.version}', "
            f"markers={n_markers}, chromosomes={n_chroms})"
        )
