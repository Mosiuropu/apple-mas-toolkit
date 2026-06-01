"""
Data Parser - Import and clean genotype data from various sources.

Handles parsing of GeneMapper output, CSV/Excel spreadsheets, and
manual data entry formats. Provides standardized DataFrames ready
for downstream analysis.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd


class DataParser:
    """Parse and standardize apple genotype data from multiple sources.

    Supports GeneMapper export files, CSV/Excel spreadsheets, and
    structured dictionaries. All outputs are standardized DataFrames
    with consistent column naming.

    Examples
    --------
    >>> parser = DataParser()
    >>> df = parser.from_csv("my_genotype_data.csv", sample_col="accession", marker_col="marker", allele_col="allele")
    >>> df.head()
    """

    MARKER_NAMES = [
        "MYB10", "RED_TE", "MA_INDEL", "BP16", "BP13",
        "ACS", "ACO", "MD_PG1",
    ]

    TRAIT_MAP = {
        "MYB10": "Flesh Color",
        "RED_TE": "Skin Color",
        "MA_INDEL": "Fruit Acidity",
        "BP16": "Bitter Pit Susceptibility",
        "BP13": "Bitter Pit Susceptibility",
        "ACS": "Fruit Texture",
        "ACO": "Fruit Texture",
        "MD_PG1": "Fruit Texture",
    }

    def __init__(self):
        """Initialize the DataParser."""
        self._validation_errors: List[str] = []

    def from_csv(
        self,
        filepath: str,
        sample_col: str = "sample_id",
        marker_col: str = "marker",
        allele_col: str = "genotype",
        delimiter: str = ",",
    ) -> pd.DataFrame:
        """Parse genotype data from a CSV file.

        The CSV should have columns for sample identifier, marker name,
        and genotype/allele call.

        Parameters
        ----------
        filepath : str
            Path to CSV file.
        sample_col : str
            Column name for sample/accession identifier.
        marker_col : str
            Column name for marker name.
        allele_col : str
            Column name for genotype or allele call.
        delimiter : str
            CSV delimiter character.

        Returns
        -------
        pd.DataFrame
            Standardized genotype data with columns:
            ['sample_id', 'marker', 'genotype', 'trait', 'no_amplification'].
        """
        df = pd.read_csv(filepath, delimiter=delimiter)
        return self._standardize_dataframe(df, sample_col, marker_col, allele_col)

    def from_excel(
        self,
        filepath: str,
        sheet_name: str = 0,
        sample_col: str = "sample_id",
        marker_col: str = "marker",
        allele_col: str = "genotype",
    ) -> pd.DataFrame:
        """Parse genotype data from an Excel file.

        Parameters
        ----------
        filepath : str
            Path to Excel file (.xlsx or .xls).
        sheet_name : str or int
            Sheet name or index.
        sample_col : str
            Column name for sample identifier.
        marker_col : str
            Column name for marker name.
        allele_col : str
            Column name for genotype call.

        Returns
        -------
        pd.DataFrame
            Standardized genotype data.
        """
        df = pd.read_excel(filepath, sheet_name=sheet_name)
        return self._standardize_dataframe(df, sample_col, marker_col, allele_col)

    def from_wide_csv(
        self,
        filepath: str,
        sample_col: str = "sample_id",
        delimiter: str = ",",
    ) -> pd.DataFrame:
        """Parse wide-format CSV where each column is a marker.

        This format has one row per sample and a column for each marker
        containing the genotype call.

        Parameters
        ----------
        filepath : str
            Path to CSV file.
        sample_col : str
            Column name for sample identifier.
        delimiter : str
            CSV delimiter.

        Returns
        -------
        pd.DataFrame
            Standardized genotype data in long format.
        """
        df = pd.read_csv(filepath, delimiter=delimiter)
        if sample_col not in df.columns:
            raise ValueError(
                f"Sample column '{sample_col}' not found. "
                f"Available: {list(df.columns)}"
            )

        marker_cols = [c for c in df.columns if c != sample_col]
        melted = df.melt(
            id_vars=[sample_col],
            value_vars=marker_cols,
            var_name="marker",
            value_name="genotype",
        )
        melted.rename(columns={sample_col: "sample_id"}, inplace=True)

        melted["trait"] = melted["marker"].map(self.TRAIT_MAP)
        melted["no_amplification"] = melted["genotype"].apply(self._is_no_amp)
        return melted

    def from_dict(
        self,
        data: Dict[str, Dict[str, str]],
    ) -> pd.DataFrame:
        """Parse genotype data from a nested dictionary.

        Expected format: {sample_id: {marker: genotype, ...}, ...}

        Parameters
        ----------
        data : dict
            Nested dictionary of genotype data.

        Returns
        -------
        pd.DataFrame
            Standardized genotype data in long format.
        """
        rows = []
        for sample_id, markers in data.items():
            for marker, genotype in markers.items():
                rows.append({
                    "sample_id": sample_id,
                    "marker": marker,
                    "genotype": str(genotype),
                })
        df = pd.DataFrame(rows)
        df["trait"] = df["marker"].map(self.TRAIT_MAP)
        df["no_amplification"] = df["genotype"].apply(self._is_no_amp)
        return df

    def from_genemapper(
        self,
        filepath: str,
        sample_col: str = "Sample Name",
        marker_col: str = "Marker",
        allele_col: str = "Allele",
        size_col: str = "Size",
    ) -> pd.DataFrame:
        """Parse GeneMapper export files.

        GeneMapper outputs tab-separated files with fragment analysis
        results. This parser handles the standard export format.

        Parameters
        ----------
        filepath : str
            Path to GeneMapper export file (tab-separated).
        sample_col : str
            Column for sample name.
        marker_col : str
            Column for marker name.
        allele_col : str
            Column for allele name.
        size_col : str
            Column for fragment size in bp.

        Returns
        -------
        pd.DataFrame
            Parsed fragment sizes per sample per marker.
        """
        df = pd.read_csv(filepath, sep="\t")

        # Try to identify the right columns
        available = list(df.columns)
        sample_col = self._find_column(available, sample_col, ["Sample", "Sample Name", "SampleID"])
        marker_col = self._find_column(available, marker_col, ["Marker", "Locus", "Marker Name"])
        size_col = self._find_column(available, size_col, ["Size", "Fragment Size", "bp", "Allele Size"])

        result = df[[sample_col, marker_col, size_col]].copy()
        result.columns = ["sample_id", "marker", "fragment_size_bp"]
        result["genotype"] = result["fragment_size_bp"].round(0).astype(int).astype(str)
        result["trait"] = result["marker"].map(self.TRAIT_MAP)
        result["no_amplification"] = False
        return result

    def create_sample_dataset(self) -> Tuple[pd.DataFrame, Dict[str, Dict[str, str]]]:
        """Create a sample dataset for demonstration and testing.

        Returns a realistic dataset based on ANABP germplasm patterns,
        with genotype distributions similar to those reported in the
        Australian National Apple Breeding Program.

        Returns
        -------
        tuple of (pd.DataFrame, dict)
            (genotype_df, raw_dict) - standardized DataFrame and raw dictionary.
        """
        import random
        random.seed(42)

        n_samples = 50
        sample_ids = [f"ANABP_{i:03d}" for i in range(1, n_samples + 1)]

        raw_data = {}
        for sid in sample_ids:
            raw_data[sid] = {
                "MYB10": "390:390" if random.random() < 0.95 else "390:490",
                "RED_TE": "750" if random.random() < 0.64 else "NB",
                "MA_INDEL": random.choices(
                    ["409:456", "409:409", "456:456", "409:452"],
                    weights=[46, 8, 10, 5],
                    k=1,
                )[0],
                "BP16": random.choices(
                    ["203:368", "368:368", "203:203", "203:206", "206:368"],
                    weights=[34, 31, 12, 4, 7],
                    k=1,
                )[0],
                "BP13": random.choices(
                    ["232:250", "224:232", "232:252", "250:252", "232:232"],
                    weights=[30, 10, 6, 5, 8],
                    k=1,
                )[0],
                "ACS": random.choices(
                    ["200:341", "200:200", "341:341"],
                    weights=[58, 20, 20],
                    k=1,
                )[0],
                "ACO": random.choices(
                    ["300:300", "237:300", "237:237"],
                    weights=[52, 30, 16],
                    k=1,
                )[0],
                "MD_PG1": random.choices(
                    ["289:292", "289:298", "292:292", "292:298", "298:298"],
                    weights=[21, 19, 14, 15, 16],
                    k=1,
                )[0],
            }

        df = self.from_dict(raw_data)
        return df, raw_data

    def _standardize_dataframe(
        self,
        df: pd.DataFrame,
        sample_col: str,
        marker_col: str,
        allele_col: str,
    ) -> pd.DataFrame:
        """Standardize a raw DataFrame to the canonical format."""
        available = list(df.columns)
        sample_col = self._find_column(available, sample_col, ["sample_id", "Sample", "accession", "ID"])
        marker_col = self._find_column(available, marker_col, ["marker", "Marker", "Locus"])
        allele_col = self._find_column(available, allele_col, ["genotype", "Genotype", "allele", "Allele"])

        result = df[[sample_col, marker_col, allele_col]].copy()
        result.columns = ["sample_id", "marker", "genotype"]
        result["trait"] = result["marker"].map(self.TRAIT_MAP)
        result["no_amplification"] = result["genotype"].apply(self._is_no_amp)
        return result

    @staticmethod
    def _find_column(
        available: List[str],
        preferred: str,
        alternatives: List[str],
    ) -> str:
        """Find a column by preferred name or alternatives."""
        if preferred in available:
            return preferred
        for alt in alternatives:
            if alt in available:
                return alt
        # Fuzzy match
        for col in available:
            for alt in [preferred] + alternatives:
                if alt.lower() in col.lower():
                    return col
        raise ValueError(
            f"Column '{preferred}' not found. Available: {available}"
        )

    @staticmethod
    def _is_no_amp(genotype: str) -> bool:
        """Check if a genotype string represents no amplification."""
        if not isinstance(genotype, str):
            return True
        genotype = genotype.strip()
        if genotype in ("", "*", "NA", "N/A", "no amp", "no_amplification"):
            return True
        if ":" in genotype and genotype.replace(":", "").replace(".", "").strip() == "":
            return True
        return False

    def validate(
        self, df: pd.DataFrame, known_markers: Optional[List[str]] = None
    ) -> Tuple[bool, List[str]]:
        """Validate a genotype DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data to validate.
        known_markers : list of str, optional
            Expected marker names. Defaults to built-in marker list.

        Returns
        -------
        tuple of (bool, list of str)
            (is_valid, list_of_issues).
        """
        issues = []

        required_cols = {"sample_id", "marker", "genotype"}
        missing = required_cols - set(df.columns)
        if missing:
            issues.append(f"Missing required columns: {missing}")
            return False, issues

        if known_markers is None:
            known_markers = self.MARKER_NAMES

        unknown_markers = set(df["marker"].unique()) - set(known_markers)
        if unknown_markers:
            issues.append(f"Unknown markers found: {unknown_markers}")

        empty_genotypes = df["genotype"].isna().sum()
        if empty_genotypes > 0:
            issues.append(f"{empty_genotypes} rows have missing genotype calls")

        no_amp = df.get("no_amplification", pd.Series(dtype=bool))
        if no_amp.any():
            n = no_amp.sum()
            issues.append(f"{n} rows represent no-amplification results")

        return len(issues) == 0, issues
