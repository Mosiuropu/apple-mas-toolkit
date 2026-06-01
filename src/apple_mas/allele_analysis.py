"""
Allele Analysis - Statistical analysis of marker allele frequencies and genotype distributions.

Provides methods for calculating allele frequencies, genotype counts,
Hardy-Weinberg equilibrium tests, polymorphism information content (PIC),
and diversity metrics across multiple markers.
"""

from collections import Counter
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats


class AlleleAnalyzer:
    """Analyze allele frequencies and genotype distributions for apple markers.

    Takes standardized genotype DataFrames (from DataParser) and computes
    frequency statistics, diversity metrics, and significance tests.

    Examples
    --------
    >>> from apple_mas import DataParser, AlleleAnalyzer
    >>> parser = DataParser()
    >>> df, _ = parser.create_sample_dataset()
    >>> analyzer = AlleleAnalyzer()
    >>> freq = analyzer.allele_frequencies(df, marker="MYB10")
    >>> print(freq)
    """

    def __init__(self, exclude_no_amplification: bool = True):
        """Initialize the analyzer.

        Parameters
        ----------
        exclude_no_amplification : bool
            If True, exclude no-amplification entries from frequency calculations.
        """
        self._exclude_no_amp = exclude_no_amplification

    def allele_frequencies(
        self, df: pd.DataFrame, marker: str, normalize: bool = True
    ) -> pd.DataFrame:
        """Calculate allele frequencies for a single marker.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data with columns ['sample_id', 'marker', 'genotype', 'no_amplification'].
        marker : str
            Marker name to analyze.
        normalize : bool
            If True, frequencies sum to 1.0.

        Returns
        -------
        pd.DataFrame
            Columns: ['allele', 'count', 'frequency'].
        """
        marker_df = self._filter_marker(df, marker)
        alleles = self._expand_alleles(marker_df["genotype"])

        counter = Counter(alleles)
        total = sum(counter.values()) if normalize else len(alleles)

        result = pd.DataFrame([
            {"allele": allele, "count": count, "frequency": count / total if total > 0 else 0}
            for allele, count in sorted(counter.items())
        ])
        return result

    def genotype_frequencies(self, df: pd.DataFrame, marker: str) -> pd.DataFrame:
        """Calculate genotype frequencies for a single marker.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        marker : str
            Marker name.

        Returns
        -------
        pd.DataFrame
            Columns: ['genotype', 'count', 'frequency', 'percentage'].
        """
        marker_df = self._filter_marker(df, marker)
        geno_counts = marker_df["genotype"].value_counts().reset_index()
        geno_counts.columns = ["genotype", "count"]
        total = geno_counts["count"].sum()
        geno_counts["frequency"] = geno_counts["count"] / total
        geno_counts["percentage"] = geno_counts["frequency"] * 100
        return geno_counts.sort_values("count", ascending=False).reset_index(drop=True)

    def all_markers_summary(self, df: pd.DataFrame) -> pd.DataFrame:
        """Get a summary table for all markers in the dataset.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.

        Returns
        -------
        pd.DataFrame
            Summary with columns: ['marker', 'trait', 'n_samples', 'n_genotypes',
            'n_unique_alleles', 'no_amp_count', 'most_common_geno', 'most_common_freq'].
        """
        summaries = []
        for marker in df["marker"].unique():
            marker_df = self._filter_marker(df, marker, include_no_amp=True)
            marker_df_filtered = self._filter_marker(df, marker)

            genotypes = marker_df_filtered["genotype"]
            alleles = self._expand_alleles(genotypes)

            most_common = genotypes.value_counts().iloc[0] if len(genotypes) > 0 else 0
            most_common_geno = genotypes.value_counts().index[0] if len(genotypes) > 0 else "N/A"

            no_amp = marker_df["no_amplification"].sum()

            summaries.append({
                "marker": marker,
                "trait": self._get_trait(marker),
                "n_samples": len(marker_df),
                "n_genotypes": len(genotypes.unique()),
                "n_unique_alleles": len(set(alleles)),
                "no_amp_count": int(no_amp),
                "most_common_geno": most_common_geno,
                "most_common_freq": most_common / len(genotypes) if len(genotypes) > 0 else 0,
            })

        return pd.DataFrame(summaries)

    def polymorphism_information_content(
        self, df: pd.DataFrame, marker: str
    ) -> float:
        """Calculate Polymorphism Information Content (PIC) for a marker.

        PIC measures the polymorphism level of a marker based on allele
        frequencies. Higher values indicate greater informativeness.

        PIC = 1 - sum(pi^2) - sum_{i<j} (2 * pi^2 * pj^2)

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        marker : str
            Marker name.

        Returns
        -------
        float
            PIC value between 0 and 1.
        """
        freq_df = self.allele_frequencies(df, marker)
        freqs = freq_df["frequency"].values

        if len(freqs) <= 1:
            return 0.0

        pic = 1.0 - np.sum(freqs ** 2)

        # Add the pairwise term
        n = len(freqs)
        pairwise_sum = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                pairwise_sum += (2 * freqs[i] ** 2 * freqs[j] ** 2)
        pic -= pairwise_sum

        return float(max(0.0, pic))

    def heterozygosity(
        self, df: pd.DataFrame, marker: str
    ) -> Dict[str, float]:
        """Calculate observed and expected heterozygosity for a marker.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        marker : str
            Marker name.

        Returns
        -------
        dict
            {'observed_het': float, 'expected_het': float, 'het_deficit': float}.
        """
        marker_df = self._filter_marker(df, marker)
        genotypes = marker_df["genotype"]

        # Observed heterozygosity: proportion of heterozygous genotypes
        heterozygous = sum(
            1 for g in genotypes if ":" in str(g) and not self._is_homozygous(str(g))
        )
        obs_het = heterozygous / len(genotypes) if len(genotypes) > 0 else 0.0

        # Expected heterozygosity: 1 - sum(pi^2) based on allele frequencies
        freq_df = self.allele_frequencies(df, marker)
        freqs = freq_df["frequency"].values
        exp_het = 1.0 - np.sum(freqs ** 2) if len(freqs) > 1 else 0.0

        return {
            "observed_het": obs_het,
            "expected_het": exp_het,
            "het_deficit": exp_het - obs_het,
        }

    def hardy_weinberg_test(
        self, df: pd.DataFrame, marker: str, significance: float = 0.05
    ) -> Dict[str, Any]:
        """Perform Hardy-Weinberg equilibrium test for biallelic markers.

        Uses chi-squared goodness-of-fit test. Only applicable to
        markers with exactly 2 alleles.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        marker : str
            Marker name.
        significance : float
            Significance level for the test.

        Returns
        -------
        dict
            {'chi2': float, 'p_value': float, 'in_hwe': bool, 'n_genotypes': int,
             'expected_counts': dict, 'observed_counts': dict}.
        """
        freq_df = self.allele_frequencies(df, marker)
        marker_df = self._filter_marker(df, marker)
        genotypes = marker_df["genotype"].value_counts().to_dict()

        if len(freq_df) != 2:
            return {
                "chi2": None,
                "p_value": None,
                "in_hwe": None,
                "n_genotypes": len(genotypes),
                "expected_counts": {},
                "observed_counts": genotypes,
                "note": "HWE test only applicable to biallelic markers",
            }

        p = freq_df.iloc[0]["frequency"]
        q = freq_df.iloc[1]["frequency"]
        n = marker_df.shape[0]

        # Expected genotype counts under HWE
        allele1 = freq_df.iloc[0]["allele"]
        allele2 = freq_df.iloc[1]["allele"]
        expected = {
            f"{allele1}:{allele1}": n * p ** 2,
            f"{allele1}:{allele2}": 2 * n * p * q,
            f"{allele2}:{allele2}": n * q ** 2,
        }

        # Map observed to match expected categories
        obs_mapped = {}
        for geno in [f"{allele1}:{allele1}", f"{allele1}:{allele2}", f"{allele2}:{allele2}"]:
            # Also check reverse order
            alt = f"{allele2}:{allele1}"
            obs_mapped[geno] = genotypes.get(geno, 0) + genotypes.get(alt, 0)

        # Chi-squared test (only for genotypes with expected > 0)
        chi2 = 0.0
        for geno in expected:
            if expected[geno] > 0:
                obs = obs_mapped.get(geno, 0)
                chi2 += (obs - expected[geno]) ** 2 / expected[geno]

        # df = number of genotypes - number of alleles = 3 - 2 = 1
        p_value = 1 - stats.chi2.cdf(chi2, df=1) if chi2 > 0 else 1.0

        return {
            "chi2": chi2,
            "p_value": p_value,
            "in_hwe": p_value > significance,
            "n_genotypes": len(genotypes),
            "expected_counts": expected,
            "observed_counts": obs_mapped,
        }

    def cross_marker_summary(self, df: pd.DataFrame) -> pd.DataFrame:
        """Get a comprehensive cross-marker analysis summary.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data for all markers.

        Returns
        -------
        pd.DataFrame
            Summary with PIC, heterozygosity, and allele counts per marker.
        """
        rows = []
        for marker in df["marker"].unique():
            pic = self.polymorphism_information_content(df, marker)
            het = self.heterozygosity(df, marker)
            freq_df = self.allele_frequencies(df, marker)

            rows.append({
                "marker": marker,
                "trait": self._get_trait(marker),
                "n_alleles": len(freq_df),
                "PIC": round(pic, 4),
                "obs_het": round(het["observed_het"], 4),
                "exp_het": round(het["expected_het"], 4),
                "het_deficit": round(het["het_deficit"], 4),
            })

        return pd.DataFrame(rows)

    def _filter_marker(
        self, df: pd.DataFrame, marker: str, include_no_amp: bool = False
    ) -> pd.DataFrame:
        """Filter DataFrame for a specific marker."""
        marker_df = df[df["marker"] == marker].copy()
        if not include_no_amp and self._exclude_no_amp:
            marker_df = marker_df[~marker_df["no_amplification"]]
        return marker_df

    @staticmethod
    def _expand_alleles(genotypes: pd.Series) -> List[str]:
        """Expand diploid genotypes into individual alleles."""
        alleles = []
        for g in genotypes:
            g_str = str(g)
            if ":" in g_str:
                parts = g_str.split(":")
                alleles.extend([p.strip() for p in parts if p.strip()])
            elif g_str and g_str not in ("", "NA", "N/A", "*"):
                alleles.append(g_str.strip())
        return alleles

    @staticmethod
    def _is_homozygous(genotype: str) -> bool:
        """Check if a genotype string represents a homozygous call."""
        if ":" not in genotype:
            return True
        parts = genotype.split(":")
        return len(parts) == 2 and parts[0].strip() == parts[1].strip()

    @staticmethod
    def _get_trait(marker: str) -> str:
        """Get trait name for a marker."""
        trait_map = {
            "MYB10": "Flesh Color",
            "RED_TE": "Skin Color",
            "MA_INDEL": "Fruit Acidity",
            "BP16": "Bitter Pit",
            "BP13": "Bitter Pit",
            "ACS": "Fruit Texture",
            "ACO": "Fruit Texture",
            "MD_PG1": "Fruit Texture",
        }
        return trait_map.get(marker, "Unknown")
