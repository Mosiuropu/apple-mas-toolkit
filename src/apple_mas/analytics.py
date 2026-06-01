"""
Analytics - Statistical analysis and visualization utilities for MAS data.

Provides functions for:
- Segregation distortion testing (chi-square) in F1, BC1, and F2 populations
- Genotype-phenotype boxplots for evaluating trait variance by marker class
- Polymorphism Information Content (PIC) calculation for multi-allelic and
  bi-allelic markers
"""

from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from apple_mas.allele_analysis import AlleleAnalyzer


# ------------------------------------------------------------------
# Segregation Distortion
# ------------------------------------------------------------------


def segregation_distortion_test(
    observed_genotypes: Dict[str, int],
    expected_ratio: Optional[Dict[str, float]] = None,
    population_type: str = "F1",
    significance: float = 0.05,
) -> Dict[str, Any]:
    """Perform a chi-square test for segregation distortion.

    Compares observed genotype counts against the expected Mendelian ratio
    for a given population type. Returns test statistics and a summary
    interpretation.

    Parameters
    ----------
    observed_genotypes : dict
        Observed genotype counts, e.g. ``{"AA": 48, "AB": 52, "BB": 0}``.
    expected_ratio : dict, optional
        Expected genotypic ratio as proportions, e.g.
        ``{"AA": 0.25, "AB": 0.50, "BB": 0.25}``. If None, uses the
        default ratio for *population_type*.
    population_type : str
        One of: 'F1', 'BC1', 'BC2', 'F2'.
    significance : float
        Significance level for the chi-square test.

    Returns
    -------
    dict
        Keys: 'chi2', 'df', 'p_value', 'significant', 'expected_counts',
        'observed_counts', 'population_type', 'interpretation'.
    """
    total = sum(observed_genotypes.values())
    if total == 0:
        return {
            "chi2": None, "df": None, "p_value": None, "significant": None,
            "expected_counts": {}, "observed_counts": observed_genotypes,
            "population_type": population_type,
            "interpretation": "No observed individuals; cannot perform test.",
        }

    default_ratios = {
        "F1": {"AA": 0.0, "AB": 1.0, "BB": 0.0},
        "BC1": {"AA": 0.5, "AB": 0.5, "BB": 0.0},
        "BC2": {"AA": 0.25, "AB": 0.50, "BB": 0.25},
        "F2": {"AA": 0.25, "AB": 0.50, "BB": 0.25},
    }

    if expected_ratio is None:
        expected_ratio = default_ratios.get(population_type, {"AA": 0.25, "AB": 0.50, "BB": 0.25})

    genotypes = sorted(set(list(observed_genotypes.keys()) + list(expected_ratio.keys())))
    expected_counts = {g: total * expected_ratio.get(g, 0) for g in genotypes}
    obs_list = [observed_genotypes.get(g, 0) for g in genotypes]
    exp_list = [expected_counts[g] for g in genotypes]

    valid = [(o, e) for o, e in zip(obs_list, exp_list) if e > 0]
    if not valid:
        return {
            "chi2": None, "df": None, "p_value": None, "significant": None,
            "expected_counts": expected_counts, "observed_counts": observed_genotypes,
            "population_type": population_type,
            "interpretation": "All expected counts are zero; test not applicable.",
        }

    obs_valid, exp_valid = zip(*valid)
    chi2 = sum((o - e) ** 2 / e for o, e in zip(obs_valid, exp_valid))
    n_genotypes_with_expected = len(valid)
    df = max(n_genotypes_with_expected - 1, 1)
    p_value = 1 - stats.chi2.cdf(chi2, df=df) if chi2 > 0 else 1.0
    significant = p_value < significance

    if significant:
        interpretation = (
            f"Significant segregation distortion detected (p={p_value:.4f}). "
            "Observed ratios deviate from expected Mendelian proportions. "
            "Possible causes: gametic selection, zygotic lethality, "
            "or marker-trait linkage."
        )
    else:
        interpretation = (
            f"No significant segregation distortion (p={p_value:.4f}). "
            "Observed genotype proportions are consistent with expected "
            f"Mendelian ratios for a {population_type} population."
        )

    return {
        "chi2": round(chi2, 4),
        "df": df,
        "p_value": round(p_value, 6),
        "significant": significant,
        "expected_counts": {g: round(v, 2) for g, v in expected_counts.items()},
        "observed_counts": observed_genotypes,
        "population_type": population_type,
        "interpretation": interpretation,
    }


def segregation_distortion_from_df(
    df: pd.DataFrame,
    marker: str,
    population_type: str = "F1",
    significance: float = 0.05,
) -> Dict[str, Any]:
    """Run segregation distortion test directly from a genotype DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data with columns ['marker', 'genotype', 'no_amplification'].
    marker : str
        Marker to test.
    population_type : str
        Population type for expected ratios.
    significance : float
        Significance level.

    Returns
    -------
    dict
        Output of :func:`segregation_distortion_test`.
    """
    marker_df = df[(df["marker"] == marker) & (~df["no_amplification"])]
    geno_counts = marker_df["genotype"].value_counts().to_dict()
    return segregation_distortion_test(
        geno_counts, population_type=population_type, significance=significance,
    )


# ------------------------------------------------------------------
# Genotype-Phenotype Boxplots
# ------------------------------------------------------------------


def genotype_phenotype_boxplot(
    genotypes: List[str],
    phenotypic_values: List[float],
    marker_name: str,
    trait_name: str = "",
    save_path: Optional[str] = None,
    figsize: Tuple[int, int] = (8, 6),
    dpi: int = 300,
    xlabel: str = "Marker Genotype",
    ylabel: str = "Phenotypic Value",
    title: Optional[str] = None,
    statistical_test: bool = True,
) -> plt.Figure:
    """Generate a box-and-whisker plot of phenotype grouped by marker genotype.

    Useful for evaluating whether phenotypic variance is associated with
    different marker allele classes (e.g., AA vs. AB vs. BB).

    Parameters
    ----------
    genotypes : list of str
        Genotype labels corresponding to each phenotypic observation.
    phenotypic_values : list of float
        Measured phenotypic values.
    marker_name : str
        Marker name (used in the plot title and file name).
    trait_name : str, optional
        Trait description (shown in axis labels).
    save_path : str, optional
        File path to save the figure.
    figsize : tuple
        Figure size in inches (width, height).
    dpi : int
        Resolution for saved figures.
    xlabel : str
        X-axis label.
    ylabel : str
        Y-axis label.
    title : str, optional
        Plot title. Auto-generated if None.
    statistical_test : bool
        If True, perform Kruskal-Wallis H-test and annotate p-value.

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure object.
    """
    plot_df = pd.DataFrame({
        "genotype": genotypes,
        "phenotype": phenotypic_values,
    })
    plot_df = plot_df.dropna(subset=["phenotype"])

    unique_genos = sorted(plot_df["genotype"].unique())
    palette = sns.color_palette("Set2", n_colors=len(unique_genos))

    fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(
        data=plot_df, x="genotype", y="phenotype", order=unique_genos,
        palette=palette, ax=ax, width=0.6, fliersize=4,
        boxprops=dict(alpha=0.8), linewidth=1.2,
    )
    sns.stripplot(
        data=plot_df, x="genotype", y="phenotype", order=unique_genos,
        color="black", alpha=0.4, size=3, jitter=True, ax=ax,
    )

    if statistical_test and len(unique_genos) >= 2:
        groups = [
            group["phenotype"].values
            for _, group in plot_df.groupby("genotype")
        ]
        groups = [g for g in groups if len(g) > 0]
        if len(groups) >= 2:
            h_stat, p_val = stats.kruskal(*groups)
            sig_text = f"Kruskal-Wallis H = {h_stat:.3f}, p = {p_val:.4f}"
            if p_val < 0.001:
                sig_text += " ***"
            elif p_val < 0.01:
                sig_text += " **"
            elif p_val < 0.05:
                sig_text += " *"
            else:
                sig_text += " (ns)"
            ax.text(
                0.02, 0.98, sig_text, transform=ax.transAxes,
                fontsize=10, verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8),
            )

    trait_label = f" ({trait_name})" if trait_name else ""
    ax.set_title(
        title or f"Phenotype by Marker Genotype — {marker_name}{trait_label}",
        fontsize=13, fontweight="bold",
    )
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    return fig


# ------------------------------------------------------------------
# PIC Calculator
# ------------------------------------------------------------------


def calculate_pic(allele_frequencies: Dict[str, float]) -> float:
    """Calculate Polymorphism Information Content (PIC) for a marker.

    For a multi-allelic marker with allele frequencies p_1, ..., p_n:

        PIC = 1 - sum(p_i^2) - sum_{i<j} (2 * p_i^2 * p_j^2)

    For bi-allelic markers this simplifies to:

        PIC = 1 - p^2 - q^2 - 2*p^2*q^2

    Parameters
    ----------
    allele_frequencies : dict
        Allele frequencies mapping allele names to proportions.
        Frequencies should sum to approximately 1.0.

    Returns
    -------
    float
        PIC value in the range [0, 1]. Returns 0.0 for monomorphic markers.

    Examples
    --------
    >>> calculate_pic({"A": 0.5, "B": 0.5})
    0.375
    >>> calculate_pic({"A": 1.0})
    0.0
    """
    freqs = list(allele_frequencies.values())
    if len(freqs) < 2:
        return 0.0

    total = sum(freqs)
    if total <= 0:
        return 0.0
    freqs = [f / total for f in freqs]

    sum_sq = sum(f ** 2 for f in freqs)
    n = len(freqs)
    pairwise_sum = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            pairwise_sum += 2 * (freqs[i] ** 2) * (freqs[j] ** 2)

    pic = 1.0 - sum_sq - pairwise_sum
    return round(max(0.0, pic), 6)


def calculate_pic_from_df(
    df: pd.DataFrame,
    marker: str,
) -> float:
    """Calculate PIC directly from a genotype DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data with columns ['marker', 'genotype', 'no_amplification'].
    marker : str
        Marker name.

    Returns
    -------
    float
        PIC value in [0, 1].
    """
    analyzer = AlleleAnalyzer()
    freq_df = analyzer.allele_frequencies(df, marker)
    allele_col = freq_df.columns[0]
    freq_dict = dict(zip(freq_df[allele_col].astype(str), freq_df["frequency"]))
    return calculate_pic(freq_dict)


# ------------------------------------------------------------------
# Batch analysis helpers
# ------------------------------------------------------------------


def batch_segregation_test(
    df: pd.DataFrame,
    population_type: str = "F1",
    significance: float = 0.05,
) -> pd.DataFrame:
    """Run segregation distortion tests for all markers in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data.
    population_type : str
        Population type for expected ratios.
    significance : float
        Significance level.

    Returns
    -------
    pd.DataFrame
        Summary with columns: ['marker', 'n_individuals', 'n_genotypes',
        'chi2', 'df', 'p_value', 'significant', 'interpretation'].
    """
    results = []
    for marker in sorted(df["marker"].unique()):
        result = segregation_distortion_from_df(
            df, marker, population_type=population_type, significance=significance,
        )
        results.append({
            "marker": marker,
            "n_individuals": sum(result["observed_counts"].values()),
            "n_genotypes": len(result["observed_counts"]),
            "chi2": result["chi2"],
            "df": result["df"],
            "p_value": result["p_value"],
            "significant": result["significant"],
            "interpretation": result["interpretation"],
        })
    return pd.DataFrame(results)


def batch_pic(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate PIC for all markers in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data.

    Returns
    -------
    pd.DataFrame
        Columns: ['marker', 'n_alleles', 'PIC'].
    """
    results = []
    for marker in sorted(df["marker"].unique()):
        pic = calculate_pic_from_df(df, marker)
        analyzer = AlleleAnalyzer()
        freq_df = analyzer.allele_frequencies(df, marker)
        results.append({
            "marker": marker,
            "n_alleles": len(freq_df),
            "PIC": pic,
        })
    return pd.DataFrame(results)
