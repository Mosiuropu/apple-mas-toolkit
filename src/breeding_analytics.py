"""
Breeding Analytics - Statistical analysis and visualization suite for apple
molecular marker data.

Provides functions for:
- Segregation distortion testing (chi-square) in F1, BC1, and F2 populations
- Genotype-phenotype mapping (box-and-whisker plots with Kruskal-Wallis tests)
- Polymorphism Information Content (PIC) calculation for multi-allelic markers
- Allele frequency analysis and diversity metrics

All functions operate on pandas DataFrames with a standardised schema:
columns ``sample_id``, ``marker``, ``genotype``, and ``trait``.

References:
    - Botstein et al., 1980, Am J Hum Genet (PIC definition)
    - Liu et al., 1996, Genetics (PIC for multi-allelic markers)
"""

from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


# ======================================================================
# Segregation Distortion Module
# ======================================================================


def segregation_distortion_test(
    observed_genotypes: Dict[str, int],
    population_type: str = "F1",
    expected_ratio: Optional[Dict[str, float]] = None,
    significance: float = 0.05,
) -> Dict[str, Any]:
    """Perform a chi-square test for segregation distortion.

    Compares observed genotype counts against the expected Mendelian ratio
    for a given population type.

    Parameters
    ----------
    observed_genotypes : dict
        Observed counts, e.g. ``{"AA": 48, "AB": 52, "BB": 0}``.
    population_type : str
        One of ``F1``, ``BC1``, ``BC2``, ``F2``.
    expected_ratio : dict, optional
        Custom expected proportions.  If ``None``, defaults for
        *population_type* are used.
    significance : float
        Alpha level for the chi-square test.

    Returns
    -------
    dict
        Keys: ``chi2``, ``df``, ``p_value``, ``significant``,
        ``expected_counts``, ``observed_counts``, ``population_type``,
        ``interpretation``.
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
        expected_ratio = default_ratios.get(
            population_type, {"AA": 0.25, "AB": 0.50, "BB": 0.25}
        )

    genotypes = sorted(
        set(list(observed_genotypes.keys()) + list(expected_ratio.keys()))
    )
    expected_counts = {g: total * expected_ratio.get(g, 0) for g in genotypes}

    valid = [
        (observed_genotypes.get(g, 0), expected_counts[g])
        for g in genotypes
        if expected_counts[g] > 0
    ]
    if not valid:
        return {
            "chi2": None, "df": None, "p_value": None, "significant": None,
            "expected_counts": expected_counts,
            "observed_counts": observed_genotypes,
            "population_type": population_type,
            "interpretation": "All expected counts are zero; test not applicable.",
        }

    obs_vals, exp_vals = zip(*valid)
    chi2 = sum((o - e) ** 2 / e for o, e in zip(obs_vals, exp_vals))
    df = max(len(valid) - 1, 1)
    p_value = 1 - stats.chi2.cdf(chi2, df=df) if chi2 > 0 else 1.0
    significant = p_value < significance

    if significant:
        interp = (
            f"Significant segregation distortion detected (p={p_value:.4f}). "
            "Observed ratios deviate from expected Mendelian proportions. "
            "Possible causes: gametic selection, zygotic lethality, "
            "or marker-trait linkage."
        )
    else:
        interp = (
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
        "interpretation": interp,
    }


def segregation_distortion_from_df(
    df: pd.DataFrame,
    marker: str,
    population_type: str = "F1",
    significance: float = 0.05,
) -> Dict[str, Any]:
    """Run segregation distortion test from a genotype DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Columns include ``marker``, ``genotype``, ``no_amplification``.
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
    mask = df["marker"] == marker
    if "no_amplification" in df.columns:
        mask = mask & ~df["no_amplification"]
    marker_df = df[mask]
    geno_counts = marker_df["genotype"].value_counts().to_dict()
    return segregation_distortion_test(
        geno_counts, population_type=population_type, significance=significance,
    )


def batch_segregation_test(
    df: pd.DataFrame,
    population_type: str = "F1",
    significance: float = 0.05,
) -> pd.DataFrame:
    """Run segregation distortion tests for all markers.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data.
    population_type : str
        Population type.
    significance : float
        Significance level.

    Returns
    -------
    pd.DataFrame
        Summary with columns: ``marker``, ``n_individuals``,
        ``n_genotypes``, ``chi2``, ``df``, ``p_value``,
        ``significant``, ``interpretation``.
    """
    results = []
    for marker in sorted(df["marker"].unique()):
        result = segregation_distortion_from_df(
            df, marker, population_type=population_type,
            significance=significance,
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


# ======================================================================
# Genotype–Phenotype Mapping
# ======================================================================


def genotype_phenotype_boxplot(
    genotypes: List[str],
    phenotypic_values: List[float],
    marker_name: str,
    trait_name: str = "",
    save_path: Optional[str] = None,
    figsize: Tuple[int, int] = (8, 6),
    dpi: int = 300,
    title: Optional[str] = None,
    statistical_test: bool = True,
) -> plt.Figure:
    """Generate a box-and-whisker plot of phenotype grouped by genotype.

    Useful for evaluating whether phenotypic variance is associated with
    different marker allele classes (e.g. Homozygous Resistant vs.
    Heterozygous).

    Parameters
    ----------
    genotypes : list of str
        Genotype labels for each observation.
    phenotypic_values : list of float
        Measured phenotypic values.
    marker_name : str
        Marker name (used in title and file name).
    trait_name : str, optional
        Trait description for axis labels.
    save_path : str, optional
        File path to save the figure.
    figsize : tuple
        Figure size (width, height) in inches.
    dpi : int
        Resolution for saved figures.
    title : str, optional
        Plot title.  Auto-generated if ``None``.
    statistical_test : bool
        If ``True``, perform Kruskal-Wallis H-test and annotate p-value.

    Returns
    -------
    matplotlib.figure.Figure
    """
    plot_df = pd.DataFrame({
        "genotype": genotypes,
        "phenotype": phenotypic_values,
    }).dropna(subset=["phenotype"])

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
            grp["phenotype"].values
            for _, grp in plot_df.groupby("genotype")
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
    ax.set_xlabel("Marker Genotype", fontsize=11)
    ax.set_ylabel("Phenotypic Value", fontsize=11)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    return fig


# ======================================================================
# Polymorphism Information Content (PIC)
# ======================================================================


def calculate_pic(allele_frequencies: Dict[str, float]) -> float:
    """Calculate PIC for a marker from allele frequencies.

    For a multi-allelic marker with frequencies p_1, ..., p_n::

        PIC = 1 - sum(p_i^2) - sum_{i<j} (2 * p_i^2 * p_j^2)

    Parameters
    ----------
    allele_frequencies : dict
        Allele name → proportion (should sum to ~1.0).

    Returns
    -------
    float
        PIC value in [0, 1].  Returns 0.0 for monomorphic markers.

    References
    ----------
    Botstein et al., 1980, Am J Hum Genet.

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
    pairwise_sum = sum(
        2 * (freqs[i] ** 2) * (freqs[j] ** 2)
        for i in range(n)
        for j in range(i + 1, n)
    )

    pic = 1.0 - sum_sq - pairwise_sum
    return round(max(0.0, pic), 6)


def calculate_pic_from_genotypes(genotypes: List[str]) -> float:
    """Calculate PIC from a list of diploid genotype strings.

    Parameters
    ----------
    genotypes : list of str
        Genotype calls, e.g. ``["390:390", "390:490", "490:490"]``.

    Returns
    -------
    float
        PIC value in [0, 1].
    """
    allele_counts: Dict[str, int] = {}
    total_alleles = 0
    for g in genotypes:
        if not isinstance(g, str) or g.strip() in ("", "NA", "N/A", "*"):
            continue
        parts = [a.strip() for a in g.split(":") if a.strip()]
        for a in parts:
            allele_counts[a] = allele_counts.get(a, 0) + 1
            total_alleles += 1

    if total_alleles == 0 or len(allele_counts) < 2:
        return 0.0

    freqs = {a: c / total_alleles for a, c in allele_counts.items()}
    return calculate_pic(freqs)


def batch_pic(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate PIC for all markers in a genotype DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data with ``marker`` and ``genotype`` columns.

    Returns
    -------
    pd.DataFrame
        Columns: ``marker``, ``n_alleles``, ``PIC``.
    """
    results = []
    for marker in sorted(df["marker"].unique()):
        marker_df = df[df["marker"] == marker]
        genotypes = marker_df["genotype"].dropna().astype(str).tolist()
        pic = calculate_pic_from_genotypes(genotypes)

        # Count unique alleles
        alleles = set()
        for g in genotypes:
            if g.strip() not in ("", "NA", "N/A", "*"):
                for a in g.split(":"):
                    a = a.strip()
                    if a:
                        alleles.add(a)

        results.append({
            "marker": marker,
            "n_alleles": len(alleles),
            "PIC": pic,
        })
    return pd.DataFrame(results)


# ======================================================================
# Allele Frequency Analysis
# ======================================================================


def allele_frequencies(
    df: pd.DataFrame,
    marker: str,
    normalize: bool = True,
) -> pd.DataFrame:
    """Calculate allele frequencies for a single marker.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data with ``marker`` and ``genotype`` columns.
    marker : str
        Marker name.
    normalize : bool
        If ``True``, frequencies sum to 1.0.

    Returns
    -------
    pd.DataFrame
        Columns: ``allele``, ``count``, ``frequency``.
    """
    from collections import Counter

    marker_df = df[df["marker"] == marker]
    alleles: List[str] = []
    for g in marker_df["genotype"].dropna():
        g_str = str(g)
        if ":" in g_str:
            alleles.extend(a.strip() for a in g_str.split(":") if a.strip())
        elif g_str.strip() and g_str.strip() not in ("NA", "N/A", "*"):
            alleles.append(g_str.strip())

    counter = Counter(alleles)
    total = sum(counter.values()) if normalize else len(alleles)

    return pd.DataFrame([
        {"allele": a, "count": c, "frequency": c / total if total > 0 else 0}
        for a, c in sorted(counter.items())
    ])


def genotype_frequencies(df: pd.DataFrame, marker: str) -> pd.DataFrame:
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
        Columns: ``genotype``, ``count``, ``frequency``, ``percentage``.
    """
    marker_df = df[df["marker"] == marker]
    geno_counts = marker_df["genotype"].value_counts().reset_index()
    geno_counts.columns = ["genotype", "count"]
    total = geno_counts["count"].sum()
    geno_counts["frequency"] = geno_counts["count"] / total
    geno_counts["percentage"] = geno_counts["frequency"] * 100
    return geno_counts.sort_values("count", ascending=False).reset_index(drop=True)


# ======================================================================
# Diversity Metrics
# ======================================================================


def heterozygosity(
    df: pd.DataFrame,
    marker: str,
) -> Dict[str, float]:
    """Calculate observed and expected heterozygosity.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data.
    marker : str
        Marker name.

    Returns
    -------
    dict
        ``observed_het``, ``expected_het``, ``het_deficit``.
    """
    marker_df = df[df["marker"] == marker]
    genotypes = marker_df["genotype"]

    heterozygous = 0
    for g in genotypes:
        g_str = str(g)
        if ":" in g_str:
            parts = [a.strip() for a in g_str.split(":")]
            if len(parts) == 2 and parts[0] != parts[1]:
                heterozygous += 1
    obs_het = heterozygous / len(genotypes) if len(genotypes) > 0 else 0.0

    freq_df = allele_frequencies(df, marker)
    freqs = freq_df["frequency"].values
    exp_het = 1.0 - np.sum(freqs ** 2) if len(freqs) > 1 else 0.0

    return {
        "observed_het": obs_het,
        "expected_het": float(exp_het),
        "het_deficit": float(exp_het) - obs_het,
    }


def cross_marker_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Comprehensive diversity summary for all markers.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data for all markers.

    Returns
    -------
    pd.DataFrame
        Columns: ``marker``, ``n_alleles``, ``PIC``,
        ``obs_het``, ``exp_het``, ``het_deficit``.
    """
    rows = []
    for marker in sorted(df["marker"].unique()):
        pic = calculate_pic_from_genotypes(
            df[df["marker"] == marker]["genotype"].dropna().astype(str).tolist()
        )
        het = heterozygosity(df, marker)
        freq_df = allele_frequencies(df, marker)
        rows.append({
            "marker": marker,
            "n_alleles": len(freq_df),
            "PIC": round(pic, 4),
            "obs_het": round(het["observed_het"], 4),
            "exp_het": round(het["expected_het"], 4),
            "het_deficit": round(het["het_deficit"], 4),
        })
    return pd.DataFrame(rows)


# ======================================================================
# Summary Report
# ======================================================================


def generate_report(
    df: pd.DataFrame,
    population_type: str = "F1",
    output_dir: str = "report",
) -> str:
    """Generate a comprehensive Markdown analysis report.

    Parameters
    ----------
    df : pd.DataFrame
        Genotype data.
    population_type : str
        Population type for segregation tests.
    output_dir : str
        Directory to write the report file.

    Returns
    -------
    str
        Path to the generated report.
    """
    from pathlib import Path

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    n_samples = df["sample_id"].nunique() if "sample_id" in df.columns else len(df)
    n_markers = df["marker"].nunique()

    pic_df = batch_pic(df)
    seg_df = batch_segregation_test(df, population_type=population_type)
    cross_df = cross_marker_summary(df)

    lines = [
        "# Apple MAS Analysis Report",
        "",
        f"**Samples analysed:** {n_samples}",
        f"**Markers analysed:** {n_markers}",
        f"**Population type:** {population_type}",
        "",
        "## Polymorphism Information Content (PIC)",
        "",
        pic_df.to_string(index=False),
        "",
        "## Segregation Distortion Tests",
        "",
        seg_df[["marker", "chi2", "p_value", "significant"]].to_string(index=False),
        "",
        "## Cross-Marker Diversity",
        "",
        cross_df.to_string(index=False),
        "",
    ]

    report_path = Path(output_dir) / "analysis_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return str(report_path.resolve())
