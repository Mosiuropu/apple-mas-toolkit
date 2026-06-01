"""
Visualization - Publication-quality plots for apple marker analysis.

Generates allele frequency bar plots, genotype distribution pie charts,
correlation heatmaps, dendrograms, and multi-marker comparison figures.
All plots use a consistent color scheme suitable for publication.
"""

from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# Consistent color palette for marker categories
TRAIT_COLORS = {
    "Flesh Color": "#C0392B",
    "Skin Color": "#E74C3C",
    "Fruit Acidity": "#2ECC71",
    "Bitter Pit": "#F39C12",
    "Bitter Pit Susceptibility": "#F39C12",
    "Fruit Texture": "#3498DB",
    "Fruit Texture (Firmness)": "#3498DB",
}

ALLELE_PALETTE = sns.color_palette("Set2", 12)


class MarkerVisualizer:
    """Generate publication-quality visualizations for marker analysis.

    Provides methods for creating allele frequency plots, genotype
    distribution charts, correlation heatmaps, and multi-panel figures
    suitable for journal publication.

    Examples
    --------
    >>> from apple_mas import DataParser, MarkerVisualizer
    >>> parser = DataParser()
    >>> df, _ = parser.create_sample_dataset()
    >>> viz = MarkerVisualizer()
    >>> fig = viz.allele_frequency_plot(df, marker="MYB10", save_path="myb10_freq.png")
    """

    def __init__(self, style: str = "seaborn-v0_8-whitegrid", dpi: int = 300):
        """Initialize the visualizer.

        Parameters
        ----------
        style : str
            Matplotlib/seaborn style name.
        dpi : int
            Resolution for saved figures.
        """
        self.dpi = dpi
        try:
            plt.style.use(style)
        except OSError:
            plt.style.use("seaborn-v0_8" if "seaborn-v0_8" in plt.style.available else "default")
        sns.set_palette("Set2")

    def allele_frequency_plot(
        self,
        df: pd.DataFrame,
        marker: str,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (8, 5),
        title: Optional[str] = None,
        color: Optional[str] = None,
    ) -> plt.Figure:
        """Create a bar plot of allele frequencies for a marker.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data with columns ['marker', 'genotype', 'no_amplification'].
        marker : str
            Marker name to plot.
        save_path : str, optional
            Path to save the figure. If None, figure is not saved.
        figsize : tuple
            Figure size (width, height) in inches.
        title : str, optional
            Plot title. If None, auto-generates from marker name.
        color : str, optional
            Bar color. If None, uses trait-based color.

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure object.
        """
        from apple_mas.allele_analysis import AlleleAnalyzer

        analyzer = AlleleAnalyzer()
        freq_df = analyzer.allele_frequencies(df, marker)
        allele_col = freq_df.columns[0] if len(freq_df.columns) > 0 else "allele"
        count_col = "count"
        freq_col = "frequency"

        fig, ax1 = plt.subplots(figsize=figsize)

        if color is None:
            trait = self._get_trait(marker)
            color = TRAIT_COLORS.get(trait, "#2C3E50")

        bars = ax1.bar(
            range(len(freq_df)),
            freq_df[freq_col].values,
            color=color,
            edgecolor="white",
            linewidth=0.5,
            alpha=0.85,
        )

        # Add percentage labels on bars
        for bar, val in zip(bars, freq_df[freq_col].values):
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.01,
                f"{val * 100:.1f}%",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

        ax1.set_xticks(range(len(freq_df)))
        ax1.set_xticklabels([str(a) for a in freq_df[allele_col].values], fontsize=11)
        ax1.set_ylabel("Allele Frequency", fontsize=12)
        ax1.set_xlabel("Allele", fontsize=12)
        ax1.set_title(title or f"Allele Frequency Distribution — {marker}", fontsize=14, fontweight="bold")
        ax1.set_ylim(0, max(freq_df[freq_col].values) * 1.2)
        ax1.grid(axis="y", alpha=0.3)

        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=self.dpi, bbox_inches="tight")
        return fig

    def genotype_distribution(
        self,
        df: pd.DataFrame,
        marker: str,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (8, 8),
        title: Optional[str] = None,
    ) -> plt.Figure:
        """Create a pie chart of genotype distribution for a marker.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        marker : str
            Marker name.
        save_path : str, optional
            Path to save figure.
        figsize : tuple
            Figure size.
        title : str, optional
            Plot title.

        Returns
        -------
        matplotlib.figure.Figure
        """
        from apple_mas.allele_analysis import AlleleAnalyzer

        analyzer = AlleleAnalyzer()
        freq_df = analyzer.genotype_frequencies(df, marker)

        fig, ax = plt.subplots(figsize=figsize)

        colors = ALLELE_PALETTE[:len(freq_df)]

        wedges, texts, autotexts = ax.pie(
            freq_df["count"].values,
            labels=freq_df["genotype"].values,
            autopct=lambda pct: f"{pct:.1f}%\n({int(round(pct / 100. * freq_df['count'].sum()))})",
            colors=colors,
            startangle=90,
            pctdistance=0.75,
            textprops={"fontsize": 10},
        )
        for autotext in autotexts:
            autotext.set_fontsize(9)
            autotext.set_fontweight("bold")

        ax.set_title(
            title or f"Genotype Distribution — {marker}",
            fontsize=14,
            fontweight="bold",
            pad=20,
        )

        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=self.dpi, bbox_inches="tight")
        return fig

    def multi_marker_frequency(
        self,
        df: pd.DataFrame,
        markers: Optional[List[str]] = None,
        save_path: Optional[str] = None,
        ncols: int = 3,
        figsize: Optional[Tuple[int, int]] = None,
    ) -> plt.Figure:
        """Create a multi-panel figure showing allele frequencies for multiple markers.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        markers : list of str, optional
            Markers to include. If None, uses all markers in the data.
        save_path : str, optional
            Path to save figure.
        ncols : int
            Number of columns in the grid.
        figsize : tuple, optional
            Figure size. Auto-calculated if None.

        Returns
        -------
        matplotlib.figure.Figure
        """
        from apple_mas.allele_analysis import AlleleAnalyzer

        analyzer = AlleleAnalyzer()

        if markers is None:
            markers = sorted(df["marker"].unique())

        nrows = (len(markers) + ncols - 1) // ncols
        if figsize is None:
            figsize = (6 * ncols, 5 * nrows)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows * ncols == 1:
            axes = np.array([axes])
        axes = axes.flatten()

        for idx, marker in enumerate(markers):
            ax = axes[idx]
            freq_df = analyzer.allele_frequencies(df, marker)
            allele_col = freq_df.columns[0]
            freq_col = "frequency"
            trait = self._get_trait(marker)
            color = TRAIT_COLORS.get(trait, "#2C3E50")

            ax.bar(
                range(len(freq_df)),
                freq_df[freq_col].values,
                color=color,
                edgecolor="white",
                alpha=0.85,
            )
            for i, (bar_x, val) in enumerate(
                zip(range(len(freq_df)), freq_df[freq_col].values)
            ):
                ax.text(
                    bar_x,
                    val + 0.01,
                    f"{val * 100:.0f}%",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )
            ax.set_xticks(range(len(freq_df)))
            ax.set_xticklabels(
                [str(a) for a in freq_df[allele_col].values], fontsize=9
            )
            ax.set_title(f"{marker}\n({trait})", fontsize=11, fontweight="bold")
            ax.set_ylabel("Freq", fontsize=9)
            ax.set_ylim(0, max(freq_df[freq_col].values) * 1.3)
            ax.grid(axis="y", alpha=0.3)

        # Hide unused axes
        for idx in range(len(markers), len(axes)):
            axes[idx].set_visible(False)

        fig.suptitle(
            "Allele Frequency Comparison Across Markers",
            fontsize=16,
            fontweight="bold",
            y=1.02,
        )
        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=self.dpi, bbox_inches="tight")
        return fig

    def genotype_heatmap(
        self,
        df: pd.DataFrame,
        markers: Optional[List[str]] = None,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (12, 8),
    ) -> plt.Figure:
        """Create a heatmap of genotype frequencies across all markers.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        markers : list of str, optional
            Markers to include.
        save_path : str, optional
            Path to save figure.
        figsize : tuple
            Figure size.

        Returns
        -------
        matplotlib.figure.Figure
        """
        from apple_mas.allele_analysis import AlleleAnalyzer

        analyzer = AlleleAnalyzer()

        if markers is None:
            markers = sorted(df["marker"].unique())

        # Build frequency matrix
        all_genotypes = set()
        marker_genos = {}
        for marker in markers:
            freq_df = analyzer.genotype_frequencies(df, marker)
            marker_genos[marker] = dict(zip(freq_df["genotype"], freq_df["frequency"]))
            all_genotypes.update(freq_df["genotype"].tolist())

        all_genotypes = sorted(all_genotypes)
        matrix = np.zeros((len(markers), len(all_genotypes)))
        for i, marker in enumerate(markers):
            for j, geno in enumerate(all_genotypes):
                matrix[i, j] = marker_genos.get(marker, {}).get(geno, 0)

        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            matrix,
            annot=True,
            fmt=".2f",
            xticklabels=all_genotypes,
            yticklabels=[f"{m}\n({self._get_trait(m)})" for m in markers],
            cmap="YlOrRd",
            linewidths=0.5,
            ax=ax,
            cbar_kws={"label": "Frequency"},
        )
        ax.set_title(
            "Genotype Frequency Heatmap Across Markers",
            fontsize=14,
            fontweight="bold",
        )
        plt.xticks(rotation=45, ha="right", fontsize=9)
        plt.yticks(fontsize=10)

        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=self.dpi, bbox_inches="tight")
        return fig

    def sample_composition_chart(
        self,
        df: pd.DataFrame,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (10, 6),
    ) -> plt.Figure:
        """Create a stacked bar chart showing sample composition across markers.

        Displays the proportion of each genotype per marker as a stacked bar.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        save_path : str, optional
            Path to save figure.
        figsize : tuple
            Figure size.

        Returns
        -------
        matplotlib.figure.Figure
        """
        from apple_mas.allele_analysis import AlleleAnalyzer

        analyzer = AlleleAnalyzer()
        markers = sorted(df["marker"].unique())

        # Build composition data
        composition = {}
        for marker in markers:
            freq_df = analyzer.genotype_frequencies(df, marker)
            composition[marker] = dict(zip(freq_df["genotype"], freq_df["percentage"]))

        all_genos = sorted(set(g for comp in composition.values() for g in comp))
        data = []
        for marker in markers:
            row = [composition[marker].get(g, 0) for g in all_genos]
            data.append(row)

        plot_df = pd.DataFrame(data, index=markers, columns=all_genos)

        fig, ax = plt.subplots(figsize=figsize)
        plot_df.plot(
            kind="barh",
            stacked=True,
            ax=ax,
            colormap="Set2",
            edgecolor="white",
            linewidth=0.5,
        )
        ax.set_xlabel("Percentage (%)", fontsize=12)
        ax.set_title(
            "Genotype Composition Across Markers",
            fontsize=14,
            fontweight="bold",
        )
        ax.legend(
            title="Genotype",
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            fontsize=8,
        )
        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=self.dpi, bbox_inches="tight")
        return fig

    def ranking_barplot(
        self,
        ranked_df: pd.DataFrame,
        top_n: int = 20,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (10, 8),
        title: str = "Parent Ranking by Total Score",
    ) -> plt.Figure:
        """Create a horizontal bar chart of parent rankings.

        Parameters
        ----------
        ranked_df : pd.DataFrame
            Output from ParentScorer.rank_accessions().
        top_n : int
            Number of top accessions to show.
        save_path : str, optional
            Path to save figure.
        figsize : tuple
            Figure size.
        title : str
            Plot title.

        Returns
        -------
        matplotlib.figure.Figure
        """
        top = ranked_df.head(top_n).copy()
        top = top.sort_values("total_score", ascending=True)

        fig, ax = plt.subplots(figsize=figsize)
        bars = ax.barh(
            range(len(top)),
            top["total_score"].values,
            color=sns.color_palette("RdYlGn", len(top)),
            edgecolor="white",
        )

        # Add score labels
        for bar, score in zip(bars, top["total_score"].values):
            ax.text(
                bar.get_width() + 0.01,
                bar.get_y() + bar.get_height() / 2,
                f"{score:.2f}",
                va="center",
                fontsize=9,
                fontweight="bold",
            )

        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["sample_id"].values, fontsize=10)
        ax.set_xlabel("Total Score", fontsize=12)
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)

        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=self.dpi, bbox_inches="tight")
        return fig

    def trait_score_radar(
        self,
        df: pd.DataFrame,
        sample_ids: List[str],
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (10, 10),
    ) -> plt.Figure:
        """Create a radar/spider chart comparing accessions across traits.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        sample_ids : list of str
            Accessions to compare (2-6 recommended).
        save_path : str, optional
            Path to save figure.
        figsize : tuple
            Figure size.

        Returns
        -------
        matplotlib.figure.Figure
        """
        from apple_mas.parent_scorer import ParentScorer

        scorer = ParentScorer()

        # Collect scores for each sample
        all_traits = set()
        sample_trait_scores = {}
        for sid in sample_ids:
            result = scorer.score_accession(df, sid)
            sample_trait_scores[sid] = {}
            for trait_name, ts in result["trait_scores"].items():
                all_traits.add(trait_name)
                sample_trait_scores[sid][trait_name] = ts["avg_score"]

        traits = sorted(all_traits)
        n_traits = len(traits)
        angles = np.linspace(0, 2 * np.pi, n_traits, endpoint=False).tolist()
        angles += angles[:1]  # close the polygon

        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))

        colors = sns.color_palette("husl", len(sample_ids))

        for idx, sid in enumerate(sample_ids):
            values = [sample_trait_scores[sid].get(t, 0) for t in traits]
            values += values[:1]  # close polygon
            ax.plot(angles, values, "o-", linewidth=2, label=sid, color=colors[idx])
            ax.fill(angles, values, alpha=0.1, color=colors[idx])

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(traits, fontsize=10)
        ax.set_title(
            "Trait Score Comparison",
            fontsize=14,
            fontweight="bold",
            pad=30,
        )
        ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))

        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=self.dpi, bbox_inches="tight")
        return fig

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
