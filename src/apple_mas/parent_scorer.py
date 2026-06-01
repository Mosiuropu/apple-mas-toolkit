"""
Parent Scorer - Rank and score apple accessions for breeding suitability.

Provides methods to score accessions based on desired allele combinations,
rank parents by trait suitability, and generate crossing recommendations.
"""

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# Default scoring weights for each marker-trait combination
DEFAULT_TRAIT_WEIGHTS = {
    "Flesh Color": {
        "MYB10": {
            "490:490": 10,
            "390:490": 7,
            "390:390": 3,
        },
    },
    "Skin Color": {
        "RED_TE": {
            "750": 8,
            "NB": 3,
        },
    },
    "Fruit Acidity": {
        "MA_INDEL": {
            "409:409": 7,
            "409:456": 6,
            "409:452": 5,
            "387:402": 4,
            "456:456": 5,
        },
    },
    "Bitter Pit Resistance": {
        "BP16": {
            "203:203": 10,
            "203:206": 8,
            "203:368": 5,
            "206:206": 4,
            "206:368": 3,
            "368:368": 1,
        },
        "BP13": {
            "226:236": 10,
            "224:232": 6,
            "232:250": 4,
            "232:252": 4,
            "250:252": 5,
            "232:232": 3,
        },
    },
    "Fruit Texture (Firmness)": {
        "ACS": {
            "341:341": 10,
            "200:341": 7,
            "200:200": 5,
        },
        "ACO": {
            "237:237": 10,
            "237:300": 7,
            "300:300": 4,
        },
        "MD_PG1": {
            "289:292": 9,
            "289:289": 9,
            "292:292": 8,
            "289:298": 6,
            "292:298": 6,
            "298:298": 3,
        },
    },
}

# Composite trait profiles for common breeding objectives
BREEDING_PROFILES = {
    "premium_table": {
        "description": "Premium table apple with red skin, balanced acidity, firm texture, and low bitter pit risk",
        "traits": {
            "Skin Color": 1.0,
            "Fruit Acidity": 1.5,
            "Fruit Texture (Firmness)": 2.0,
            "Bitter Pit Resistance": 1.5,
        },
    },
    "juice_cider": {
        "description": "Apple for juice/cider production with high acidity and good color",
        "traits": {
            "Skin Color": 0.5,
            "Fruit Acidity": 2.5,
            "Fruit Texture (Firmness)": 0.5,
            "Bitter Pit Resistance": 0.5,
        },
    },
    "red_flesh_novelty": {
        "description": "Novel red-fleshed apple for specialty market with good overall quality",
        "traits": {
            "Flesh Color": 3.0,
            "Skin Color": 1.0,
            "Fruit Acidity": 1.0,
            "Fruit Texture (Firmness)": 1.5,
            "Bitter Pit Resistance": 1.0,
        },
    },
    "storage_shelf_life": {
        "description": "Apple optimized for long storage and shelf life",
        "traits": {
            "Skin Color": 1.0,
            "Fruit Acidity": 1.0,
            "Fruit Texture (Firmness)": 3.0,
            "Bitter Pit Resistance": 2.0,
        },
    },
    "balanced": {
        "description": "Balanced selection across all traits with equal weighting",
        "traits": {
            "Skin Color": 1.0,
            "Fruit Acidity": 1.0,
            "Fruit Texture (Firmness)": 1.0,
            "Bitter Pit Resistance": 1.0,
        },
    },
}


class ParentScorer:
    """Score and rank apple accessions for breeding suitability.

    Provides scoring based on desired allele combinations, composite trait
    profiles, and customizable weights. Generates ranked tables and crossing
    recommendations.

    Examples
    --------
    >>> from apple_mas import DataParser, ParentScorer
    >>> parser = DataParser()
    >>> df, _ = parser.create_sample_dataset()
    >>> scorer = ParentScorer()
    >>> ranked = scorer.rank_accessions(df, profile="premium_table")
    >>> print(ranked.head(10))
    """

    def __init__(
        self,
        trait_weights: Optional[Dict] = None,
        no_amplification_penalty: float = 0.5,
    ):
        """Initialize the ParentScorer.

        Parameters
        ----------
        trait_weights : dict, optional
            Custom trait weights. If None, uses DEFAULT_TRAIT_WEIGHTS.
        no_amplification_penalty : float
            Penalty multiplier (0-1) for accessions with no amplification at any marker.
            0 = no penalty, 1 = complete exclusion.
        """
        self.trait_weights = trait_weights or DEFAULT_TRAIT_WEIGHTS
        self.no_amp_penalty = no_amplification_penalty

    def score_accession(
        self,
        df: pd.DataFrame,
        sample_id: str,
        trait_weights: Optional[Dict] = None,
    ) -> Dict[str, Any]:
        """Score a single accession across all markers.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        sample_id : str
            Accession identifier to score.
        trait_weights : dict, optional
            Override default weights for this scoring.

        Returns
        -------
        dict
            Scoring details including per-marker scores, trait scores, and total.
        """
        weights = trait_weights or self.trait_weights
        sample_df = df[df["sample_id"] == sample_id]

        marker_scores = {}
        trait_scores = {}
        no_amp_markers = []

        for trait_name, marker_weights in weights.items():
            trait_total = 0.0
            trait_count = 0

            for marker, geno_scores in marker_weights.items():
                marker_data = sample_df[sample_df["marker"] == marker]
                if marker_data.empty:
                    marker_scores[marker] = {"genotype": "N/A", "score": 0, "max": 10}
                    no_amp_markers.append(marker)
                    continue

                genotype = marker_data.iloc[0]["genotype"]
                no_amp = marker_data.iloc[0].get("no_amplification", False)

                if no_amp:
                    score = 0
                    no_amp_markers.append(marker)
                else:
                    # Try exact match, then partial match
                    score = geno_scores.get(genotype, 0)
                    if score == 0:
                        # Try reverse genotype (e.g., 341:200 for 200:341)
                        if ":" in str(genotype):
                            parts = str(genotype).split(":")
                            reverse = f"{parts[1].strip()}:{parts[0].strip()}"
                            score = geno_scores.get(reverse, 0)

                max_score = max(geno_scores.values()) if geno_scores else 10
                marker_scores[marker] = {
                    "genotype": genotype,
                    "score": score,
                    "max": max_score,
                    "normalized": score / max_score if max_score > 0 else 0,
                }
                trait_total += score
                trait_count += 1

            trait_scores[trait_name] = {
                "raw_score": trait_total,
                "n_markers": trait_count,
                "avg_score": trait_total / trait_count if trait_count > 0 else 0,
            }

        # Calculate total composite score
        total_raw = 0.0
        total_weight = 0.0
        for trait_name, ts in trait_scores.items():
            w = 1.0
            total_raw += ts["raw_score"] * w
            total_weight += w

        total_score = total_raw / total_weight if total_weight > 0 else 0

        # Apply no-amplification penalty
        if no_amp_markers and self.no_amp_penalty > 0:
            penalty_factor = 1.0 - (len(no_amp_markers) / len(marker_scores) * self.no_amp_penalty)
            total_score *= max(0, penalty_factor)

        return {
            "sample_id": sample_id,
            "marker_scores": marker_scores,
            "trait_scores": trait_scores,
            "total_score": round(total_score, 4),
            "no_amp_markers": no_amp_markers,
            "n_markers_scored": len(marker_scores) - len(no_amp_markers),
        }

    def rank_accessions(
        self,
        df: pd.DataFrame,
        profile: Optional[str] = None,
        trait_weights: Optional[Dict] = None,
        top_n: Optional[int] = None,
    ) -> pd.DataFrame:
        """Rank all accessions based on a breeding profile.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        profile : str, optional
            Name of a predefined breeding profile (see BREEDING_PROFILES).
            If None and trait_weights is None, uses 'balanced'.
        trait_weights : dict, optional
            Custom trait weights dictionary.
        top_n : int, optional
            Return only the top N accessions.

        Returns
        -------
        pd.DataFrame
            Ranked accessions with scores per trait.
        """
        # Resolve weights
        if trait_weights is None and profile is not None:
            if profile not in BREEDING_PROFILES:
                available = ", ".join(BREEDING_PROFILES.keys())
                raise ValueError(
                    f"Profile '{profile}' not found. Available: {available}"
                )
            profile_info = BREEDING_PROFILES[profile]
            trait_weights = self._profile_to_weights(profile_info["traits"])

        if trait_weights is None:
            trait_weights = self._profile_to_weights(
                {t: 1.0 for t in self.trait_weights.keys()}
            )

        # Score each accession
        results = []
        for sample_id in df["sample_id"].unique():
            score = self.score_accession(df, sample_id, trait_weights)
            row = {"sample_id": sample_id, "total_score": score["total_score"]}
            for trait_name, ts in score["trait_scores"].items():
                row[f"{trait_name}_score"] = round(ts["raw_score"], 2)
                row[f"{trait_name}_avg"] = round(ts["avg_score"], 2)
            row["n_no_amp"] = len(score["no_amp_markers"])
            results.append(row)

        ranked = pd.DataFrame(results)
        ranked = ranked.sort_values("total_score", ascending=False).reset_index(drop=True)
        ranked.index += 1
        ranked.index.name = "rank"

        if top_n:
            ranked = ranked.head(top_n)

        return ranked

    def find_complementary_parents(
        self,
        df: pd.DataFrame,
        n_pairs: int = 5,
        profile: Optional[str] = None,
    ) -> pd.DataFrame:
        """Find complementary parent pairs for crossing.

        Identifies pairs of accessions that complement each other's
        weaknesses, maximizing the probability of desirable offspring.

        Parameters
        ----------
        df : pd.DataFrame
            Genotype data.
        n_pairs : int
            Number of complementary pairs to return.
        profile : str, optional
            Breeding profile to use for scoring.

        Returns
        -------
        pd.DataFrame
            Complementary parent pairs with combined scores and reasoning.
        """
        ranked = self.rank_accessions(df, profile=profile)
        sample_ids = ranked["sample_id"].tolist()

        pairs = []
        for i in range(len(sample_ids)):
            for j in range(i + 1, min(i + 20, len(sample_ids))):
                s1, s2 = sample_ids[i], sample_ids[j]
                row1 = ranked[ranked["sample_id"] == s1].iloc[0]
                row2 = ranked[ranked["sample_id"] == s2].iloc[0]

                # Calculate complementarity: where one is weak, the other is strong
                complementarity = 0
                n_traits = 0
                for col in ranked.columns:
                    if col.endswith("_avg"):
                        trait_name = col.replace("_avg", "")
                        v1 = row1[col]
                        v2 = row2[col]
                        if v1 + v2 > 0:
                            complementarity += abs(v1 - v2) / max(v1, v2)
                        n_traits += 1

                combined_score = (row1["total_score"] + row2["total_score"]) / 2
                comp_ratio = complementarity / n_traits if n_traits > 0 else 0

                pairs.append({
                    "parent_1": s1,
                    "parent_2": s2,
                    "parent_1_score": round(row1["total_score"], 4),
                    "parent_2_score": round(row2["total_score"], 4),
                    "combined_score": round(combined_score, 4),
                    "complementarity": round(comp_ratio, 4),
                    "fitness_score": round(combined_score * (1 + comp_ratio * 0.3), 4),
                })

        pairs_df = pd.DataFrame(pairs)
        pairs_df = pairs_df.sort_values("fitness_score", ascending=False)
        pairs_df.index += 1
        pairs_df.index.name = "rank"
        return pairs_df.head(n_pairs)

    def list_profiles(self) -> pd.DataFrame:
        """List all available breeding profiles.

        Returns
        -------
        pd.DataFrame
            Profile names, descriptions, and trait weightings.
        """
        rows = []
        for name, info in BREEDING_PROFILES.items():
            rows.append({
                "profile": name,
                "description": info["description"],
                "traits": ", ".join(f"{k}: {v}" for k, v in info["traits"].items()),
            })
        return pd.DataFrame(rows)

    def _profile_to_weights(
        self, trait_multipliers: Dict[str, float]
    ) -> Dict[str, Any]:
        """Convert a trait multiplier dict into full scoring weights."""
        weights = {}
        for trait_name, multiplier in trait_multipliers.items():
            if trait_name in self.trait_weights:
                marker_weights = {}
                for marker, geno_scores in self.trait_weights[trait_name].items():
                    marker_weights[marker] = {
                        g: s * multiplier for g, s in geno_scores.items()
                    }
                weights[trait_name] = marker_weights
        return weights

    @staticmethod
    def available_profiles() -> List[str]:
        """List available predefined breeding profile names."""
        return list(BREEDING_PROFILES.keys())
