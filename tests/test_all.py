"""
Tests for Apple MAS Toolkit.

Run with: pytest tests/ -v
"""

import json
import os
import sys
import tempfile
import pytest
import pandas as pd
import numpy as np

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from apple_mas.marker_db import MarkerDatabase
from apple_mas.data_parser import DataParser
from apple_mas.allele_analysis import AlleleAnalyzer
from apple_mas.parent_scorer import ParentScorer, BREEDING_PROFILES
from apple_mas.visualization import MarkerVisualizer
from apple_mas.selector import MarkerSelector
from apple_mas.analytics import (
    segregation_distortion_test,
    segregation_distortion_from_df,
    genotype_phenotype_boxplot,
    calculate_pic,
    calculate_pic_from_df,
    batch_segregation_test,
    batch_pic,
)


# ==========================================
# MarkerDatabase Tests
# ==========================================

class TestMarkerDatabase:
    """Tests for the MarkerDatabase class."""

    def setup_method(self):
        self.db = MarkerDatabase()

    def test_load_database(self):
        """Database loads successfully."""
        assert self.db.version is not None
        assert len(self.db.list_markers()) > 0

    def test_list_markers(self):
        """All expected markers are present."""
        markers = self.db.list_markers()
        expected = [
            "MYB10", "RED_TE", "MA_INDEL", "BP16", "BP13", "ACS", "ACO", "MD_PG1",
            "Vf_SCAR", "AL07", "PL2_SSR", "CO_INDEL", "CRISP_SSR",
        ]
        assert markers == expected

    def test_get_marker(self):
        """Retrieve a specific marker."""
        info = self.db.get_marker("MYB10")
        assert info["trait"] == "Flesh Color"
        assert info["chromosome"] == 9
        assert info["marker_type"] == "SSR"

    def test_get_marker_not_found(self):
        """KeyError raised for unknown marker."""
        with pytest.raises(KeyError):
            self.db.get_marker("NONEXISTENT")

    def test_get_markers_by_trait(self):
        """Retrieve markers by trait category."""
        color_markers = self.db.get_markers_by_trait("color")
        assert "MYB10" in color_markers
        assert "RED_TE" in color_markers

    def test_get_markers_by_trait_growth_habit(self):
        """Retrieve growth habit markers."""
        gh_markers = self.db.get_markers_by_trait("growth_habit")
        assert "CO_INDEL" in gh_markers

    def test_get_markers_by_chromosome(self):
        """Retrieve markers by chromosome."""
        chr9 = self.db.get_markers_by_chromosome("9")
        assert "MYB10" in chr9
        assert "RED_TE" in chr9

    def test_new_markers_present(self):
        """New disease resistance and growth habit markers load correctly."""
        for name in ["Vf_SCAR", "AL07", "PL2_SSR", "CO_INDEL", "CRISP_SSR"]:
            info = self.db.get_marker(name)
            assert "trait" in info
            assert "chromosome" in info

    def test_chromosomal_locations_updated(self):
        """New chromosomal locations are indexed."""
        chr1 = self.db.get_markers_by_chromosome("1")
        assert "Vf_SCAR" in chr1
        assert "AL07" in chr1
        chr11 = self.db.get_markers_by_chromosome("11")
        assert "PL2_SSR" in chr11

    def test_get_allele_info(self):
        """Retrieve allele information."""
        info = self.db.get_allele_info("MYB10", "390")
        assert "phenotype" in info
        assert "white" in info["phenotype"].lower()

    def test_get_genotype_phenotype(self):
        """Retrieve genotype-phenotype mapping."""
        result = self.db.get_genotype_phenotype("MYB10", "390:390")
        assert "phenotype" in result

    def test_get_breeding_notes(self):
        """Retrieve breeding notes."""
        notes = self.db.get_breeding_notes("MYB10")
        assert len(notes) > 0

    def test_search_markers(self):
        """Search markers by keyword."""
        results = self.db.search_markers("texture")
        assert len(results) > 0
        assert "ACS" in results or "ACO" in results or "MD_PG1" in results

    def test_search_scab_markers(self):
        """Search markers finds apple scab resistance markers."""
        results = self.db.search_markers("scab")
        assert len(results) > 0
        assert "Vf_SCAR" in results

    def test_get_pcr_methodology(self):
        """Retrieve PCR methodology."""
        method = self.db.get_pcr_methodology()
        assert "PCR_conditions" in method

    def test_summary(self):
        """Summary returns a non-empty string."""
        summary = self.db.summary()
        assert isinstance(summary, str)
        assert len(summary) > 0

    def test_repr(self):
        """String representation works."""
        r = repr(self.db)
        assert "MarkerDatabase" in r


# ==========================================
# DataParser Tests
# ==========================================

class TestDataParser:
    """Tests for the DataParser class."""

    def setup_method(self):
        self.parser = DataParser()

    def test_create_sample_dataset(self):
        """Sample dataset creation works."""
        df, raw = self.parser.create_sample_dataset()
        assert len(df) > 0
        assert "sample_id" in df.columns
        assert "marker" in df.columns
        assert "genotype" in df.columns
        assert df["marker"].nunique() == 8

    def test_from_dict(self):
        """Parse data from dictionary."""
        data = {
            "S1": {"MYB10": "390:390", "RED_TE": "750"},
            "S2": {"MYB10": "390:490", "RED_TE": "NB"},
        }
        df = self.parser.from_dict(data)
        assert len(df) == 4
        assert df["sample_id"].nunique() == 2

    def test_from_csv(self):
        """Parse data from CSV file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("sample_id,marker,genotype\n")
            f.write("S1,MYB10,390:390\n")
            f.write("S1,RED_TE,750\n")
            f.write("S2,MYB10,390:490\n")
            f.write("S2,RED_TE,NB\n")
            f.flush()
            fpath = f.name

        try:
            df = self.parser.from_csv(fpath)
            assert len(df) == 4
        finally:
            os.unlink(fpath)

    def test_from_wide_csv(self):
        """Parse wide-format CSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("sample_id,MYB10,RED_TE\n")
            f.write("S1,390:390,750\n")
            f.write("S2,390:490,NB\n")
            f.flush()
            fpath = f.name

        try:
            df = self.parser.from_wide_csv(fpath)
            assert len(df) == 4
        finally:
            os.unlink(fpath)

    def test_is_no_amp(self):
        """No-amplification detection."""
        assert DataParser._is_no_amp("*")
        assert DataParser._is_no_amp("NA")
        assert DataParser._is_no_amp("")
        assert not DataParser._is_no_amp("390:390")
        assert not DataParser._is_no_amp("750")

    def test_validate(self):
        """Data validation works."""
        df, _ = self.parser.create_sample_dataset()
        is_valid, issues = self.parser.validate(df)
        assert is_valid

    def test_validate_missing_columns(self):
        """Validation catches missing columns."""
        df = pd.DataFrame({"a": [1], "b": [2]})
        is_valid, issues = self.parser.validate(df)
        assert not is_valid


# ==========================================
# AlleleAnalyzer Tests
# ==========================================

class TestAlleleAnalyzer:
    """Tests for the AlleleAnalyzer class."""

    def setup_method(self):
        parser = DataParser()
        self.df, _ = parser.create_sample_dataset()
        self.analyzer = AlleleAnalyzer()

    def test_allele_frequencies(self):
        """Allele frequency calculation works."""
        freq = self.analyzer.allele_frequencies(self.df, "MYB10")
        assert len(freq) > 0
        assert "frequency" in freq.columns
        assert abs(freq["frequency"].sum() - 1.0) < 0.01

    def test_genotype_frequencies(self):
        """Genotype frequency calculation works."""
        geno = self.analyzer.genotype_frequencies(self.df, "MYB10")
        assert len(geno) > 0
        assert "count" in geno.columns
        assert "percentage" in geno.columns

    def test_all_markers_summary(self):
        """Cross-marker summary works."""
        summary = self.analyzer.all_markers_summary(self.df)
        assert len(summary) == 8
        assert "marker" in summary.columns
        assert "n_samples" in summary.columns

    def test_pic(self):
        """PIC calculation returns valid value."""
        pic = self.analyzer.polymorphism_information_content(self.df, "ACS")
        assert 0 <= pic <= 1

    def test_heterozygosity(self):
        """Heterozygosity calculation works."""
        het = self.analyzer.heterozygosity(self.df, "ACS")
        assert "observed_het" in het
        assert "expected_het" in het
        assert 0 <= het["observed_het"] <= 1

    def test_hwe_test(self):
        """HWE test works for biallelic markers."""
        result = self.analyzer.hardy_weinberg_test(self.df, "RED_TE")
        assert "chi2" in result
        assert "p_value" in result

    def test_cross_marker_summary(self):
        """Cross-marker diversity summary works."""
        summary = self.analyzer.cross_marker_summary(self.df)
        assert "PIC" in summary.columns
        assert "obs_het" in summary.columns
        assert len(summary) == 8


# ==========================================
# ParentScorer Tests
# ==========================================

class TestParentScorer:
    """Tests for the ParentScorer class."""

    def setup_method(self):
        parser = DataParser()
        self.df, _ = parser.create_sample_dataset()
        self.scorer = ParentScorer()

    def test_score_accession(self):
        """Single accession scoring works."""
        sample_id = self.df["sample_id"].iloc[0]
        result = self.scorer.score_accession(self.df, sample_id)
        assert "total_score" in result
        assert result["total_score"] >= 0
        assert "marker_scores" in result

    def test_rank_accessions(self):
        """Accession ranking works."""
        ranked = self.scorer.rank_accessions(self.df, profile="balanced")
        assert len(ranked) > 0
        assert "total_score" in ranked.columns
        assert ranked.iloc[0]["total_score"] >= ranked.iloc[-1]["total_score"]

    def test_rank_with_profile(self):
        """Ranking with different profiles works."""
        for profile in BREEDING_PROFILES:
            ranked = self.scorer.rank_accessions(self.df, profile=profile)
            assert len(ranked) > 0

    def test_find_complementary_pairs(self):
        """Complementary parent finding works."""
        pairs = self.scorer.find_complementary_parents(self.df, n_pairs=3)
        assert len(pairs) == 3
        assert "parent_1" in pairs.columns
        assert "parent_2" in pairs.columns
        assert "complementarity" in pairs.columns

    def test_list_profiles(self):
        """Profile listing works."""
        profiles = self.scorer.list_profiles()
        assert len(profiles) == len(BREEDING_PROFILES)

    def test_available_profiles(self):
        """Static profile list works."""
        profiles = ParentScorer.available_profiles()
        assert "premium_table" in profiles
        assert "balanced" in profiles

    def test_invalid_profile(self):
        """Invalid profile raises ValueError."""
        with pytest.raises(ValueError):
            self.scorer.rank_accessions(self.df, profile="nonexistent")


# ==========================================
# MarkerVisualizer Tests
# ==========================================

class TestMarkerVisualizer:
    """Tests for the MarkerVisualizer class."""

    def setup_method(self):
        parser = DataParser()
        self.df, _ = parser.create_sample_dataset()
        self.viz = MarkerVisualizer()

    def test_allele_frequency_plot(self):
        """Allele frequency plot generates without error."""
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            path = f.name
        try:
            fig = self.viz.allele_frequency_plot(self.df, "MYB10", save_path=path)
            assert fig is not None
            assert os.path.exists(path)
        finally:
            os.unlink(path)

    def test_genotype_distribution(self):
        """Genotype distribution plot generates."""
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            path = f.name
        try:
            fig = self.viz.genotype_distribution(self.df, "MYB10", save_path=path)
            assert fig is not None
            assert os.path.exists(path)
        finally:
            os.unlink(path)

    def test_multi_marker_frequency(self):
        """Multi-marker frequency plot generates."""
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            path = f.name
        try:
            fig = self.viz.multi_marker_frequency(self.df, save_path=path)
            assert fig is not None
            assert os.path.exists(path)
        finally:
            os.unlink(path)

    def test_genotype_heatmap(self):
        """Genotype heatmap generates."""
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            path = f.name
        try:
            fig = self.viz.genotype_heatmap(self.df, save_path=path)
            assert fig is not None
            assert os.path.exists(path)
        finally:
            os.unlink(path)

    def test_ranking_barplot(self):
        """Ranking barplot generates."""
        scorer = ParentScorer()
        ranked = scorer.rank_accessions(self.df, profile="balanced")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            path = f.name
        try:
            fig = self.viz.ranking_barplot(ranked, save_path=path)
            assert fig is not None
            assert os.path.exists(path)
        finally:
            os.unlink(path)


# ==========================================
# MarkerSelector Tests
# ==========================================

class TestMarkerSelector:
    """Tests for the MarkerSelector class."""

    def setup_method(self):
        self.selector = MarkerSelector()

    def test_list_available_traits(self):
        """Trait listing works."""
        traits = self.selector.list_available_traits()
        assert len(traits) > 0
        assert "category" in traits.columns
        assert "marker_count" in traits.columns

    def test_list_available_markers(self):
        """Marker listing works."""
        markers = self.selector.list_available_markers()
        assert len(markers) >= 13
        assert "marker" in markers.columns
        assert "chromosome" in markers.columns

    def test_search_markers(self):
        """Marker search works."""
        results = self.selector.search_markers("firmness")
        assert len(results) > 0

    def test_get_markers_for_trait(self):
        """Trait-to-marker retrieval works."""
        markers = self.selector.get_markers_for_trait("disease_resistance")
        assert "Vf_SCAR" in markers
        assert "PL2_SSR" in markers

    def test_get_marker_primer_info(self):
        """Primer information retrieval works."""
        info = self.selector.get_marker_primer_info("Vf_SCAR")
        assert "primer_sequences" in info
        assert "forward" in info["primer_sequences"]

    def test_get_expected_alleles(self):
        """Allele definition retrieval works."""
        alleles = self.selector.get_expected_alleles("MA_INDEL")
        assert "409" in alleles
        assert "456" in alleles

    def test_get_marker_plan_table(self):
        """Genotyping plan table generation works."""
        table = self.selector.get_marker_plan_table(["Vf_SCAR", "ACS"])
        assert len(table) == 2
        assert "primer_forward" in table.columns

    def test_build_genotyping_plan(self):
        """Genotyping plan construction works."""
        plan = self.selector.build_genotyping_plan(
            target_traits=["fruit_quality", "disease_resistance"],
            population_type="BC1",
            population_size=100,
        )
        assert plan["marker_count"] > 0
        assert plan["population"]["type"] == "BC1"
        assert plan["population"]["size"] == 100
        assert len(plan["markers"]) > 0

    def test_build_plan_with_extra_markers(self):
        """Plan with additional markers works."""
        plan = self.selector.build_genotyping_plan(
            target_traits=["growth_habit"],
            additional_markers=["MYB10"],
        )
        assert "CO_INDEL" in plan["markers"]
        assert "MYB10" in plan["markers"]

    def test_build_plan_invalid_traits(self):
        """Invalid trait raises ValueError."""
        with pytest.raises(ValueError):
            self.selector.build_genotyping_plan(target_traits=["nonexistent_trait"])

    def test_export_plan_json(self):
        """Plan export to JSON works."""
        plan = self.selector.build_genotyping_plan(target_traits=["color"])
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            path = f.name
        try:
            self.selector.export_plan(plan, path, fmt="json")
            assert os.path.exists(path)
            with open(path) as fh:
                loaded = json.load(fh)
            assert loaded["marker_count"] > 0
        finally:
            os.unlink(path)

    def test_export_plan_csv(self):
        """Plan export to CSV works."""
        plan = self.selector.build_genotyping_plan(target_traits=["color"])
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            path = f.name
        try:
            self.selector.export_plan(plan, path, fmt="csv")
            assert os.path.exists(path)
            df = pd.read_csv(path)
            assert len(df) > 0
        finally:
            os.unlink(path)

    def test_plan_summary(self):
        """Plan summary generation works."""
        plan = self.selector.build_genotyping_plan(target_traits=["growth_habit"])
        summary = self.selector.plan_summary(plan)
        assert "GENOTYPING PLAN SUMMARY" in summary
        assert "CO_INDEL" in summary


# ==========================================
# Analytics Tests
# ==========================================

class TestAnalytics:
    """Tests for the analytics module."""

    def setup_method(self):
        parser = DataParser()
        self.df, _ = parser.create_sample_dataset()

    def test_segregation_distortion_test_basic(self):
        """Basic segregation distortion test works."""
        result = segregation_distortion_test(
            observed_genotypes={"AA": 50, "AB": 50, "BB": 0},
            population_type="F1",
        )
        assert "chi2" in result
        assert "p_value" in result
        assert "interpretation" in result

    def test_segregation_distortion_no_distortion(self):
        """Balanced F2 population shows no distortion."""
        result = segregation_distortion_test(
            observed_genotypes={"AA": 25, "AB": 50, "BB": 25},
            population_type="F2",
        )
        assert result["significant"] == False

    def test_segregation_distortion_with_distortion(self):
        """Unbalanced population shows significant distortion."""
        result = segregation_distortion_test(
            observed_genotypes={"AA": 90, "AB": 10, "BB": 0},
            population_type="BC1",
        )
        assert result["significant"] == True

    def test_segregation_distortion_empty(self):
        """Empty data returns appropriate result."""
        result = segregation_distortion_test(
            observed_genotypes={},
            population_type="F1",
        )
        assert result["chi2"] is None

    def test_segregation_distortion_from_df(self):
        """Segregation test from DataFrame works."""
        result = segregation_distortion_from_df(
            self.df, marker="ACS", population_type="F1",
        )
        assert "chi2" in result
        assert result["population_type"] == "F1"

    def test_calculate_pic_biallelic(self):
        """PIC calculation for bi-allelic marker."""
        pic = calculate_pic({"A": 0.5, "B": 0.5})
        assert abs(pic - 0.375) < 0.001

    def test_calculate_pic_monoallelic(self):
        """PIC for monomorphic marker is 0."""
        pic = calculate_pic({"A": 1.0})
        assert pic == 0.0

    def test_calculate_pic_multiallelic(self):
        """PIC calculation for multi-allelic marker."""
        pic = calculate_pic({"A": 0.4, "B": 0.3, "C": 0.2, "D": 0.1})
        assert 0 < pic < 1

    def test_calculate_pic_from_df(self):
        """PIC from DataFrame works."""
        pic = calculate_pic_from_df(self.df, "ACS")
        assert 0 <= pic <= 1

    def test_batch_segregation_test(self):
        """Batch segregation tests work."""
        result_df = batch_segregation_test(self.df, population_type="F1")
        assert len(result_df) == 8
        assert "chi2" in result_df.columns
        assert "significant" in result_df.columns

    def test_batch_pic(self):
        """Batch PIC calculation works."""
        result_df = batch_pic(self.df)
        assert len(result_df) == 8
        assert "PIC" in result_df.columns
        assert "n_alleles" in result_df.columns

    def test_genotype_phenotype_boxplot(self):
        """Boxplot generation works."""
        genotypes = ["AA", "AA", "AB", "AB", "BB", "BB"] * 5
        phenotypes = [10, 11, 15, 14, 20, 19] * 5
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            path = f.name
        try:
            fig = genotype_phenotype_boxplot(
                genotypes=genotypes,
                phenotypic_values=phenotypes,
                marker_name="TEST",
                trait_name="Test Trait",
                save_path=path,
            )
            assert fig is not None
            assert os.path.exists(path)
        finally:
            os.unlink(path)
