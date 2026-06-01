"""
Tests for Apple MAS Toolkit.

Run with: pytest tests/ -v
"""

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
        expected = ["MYB10", "RED_TE", "MA_INDEL", "BP16", "BP13", "ACS", "ACO", "MD_PG1"]
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

    def test_get_markers_by_chromosome(self):
        """Retrieve markers by chromosome."""
        chr9 = self.db.get_markers_by_chromosome("9")
        assert "MYB10" in chr9
        assert "RED_TE" in chr9

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
