"""
Command-line interface for Apple MAS Toolkit.

Provides commands for analyzing genotype data, scoring accessions,
generating reports, and visualizing results.
"""

import argparse
import os
import sys


def main():
    """Main entry point for the apple-mas CLI."""
    parser = argparse.ArgumentParser(
        prog="apple-mas",
        description="Apple Marker-Assisted Selection (MAS) Toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  apple-mas info                        Show marker database summary
  apple-mas info --marker MYB10         Show details for a specific marker
  apple-mas analyze --input data.csv    Analyze genotype data
  apple-mas rank --input data.csv --profile premium_table
  apple-mas visualize --input data.csv  Generate all plots
  apple-mas sample-data                 Create sample dataset for testing
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # --- info command ---
    info_parser = subparsers.add_parser("info", help="Show marker database information")
    info_parser.add_argument("--marker", "-m", help="Show details for a specific marker")
    info_parser.add_argument("--trait", "-t", help="Filter markers by trait category")
    info_parser.add_argument("--search", "-s", help="Search markers by keyword")
    info_parser.add_argument("--list", "-l", action="store_true", help="List all markers")
    info_parser.add_argument("--pcr", action="store_true", help="Show PCR methodology")
    info_parser.add_argument("--json", action="store_true", help="Output as JSON")

    # --- analyze command ---
    analyze_parser = subparsers.add_parser("analyze", help="Analyze genotype data")
    analyze_parser.add_argument("--input", "-i", required=True, help="Input CSV/Excel file")
    analyze_parser.add_argument("--format", "-f", choices=["csv", "excel", "wide", "dict"], default="csv", help="Input format")
    analyze_parser.add_argument("--sample-col", default="sample_id", help="Sample ID column name")
    analyze_parser.add_argument("--marker-col", default="marker", help="Marker column name")
    analyze_parser.add_argument("--genotype-col", default="genotype", help="Genotype column name")
    analyze_parser.add_argument("--marker", "-m", help="Analyze a specific marker only")
    analyze_parser.add_argument("--output", "-o", help="Output directory for results")
    analyze_parser.add_argument("--json", action="store_true", help="Output as JSON")

    # --- rank command ---
    rank_parser = subparsers.add_parser("rank", help="Rank accessions for breeding")
    rank_parser.add_argument("--input", "-i", required=True, help="Input CSV/Excel file")
    rank_parser.add_argument("--profile", "-p", choices=["premium_table", "juice_cider", "red_flesh_novelty", "storage_shelf_life", "balanced"], default="balanced", help="Breeding profile")
    rank_parser.add_argument("--top", "-n", type=int, default=20, help="Show top N accessions")
    rank_parser.add_argument("--pairs", type=int, help="Find complementary parent pairs")
    rank_parser.add_argument("--output", "-o", help="Output CSV file")
    rank_parser.add_argument("--format", "-f", choices=["csv", "excel"], default="csv", help="Output format")

    # --- visualize command ---
    viz_parser = subparsers.add_parser("visualize", help="Generate visualization plots")
    viz_parser.add_argument("--input", "-i", required=True, help="Input CSV/Excel file")
    viz_parser.add_argument("--output", "-o", default="plots", help="Output directory for plots")
    viz_parser.add_argument("--format", "-f", choices=["png", "pdf", "svg"], default="png", help="Image format")
    viz_parser.add_argument("--dpi", type=int, default=300, help="Image resolution")
    viz_parser.add_argument("--marker", "-m", help="Generate plots for a specific marker only")
    viz_parser.add_argument("--type", "-t", choices=["all", "frequency", "genotype", "heatmap", "composition", "ranking"], default="all", help="Plot type")

    # --- sample-data command ---
    sample_parser = subparsers.add_parser("sample-data", help="Create sample dataset for testing")
    sample_parser.add_argument("--output", "-o", default="sample_data", help="Output directory")
    sample_parser.add_argument("--format", "-f", choices=["csv", "excel", "both"], default="both", help="Output format")

    # --- report command ---
    report_parser = subparsers.add_parser("report", help="Generate comprehensive analysis report")
    report_parser.add_argument("--input", "-i", required=True, help="Input CSV/Excel file")
    report_parser.add_argument("--profile", "-p", default="balanced", help="Breeding profile for ranking")
    report_parser.add_argument("--output", "-o", default="report", help="Output directory")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    # Dispatch commands
    if args.command == "info":
        _cmd_info(args)
    elif args.command == "analyze":
        _cmd_analyze(args)
    elif args.command == "rank":
        _cmd_rank(args)
    elif args.command == "visualize":
        _cmd_visualize(args)
    elif args.command == "sample-data":
        _cmd_sample_data(args)
    elif args.command == "report":
        _cmd_report(args)


def _cmd_info(args):
    """Show marker database information."""
    from apple_mas.marker_db import MarkerDatabase
    import json

    db = MarkerDatabase()

    if args.marker:
        try:
            marker = db.get_marker(args.marker)
            if args.json:
                print(json.dumps(marker, indent=2, default=str))
            else:
                _print_marker_info(args.marker, marker)
        except KeyError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.trait:
        try:
            markers = db.get_markers_by_trait(args.trait)
            desc = db.get_trait_description(args.trait)
            print(f"\nTrait: {args.trait}")
            print(f"Description: {desc}")
            print(f"Markers: {', '.join(markers)}")
        except KeyError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.search:
        results = db.search_markers(args.search)
        if results:
            print(f"\nMarkers matching '{args.search}':")
            for m in results:
                info = db.get_marker(m)
                print(f"  {m:12s} -> {info['trait']}")
        else:
            print(f"No markers found matching '{args.search}'")
    elif args.pcr:
        methodology = db.get_pcr_methodology()
        if args.json:
            print(json.dumps(methodology, indent=2))
        else:
            _print_pcr_methodology(methodology)
    elif args.list:
        for marker in db.list_markers():
            info = db.get_marker(marker)
            print(f"  {marker:12s} | Chr {info.get('chromosome', '?'):2s} | {info['marker_type']:4s} | {info['trait']}")
    else:
        print(db.summary())


def _cmd_analyze(args):
    """Analyze genotype data."""
    from apple_mas.data_parser import DataParser
    from apple_mas.allele_analysis import AlleleAnalyzer
    import json

    parser = DataParser()
    analyzer = AlleleAnalyzer()

    # Load data
    try:
        if args.format == "csv":
            df = parser.from_csv(args.input, args.sample_col, args.marker_col, args.genotype_col)
        elif args.format == "wide":
            df = parser.from_wide_csv(args.input, args.sample_col)
        elif args.format == "excel":
            df = parser.from_excel(args.input, sample_col=args.sample_col, marker_col=args.marker_col, allele_col=args.genotype_col)
        else:
            print(f"Format '{args.format}' not yet implemented.", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(df)} genotype records for {df['sample_id'].nunique()} samples across {df['marker'].nunique()} markers\n")

    if args.marker:
        markers = [args.marker]
    else:
        markers = sorted(df["marker"].unique())

    for marker in markers:
        try:
            freq = analyzer.allele_frequencies(df, marker)
            geno = analyzer.genotype_frequencies(df, marker)
            het = analyzer.heterozygosity(df, marker)
            pic = analyzer.polymorphism_information_content(df, marker)

            print(f"{'='*50}")
            print(f"  Marker: {marker}")
            print(f"{'='*50}")
            print(f"\n  Allele Frequencies:")
            for _, row in freq.iterrows():
                allele_col = freq.columns[0]
                print(f"    {str(row[allele_col]):>8s}: {row['frequency']:.4f} ({row['count']})")

            print(f"\n  Genotype Distribution:")
            for _, row in geno.iterrows():
                print(f"    {str(row['genotype']):>12s}: {row['count']:4d} ({row['percentage']:.1f}%)")

            print(f"\n  Diversity Metrics:")
            print(f"    PIC:                     {pic:.4f}")
            print(f"    Observed Heterozygosity:  {het['observed_het']:.4f}")
            print(f"    Expected Heterozygosity:  {het['expected_het']:.4f}")
            print()

        except Exception as e:
            print(f"  Error analyzing {marker}: {e}\n")


def _cmd_rank(args):
    """Rank accessions for breeding."""
    from apple_mas.data_parser import DataParser
    from apple_mas.parent_scorer import ParentScorer

    parser = DataParser()
    scorer = ParentScorer()

    try:
        if args.format == "csv":
            df = parser.from_csv(args.input)
        else:
            df = parser.from_excel(args.input)
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(df)} records for {df['sample_id'].nunique()} samples\n")

    if args.pairs:
        print(f"Finding top {args.pairs} complementary parent pairs...\n")
        pairs = scorer.find_complementary_parents(df, n_pairs=args.pairs, profile=args.profile)
        print(pairs.to_string())
        if args.output:
            pairs.to_csv(args.output, index=True)
            print(f"\nSaved to {args.output}")
    else:
        ranked = scorer.rank_accessions(df, profile=args.profile, top_n=args.top)
        print(f"Top {args.top} accessions (profile: {args.profile}):\n")
        print(ranked.to_string())
        if args.output:
            if args.format == "excel":
                ranked.to_excel(args.output)
            else:
                ranked.to_csv(args.output)
            print(f"\nSaved to {args.output}")


def _cmd_visualize(args):
    """Generate visualization plots."""
    from apple_mas.data_parser import DataParser
    from apple_mas.visualization import MarkerVisualizer
    from apple_mas.allele_analysis import AlleleAnalyzer

    parser = DataParser()
    viz = MarkerVisualizer(dpi=args.dpi)
    analyzer = AlleleAnalyzer()

    try:
        df = parser.from_csv(args.input)
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)
    fmt = args.format
    ext = f".{fmt}"

    markers = [args.marker] if args.marker else sorted(df["marker"].unique())

    print(f"Generating plots for {len(markers)} marker(s)...\n")

    for marker in markers:
        if args.type in ("all", "frequency"):
            path = os.path.join(args.output, f"{marker}_allele_frequency{ext}")
            viz.allele_frequency_plot(df, marker, save_path=path)
            print(f"  Saved: {path}")

        if args.type in ("all", "genotype"):
            path = os.path.join(args.output, f"{marker}_genotype_distribution{ext}")
            viz.genotype_distribution(df, marker, save_path=path)
            print(f"  Saved: {path}")

    if args.type in ("all", "heatmap") and not args.marker:
        path = os.path.join(args.output, f"genotype_heatmap{ext}")
        viz.genotype_heatmap(df, markers, save_path=path)
        print(f"  Saved: {path}")

    if args.type in ("all", "composition") and not args.marker:
        path = os.path.join(args.output, f"sample_composition{ext}")
        viz.sample_composition_chart(df, save_path=path)
        print(f"  Saved: {path}")

    if args.type in ("all", "frequency") and not args.marker:
        path = os.path.join(args.output, f"multi_marker_frequency{ext}")
        viz.multi_marker_frequency(df, markers, save_path=path)
        print(f"  Saved: {path}")

    print(f"\nAll plots saved to {args.output}/")


def _cmd_sample_data(args):
    """Create sample dataset."""
    from apple_mas.data_parser import DataParser

    parser = DataParser()
    df, raw = parser.create_sample_dataset()
    os.makedirs(args.output, exist_ok=True)

    if args.format in ("csv", "both"):
        path = os.path.join(args.output, "sample_genotypes.csv")
        df.to_csv(path, index=False)
        print(f"  Saved: {path}")

    if args.format in ("excel", "both"):
        path = os.path.join(args.output, "sample_genotypes.xlsx")
        df.to_excel(path, index=False)
        print(f"  Saved: {path}")

    print(f"\nSample dataset created with {df['sample_id'].nunique()} accessions and {df['marker'].nunique()} markers")
    print(f"Location: {args.output}/")


def _cmd_report(args):
    """Generate a comprehensive analysis report."""
    from apple_mas.data_parser import DataParser
    from apple_mas.allele_analysis import AlleleAnalyzer
    from apple_mas.parent_scorer import ParentScorer
    from apple_mas.visualization import MarkerVisualizer

    parser = DataParser()
    analyzer = AlleleAnalyzer()
    scorer = ParentScorer()
    viz = MarkerVisualizer()

    try:
        df = parser.from_csv(args.input)
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)
    plots_dir = os.path.join(args.output, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print(f"Generating comprehensive report for {df['sample_id'].nunique()} samples...\n")

    # 1. Summary table
    summary = analyzer.all_markers_summary(df)
    summary_path = os.path.join(args.output, "marker_summary.csv")
    summary.to_csv(summary_path, index=False)
    print(f"  Saved: {summary_path}")

    # 2. Cross-marker summary
    cross = analyzer.cross_marker_summary(df)
    cross_path = os.path.join(args.output, "cross_marker_summary.csv")
    cross.to_csv(cross_path, index=False)
    print(f"  Saved: {cross_path}")

    # 3. Rankings
    ranked = scorer.rank_accessions(df, profile=args.profile)
    rank_path = os.path.join(args.output, "parent_rankings.csv")
    ranked.to_csv(rank_path, index=True)
    print(f"  Saved: {rank_path}")

    # 4. Complementary pairs
    pairs = scorer.find_complementary_parents(df, n_pairs=10, profile=args.profile)
    pairs_path = os.path.join(args.output, "complementary_pairs.csv")
    pairs.to_csv(pairs_path, index=True)
    print(f"  Saved: {pairs_path}")

    # 5. Plots
    markers = sorted(df["marker"].unique())
    for marker in markers:
        viz.allele_frequency_plot(df, marker, save_path=os.path.join(plots_dir, f"{marker}_frequency.png"))
        viz.genotype_distribution(df, marker, save_path=os.path.join(plots_dir, f"{marker}_genotype.png"))

    viz.multi_marker_frequency(df, markers, save_path=os.path.join(plots_dir, "all_markers_frequency.png"))
    viz.genotype_heatmap(df, markers, save_path=os.path.join(plots_dir, "genotype_heatmap.png"))
    viz.sample_composition_chart(df, save_path=os.path.join(plots_dir, "composition.png"))
    viz.ranking_barplot(ranked, save_path=os.path.join(plots_dir, "ranking.png"))
    print(f"  Plots saved to: {plots_dir}/")

    # 6. Report markdown
    md_lines = [
        "# Apple MAS Analysis Report\n",
        f"**Samples analyzed:** {df['sample_id'].nunique()}",
        f"**Markers analyzed:** {df['marker'].nunique()}",
        f"**Breeding profile:** {args.profile}\n",
        "## Marker Summary\n",
        summary.to_string(index=False),
        "\n## Cross-Marker Diversity\n",
        cross.to_string(index=False),
        "\n## Top 20 Parent Rankings\n",
        ranked.head(20).to_string(index=True),
        "\n## Top 10 Complementary Pairs\n",
        pairs.to_string(index=True),
    ]
    report_path = os.path.join(args.output, "report.md")
    with open(report_path, "w") as f:
        f.write("\n".join(md_lines))
    print(f"  Report saved to: {report_path}")

    print(f"\nReport complete! All files in {args.output}/")


def _print_marker_info(marker_name, info):
    """Pretty-print marker information."""
    print(f"\n{'='*60}")
    print(f"  Marker: {marker_name}")
    print(f"  Full Name: {info.get('full_name', 'N/A')}")
    print(f"  Trait: {info.get('trait', 'N/A')}")
    print(f"  Chromosome: {info.get('chromosome', 'N/A')}")
    print(f"  Type: {info.get('marker_type', 'N/A')}")
    print(f"{'='*60}\n")

    if "primer_sequences" in info:
        print("  Primer Sequences:")
        for direction, seq in info["primer_sequences"].items():
            print(f"    {direction}: {seq}")
        print()

    if "allele_definitions" in info:
        print("  Allele Definitions:")
        for allele, defn in info["allele_definitions"].items():
            print(f"    {allele}: {defn.get('phenotype', 'N/A')}")
        print()

    if "notes" in info:
        print(f"  Notes: {info['notes']}\n")

    if "breeding_notes" in info:
        print(f"  Breeding Notes: {info['breeding_notes']}\n")

    if "references" in info:
        print("  References:")
        for ref in info["references"]:
            print(f"    - {ref}")
        print()


def _print_pcr_methodology(methodology):
    """Pretty-print PCR methodology."""
    print("\n  Standard PCR Methodology")
    print("  " + "="*40)

    if "PCR_conditions" in methodology:
        pcr = methodology["PCR_conditions"]
        print(f"\n  Reaction Volume: {pcr.get('total_volume_ul', 'N/A')} uL")
        if "components" in pcr:
            print("  Components:")
            for k, v in pcr["components"].items():
                print(f"    {k}: {v}")
        if "thermal_cycling" in pcr:
            tc = pcr["thermal_cycling"]
            print(f"\n  Thermal Cycling:")
            print(f"    Initial denaturation: {tc.get('initial_denaturation', {}).get('temp_C')}C / {tc.get('initial_denaturation', {}).get('time_min')} min")
            print(f"    Cycles: {tc.get('cycles')}")
            print(f"    Denaturation: {tc.get('denaturation', {}).get('temp_C')}C / {tc.get('denaturation', {}).get('time_sec')} sec")
            print(f"    Annealing: {tc.get('annealing', {}).get('temp_C')}C / {tc.get('annealing', {}).get('time_sec')} sec")
            print(f"    Extension: {tc.get('extension', {}).get('temp_C')}C / {tc.get('extension', {}).get('time_sec')} sec")
            print(f"    Final extension: {tc.get('final_extension', {}).get('temp_C')}C / {tc.get('final_extension', {}).get('time_min')} min")

    if "fragment_analysis" in methodology:
        fa = methodology["fragment_analysis"]
        print(f"\n  Fragment Analysis:")
        print(f"    Platform: {fa.get('platform', 'N/A')}")
        print(f"    Software: {fa.get('software', 'N/A')}")
        print(f"    Size Standard: {fa.get('size_standard', 'N/A')}")

    print()


if __name__ == "__main__":
    main()
