"""Figure 3: CARD ARO mechanisms associated with filtered ABRicate gene-locus hits."""

from collections import Counter, defaultdict

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from common import (ROOT, RULE, configure, export_figure,
                    parser, read_csv, read_hits, verify_frozen_sources, write_csv)

ORDER = ["antibiotic efflux", "antibiotic inactivation", "antibiotic target alteration",
         "antibiotic target replacement", "antibiotic target protection",
         "reduced permeability to antibiotic"]
LABELS = {"antibiotic efflux": "Antibiotic efflux",
          "antibiotic inactivation": "Antibiotic inactivation",
          "antibiotic target alteration": "Antibiotic target alteration",
          "antibiotic target replacement": "Antibiotic target replacement",
          "antibiotic target protection": "Antibiotic target protection",
          "reduced permeability to antibiotic": "Reduced permeability\nto antibiotic"}
COLORS = {"antibiotic efflux": "#2563EB",
          "antibiotic inactivation": "#D97706",
          "antibiotic target alteration": "#059669",
          "antibiotic target replacement": "#059669",
          "antibiotic target protection": "#059669",
          "reduced permeability to antibiotic": "#DC2626"}


def main():
    args = parser(__doc__).parse_args()
    configure(args)
    verify_frozen_sources()
    rows, audit = read_hits(ROOT / "data/raw/card_abricate.tsv", "card", args.min_coverage, args.min_identity)
    mappings = defaultdict(list)
    for entry in read_csv(args.mapping_dir / "card_mechanism_mapping.csv"):
        mappings[entry["abricate_gene"]].append(entry)
    assignments = []
    for row in rows:
        gene = row["GENE"]
        if gene not in mappings:
            raise ValueError(f"No supported CARD mapping for {gene}")
        observed = set()
        for mapping in mappings[gene]:
            mechanism = mapping["mechanism"]
            if mechanism not in ORDER or mechanism in observed:
                raise ValueError(f"Unknown or duplicate mechanism for {gene}: {mechanism}")
            observed.add(mechanism)
            assignments.append({"gene": gene, "contig": row["SEQUENCE"], "start": row["START"], "end": row["END"],
                                "coverage_percent": row["%COVERAGE"], "identity_percent": row["%IDENTITY"],
                                **{k: mapping[k] for k in ["aro_accession", "mechanism", "mechanism_aro_accession",
                                                          "efflux_regulator", "ontology_path", "source_url"]}})
    counts = Counter(a["mechanism"] for a in assignments)
    mechanisms = [m for m in ORDER if counts[m]]
    denominator = len(rows)
    summary = [{"mechanism": m, "count_gene_locus_hits": counts[m], "denominator_hits": denominator,
                "percent_of_hits": f"{100 * counts[m] / denominator:.1f}"} for m in mechanisms]
    genes_multilabel = sorted({r["GENE"] for r in rows if len(mappings[r["GENE"]]) > 1})
    regulators = sorted({r["GENE"] for r in rows if any(m["efflux_regulator"] == "true" for m in mappings[r["GENE"]])})

    fig = plt.figure(figsize=(178 / 25.4, 90 / 25.4), dpi=120)
    # Plot labels/values only; explanatory text is retained in the separate caption file.
    ax = fig.add_axes([0.355, 0.190, 0.595, 0.770])
    values = [counts[m] for m in mechanisms]
    ax.barh(range(len(mechanisms)), values, height=0.55,
            color=[COLORS[m] for m in mechanisms], edgecolor="none", zorder=3)
    ax.set_yticks(range(len(mechanisms)), [LABELS[m] for m in mechanisms], fontsize=8.5)
    ax.invert_yaxis()
    maximum = max(values)
    ax.set_xlim(0, max(10, maximum * 1.29))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True, min_n_ticks=3))
    ax.tick_params(axis="y", length=0, pad=8)
    ax.tick_params(axis="x", length=3, color=RULE, labelsize=8.25)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(RULE)
    ax.set_xlabel("Number of gene-locus hits", fontsize=8.5, labelpad=8)
    for position, mechanism in enumerate(mechanisms):
        value = counts[mechanism]
        ax.text(value + maximum * 0.025, position,
                f"{value} ({100 * value / denominator:.1f}%)", ha="left", va="center", fontsize=8.5)
    write_csv(args.outdir / "tables/figure3_gene_mechanism_assignments.csv", assignments)
    write_csv(args.outdir / "tables/figure3_mechanism_counts.csv", summary)
    write_csv(args.outdir / "tables/figure3_filter_audit.csv", audit)
    export_figure(fig, "Figure_3_S11_card_mechanisms", args, {
        "description": "CARD ontology resistance mechanisms associated with ABRicate hits in E. coli S11",
        "input_files": ["data/raw/card_abricate.tsv", "data/mappings/card_mechanism_mapping.csv"],
        "retained_gene_locus_hits": denominator, "gene_mechanism_pairs": len(assignments),
        "mechanism_counts": dict(counts), "multimechanism_genes": genes_multilabel,
        "efflux_regulatory_genes": regulators, "efflux_regulatory_gene_count": len(regulators),
        "counting_unit": "One distinct gene-locus hit per associated CARD mechanism; categories may overlap",
        "interpretation": "Homology-associated mechanisms, including regulatory determinants; expression or resistance phenotype is not measured",
        "uncertainty": "Descriptive counts from one genome; no sampling interval",
        "presentation": {"plot_only": True, "title": False, "subtitle": False, "footnotes": False,
                         "colors_matched_to_original_manuscript": {m: COLORS[m] for m in mechanisms}},
    })
    print(f"Figure 3: {denominator} hits, {len(assignments)} gene-mechanism pairs; {dict(counts)}")


if __name__ == "__main__":
    main()
