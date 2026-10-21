"""Figure 2: data-derived ResFinder gene-class Sankey, at >=300 dpi."""

from collections import Counter, defaultdict

import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle

from common import (ROOT, configure, export_figure,
                    parser, read_csv, read_hits, verify_frozen_sources, write_csv)

CLASS_ORDER = ["Beta-lactam", "Aminoglycoside", "Quinolone"]
DISPLAY = {"Beta-lactam": "Beta-lactams", "Aminoglycoside": "Aminoglycosides", "Quinolone": "Quinolones"}
COLORS = {"Beta-lactam": "#3B82F6", "Aminoglycoside": "#F97316", "Quinolone": "#10B981"}
GENE_NODE_COLOR = "#334155"


def derive_associations(rows, mapping):
    drug_classes = {r["resfinder_resistance_term"]: r["antimicrobial_class"] for r in mapping}
    evidence = defaultdict(set)
    loci = defaultdict(set)
    for row in rows:
        gene = row["PRODUCT"].strip()
        if not gene:
            raise ValueError("A ResFinder hit has no product/display gene")
        for drug in row["RESISTANCE"].split(";"):
            if drug not in drug_classes:
                raise ValueError(f"Unmapped resistance term: {drug}")
            key = (gene, drug_classes[drug])
            evidence[key].add(drug)
            loci[key].add(f"{row['SEQUENCE']}:{row['START']}-{row['END']}")
    unknown = {g for _, g in evidence} - set(CLASS_ORDER)
    if unknown:
        raise ValueError(f"Add explicit class ordering and colors for {unknown}")
    return [{"gene": gene, "antimicrobial_class": group, "association_weight": 1,
             "supporting_resfinder_terms": ";".join(sorted(drugs)),
             "supporting_loci": ";".join(sorted(loci[(gene, group)]))}
            for (gene, group), drugs in sorted(evidence.items(), key=lambda p: (CLASS_ORDER.index(p[0][1]), p[0][0]))]


def ribbon(ax, x0, x1, source_bottom, target_bottom, thickness, color):
    control = (x1 - x0) * 0.42
    vertices = [(x0, source_bottom), (x0 + control, source_bottom),
                (x1 - control, target_bottom), (x1, target_bottom),
                (x1, target_bottom + thickness), (x1 - control, target_bottom + thickness),
                (x0 + control, source_bottom + thickness), (x0, source_bottom + thickness),
                (x0, source_bottom)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(vertices, codes), facecolor=color, alpha=0.42,
                           edgecolor="none", linewidth=0, zorder=1))


def main():
    args = parser(__doc__).parse_args()
    configure(args)
    verify_frozen_sources()
    raw = ROOT / "data/raw/resfinder_abricate.tsv"
    rows, audit = read_hits(raw, "resfinder", args.min_coverage, args.min_identity)
    associations = derive_associations(rows, read_csv(args.mapping_dir / "resfinder_drug_class_mapping.csv"))
    class_counts = Counter(row["antimicrobial_class"] for row in associations)
    gene_classes = defaultdict(list)
    for row in associations:
        gene_classes[row["gene"]].append(row["antimicrobial_class"])
    genes = sorted(gene_classes, key=lambda g: (min(CLASS_ORDER.index(c) for c in gene_classes[g]), g))
    classes = [c for c in CLASS_ORDER if class_counts[c]]
    if len(genes) > 10:
        raise ValueError("More than ten genes need a revised manuscript layout; do not shrink labels")

    fig = plt.figure(figsize=(178 / 25.4, 100 / 25.4), dpi=120)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    # Only node labels and ribbons are drawn. Captions and methods stay outside the image.
    unit, gap, middle = 0.090, 0.035, 0.500
    left, right, node_width = 0.200, 0.760, 0.011
    total = len(associations) * unit + (len(genes) - 1) * gap
    cursor = middle + total / 2
    gene_top = {}
    for gene in genes:
        height = len(gene_classes[gene]) * unit
        gene_top[gene] = cursor
        ax.add_patch(Rectangle((left, cursor - height), node_width, height, facecolor=GENE_NODE_COLOR, edgecolor="none", zorder=3))
        ax.text(left - 0.014, cursor - height / 2, gene, ha="right", va="center", fontsize=8.5)
        cursor -= height + gap
    total = len(associations) * unit + (len(classes) - 1) * gap
    cursor = middle + total / 2
    class_top = {}
    for group in classes:
        height = class_counts[group] * unit
        class_top[group] = cursor
        ax.add_patch(Rectangle((right, cursor - height), node_width, height, facecolor=COLORS[group], edgecolor="none", zorder=3))
        ax.text(right + node_width + 0.014, cursor - height / 2,
                f"{DISPLAY[group]} ({class_counts[group]})", ha="left", va="center", fontsize=8.5, fontweight="semibold")
        cursor -= height + gap
    gene_used, class_used = Counter(), Counter()
    for gene in genes:
        for group in sorted(gene_classes[gene], key=CLASS_ORDER.index):
            source = gene_top[gene] - (gene_used[gene] + 1) * unit
            target = class_top[group] - (class_used[group] + 1) * unit
            ribbon(ax, left + node_width, right, source, target, unit, COLORS[group])
            gene_used[gene] += 1
            class_used[group] += 1
    assert sum(gene_used.values()) == sum(class_used.values()) == len(associations)
    write_csv(args.outdir / "tables/figure2_gene_class_associations.csv", associations)
    write_csv(args.outdir / "tables/figure2_filter_audit.csv", audit)
    export_figure(fig, "Figure_2_S11_resfinder_sankey", args, {
        "description": "Acquired resistance genes and antimicrobial class associations in E. coli S11",
        "input_files": ["data/raw/resfinder_abricate.tsv", "data/mappings/resfinder_drug_class_mapping.csv"],
        "unique_genes": len(genes), "gene_class_pairs": len(associations), "class_counts": dict(class_counts),
        "association_unit": "One unique gene-class pair; multiple drug terms within a class count once",
        "uncertainty": "Descriptive counts from one genome; no sampling interval",
        "sankey_flow_conserved": True,
        "presentation": {"plot_only": True, "title": False, "subtitle": False, "footnotes": False,
                         "colors_matched_to_original_manuscript": COLORS, "gene_node_color": GENE_NODE_COLOR,
                         "ribbon_alpha": 0.42},
    })
    print(f"Figure 2: {len(genes)} genes, {len(associations)} associations; {dict(class_counts)}")


if __name__ == "__main__":
    main()
