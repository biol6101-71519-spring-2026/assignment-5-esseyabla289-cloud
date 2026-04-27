#!/usr/bin/env python3
"""
Calculate Gene Concordance Factors (gCF) for each internal node of the
species tree using individual gene trees.

For each bipartition in the species tree:
  concordant    - gene trees that support the bipartition
  conflicting   - gene trees that contradict the bipartition
  uninformative - gene trees that lack taxa needed to evaluate it

  gCF = concordant / (concordant + conflicting)

Called by Snakemake; receives input/output/params via the snakemake object.
"""

import logging
from ete3 import Tree

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger()


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_tree(filepath):
    with open(filepath) as fh:
        content = fh.read().strip()
    return Tree(content.split("\n")[0])


def get_bipartitions(tree, all_taxa):
    """
    Return the set of unrooted bipartitions for *tree* as frozensets of
    frozenset pairs, restricted to *all_taxa*.
    """
    bps = set()
    tree = tree.copy("newick")
    tree.unroot()
    for node in tree.traverse():
        if node.is_leaf() or node.is_root():
            continue
        side = frozenset(node.get_leaf_names()) & all_taxa
        other = all_taxa - side
        if len(side) >= 1 and len(other) >= 1:
            bps.add(frozenset([side, other]))
    return bps


def evaluate_bipart(target_bp, gene_tree, all_taxa):
    """
    Return 'concordant', 'conflicting', or 'uninformative' for *target_bp*
    in the context of *gene_tree*.
    """
    side1, side2 = list(target_bp)
    gene_taxa = frozenset(gene_tree.get_leaf_names())

    s1_in = side1 & gene_taxa
    s2_in = side2 & gene_taxa

    # Need at least 2 taxa on each side to be informative
    if len(s1_in) < 2 or len(s2_in) < 2:
        return "uninformative"

    # Get bipartitions restricted to species tree taxa
    gene_bps = get_bipartitions(gene_tree, all_taxa)

    # Restricted target bipartition (only taxa present in gene tree)
    restricted = frozenset([s1_in, s2_in])

    if restricted in gene_bps:
        return "concordant"

    # Check for conflict: any bipartition that splits s1_in and s2_in
    for bp in gene_bps:
        bp_sides = list(bp)
        for bs in bp_sides:
            if bs & s1_in and bs & s2_in:
                return "conflicting"

    return "uninformative"


# ── Load inputs ───────────────────────────────────────────────────────────────

gene_tree_files = snakemake.input.gene_trees
species_tree_file = snakemake.input.species_tree
genes = list(snakemake.params.genes)

log.info("Loading species tree: %s", species_tree_file)
species_tree = load_tree(species_tree_file)
all_taxa = frozenset(species_tree.get_leaf_names())

log.info("Loading %d gene trees", len(genes))
gene_trees = []
for gene, fp in zip(genes, gene_tree_files):
    gene_trees.append(load_tree(fp))
    log.info("  Loaded %s from %s", gene, fp)


# ── Calculate gCF ────────────────────────────────────────────────────────────

results = []

ref = species_tree.copy("newick")
ref.unroot()

for node in ref.traverse():
    if node.is_leaf() or node.is_root():
        continue

    node_taxa    = frozenset(node.get_leaf_names())
    complement   = all_taxa - node_taxa
    target_bp    = frozenset([node_taxa, complement])

    concordant_genes  = []
    conflicting_genes = []
    uninformative_ct  = 0

    for gene, gtree in zip(genes, gene_trees):
        verdict = evaluate_bipart(target_bp, gtree, all_taxa)
        if verdict == "concordant":
            concordant_genes.append(gene)
        elif verdict == "conflicting":
            conflicting_genes.append(gene)
        else:
            uninformative_ct += 1

    informative = len(concordant_genes) + len(conflicting_genes)
    gcf = len(concordant_genes) / informative if informative > 0 else float("nan")

    results.append({
        "clade":             sorted(node_taxa),
        "concordant":        len(concordant_genes),
        "conflicting":       len(conflicting_genes),
        "uninformative":     uninformative_ct,
        "gCF":               gcf,
        "concordant_genes":  concordant_genes,
        "conflicting_genes": conflicting_genes,
    })

    log.info("Clade %s: gCF=%.3f (%d conc / %d conf / %d uninf)",
             "+".join(sorted(node_taxa)[:3]), gcf,
             len(concordant_genes), len(conflicting_genes), uninformative_ct)


# ── Write output ──────────────────────────────────────────────────────────────

with open(snakemake.output.concordance_analysis, "w") as fh:
    def w(line=""):
        fh.write(line + "\n")

    w("Gene Concordance Factor (gCF) Analysis")
    w("=" * 65)
    w()
    w(f"Reference tree  : {species_tree_file}")
    w(f"Gene trees      : {len(gene_trees)}")
    w(f"Genes           : {', '.join(genes)}")
    w(f"Total taxa      : {len(all_taxa)}")
    w()

    # Table header
    w(f"{'Clade (representative taxa)':<38} {'gCF':>6} {'Conc':>6} {'Conf':>6} {'Uninf':>6}")
    w("-" * 65)
    for r in results:
        label = ", ".join(r["clade"][:3])
        if len(r["clade"]) > 3:
            label += f" +{len(r['clade']) - 3}"
        gcf_str = f"{r['gCF']:.3f}" if r["gCF"] == r["gCF"] else "NA"
        w(f"{label:<38} {gcf_str:>6} {r['concordant']:>6} {r['conflicting']:>6} {r['uninformative']:>6}")

    w()
    w("Detailed Analysis")
    w("-" * 65)
    for r in results:
        w()
        w(f"Clade : {', '.join(r['clade'])}")
        gcf_str = f"{r['gCF']:.3f}" if r["gCF"] == r["gCF"] else "NA"
        w(f"  gCF               : {gcf_str}")
        w(f"  Concordant ({r['concordant']})   : {', '.join(r['concordant_genes']) or 'none'}")
        w(f"  Conflicting ({r['conflicting']})  : {', '.join(r['conflicting_genes']) or 'none'}")
        w(f"  Uninformative ({r['uninformative']}): —")

    w()
    w("Interpretation")
    w("-" * 65)
    w("  gCF = 1.0  → All informative gene trees support this bipartition")
    w("  gCF = 0.5  → Equal concordance and conflict")
    w("  gCF < 0.5  → More conflicting than concordant gene trees")
    w("  Low gCF values suggest Incomplete Lineage Sorting (ILS) or")
    w("  other sources of gene tree heterogeneity.")

log.info("Concordance analysis complete: %s", snakemake.output.concordance_analysis)
