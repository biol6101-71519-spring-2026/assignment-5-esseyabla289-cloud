#!/usr/bin/env python3
"""
Calculate Robinson-Foulds (RF) distances between all tree pairs:
  - each gene tree vs. coalescent species tree
  - each gene tree vs. concatenated (supermatrix) tree
  - all pairwise gene tree comparisons
  - coalescent tree vs. concatenated tree

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
    """Load first Newick tree from *filepath*."""
    with open(filepath) as fh:
        content = fh.read().strip()
    first_tree = content.split("\n")[0]
    return Tree(first_tree)


def calc_rf(tree1, tree2, label1="t1", label2="t2"):
    """
    Return (rf, max_rf, normalized_rf) for two unrooted trees.
    Returns (None, None, None) on failure.
    """
    try:
        rf, max_rf, *_ = tree1.robinson_foulds(tree2, unrooted_trees=True)
        norm = rf / max_rf if max_rf > 0 else 0.0
        log.info("%s vs %s: RF=%d  max=%d  norm=%.4f", label1, label2, rf, max_rf, norm)
        return rf, max_rf, norm
    except Exception as exc:
        log.warning("RF failed for %s vs %s: %s", label1, label2, exc)
        return None, None, None


# ── Load inputs ───────────────────────────────────────────────────────────────

gene_tree_files = snakemake.input.gene_trees
coalescent_file = snakemake.input.coalescent_tree
concat_file     = snakemake.input.concat_tree
genes           = list(snakemake.params.genes)

log.info("Loading %d gene trees", len(genes))
gene_trees = {}
for gene, fp in zip(genes, gene_tree_files):
    gene_trees[gene] = load_tree(fp)

log.info("Loading coalescent species tree: %s", coalescent_file)
coalescent_tree = load_tree(coalescent_file)

log.info("Loading concatenated tree: %s", concat_file)
concat_tree = load_tree(concat_file)


# ── Calculate RF distances ────────────────────────────────────────────────────

rows = []  # list of dict rows for the output table

# Gene trees vs coalescent
for gene, tree in gene_trees.items():
    rf, max_rf, norm = calc_rf(tree, coalescent_tree, gene, "coalescent")
    if rf is not None:
        rows.append({"Tree1": gene, "Tree2": "coalescent",
                     "RF_Distance": rf, "Max_RF": max_rf, "Normalized_RF": f"{norm:.4f}"})

# Gene trees vs concatenated
for gene, tree in gene_trees.items():
    rf, max_rf, norm = calc_rf(tree, concat_tree, gene, "concatenated")
    if rf is not None:
        rows.append({"Tree1": gene, "Tree2": "concatenated",
                     "RF_Distance": rf, "Max_RF": max_rf, "Normalized_RF": f"{norm:.4f}"})

# Pairwise gene trees
gene_names = list(gene_trees.keys())
for i in range(len(gene_names)):
    for j in range(i + 1, len(gene_names)):
        g1, g2 = gene_names[i], gene_names[j]
        rf, max_rf, norm = calc_rf(gene_trees[g1], gene_trees[g2], g1, g2)
        if rf is not None:
            rows.append({"Tree1": g1, "Tree2": g2,
                         "RF_Distance": rf, "Max_RF": max_rf, "Normalized_RF": f"{norm:.4f}"})

# Coalescent vs concatenated
rf, max_rf, norm = calc_rf(coalescent_tree, concat_tree, "coalescent", "concatenated")
if rf is not None:
    rows.append({"Tree1": "coalescent", "Tree2": "concatenated",
                 "RF_Distance": rf, "Max_RF": max_rf, "Normalized_RF": f"{norm:.4f}"})


# ── Write RF distance table ───────────────────────────────────────────────────

with open(snakemake.output.rf_distances, "w") as fh:
    fh.write("Tree1\tTree2\tRF_Distance\tMax_RF\tNormalized_RF\n")
    for r in rows:
        fh.write(f"{r['Tree1']}\t{r['Tree2']}\t{r['RF_Distance']}\t"
                 f"{r['Max_RF']}\t{r['Normalized_RF']}\n")

log.info("Wrote RF distances: %s", snakemake.output.rf_distances)


# ── Write topological summary ─────────────────────────────────────────────────

gene_vs_coal   = [r for r in rows if r["Tree2"] == "coalescent"   and r["Tree1"] in genes]
gene_vs_concat = [r for r in rows if r["Tree2"] == "concatenated" and r["Tree1"] in genes]
coal_vs_concat = [r for r in rows if r["Tree1"] == "coalescent"   and r["Tree2"] == "concatenated"]


def avg_norm(row_list):
    vals = [float(r["Normalized_RF"]) for r in row_list]
    return sum(vals) / len(vals) if vals else float("nan")


with open(snakemake.output.topo_summary, "w") as fh:
    def w(line=""):
        fh.write(line + "\n")

    w("Topological Summary")
    w("=" * 60)
    w()
    w(f"Genes analysed : {len(gene_trees)}")
    w(f"Gene names     : {', '.join(genes)}")
    w()

    w("Gene Trees vs. Coalescent Species Tree")
    w("-" * 40)
    for r in gene_vs_coal:
        w(f"  {r['Tree1']:20s}  RF={r['RF_Distance']:4}  Norm={r['Normalized_RF']}")
    if gene_vs_coal:
        w(f"  Mean Normalized RF : {avg_norm(gene_vs_coal):.4f}")
    w()

    w("Gene Trees vs. Concatenated Tree")
    w("-" * 40)
    for r in gene_vs_concat:
        w(f"  {r['Tree1']:20s}  RF={r['RF_Distance']:4}  Norm={r['Normalized_RF']}")
    if gene_vs_concat:
        w(f"  Mean Normalized RF : {avg_norm(gene_vs_concat):.4f}")
    w()

    w("Coalescent Tree vs. Concatenated Tree")
    w("-" * 40)
    for r in coal_vs_concat:
        w(f"  RF={r['RF_Distance']}  Max={r['Max_RF']}  Normalized={r['Normalized_RF']}")
    w()

    w("Interpretation")
    w("-" * 40)
    w("  Normalized RF = 0.0  →  identical topologies")
    w("  Normalized RF = 1.0  →  maximally different")
    w("  Normalized RF < 0.3  →  largely concordant")
    w("  Normalized RF ≥ 0.3  →  significant topological conflict")

log.info("Wrote topological summary: %s", snakemake.output.topo_summary)
log.info("Tree comparison complete.")
