# visualize_trees.R
# ─────────────────────────────────────────────────────────────────────────────
# Unified tree visualization script called by multiple Snakemake rules.
# Behaviour is controlled by snakemake@params$mode:
#   "gene_tree"  - individual gene tree with bootstrap labels
#   "species_tree" - coalescent or concatenated species tree
#   "comparison" - side-by-side coalescent vs. supermatrix comparison
# ─────────────────────────────────────────────────────────────────────────────

# ── Ensure required packages are installed ─────────────────────────────
required_cran <- c("ape", "ggplot2", "cowplot")
required_bioc <- c("ggtree", "treeio")
required_cran_extra <- c("phytools", "phangorn")

# Install CRAN packages if missing
for (pkg in required_cran) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

# Install Bioconductor manager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install Bioconductor packages
for (pkg in required_bioc) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
}

# Install extra CRAN phylo packages
for (pkg in required_cran_extra) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

# Redirect all output to the Snakemake log file
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

suppressPackageStartupMessages({
    library(ape)
    library(ggtree)
    library(treeio)
    library(ggplot2)
    library(cowplot)
})

mode  <- snakemake@params$mode
title <- snakemake@params$title

message("Mode: ", mode)
message("Title: ", title)


# ── Helper: safe tree loader ──────────────────────────────────────────────────
load_tree <- function(path) {
    message("Loading tree: ", path)
    tryCatch(
        read.tree(path),
        error = function(e) {
            message("WARNING: could not parse tree — ", e$message)
            NULL
        }
    )
}


# ── Helper: standard ggtree panel ────────────────────────────────────────────
tree_panel <- function(tree, panel_title, show_node_labels = TRUE,
                       tip_size = 3, node_size = 2.5) {
    if (is.null(tree)) {
        return(ggplot() + ggtitle(panel_title) +
               annotate("text", x = 0.5, y = 0.5, label = "Tree not available",
                        size = 5) + theme_void())
    }
    p <- ggtree(tree, branch.length = "branch.length") +
        geom_tiplab(size = tip_size, fontface = "italic") +
        geom_treescale(fontsize = 3) +
        ggtitle(panel_title) +
        theme_tree2(legend.position = "none")

    if (show_node_labels) {
        p <- p + geom_nodelab(aes(label = label),
                              size = node_size, hjust = 1.2, vjust = -0.4,
                              color = "steelblue")
    }
    p
}


# ═════════════════════════════════════════════════════════════════════════════
# MODE: gene_tree
# ═════════════════════════════════════════════════════════════════════════════
if (mode == "gene_tree") {

    tree <- load_tree(snakemake@input$tree)

    p <- tree_panel(tree, title, show_node_labels = TRUE, tip_size = 3.5)

    ggsave(snakemake@output$pdf, p,
           width = 10, height = 8, device = "pdf")
    message("Saved: ", snakemake@output$pdf)


# ═════════════════════════════════════════════════════════════════════════════
# MODE: species_tree  (coalescent OR concatenated)
# ═════════════════════════════════════════════════════════════════════════════
} else if (mode == "species_tree") {

    tree   <- load_tree(snakemake@input$tree)
    method <- snakemake@params$method   # "coalescent" | "concatenated"

    node_color <- if (method == "coalescent") "darkorange" else "steelblue"

    if (!is.null(tree)) {
        p <- ggtree(tree) +
            geom_tiplab(size = 3.8, fontface = "italic") +
            geom_nodelab(aes(label = label),
                         size = 3, hjust = 1.2, vjust = -0.5,
                         color = node_color) +
            geom_treescale(fontsize = 3) +
            ggtitle(title) +
            theme_tree2()
    } else {
        p <- ggplot() + ggtitle(title) +
             annotate("text", x = 0.5, y = 0.5, label = "Tree not available") +
             theme_void()
    }

    ggsave(snakemake@output$pdf, p,
           width = 12, height = 9, device = "pdf")
    message("Saved: ", snakemake@output$pdf)


# ═════════════════════════════════════════════════════════════════════════════
# MODE: comparison  (side-by-side panels + RF annotation)
# ═════════════════════════════════════════════════════════════════════════════
} else if (mode == "comparison") {

    coal <- load_tree(snakemake@input$coalescent)
    conc <- load_tree(snakemake@input$concatenated)

    # Read RF distances
    rf_table <- tryCatch(
        read.table(snakemake@input$rf_distances,
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE),
        error = function(e) { message("WARNING: RF table not readable"); NULL }
    )

    # Extract the coalescent-vs-concatenated row
    rf_label <- "RF distance: N/A"
    if (!is.null(rf_table)) {
        row <- rf_table[rf_table$Tree1 == "coalescent" &
                        rf_table$Tree2 == "concatenated", ]
        if (nrow(row) > 0) {
            rf_label <- sprintf("Normalized RF distance = %s", row$Normalized_RF[1])
        }
    }

    p1 <- tree_panel(coal, "Coalescent Species Tree (ASTRAL)",
                     tip_size = 3, node_size = 2.5)
    p2 <- tree_panel(conc,
                     paste0("Supermatrix Tree (IQ-TREE)\n", rf_label),
                     tip_size = 3, node_size = 2.5)

    combined <- plot_grid(p1, p2, ncol = 2, labels = c("A", "B"),
                          label_size = 14)

    ggsave(snakemake@output$comparison, combined,
           width = 20, height = 10, device = "pdf")
    message("Saved: ", snakemake@output$comparison)


} else {
    stop("Unknown mode: ", mode,
         "  (expected 'gene_tree', 'species_tree', or 'comparison')")
}

sink(type = "message")
sink()
