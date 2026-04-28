# visualize_trees.R
# ─────────────────────────────────────────────────────────────────────────────
# Unified tree visualization script called by multiple Snakemake rules.
# Behaviour is controlled by snakemake@params$mode:
#   "gene_tree"  - individual gene tree with bootstrap labels
#   "species_tree" - coalescent or concatenated species tree
#   "comparison" - side-by-side coalescent vs. supermatrix comparison
#                  (also validates that gene counts agree between methods)
# ─────────────────────────────────────────────────────────────────────────────

# ── Ensure required packages are installed ─────────────────────────────
required_cran <- c("ape", "ggplot2", "cowplot")
required_bioc <- c("ggtree", "treeio")
required_cran_extra <- c("phytools", "phangorn")

for (pkg in required_cran) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in required_bioc) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
}

for (pkg in required_cran_extra) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

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


# ── Build accession → species map ────────────────────────────────────────────
# After the rename_fasta_headers pipeline step, tree tips already use species
# names with underscores (e.g. "Homo_sapiens").  This function handles both:
#   (a) new-style tips: replaces underscores with spaces for display
#   (b) legacy accession tips: falls back to matching against FASTA headers
build_taxa_map <- function(fasta_dir) {
    taxa_map <- list()
    fasta_files <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)
    for (fpath in fasta_files) {
        lines <- readLines(fpath)
        headers <- lines[startsWith(lines, ">")]
        for (h in headers) {
            parts  <- strsplit(sub("^>", "", h), "\\s+")[[1]]
            acc    <- parts[1]
            if (length(parts) >= 3) {
                species <- paste(parts[2], parts[3])
            } else if (length(parts) == 2) {
                species <- parts[2]
            } else {
                species <- acc
            }
            if (is.null(taxa_map[[acc]])) {
                taxa_map[[acc]] <- species
            }
        }
    }
    message("Built taxa map with ", length(taxa_map), " entries from ", length(fasta_files), " FASTA files")
    taxa_map
}

# ── Rename tree tips ──────────────────────────────────────────────────────────
rename_tips <- function(tree, taxa_map) {
    if (is.null(tree)) return(NULL)
    tree$tip.label <- sapply(tree$tip.label, function(label) {
        # New-style tip: Genus_species(_subspecies) — just replace underscores
        if (grepl("^[A-Z][a-z]+_[a-z]", label)) {
            return(gsub("_", " ", label))
        }
        # Legacy accession tip: look up in taxa map
        sp <- taxa_map[[label]]
        if (!is.null(sp)) sp else label
    })
    tree
}


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


# ── Build taxa map (needed by all modes) ─────────────────────────────────────
fasta_dir <- snakemake@params$fasta_dir
taxa_map  <- build_taxa_map(fasta_dir)


# ═════════════════════════════════════════════════════════════════════════════
# MODE: gene_tree
# ═════════════════════════════════════════════════════════════════════════════
if (mode == "gene_tree") {

    tree <- load_tree(snakemake@input$tree)
    tree <- rename_tips(tree, taxa_map)

    p <- tree_panel(tree, title, show_node_labels = TRUE, tip_size = 3.5)

    ggsave(snakemake@output$pdf, p,
           width = 10, height = 8, device = "pdf")
    message("Saved: ", snakemake@output$pdf)


# ═════════════════════════════════════════════════════════════════════════════
# MODE: species_tree  (coalescent OR concatenated)
# ═════════════════════════════════════════════════════════════════════════════
} else if (mode == "species_tree") {

    tree   <- load_tree(snakemake@input$tree)
    tree   <- rename_tips(tree, taxa_map)
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
# MODE: comparison  (side-by-side panels + RF annotation + gene count check)
# ═════════════════════════════════════════════════════════════════════════════
} else if (mode == "comparison") {

    # ── Gene count validation ─────────────────────────────────────────────────
    # Count gene trees used by ASTRAL
    astral_input_path <- snakemake@input$astral_input
    astral_lines <- readLines(astral_input_path)
    # Each non-empty line is one gene tree
    n_coalescent_genes <- sum(nzchar(trimws(astral_lines)))

    # Count partitions used in concatenated alignment
    partition_path <- snakemake@input$partitions
    part_lines <- readLines(partition_path)
    n_concat_genes <- sum(nzchar(trimws(part_lines)))

    message("Gene count check:")
    message("  ASTRAL gene trees : ", n_coalescent_genes)
    message("  Concatenated partitions : ", n_concat_genes)

    if (n_coalescent_genes != n_concat_genes) {
        warning(sprintf(
            "Gene count MISMATCH: ASTRAL used %d gene tree(s) but concatenated alignment has %d partition(s). ",
            n_coalescent_genes, n_concat_genes
        ))
    } else {
        message("  Gene counts match: both methods used ", n_coalescent_genes, " gene(s). OK")
    }

    gene_count_label <- sprintf("Genes used — coalescent: %d | concatenated: %d%s",
                                n_coalescent_genes, n_concat_genes,
                                if (n_coalescent_genes != n_concat_genes) "  *** MISMATCH ***" else "  (match)")

    # ── Load and rename trees ─────────────────────────────────────────────────
    coal <- load_tree(snakemake@input$coalescent)
    conc <- load_tree(snakemake@input$concatenated)
    coal <- rename_tips(coal, taxa_map)
    conc <- rename_tips(conc, taxa_map)

    # Read RF distances
    rf_table <- tryCatch(
        read.table(snakemake@input$rf_distances,
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE),
        error = function(e) { message("WARNING: RF table not readable"); NULL }
    )

    rf_label <- "RF distance: N/A"
    if (!is.null(rf_table)) {
        row <- rf_table[rf_table$Tree1 == "coalescent" &
                        rf_table$Tree2 == "concatenated", ]
        if (nrow(row) > 0) {
            rf_label <- sprintf("Normalized RF distance = %s", row$Normalized_RF[1])
        }
    }

    p1 <- tree_panel(coal,
                     paste0("Coalescent Species Tree (ASTRAL)\n", gene_count_label),
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
