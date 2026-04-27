# ─────────────────────────────────────────────
# Snakemake R wrapper for rendering Rmarkdown
# ─────────────────────────────────────────────

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

suppressPackageStartupMessages(library(rmarkdown))

rmd_file    <- "scripts/generate_report.Rmd"
out_file    <- snakemake@output$report
out_dir     <- dirname(out_file)
project_dir <- getwd()

message("Rendering report:")
message("  Rmd        : ", rmd_file)
message("  HTML       : ", out_file)
message("  Project dir: ", project_dir)

if (!file.exists(rmd_file)) {
    stop("Rmd file not found: ", rmd_file)
}

rmarkdown::render(
    input         = rmd_file,
    output_file   = basename(out_file),
    output_dir    = out_dir,
    knit_root_dir = project_dir,
    params = list(
        rf_distances      = file.path(project_dir, snakemake@input$rf_distances),
        topo_summary      = file.path(project_dir, snakemake@input$topological_summary),
        concordance       = file.path(project_dir, snakemake@input$concordance_analysis),
        model_summary     = file.path(project_dir, snakemake@input$model_summary),
        alignment_summary = file.path(project_dir, snakemake@input$align_summary_report),
        concat_tree       = file.path(project_dir, snakemake@input$concat_tree),
        sp_tree           = file.path(project_dir, snakemake@input$sp_concatenated_tree),
        coalescent_tree   = file.path(project_dir, snakemake@input$coalescent_tree),
        genes             = snakemake@params$genes
    ),
    quiet = FALSE
)

message("Report complete: ", out_file)

sink(type = "message")
sink()