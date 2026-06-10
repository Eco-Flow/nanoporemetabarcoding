#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)
library(dplyr)
library(tidyr)

parser <- ArgumentParser(description = "Stacked bar chart of ASV counts per phylum per barcode")
parser$add_argument("asv_table",       help = "Path to joined ASV table (TSV/CSV)")
parser$add_argument("--min_fraction",  type = "double", default = 0.01,
                    help = "Phyla below this fraction of total ASVs are merged into 'Other'")
parser$add_argument("--width",         type = "double", default = 10, help = "Plot width in inches")
parser$add_argument("--height",        type = "double", default = 6,  help = "Plot height in inches")
parser$add_argument("--tax",           nargs = "+", default = "phylum",
                    help = "Taxonomic level(s) to plot")
args <- parser$parse_args()

df   <- read.csv(args$asv_table, stringsAsFactors = FALSE)
taxa <- unlist(strsplit(args$tax, ","))

for (tax in taxa) {

    print(paste("Processing taxonomic level:", tax))

    df[[tax]][is.na(df[[tax]]) | df[[tax]] == ""] <- "Unclassified"

    print(paste("Calculating ASV counts per barcode and", tax))

    counts <- df %>%
        group_by(barcode, .data[[tax]]) %>%
        summarise(n_asvs = n(), .groups = "drop")

    print(paste("Calculating total ASVs per", tax, "to identify rare taxa"))

    totals <- counts %>%
        group_by(.data[[tax]]) %>%
        summarise(total = sum(n_asvs))

    print(paste("Identifying rare taxa below", args$min_fraction * 100, "% of total ASVs"))

    rare <- totals[[tax]][totals$total / sum(totals$total) < args$min_fraction &
                          totals[[tax]] != "Unclassified"]

    print(paste("Merging", length(rare), "rare taxa into 'Other'"))

    if (length(rare) > 0) {
        counts[[tax]][counts[[tax]] %in% rare] <- "Other"
        counts <- counts %>%
            group_by(barcode, .data[[tax]]) %>%
            summarise(n_asvs = sum(n_asvs), .groups = "drop")
    }

    counts <- counts %>%
        group_by(barcode) %>%
        mutate(relative_abundance = n_asvs / sum(n_asvs)) %>%
        ungroup()

    level_order <- c(
        sort(setdiff(unique(counts[[tax]]), c("Unclassified", "Other"))),
        intersect(c("Unclassified", "Other"), unique(counts[[tax]]))
    )
    counts[[tax]] <- factor(counts[[tax]], levels = level_order)

    p <- ggplot(counts, aes(x = barcode, y = relative_abundance, fill = .data[[tax]])) +
        geom_col() +
        scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
        labs(x = "Barcode", y = "Relative Abundance (%)", fill = toupper(tax),
             title = paste("Relative Abundance per", toupper(tax), "per barcode")) +
        theme_classic() +
        theme(
            axis.text.x    = element_text(angle = 45, hjust = 1),
            legend.position = "right"
        )

    ggsave(paste0(tax, "_asv_per_barcode.png"), plot = p, width = args$width, height = args$height, dpi = 150)

}
