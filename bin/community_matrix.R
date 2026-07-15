#!/usr/bin/env Rscript

# Written by Fernando Duarte

library(argparse)
library(dplyr)
library(tidyr)
library(vegan)

parser <- ArgumentParser(description = "Build community matrix of ASV presence/absence and abundance per barcode")
parser$add_argument("asv_table", help = "Path to joined ASV table (TSV/CSV)")
parser$add_argument("--tax", default = "species",
                    help = "Taxonomic level to include in the community matrix")
parser$add_argument("--prefix", default = "community_matrix",
                    help = "Prefix for output files (default: community_matrix)")
args <- parser$parse_args()

# Load ASV table
df <- read.csv(args$asv_table)

# Replace NA or empty taxonomic assignments with "Unclassified"
df[[args$tax]][is.na(df[[args$tax]]) | df[[args$tax]] == ""] <- "Unclassified"

# Keep only relevant columns for community matrix
df <- df[c("sample_name", "read_count", args$tax)]

# Pivot to wide format to create abundance community matrix
df_abundance <- df |>
    pivot_wider(names_from = all_of(args$tax), values_from = read_count, values_fill = 0, values_fn = sum)

# Pivot to wide format to create presence/absence community matrix
df_presence_absence <- df |>
    mutate(presence = ifelse(read_count > 0, 1, 0)) |>
    select(!read_count) |>
    pivot_wider(names_from = all_of(args$tax), values_from = presence, values_fill = 0, values_fn = max)

# Save output files
write.csv(df_abundance, paste0(args$prefix, "_abundance.csv"), row.names = FALSE)
write.csv(df_presence_absence, paste0(args$prefix, "_presence_absence.csv"), row.names = FALSE)
