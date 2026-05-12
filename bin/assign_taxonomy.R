
# Written by Fernando Duarte and Jordan Cuff

library(argparse)
library(dplyr)
library(tidyr)
library(stringr)
library(taxonomizr)


# Parse command-line arguments
parser <- ArgumentParser(description = 'Assign taxonomy to blast hits')
parser$add_argument('blast_hits', type = 'character', help = 'Path to blast hits in outfmt 6')
parser$add_argument('read_counts', type = 'character', help = 'Path to read counts per ASV')
parser$add_argument('--db_type', type = 'character', help = 'Taxonomy database type: taxonomizr or custom tsv')
parser$add_argument('--sql_db', type = 'character', help = 'Path to SQL ncbi taxonomy database')
parser$add_argument('--spident', type = 'numeric', help = 'Identity threshold (in %) for taxonomy assignment at species level', default = 99)
parser$add_argument('--gpident', type = 'numeric', help = 'Identity threshold (in %) for taxonomy assignment at genus level', default = 90)
parser$add_argument('--fpident', type = 'numeric', help = 'Identity threshold (in %) for taxonomy assignment at family level', default = 80)
parser$add_argument('--opident', type = 'numeric', help = 'Identity threshold (in %) for taxonomy assignment at oder level', default = 70)


args <- parser$parse_args()

# Load blast hits
Blastout <- read.table(args$blast_hits)

# Filter blast hits to retain only the best hit for each ASV based on bitscore, evalue, and percent identity
blast_filtered <- Blastout %>%
  rename(qseqid = V1, sseqid = V2, pident = V3, length = V4,
         mismatch = V5, gapopen = V6, qstart = V7, qend = V8,
         sstart = V9, send = V10, evalue = V11, bitscore = V12) %>%
  group_by(qseqid) %>%
  filter(bitscore == max(bitscore)) %>%
  filter(evalue   == min(evalue)) %>%
  filter(pident   == max(pident)) %>%
  ungroup() %>%
  dplyr::select(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

# Load read counts
read_counts <- read.csv(args$read_counts, header=FALSE)

# Column names for read counts - assuming the first column is ASV and the second column is read count
colnames(read_counts) <- c("ASV", "read_count")

# Debugging: Check the structure of the loaded data
head(Blastout)
head(read_counts)

# Create SQLite database - can be stored centrally to
# avoid replication across projects - this seems to variably work, so the below steps avoid it
# prepareDatabase('accessionTaxa.sql')

print(accessionToTaxa(data.frame(blast_filtered)$sseqid,args$sql_db))

sseqids <- data.frame(blast_filtered) %>%
  #dplyr::select(qseqid, pident, sseqid) %>%
  mutate(seqid2 = paste0(sseqid, ".1")) %>% # Prob should remove the seqid2 column
  mutate(taxaId = accessionToTaxa(sseqid, args$sql_db)) %>%
  mutate(Taxonomic.ranks = getTaxonomy(taxaId, args$sql_db)) %>%
  rename("ASV" = "qseqid")

taxonomic.df <- as.data.frame(sseqids$Taxonomic.ranks, stringsAsFactors = FALSE)
sseqids <- cbind(sseqids, taxonomic.df)
sseqids$Taxonomic.ranks <- NULL

write.csv(sseqids, "ASV_table_pre-assigned.csv", row.names = FALSE)

# Assign taxonomic ranks to ASVs based on percent identity thresholds and resolve conflicts
ASV.ids <- sseqids %>%
  # Assign taxonomic rank based on percent identity thresholds
  mutate(Taxon = if_else(pident > args$spident, "species", if_else(pident > args$gpident, "genus", if_else(pident > args$fpident, "family", if_else(pident > args$opident, "order", "class"))))) %>%
  # Set the order of taxonomic ranks for filtering
  mutate(Taxon = factor(Taxon, levels = c("species", "genus", "family", "order", "class"), ordered = TRUE)) %>%
  # Filter unnecessary columns
  dplyr::select(!c(sseqid,seqid2)) %>%
  group_by(ASV) %>%
  # Resolve taxonomic assignment. If the ASV has multiple hits,
  # assign the most specific taxon that is consistent across all hits.
  # If there are conflicting assignments at a given rank, move up to
  # the next rank until a consistent assignment is found or assign "Unassigned"
  # if no consistent assignment is found.
  mutate(Resolved.taxon = case_when(Taxon == "species" & n_distinct(species) == 1 ~ coalesce(first(species), "Unassigned"),
                                    Taxon %in% c("species", "genus") & n_distinct(genus) == 1 ~ coalesce(first(genus), "Unassigned"),
                                    Taxon %in% c("species", "genus", "family") & n_distinct(family) == 1 ~ coalesce(first(family), "Unassigned"),
                                    Taxon %in% c("species", "genus", "family", "order") & n_distinct(order) == 1 ~ coalesce(first(order), "Unassigned"),
                                    Taxon %in% c("species", "genus", "family", "order", "class") & n_distinct(class) == 1 ~ coalesce(first(class), "Unassigned"),
                                    TRUE ~ coalesce(phylum, "Unassigned"))) %>%
  filter(Taxon == min(Taxon)) %>%
  ungroup() %>% # Might be rdundant
  group_by(ASV) %>% # Might be redundant
  summarise(
    pident         = mean(pident), # Should be the same across all hits for a given ASV, but we take the mean just in case
    length         = mean(length),
    mismatch       = mean(mismatch),
    #gapopen        = mean(gapopen),
    #qstart         = mean(qstart),
    #qend           = mean(qend),
    #sstart         = mean(sstart),
    #send           = mean(send),
    evalue         = mean(evalue), # Should be the same across all hits for a given ASV, but we take the mean just in case
    bitscore       = mean(bitscore), # Should be the same across all hits for a given ASV, but we take the mean just in case
    taxaId         = getId(Resolved.taxon, args$sql_db),
    phylum         = if_else(n_distinct(phylum)  == 1, first(phylum),  NA_character_),
    class          = if_else(n_distinct(class)   == 1, first(class),   NA_character_),
    order          = if_else(n_distinct(order)   == 1, first(order),   NA_character_),
    family         = if_else(n_distinct(family)  == 1, first(family),  NA_character_),
    genus          = if_else(n_distinct(genus)   == 1, first(genus),   NA_character_),
    species        = if_else(n_distinct(species) == 1, first(species), NA_character_),
    Resolved.taxon = first(Resolved.taxon),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(sample_name = str_remove(ASV, "_\\d+_\\d+$")) %>%
  relocate(sample_name, .before = 1)

# Write the ASV table with assigned taxonomy to a CSV file
#write.csv(ASV.ids, "ASV_table_assigned.csv", row.names = FALSE)

# Merge reads counts bease on the ASV column
ASV.ids.read_counts <- merge(ASV.ids, read_counts, by = "ASV", all.x = TRUE) %>%
                        relocate(read_count, .before = 3)

# Write the final ASV table with assigned taxonomy and read counts to a CSV file
write.csv(ASV.ids.read_counts, "ASV_table_final.csv", row.names = FALSE)



