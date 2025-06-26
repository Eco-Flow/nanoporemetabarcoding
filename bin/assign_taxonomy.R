
library(argparse)
library(dplyr)
library(tidyr)
library(stringr)
library(taxonomizr)


# Parse command-line arguments
parser <- ArgumentParser(description = 'Assign taxonomy to blast hits')
parser$add_argument('blast_hits', type = 'character', help = 'Path to blast hits in outfmt 6')
parser$add_argument('--db_type', type = 'character', help = 'Taxonomy database type: taxonomizr or custom tsv')
parser$add_argument('--sql_db', type = 'character', help = 'Path to SQL ncbi taxonomy database')

args <- parser$parse_args()

Blastout <- read.table(args$blast_hits)

# Load hits table
blast_filtered <- Blastout %>%
  mutate(qseqid = V1,
         sseqid = V2,
         pident = V3,
         length = V4,
         mismatch = V5,
         gapopen = V6,
         qstart = V7,
         qend = V8,
         sstart = V9,
         send = V10,
         evalue = V11,
         bitscore = V12) %>%
  dplyr::select(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

# Create SQLite database - can be stored centrally to avoid replication across projects - this seems to variably work, so the below steps avoid it
#prepareDatabase('accessionTaxa.sql')

print(data.frame(blast_filtered)$sseqid)

print(accessionToTaxa(data.frame(blast_filtered)$sseqid,args$sql_db))

sseqids <- data.frame(blast_filtered) %>%
  dplyr::select(qseqid, pident, sseqid) %>%
  mutate(seqid2 = paste0(sseqid, ".1")) %>%
  mutate(taxaId = accessionToTaxa(sseqid, args$sql_db)) %>%
  mutate(Taxonomic.ranks = getTaxonomy(taxaId, args$sql_db)) %>%
  rename("ASV" = "qseqid")

taxonomic.df <- as.data.frame(sseqids$Taxonomic.ranks, stringsAsFactors = FALSE)
sseqids <- cbind(sseqids, taxonomic.df)
sseqids$Taxonomic.ranks <- NULL

write.csv(sseqids, "ASV_taxa.csv")

print(sseqids)

ASV.ids <- sseqids %>%
  mutate(Taxon = if_else(order == "Araneae",
                         if_else(pident > 90, genus,
                                         if_else(pident > 80, family,
                                                 if_else(pident > 70, order,
                                                         phylum))),
         if_else(pident > 80, family,
                 if_else(pident > 70, order,
                         phylum)))) %>%
  dplyr::select(ASV, Taxon) %>%
  mutate(ASV = str_remove(ASV, "_"))

print("debug")

print(ASV.ids)

Plate.metabar <- merge(ASV.ids, asv_tab2, by = "ASV") %>%
  dplyr::select(-ASV) %>%
  pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "Reads") %>%
  group_by(Taxon, Sample) %>%
  summarise(Reads = sum(Reads, na.rm = TRUE)) %>%
  pivot_wider(names_from = "Sample", values_from = "Reads")

write.csv(Plate.metabar, "aasign_tax_output.csv")



