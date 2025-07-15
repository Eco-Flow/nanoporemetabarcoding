
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
parser$add_argument('--spident', type = 'numeric', help = 'Identity threshold (in %) for taxonomy assignment at species level', default = 99)
parser$add_argument('--gpident', type = 'numeric', help = 'Identity threshold (in %) for taxonomy assignment at genus level', default = 90)
parser$add_argument('--fpident', type = 'numeric', help = 'Identity threshold (in %) for taxonomy assignment at family level', default = 80)
parser$add_argument('--opident', type = 'numeric', help = 'Identity threshold (in %) for taxonomy assignment at oder level', default = 70)


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
  mutate(seqid2 = paste0(sseqid, ".1")) %>% # Prob should remove the seqid2 column
  mutate(taxaId = accessionToTaxa(sseqid, args$sql_db)) %>%
  mutate(Taxonomic.ranks = getTaxonomy(taxaId, args$sql_db)) %>%
  rename("ASV" = "qseqid")

taxonomic.df <- as.data.frame(sseqids$Taxonomic.ranks, stringsAsFactors = FALSE)
sseqids <- cbind(sseqids, taxonomic.df)
sseqids$Taxonomic.ranks <- NULL

write.csv(sseqids, "ASV_taxa.csv", row.names = FALSE)

print(sseqids)

ASV.ids <- sseqids %>%
  mutate(assigned_taxon = if_else(!is.na(order) & order == "Araneae",
                                if_else(pident > args$gpident, genus,
                                         if_else(pident > args$fpident, family,
                                                 if_else(pident > args$opident, order,
                                                         phylum))),
        if_else(pident > args$spident, species, # Species and genus level assignment where not in Jordan's original code. Is there a reason for this?
            if_else(pident > args$gpident, genus,
                if_else(pident > args$fpident, family,
                     if_else(pident > args$opident, order,
                         phylum)))))) %>%
  #dplyr::select(taxaId, ASV, Taxon) %>%
  dplyr::select(!c(sseqid,seqid2)) %>%
  #mutate(ASV = str_remove(ASV, "_")) %>%
  mutate(sample_name = str_remove(ASV, "_\\d+_\\d+$")) %>%
  relocate(sample_name, .before = 1)

print("debug")

print(ASV.ids)

write.csv(ASV.ids, "ASV_filtered.csv", row.names = FALSE)

#Plate.metabar <- merge(ASV.ids, asv_tab2, by = "ASV") %>%
#  dplyr::select(-ASV) %>%
#  pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "Reads") %>%
#  group_by(Taxon, Sample) %>%
#  summarise(Reads = sum(Reads, na.rm = TRUE)) %>%
#  pivot_wider(names_from = "Sample", values_from = "Reads")

#write.csv(Plate.metabar, "asign_tax_output.csv")



