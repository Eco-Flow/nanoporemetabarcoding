
library(argparse)
library(dplyr)


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

print(blast_filtered)
