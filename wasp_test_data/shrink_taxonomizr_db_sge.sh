#!/bin/bash -l
#$ -l h_rt=1:00:0
#$ -l mem=16G
#$ -N shrink_taxonomizr_db

conda activate taxonmizr

python3 shrink_taxonomizr_db.py \
    --custom-db custom_db.fasta \
    --source-db accessionTaxa.sql \
    --output-db nameNode.small.sqlite
