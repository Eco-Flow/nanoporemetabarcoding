# nf-core/nanoporemetabarcoding: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.1 - [date]

### `Added`

- Optional chimera detection before the final BLAST, enabled with `--chimera_detection` and selected via `--chimera_method`:
  - `blast` (default): reference-based detection that splits each polished consensus into 5′/3′ halves, BLASTs them against the reference database, and flags a sequence as chimeric only when the two halves share no reference subject in common among their top hits (tolerating redundant same-organism database entries rather than requiring an exact accession match). Configurable with `--chimera_min_identity`, `--chimera_min_coverage`, `--chimera_top_n` and `--chimera_bitscore_margin`. Emits a non-chimeric FASTA and a per-sequence report to `${outdir}/chimera_detection`
  - `vsearch`: de novo detection with [`vsearch`](https://github.com/torognes/vsearch) `--uchime3_denovo` (abundance-based, no reference required). Per-ASV read counts from `amplicon_sorter` are written onto the consensus headers as `;size=N` annotations so vsearch uses real abundances. Emits non-chimeric/chimeric FASTAs and a UCHIME report to `${outdir}/chimera_detection`
- Optional [`Porechop`](https://github.com/rrwick/Porechop) step on raw reads, enabled with `--run_porechop`, to trim residual adapters and split reads containing internal adapters (ligation/split chimeras) before filtering and demultiplexing

### `Changed`

- nf-schema plugin bumped to `2.7.3` for compatibility with Nextflow 26.04, restoring the `validation` config scope (and coercion of bare boolean CLI flags such as `--skip_nanoplot`)

---

## v1.0.0dev - [date]

### `Added`

- Community matrix output from taxonomy assignment results, with both abundance and presence/absence matrices ([`vegan`](https://github.com/vegandevs/vegan), [`tidyr`](https://tidyr.tidyverse.org/))
- Configurable taxonomic level for community matrix construction (`--community_tax`)
- Taxonomy proportion plots per barcode and per sample using [`ggplot2`](https://ggplot2.tidyverse.org/), with configurable taxonomic levels (`--tax_list`) and minimum fraction threshold for rare taxa grouping
- Sample name validation now rejects underscores in addition to spaces
- Synthetic wasp test dataset (`wasp_test_data/`) built from real Nanopore reads, replacing the earlier fully-synthetic version that was too sparse for `amplicon_sorter` to cluster reliably
- New `test_synth` profile (`conf/test_synth.config`) to run the pipeline against the synthetic wasp dataset

### `Changed`

- Local modules (`AMPLICON_SORTER`, `ASSIGN_TAXONOMY`, `PLOT_TAXONOMY`, `COMMUNITY_MATRIX`) migrated to nf-core module structure with `environment.yml`, `meta.yml`, and test stubs
- `AMPLICON_SORTER` and `MEDAKA` intermediate outputs disabled from `publishDir` by default to reduce storage usage
- `FIND_CONCATENATE` intermediate outputs disabled from `publishDir` by default
- `BLAST_MAKEBLASTDB` `publishDir` `enable` typo corrected to `enabled`
- `MEDAKA` input staging mode set to `copy` to avoid symlink issues
- `MEDAKA` output prefix now includes the barcode identifier (`old_id`) to ensure unique filenames across barcodes

### `Fixed`

- `BLAST_MAKEBLASTDB` publish directory config used incorrect key `enable` instead of `enabled`
- Community matrix `pivot_wider` now correctly aggregates multiple ASVs mapping to the same taxon per sample (`values_fn = sum` for abundance, `values_fn = max` for presence/absence)

---

## v1.0.0dev - [date] (initial)

Initial release of nf-core/nanoporemetabarcoding, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Read filtering and trimming with [`NanoFilt`](https://github.com/wdecoster/nanofilt)
- Read QC with [`NanoPlot`](https://github.com/wdecoster/NanoPlot) on raw, filtered, and demultiplexed reads
- MultiQC report aggregating NanoPlot QC results
- Two-step tag+primer demultiplexing with [`Cutadapt`](https://github.com/marcelm/cutadapt): first on forward tags/primers, then on matching reverse tags/primers
- Optional metadata file to map primer-tag combinations to sample names
- Samplesheet and metadata input validation
- Amplicon grouping into consensus sequences with [`amplicon_sorter`](https://github.com/avierstr/amplicon_sorter)
- Consensus sequence polishing with [`Medaka`](https://github.com/nanoporetech/medaka)
- Custom BLAST database creation with [`makeblastdb`](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- Consensus sequence annotation with [`blastn`](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- Taxonomy assignment from BLAST hits using [`taxonomizr`](https://github.com/sherrillmix/taxonomizr) (NCBI GenBank and RefSeq accessions only)
- Configurable percent identity thresholds for taxonomy assignment at species, genus, family, and order level (`--spident`, `--gpident`, `--fpident`, `--opident`)
- Lowest Common Ancestor (LCA) consensus logic: when an ASV has multiple equally good BLAST hits, the most specific taxonomic rank consistent across all hits is assigned; falls back to coarser ranks if hits disagree, labelling as `Unassigned` if no consensus is found
- Read count per ASV included in the final taxonomy output table
- Version tracking for all pipeline tools including R packages (`argparse`, `dplyr`, `tidyr`, `stringr`, `taxonomizr`)

### `Fixed`

- Medaka consensus polishing not running correctly
- NanoPlot results not being picked up by MultiQC
- Tool versions not being reported correctly in pipeline output
- Sample renaming using incorrect sample identifier

### `Dependencies`

- [`NanoFilt`](https://github.com/wdecoster/nanofilt)
- [`NanoPlot`](https://github.com/wdecoster/NanoPlot)
- [`Cutadapt`](https://github.com/marcelm/cutadapt)
- [`amplicon_sorter`](https://github.com/avierstr/amplicon_sorter)
- [`Medaka`](https://github.com/nanoporetech/medaka)
- BLAST+ (`makeblastdb`, `blastn`)
- R packages: `argparse`, `dplyr`, `tidyr`, `stringr`, `taxonomizr`

### `Deprecated`
