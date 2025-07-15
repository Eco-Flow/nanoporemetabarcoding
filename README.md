<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-nanoporemetabarcoding_logo_dark.png">
    <img alt="nf-core/nanoporemetabarcoding" src="docs/images/nf-core-nanoporemetabarcoding_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/nanoporemetabarcoding/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/nanoporemetabarcoding/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/nanoporemetabarcoding/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/nanoporemetabarcoding/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/nanoporemetabarcoding/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/nanoporemetabarcoding)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23nanoporemetabarcoding-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/nanoporemetabarcoding)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/nanoporemetabarcoding** is a bioinformatics pipeline for processing nanopore metabarcoding data.

Overview:

![pipeline_diagram](docs/images/pipeline_overview.png)

Steps:

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Filtering and trimming ([`NanoFilt`](https://github.com/wdecoster/nanofilt))
2. Read QC ([`NanoPlot`](https://github.com/wdecoster/NanoPlot))
    - Ran on both raw and filtered and trimmed reads
3. Tags+primer based demultiplexing and trimming ([`Cutadapt`](https://github.com/marcelm/cutadapt)). Divided in two steps:
    1. First, demultiplexing based on the forward tags+primers
    2. Then, demultiplexing is based on the combination of forward and reverse tags+primers
4. Group amplicons reads into "species" (consensus sequences) ([`amplicon_sorter`](https://github.com/avierstr/amplicon_sorter))
5. Consensus sequence correction ([`Medaka`](https://github.com/nanoporetech/medaka))
   - Correction using the consensus sequence from aplicon_sorter as reference and the grouped amplicon reads as the based called data
6. Create custom BLAST database ([`makeblastdb`](https://www.ncbi.nlm.nih.gov/books/NBK279690/))
7. Consensus sequence annotation based on database ([`blastn`](https://www.ncbi.nlm.nih.gov/books/NBK279690/))
   - Annotation is based on the best blast hit per consensus. And best blast hit is based on:
     1. First on the bit score
     2. Second on the e-value
8. Assign taxonomy to blast hits using taxonomizr ([`taxonomizr`](https://github.com/sherrillmix/taxonomizr)). It only works with NCBI accessions (GenBank and RefSeq)

<!-- 1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/)) -->


## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

To use this pipeline, first clone this repo:

```bash
git clone https://github.com/Eco-Flow/nanoporemetabarcoding.git
```

Before running on your data, you can test if the pipeline is suited to your setup by running:

```bash
nextflow run main.nf \
   -profile test,<docker/singularity/conda/.../institute> \
   --outdir <OUTDIR>
```

The running time should be around 8 minutes on a machine with 12 cpus/threads.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):
-->

For running the pipleine on your data, prepare a samplesheet with your input that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq
sample_name,path/to/sample.fastq.gz
```

The first column represents the sample name or sample id, and the second the location to the corresponding FASTQ file.


Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run main.nf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/nanoporemetabarcoding/usage) and the [parameter documentation](https://nf-co.re/nanoporemetabarcoding/parameters).

## Params

Check the config file `nextflow.config` in the params block and fill parameters with the proper values:

**Input options:**

Optionally supply a metadata file in csv format with the primer-tag combination (first column) and the correspoding sample (second column):

```
    // Metadata with primer combinations
    metadata                   = path/to/metadata.csv
```

Check `./test_data/metadata.csv` for an example. If `null` (value by default), the pipeline will use the primer-tag combination as ID.

**Demultiplex options:**

List of forward primer-tag combinations in fasta format:

```
    // Cutadapt options
    tags_f                      = 'path/to/forward/primer-tag.fasta' // List of forward primer-tag combinations in fasta format
```

Check `./test_data/primers_f.fasta` and `./test_data/primers_r.fasta` for examples of `tags_f` and `tags_r` files, respectively.

List of reverse primer-tag combinations in fasta format:

```
    tags_r                      = 'path/to/reverse/primer-tag.fasta' // List of forward primer-tag combinations in fasta format
```

Error rate for adapter removal. It can be a value between 0 and 1, if a maximum error rate wants to be applied to all adapters, or it can be > 1, in which case it will be converted to a maximum error rate depending on the adapter length. Check [cutadatpt documentation](https://cutadapt.readthedocs.io/en/stable/guide.html) for more information:

```
    error_rate                  = 2 // Error rate for adapter removal
```

Filter FASTQs with low number of reads:

```
    filt_fastq                  = 100
```

**Blast options:**

Path to an already built blast database:

```
    blast_db                   = 'path/to/blast.db' // Path to already built database
```

Build local blast database based on a collections of sequences:

```
    custom_db                  = 'path/to/local/database.fasta' // Path to database to be built
```

You should supply either the path to the alraedy built database or provide a custom fasta with sequences to built a custom one. If both are supplied the pipeline will fail.

Output format of blast results. This should probably be hard coded. Do not change:

```
    outfmt                     = 6 // Output format of blast results
```

The E value describes the number of one can “expect” to see by chance when searching a database of a particular size. The lower the value, the more significant a match is:

```
    evalue                     = 0.001 // evalue cutoff
```

Maximum number of hits per query (in this case, per consensus sequences):

```
    max_target_seqs            = 5 // Maximum number of hits per query
```

Maximum number of High-scoring Segment Pair per subject (subjects in the database):

```
    max_hsps                   = null // Maximum number of HSPs per subject (subjects in the database)
```

Minimum query (consensus) coverage per High-scoring Segment Pair:

```
    qcov_hsp_perc              = 90 // Minimum query coverage per HSP
```

For more information, check the appendix section of the [blast documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

**Nanofilt options:**

```
    // Nanofilt options
    nano_quality                = null
    nano_read_length            = 250
```
**Assign taxonomy options:**

Assign taxonomy according to an identity treshold:

```
    // Assign taxonomy options
    spident                     = null // Identity threshold (in %) for taxonomy assignment at species level. > 100 if not wanting to assign species
    gpident                     = null // Identity threshold (in %) for taxonomy assignment at genus level
    fpident                     = null // Identity threshold (in %) for taxonomy assignment at family level
    opident                     = null // Identity threshold (in %) for taxonomy assignment at order level
```

Path to taxonomizr SQL database:

```
    sql_db                      = 'path/to/taxonomizr/database.sqlite' // Path to the taxonomizr SQLite database
```

Alternatively, you can provide/modify pipeline options calling them as command line arguments. For example:

```bash
nextflow run main.nf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
   --nano_quality 10
   --evalue 1e-10
   --tags_f path/to/forward/tags+primers.fasta
   --tags_r path/to/reverse/tags+primers.fasta
   ...
```

## Pipeline output

<!-- To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/nanoporemetabarcoding/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/nanoporemetabarcoding/output). -->

<!-- For an example of some of the relevant results check the `./results` folder in this repo. -->

## Credits

nf-core/nanoporemetabarcoding was originally written by Fernando Duarte, Chris Wyatt.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#nanoporemetabarcoding` channel](https://nfcore.slack.com/channels/nanoporemetabarcoding) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/nanoporemetabarcoding for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
