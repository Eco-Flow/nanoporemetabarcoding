#!/usr/bin/env python3
"""
Build the wasp_test_data test dataset from REAL Nanopore reads, replacing the
earlier fully-synthetic version (generate_test_data.py) which used independent
per-read random error simulation and too few reads/well -- amplicon_sorter
couldn't cluster them (its default similar_species=85%/similar_consensus=96%
thresholds need the correlated, context-driven error profile real ONT reads
have, plus real per-well redundancy).

Real reads come from Jordan's actual wasp gut-content Nanopore run
(NanoporeJordan_romauld/barcode06, barcode07 -- the "WaspEx" plates; NOT
barcode01/02, which use a different, unrelated primer scheme -- see that
project's Readme.md). Provenance is documented in REAL_WELLS below: each
maps a role in our training scenario to the exact real primer_comb it came
from and the real BLAST hit that run's full pipeline assigned it (from
final_results_fixed/assign_taxa/*/ASV_taxa.csv), which is also how the
reference species in custom_db.fasta were chosen -- they are the actual
real top BLAST hits for these wells, fetched from NCBI by accession, so
this small custom_db gives correct hits instead of "no hit everywhere".

Real reads were demultiplexed once, upstream of this script, with the
pipeline's own cutadapt container and exact conf/modules.config args:

    docker run --rm -v $PWD:/data quay.io/biocontainers/cutadapt:5.2--py311haab0aaa_0 \\
        cutadapt --cores 4 --rc -e 2 \\
        --untrimmed-output=/data/unknown.trim.fastq.gz \\
        -g file:/data/tags_f_real.fasta -o /data/{name}.trim.fastq.gz /data/merged_barcodeNN.fastq.gz

    docker run --rm -v $PWD:/data quay.io/biocontainers/cutadapt:5.2--py311haab0aaa_0 \\
        cutadapt --cores 4 --rc --discard-untrimmed -e 2 \\
        -g file:/data/tags_r_real.fasta -o /data/{fwd_name}_{name}.trim.fastq.gz /data/{fwd_name}.trim.fastq.gz

(tags_f_real.fasta / tags_r_real.fasta = the subset of the real project's
primers_f.fasta / primers_r.fasta needed for the combos in REAL_WELLS.)

The FASTQs cutadapt produced are pure, already tag/primer-trimmed inserts,
carrying real ONT reads' actual (correlated, not independently random)
error profile -- so unlike the old synthetic reads, subsamples of these
reliably cluster in amplicon_sorter, because that's the exact same data
(same barcodes, same tags, same cutadapt version/args) the real pipeline
run already clustered successfully to produce final_results_fixed.

This script then re-flanks a random subsample of each real well's reads
with OUR training scenario's fixed tag+primer sequences (F1-F3/R1-R2, the
ones documented in docs/nanopore_metabarcoding.md Step 4) so the dataset
keeps that doc's simple 3x2 design regardless of which real tag combo the
reads actually originated from. POS_CON and BLANK still have no usable
real signal anywhere in this experiment (every real POS_CON well across
all 4 WaspEx plates had <=4 reads total -- the real positive control
failed to amplify), so those two wells remain synthetic, clearly marked
as such below.

Usage:
    python3 build_real_test_data.py --real-reads-dir /path/to/cutadapt/output --outdir .
"""

import argparse
import gzip
import random
from pathlib import Path

IUPAC = {
    "W": "AT", "S": "CG", "R": "AG", "Y": "CT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Training scenario's fixed tags (docs/nanopore_metabarcoding.md Step 4).
PRIMERS_F = {
    "F1_WaspExF_Tab1": "AACAAGCCCCTTTATCWTSWRRWWTTGS",
    "F2_WaspExF_Tab2": "GGAATGAGTCCTTTATCWTSWRRWWTTGS",
    "F3_WaspExF_Tab3": "AATTGCCGGTCCTTTATCWTSWRRWWTTGS",
}
PRIMERS_R = {
    "R1_LuthienR_Tab29": "GAGTAACCACTTCWGGRTGWCCAAARAAYCA",
    "R2_LuthienR_Tab54": "CGATGAGTTACTTCWGGRTGWCCAAARAAYCA",
}

# Reference species for --custom_db: the ACTUAL real top BLAST hits for the
# real wells below, fetched from NCBI by the exact accession the full
# pipeline run found (final_results_fixed/assign_taxa/*/ASV_taxa.csv), plus
# one synthetic spike-in for the synthetic POS_CON well.
REFERENCES = {
    # Real hit for CM4_1_4, CM9_4_1, CM4_1_3 (plate06) and the wasp itself
    # showing up in EXT_NEG as contamination.
    "Belonogaster_juncea": (
        "LC510521.1",
        "AACTTTATATTTTATATTTAGCCTTTGAGCTGGAACATTAGGAGCTAGCTTAAGTATAATCATTCGTTTA"
        "GAGTTAAGAACACCAGGAATATTAATTAGAGATGACCAATTATTTAATACTATTGTAACATCTCATGCTT"
        "TAATTATAATTTTTTTTATAGTTATACCTTTTATAATTGGAGGATTTGGTAACTGACTTATTCCTTTAAT"
        "ACTAGGAGTACCTGATATAGCTTTTCCACGAATAAATAATATAAGATTCTGACTTTTACCACCATCATTA"
        "ATTATACTTATATTAAGAAGAATTATTGGAATAGGTGTTGGAACAGGATGAACACTTTATCCACCTTTAT"
        "CATCAATCTTAGGACATAACTCAATTTCAGTAGATTTAAGAATTTTTTCATTACATATTGCAGGAATTTC"
        "TTCAATTATAGGAGCAATTAATTTTATTGTAACTATTCTAAATATACATATTAAAACATCTTCACTTAAT"
        "TTTATTCCTTTATTTGCTTGATCAATCTTAATTACTACAATTTTACTTTTACTTTCTTTACCTGTATTAG"
        "CTGGCGCAATTACTATACTTTTAACAGATCGTAATTTAAATACAACTTTTTTTGATCCAACTGGTGGAGG"
        "AGACCCAATTTTATTTCAACATTTATTT",
    ),
    # Real hit for CM11_4_4 (plate07).
    "Noctuidae_sp_BOLD_AAL8777": (
        "HM893567.1",
        "AACTTTATATTTTATTTTTGGAATTTGAGCAGGAATAGTAGGAACCTCTTTAAGTTTATTAATTCGAGCT"
        "GAATTAGGTAACCCTGGATCATTGATTGGTGATGATCAAATTTATAATACTATTGTTACAGCTCATGCAT"
        "TTATTATAATTTTTTTTATAGTAATACCTATTATAATTGGTGGATTTGGAAATTGATTAGTCCCTTTAAT"
        "ATTAGGAGCTCCTGATATAGCTTTCCCTCGAATAAATAATATAAGATTTTGACTTCTTCCTCCTTCATTA"
        "ACTCTTTTAATTTCAAGAAGAATTGTAGAAAATGGTGCAGGAACTGGATGAACTGTATACCCCCCTTTAT"
        "CTTCTAATATTGCCCATAGAGGAAGATCTGTAGATTTAGCAATTTTTTCCTTACATTTAGCGGGAATTTC"
        "CTCAATTCTAGGAGCTATTAATTTTATTACAACAATTATTAATATACGATTAAATAGTTTAATATTTGAT"
        "CAAATACCTTTATTTATTTGAGCAGTAGGAATTACAGCATTTTTACTTTTACTTTCTTTACCTGTATTAG"
        "CTGGTGCTATTACTATACTTTTAACCGATCGTAATTTAAATACATCATTTTTTGACCCTGCTGGTGGAGG"
        "AGACCCAATTCTTTATCAACATTTATTT",
    ),
    # Real hit for CM14_2_6/CM14_1_1 (plate07) and EXT_NEG_10/11 contamination.
    # Just the ~425bp region cutadapt/BLAST actually matched in the genome,
    # not the whole chromosome.
    "Providencia_rettgeri": (
        "CP157876.1:2283054-2283478",
        "TGTTCTTCAGTTAAATAGTTGAGTTTCAATGCAGACTCTTTCAAAGTCAGCCCTTCTTTGTGCGCTTTTT"
        "TAGCGATTTCTGCCGCTTTATCGTACCCAATGTGAGTATTTAAAGCTGTTACTAGCATTAATGATTCATG"
        "AAGTAATTTTTCAATACGCTCACGGTTTGGTTCAATTCCAATCGCACAGTGCTCATTAAAGCTACGCATA"
        "CCATCCGCTAACAAACGAACGGATTGTAAGAAGTTATCAATTAGCATTGGGCGGAATACATTCAGCTCAA"
        "AATTACCGGATGCACCACCAATATTAACAGCAACGTCATTTCCCATCACTTGTGCGCATAACATGGTTAA"
        "TGCTTCACACTGCGTAGGGTTAACTTTACCTGGCATGATTGAGCTGCCTGGCTCATTTTCTGGAATCGAA"
        "ATTTC",
    ),
    # Second real hit found alongside Providencia in CM14_2_6/CM14_1_1 -- a
    # relative of the wasp itself, not the same species as Belonogaster_juncea.
    "Belonogaster_sp2": (
        "MT444667.1",
        "ACTTTATATTTTATATTTAGCCTTTGAGCTGGAACATTAGGAGCTAGATTAAGTATAATTATTCGTTTAG"
        "AACTAAGAACACCAGGAATATTAATTGGAAATGATCAATTATTCAATACTATTGTAACATCACATGCTTT"
        "AATTATAATTTTTTTTATAGTTATACCTTTTATAATTGGAGGATTTGGTAATTGACTTATTCCTTTAATA"
        "TTAGGAGTACCTGATATAGCTTTTCCACGAATAAATAATATAAGATTCTGACTTTTACCTCCATCATTAA"
        "TTATACTTATATTAAGAAGAATTATTGGAATAGGTGTTGGAACAGGATGAACACTTTACCCACCTTTATC"
        "ATCAATTTTAGGACATAATTCAATTTCAGTAGATTTAAGAATTTTTTCATTACATATCGCAGGAATCTCT"
        "TCAATTATAGGAGCAATTAATTTTATTGTAACTATTTTAAATATACATATTAAAACATCCTCACTTAATT"
        "TTATTCCTTTATTTGCTTGATCAATTTTAATTACTACAATTTTACTTTTACTTTCTTTACCTGTCTTAGC"
        "TGGTGCAATTACTATGCTTTTAACAGATCGTAATTTAAATACAACTTTTTTTGATCCAACTGGTGGAGGA"
        "GATCCAATTTTATTTCAACATTTATTTTGATTTTTTGGACATCCAGAAGTATATATTTTAATTTTACCTG"
        "GATTTGGAATTATTTCCCATATTGTTACTAATGAAACTGGTAAAAAAGAAATTTTTGGAACTCTTGGAAT"
        "AATTTATGCAATAATTGCTATTGGAGCACTAGGCTTTATTGTTTGAGCTCATCATATATTTACTGTTGGA"
        "TTAGATATTGATACCCGAGCATATTTCACTTCAGCAACTATAGTAATTGCTATCCCAACTGGAATTAAAG"
        "TATTTAGATGAATAGCTTCAATTTATGGATCAAAAATAATTTTTTCTCCTGCAATAATATGAAGAATGGG"
        "TTTTATTTTTTTATTTACAATTGGTGGATTAACGGGAATTATTTTATCAAATTCATCATTAGATATTATA"
        "CTTCACGATACTTATTATGTTATTGGACACTTTCATTACGTATTATCAATAGGTGCAGTATTTGCAATTA"
        "TTGCAGGATTTATTCATTGATTTCCTTTATTTTTTGGTATATCTTTAAATAAAATTTGATTAAAAATTCA"
        "ATTCTCAATTATATTTATTGGAGTTAATATAACTTTTTTCCCTCATCATTTCTTAGGCCTTTCTGGTTTC"
        "CCTCGACGATACTCTGATTATCCAGATATTTTTTTAACATGAAATGTAATTTCATCATTAGGATCAATTA"
        "TATCCTTCTTATCAATAATTTTATTTATTTATATTATTTGAGAAAGATTTTATATTCAACGATTTATTAT"
        "TTATAAATTTTATATACCTTCAAATATAGAATGAATTCACTCAATACCTCCAACTACTCACACTTATTCT"
        "CAAATTCCTTATACTACAAATTTA",
    ),
    # Synthetic spike-in for the synthetic POS_CON well (see REAL_WELLS).
    "Drosophila_melanogaster": (
        "OZ481787.1",
        "AATAATTTCTCATATTATTAGACAAGAATCAGGAAAAAAGGAAACTTTTGGTTCTCTAGGAATAATTTAT"
        "GCTATATTAGCTATTGGATTATTAGGATTTATTGTATGAGCTCATCATATATTTACCGTTGGAATAGATG"
        "TAGATACTCGAGCTTATTTTACCTCAGCTACTATAATTATTGCAGTTCCTACTGGAATTAAAATTTTTAG"
        "TTGATTAGCTACTTTACATGGAACTCAACTTTCTTATTCTCCAGCTATTTTATGAGCTTTAGGATTTGTT"
        "TTTTTATTTACAGTAGGAGGATTAACAGGAGTTGTTTTAGCTAATTCATCAGTAGATATTATTTTACATG"
        "ATACTTATTATGTAGTAGCTCATTTTCATTATGTTTTATCTATAGGAGCTGTATTTGCTATTATAGCAGG"
        "TTTTATTCACTGATACCCTTTATTTACTGGATTAACGTTAAATAATAAATGATTAAAAAGTCATTTCATT"
        "ATTATATTTATTGGAGTTAATTTAACATTTTTTCCTCAACATTTTTTAGGATTGGCTGGAATACCTCGAC"
        "GTTATTCAGATTACCCAGATGCTTACACAACATGAAATATTGTATCAACTATTGGATCAACTATTTCATT"
        "ATTAGGAATTTTATTCTTTTTTTTTATTATTTGAGAAAGTTTAGTATCACAACGACAAGTAATTTACCCA"
        "ATTCAACTAAATTCATCAATTGAATGATACCAAAATACTCCACC",
    ),
}

# real_reads_dir-relative path to the cutadapt output for each real well,
# plus its provenance: which real barcode/plate/primer_comb/sample it is,
# and what the real full-pipeline run found there.
REAL_WELLS = {
    "plate01": {  # <- real barcode06 (WaspEx Plate1, "Woodland")
        "F1_WaspExF_Tab1_R1_LuthienR_Tab29": (
            "EXT_NEG_F1_R1",
            "b06_r/WaspExF_Tab1_LuthienR_Tab29.trim.fastq.gz",
            "real EXT_NEG well; real hit = Belonogaster juncea (contamination)",
        ),
        "F2_WaspExF_Tab2_R1_LuthienR_Tab29": (
            "WD_wasp01",
            "b06_r/WaspExF_Tab8_LuthienR_Tab31.trim.fastq.gz",
            "real CM9_4_1; real hit = Belonogaster juncea",
        ),
        "F3_WaspExF_Tab3_R1_LuthienR_Tab29": (
            "WD_wasp02",
            "b06_r/WaspExF_Tab1_LuthienR_Tab31.trim.fastq.gz",
            "real CM4_1_3; real hit = Belonogaster juncea",
        ),
        "F1_WaspExF_Tab1_R2_LuthienR_Tab54": (
            "WD_wasp03",
            "b06_r/WaspExF_Tab1_LuthienR_Tab54.trim.fastq.gz",
            "real CM4_1_4; real hit = Belonogaster juncea",
        ),
    },
    "plate02": {  # <- real barcode07 (WaspEx Plate2, "Meadow")
        "F1_WaspExF_Tab1_R1_LuthienR_Tab29": (
            "EXT_NEG_F1_R1",
            "b07_r/WaspExF_Tab6_LuthienR_Tab38.trim.fastq.gz",
            "real EXT_NEG_10 well; real hit = Providencia rettgeri (contamination)",
        ),
        "F2_WaspExF_Tab2_R1_LuthienR_Tab29": (
            "MW_wasp01",
            "b07_r/WaspExF_Tab1_LuthienR_Tab39.trim.fastq.gz",
            "real CM11_4_4; real hit = Noctuidae sp. BOLD:AAL8777",
        ),
        "F3_WaspExF_Tab3_R1_LuthienR_Tab29": (
            "MW_wasp02",
            "b07_r/WaspExF_Tab6_LuthienR_Tab39.trim.fastq.gz",
            "real CM14_2_6; real hits = Providencia rettgeri + Belonogaster sp. 2 (mixed well)",
        ),
        "F1_WaspExF_Tab1_R2_LuthienR_Tab54": (
            "MW_wasp03",
            "b07_r/WaspExF_Tab6_LuthienR_Tab29.trim.fastq.gz",
            "real CM14_1_1; real hits = Belonogaster sp. 2 + Providencia rettgeri (mixed well)",
        ),
    },
}

# Wells with no usable real signal: every real POS_CON well across all 4
# WaspEx plates (06-09) had <=4 reads total in the full run -- the real
# positive control failed to amplify. BLANK legitimately having ~0 real
# template is also just correct PCR-NTC behaviour, so both stay synthetic.
SYNTHETIC_WELLS = {
    "F2_WaspExF_Tab2_R2_LuthienR_Tab54": ("POS_CON_F2_R2", "Drosophila_melanogaster", 8),
    "F3_WaspExF_Tab3_R2_LuthienR_Tab54": ("BLANK_F3_R2", None, 2),
}

PLATES = ["plate01", "plate02"]


def resolve_iupac(seq, rng):
    return "".join(rng.choice(IUPAC[b]) if b in IUPAC else b for b in seq)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_dna(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def apply_ont_errors(seq, rng, sub_rate=0.05, ins_rate=0.015, del_rate=0.015):
    out = []
    for base in seq:
        r = rng.random()
        if r < del_rate:
            continue
        elif r < del_rate + ins_rate:
            out.append(rng.choice("ACGT"))
            out.append(base)
        elif r < del_rate + ins_rate + sub_rate:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        else:
            out.append(base)
    return "".join(out)


def flank_quality(n, rng, mean_q=14):
    return "".join(chr(33 + max(2, min(40, int(rng.gauss(mean_q, 4))))) for _ in range(n))


def read_fastq_records(path):
    records = []
    with gzip.open(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()
            records.append((seq, qual))
    return records


def write_fastq_record(fh, read_id, seq, qual):
    fh.write(f"@{read_id}\n{seq}\n+\n{qual}\n".encode())


def build_real_read(insert_seq, insert_qual, fwd_resolved, rev_resolved, rng):
    seq = fwd_resolved + insert_seq + revcomp(rev_resolved)
    qual = flank_quality(len(fwd_resolved), rng) + insert_qual + flank_quality(len(rev_resolved), rng)
    return seq, qual


def build_synthetic_read(insert_source, fwd_resolved, rev_resolved, rng):
    if insert_source is None:
        insert = random_dna(rng.randint(0, 6), rng)
    else:
        insert = REFERENCES[insert_source][1]
    raw = fwd_resolved + insert + revcomp(rev_resolved)
    noisy = apply_ont_errors(raw, rng)
    qual = flank_quality(len(noisy), rng)
    return noisy, qual


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--real-reads-dir", required=True, help="Dir containing b06_r/, b07_r/ cutadapt output")
    parser.add_argument("--outdir", default=".", help="Output directory (default: current dir)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--reads-per-well", type=int, default=100)
    args = parser.parse_args()

    real_dir = Path(args.real_reads_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    with open(outdir / "primers_f.fasta", "w") as f:
        for header, seq in PRIMERS_F.items():
            f.write(f">{header}\n{seq}\n")
    with open(outdir / "primers_r.fasta", "w") as f:
        for header, seq in PRIMERS_R.items():
            f.write(f">{header}\n{seq}\n")

    with open(outdir / "custom_db.fasta", "w") as f:
        for name, (accession, seq) in REFERENCES.items():
            f.write(f">{name} {accession}\n{seq}\n")

    metadata_rows = ["id,primer_comb,sample"]
    samplesheet_rows = ["id,fastq"]
    barcode_dirs = {"plate01": "barcode01", "plate02": "barcode02"}
    fastq_names = {"plate01": "plate1_combined.fastq.gz", "plate02": "plate2_combined.fastq.gz"}

    for plate_id in PLATES:
        barcode_dir = barcode_dirs[plate_id]
        fastq_name = fastq_names[plate_id]
        plate_dir = outdir / barcode_dir
        plate_dir.mkdir(exist_ok=True)
        samplesheet_rows.append(f"{plate_id},wasp_test_data/{barcode_dir}/{fastq_name}")

        wells = dict(REAL_WELLS[plate_id])
        wells.update(SYNTHETIC_WELLS)
        read_counter = 0

        with gzip.open(plate_dir / fastq_name, "wb") as fq:
            for primer_comb, well in wells.items():
                fwd_header, rev_header = primer_comb.split("_R", 1)
                rev_header = "R" + rev_header
                sample = well[0]
                metadata_rows.append(f"{plate_id},{primer_comb},{sample}")

                fwd_resolved = resolve_iupac(PRIMERS_F[fwd_header], rng)
                rev_resolved = resolve_iupac(PRIMERS_R[rev_header], rng)

                if primer_comb in REAL_WELLS[plate_id]:
                    _, rel_path, _provenance = well
                    records = read_fastq_records(real_dir / rel_path)
                    rng.shuffle(records)
                    chosen = records[: args.reads_per_well]
                    for insert_seq, insert_qual in chosen:
                        seq, qual = build_real_read(insert_seq, insert_qual, fwd_resolved, rev_resolved, rng)
                        read_counter += 1
                        read_id = f"{plate_id}-{read_counter:08x}-real runid=wasp_test"
                        write_fastq_record(fq, read_id, seq, qual)
                else:
                    _, insert_source, n_reads = well
                    for _ in range(n_reads):
                        seq, qual = build_synthetic_read(insert_source, fwd_resolved, rev_resolved, rng)
                        read_counter += 1
                        read_id = f"{plate_id}-{read_counter:08x}-sim runid=wasp_test"
                        write_fastq_record(fq, read_id, seq, qual)

        print(f"{plate_id}: wrote {read_counter} reads -> {plate_dir / fastq_name}")

    with open(outdir / "metadata.csv", "w") as f:
        f.write("\n".join(metadata_rows) + "\n")
    with open(outdir / "samplesheet.csv", "w") as f:
        f.write("\n".join(samplesheet_rows) + "\n")

    print(f"\nDone. Test dataset written to {outdir.resolve()}")


if __name__ == "__main__":
    main()
