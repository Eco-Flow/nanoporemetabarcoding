#!/usr/bin/env python3
"""
Generate a small, comprehensive raw-read test dataset for the nanoporemetabarcoding
pipeline, matching the wasp/caterpillar experiment in docs/nanopore_metabarcoding.md.

Produces, per plate, RAW (pre-demultiplex) Nanopore-style reads: every read still
carries its well's forward tag+primer at the 5' end and the reverse-complement of
its reverse tag+primer at the 3' end, exactly as real pooled ONT reads would before
Cutadapt splits them. Reads also carry simulated ONT-like sequencing error, so this
is NOT the same as the pipeline's own bundled test_data/ (which the doc notes isn't
representative of this experiment).

Woodland (plate01) and Meadow (plate02) are deliberately given different
species diversity, not just different sample names: Woodland mixes 5 real
prey species across its wasp wells (including a 3-way split in one well),
Meadow is dominated by 1 species with two others only as minor traces --
see WELLS_BY_PLATE.

Usage:
    python3 generate_test_data.py [--outdir .] [--seed 42]

Outputs (relative to --outdir):
    primers_f.fasta, primers_r.fasta   - tag+primer FASTAs (from Step 4 of the doc)
    metadata.csv                       - one row per well, both plates
    samplesheet.csv                    - one row per plate
    custom_db.fasta                    - reference DB for --custom_db / BLAST
    barcode01/plate1_combined.fastq.gz - plate01 (Woodland) raw reads
    barcode02/plate2_combined.fastq.gz - plate02 (Meadow) raw reads
"""

import argparse
import gzip
import random
from pathlib import Path

# ---------------------------------------------------------------------------
# Real reference sequences, fetched from NCBI nucleotide via E-utils (COI gene,
# partial cds, mitochondrial). Embedded here so the generator runs offline.
# ---------------------------------------------------------------------------

REFERENCES = {
    # Woodland-site prey: winter moth, common on oak/hazel.
    "Operophtera_brumata": (
        "OR614726.1",
        "AACTTTATACTTTATTTTTGGAATTTGAGCCGGTATAATTGGAACTTCACTAAGTTTATTAATTCGAGCT"
        "GAATTAGGTAACCCTGGTTCTTTAATTGGGGATGACCAAATTTACAACACTATTGTTACAGCACATGCTT"
        "TTATTATAATTTTTTTTATAGTTATACCAATTATAATTGGAGGATTTGGTAATTGATTAGTACCTTTAAT"
        "ACTTGGAGCTCCTGATATAGCTTTCCCCCGAATAAATAATATAAGATTTTGATTATTACCTCCCTCTATT"
        "ACTCTTTTAATTTCTAGAAGAATTGTAGAAAATGGGGCAGGAACTGGATGAACTGTTTACCCCCCTTTAT"
        "CTTCTAATATTGCCCATGGAGGAAGATCCGTAGATCTAGCTATCTTTTCTCTTCATTTAGCTGGTATTTC"
        "CTCAATTTTAGGTGCAATTAACTTTATTACCACTATTATCAATATACGATTAAATAATATATTTTTTGAC"
        "CAATTACCATTATTTGTTTGAGCTGTAGGAATCACAGCATTTTTACTTTTATTGTCATTACCAGTATTAG"
        "CGGGAGCTATTACTATATTATTAACAGATCGAAATTTAAATACATCATTTTTCGATCCTGCTGGGGGGGG"
        "AGATCCTATTCTTTATCAACACTTATTT",
    ),
    # Meadow-site prey: large white butterfly, feeds on brassicas/forbs.
    "Pieris_brassicae": (
        "OR891167.1",
        "TAAGCCTACTAATTCGAACTGAATTAGGAAATCCAGGATTTTTAATTGGAGACGATCAAATTTACAATAC"
        "TATTGTAACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTCGGA"
        "AATTGATTAGTACCTTTAATATTAGGGGCCCCTGATATAGCTTTCCCCCGAATAAATAATATAAGATTCT"
        "GATTATTACCCCCATCTTTAACCCTTCTAATTTCAAGAAGAATCGTAGAAAATGGAGCAGGAACAGGATG"
        "AACAGTTTACCCCCCACTTTCATCTAATATTGCTCATAGAGGAGCATCAGTAGACTTAGCTATTTTTTCT"
        "CTTCATTTAGCTGGGATTTCATCAATTTTAGGTGCAATTAATTTTATCACTACTATTATCAATATACGAA"
        "TTAGAAATATATCTTTTGATCAAATACCTTTATTTGTTTGAGCAGTTGGAATTACTGCTTTATTATTACT"
        "TCTATCTCTCCCTGTTCTCGCAGGAGCAATCACTATACTTTTAACAGATCGTAATTTAAATACATCATTT"
        "TTTGATCCAGCTGGAGGAGGAGATCCTATTTTATACCAACATTTATTC",
    ),
    # Positive control: distinctive species not expected in wasp/caterpillar samples.
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
    # Second Woodland-site prey: oak leafroller, another classic oak-canopy
    # caterpillar -- gives the Woodland plate a second real species.
    "Tortrix_viridana": (
        "OQ563582.1",
        "AACATTATATTTTATTTTTGGAATTTGAGCAGGTATAATTGGAACTTCTTTAAGTCTTTTAATTCGGGCA"
        "GAATTAGGAAATCCAGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCCCATGCAT"
        "TTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGTAATTGATTAGTACCTTTAAT"
        "ATTAGGAGCTCCTGATATAGCTTTTCCACGAATAAATAATATAAGTTTCTGACTTCTCCCCCCCTCTATT"
        "ATACTTTTAATTTCAAGTAGAATTGTAGAAAACGGAGCAGGAACAGGTTGAACAGTTTATCCCCCCCTTT"
        "CTTCTAATATTGCTCATAGTGGAAGCTCAGTAGATTTAGCAATTTTTTCTTTACATTTAGCTGGAATTTC"
        "CTCAATTTTAGGTGCAGTAAATTTTATTACAACAATTATTAATATACGACCTAATAATATATCATTAGAT"
        "CAAATACCCTTATTTGTTTGAGCTGTAGGAATTACAGCTTTATTATTACTTTTATCTTTACCAGTTTTAG"
        "CAGGTGCTATTACTATACTATTAACAGATCGAAATTTAAATACCTCATTTTTTGACCCAGCTGGGGGAGG"
        "AGATCCAATTTTATATCAACATTTATTT",
    ),
    # Third Woodland-site prey: mottled umber moth, also common on oak/hazel.
    "Erannis_defoliaria": (
        "OR369281.1",
        "AACATTATACTTTATTTTTGGTATTTGAGCTGGAATAGTTGGAACTTCTTTAAGTTTATTAATTCGAGCA"
        "GAATTAGGAAATCCTGGATCTTTAATTGGGGATGATCAAATTTATAACACTATTGTAACAGCCCATGCAT"
        "TTATTATAATTTTTTTTATAGTTATACCAATTATAATTGGAGGTTTTGGAAATTGATTAGTACCTTTAAT"
        "ACTGGGTGCCCCTGATATAGCTTTCCCACGAATAAATAATATAAGATTTTGATTATTACCCCCATCTATT"
        "ACTCTTTTAATTTCAAGAAGAATTGTAGAAAATGGGGCAGGAACTGGTTGAACGGTTTACCCGCCTTTAT"
        "CCTCTAATATCGCTCATGGAGGAAGTTCAGTAGATTTAGCTATTTTTTCACTACATTTAGCTGGTATTTC"
        "TTCAATTTTAGGAGCTATTAATTTTATTACAACAATTATTAATATACGATTAAATAATTTATCATTTGAT"
        "CAAATACCTTTATTTGTTTGATCTGTAGGAATTACAGCATTCTTACTATTATTATCTTTACCAGTTTTAG"
        "CTGGGGCTATTACAATATTATTAACTGATCGAAATTTAAATACATCATTTTTCGACCCCGCAGGAGGGGG"
        "AGACCCAATTCTTTATCAACACTTATTT",
    ),
    # Minor second Meadow-site prey: silver Y moth, common on forbs/grasses --
    # kept to a single, small well so Meadow's diversity stays clearly lower
    # than Woodland's.
    "Autographa_gamma": (
        "PQ525744.1",
        "AACTTTATATTTTATTTTTGGTATTTGAGCTGGAATAGTTGGTACATCTTTAAGATTACTAATTCGAGCA"
        "GAATTAGGAACCCCTGGATCTTTAATTGGTGATGATCAAATTTATAATACTATTGTTACAGCTCATGCAT"
        "TTATTATAATTTTTTTTATAGTTATGCCTATTATAATTGGAGGATTTGGTAATTGACTCGTTCCTCTAAT"
        "ATTAGGAGCTCCTGATATAGCTTTCCCTCGTATAAATAACATAAGTTTTTGACTTTTACCCCCATCTTTA"
        "ACTCTTTTAATTTCTAGAAGAATTGTAGAAAATGGAGCTGGTACTGGATGAACAGTTTATCCCCCACTTT"
        "CATCTAATATCGCCCATGGTGGAAGATCTGTTGATTTAGCTATTTTTTCTTTACATTTAGCTGGAATTTC"
        "ATCAATTTTAGGAGCAATTAATTTTATTACAACAATTATTAATATACGATTAAATAGTTTATCTTTTGAT"
        "CAAATACCTTTATTTATCTGAGCTGTTGGAATTACAGCTTTCCTTTTATTACTTTCTTTACCTGTTTTAG"
        "CAGGGGCAATTACTATACTTTTAACAGATCGTAATTTAAATACTTCTTTTTTTGATCCTGCTGGAGGAGG"
        "AGACCCAATCTTATACCAACATTTATTT",
    ),
    # Fourth Woodland-site prey: brown hairstreak, another oak/hazel specialist.
    "Thecla_betulae": (
        "OR891377.1",
        "GGAACATATTTAAGAATTCTAATTCGTATAGAATTAGGAACTCCTGGATCTTTAATTGGAGACGATCAAA"
        "TTTATAATACTATTGTAACAGCACATGCTTTTATTATAATTTTTTTTATAGTAATACCAATTATAATTGG"
        "AGGATTTGGAAATTGATTAGTACCTTTAATATTAGGAGCCCCAGATATAGCATTTCCTCGAATAAATAAT"
        "ATAAGATTTTGATTATTACCTCCTTCATTAATATTATTAATTTCAAGAAGAATTGTAGAAAATGGAGCAG"
        "GAACAGGATGAACAGTGTACCCCCCACTTTCATCTAATATTGCACATAGAGGAGCCTCTGTAGATTTAGC"
        "TATTTTTTCATTACACTTAGCAGGTATCTCATCAATTTTAGGAGCTATTAATTTTATTACAACTATTATT"
        "AATATACGAATTAATAACTTAAATTTTGATCAAATATCCTTATTCATTTGAGCTGTAGGAATCACAGCTT"
        "TATTATTATTATTATCTCTCCCAGTATTAGCTGGTGCTATTACTATATTATTAACTGACCGAAATTTAAA"
        "TACATCTTTTTTTGACCCTGCTGGAGGAGGAGATCCAATCCTTTATCAACATTTATTT",
    ),
    # Fifth Woodland-site prey: feathered thorn moth, also oak/hazel.
    "Colotois_pennaria": (
        "OR614513.1",
        "AACATTATATTTTATTTTTGGAATTTGAGCAGGAATAGTAGGAACTTCATTAAGATTACTAATTCGAGCT"
        "GAATTAGGAAACCCTGGATCTTTAATTGGGGATGATCAAATTTACAATACTATTGTAACTGCACACGCTT"
        "TTATTATAATTTTTTTTATGGTAATACCAATTATAATTGGAGGATTTGGTAATTGATTGGTACCTTTAAT"
        "ATTAGGGGCCCCAGATATAGCTTTCCCCCGAATAAATAATATAAGATTTTGATTACTACCCCCATCTTTA"
        "ACTCTTTTAATTTCAAGAAGAGTTGTAGAAAATGGAGCTGGCACAGGATGAACAGTTTACCCCCCATTAT"
        "CCTCTAATATTGCACATGGTGGTAGATCTGTAGATTTAGCTATTTTTTCATTACATTTAGCTGGTATTTC"
        "CTCTATTTTAGGAGCTATTAATTTTATTACAACAATCATTAATATACGATTAAATAATATATCTTTTGAT"
        "CAAATACCATTATTTGTCTGAGCTGTAGGTATTACAGCTTTTTTATTATTGTTATCCTTACCTGTTTTAG"
        "CTGGAGCTATCACTATATTACTAACCGATCGAAATTTAAATACATCATTTTTTGATCCTGCAGGAGGAGG"
        "AGACCCAATTTTATATCAACATTTATTT",
    ),
    # Third Meadow-site prey: large yellow underwing, a common grassland moth.
    "Noctua_pronuba": (
        "OQ563874.1",
        "AACATTATATTTTATTTTTGGAATTTGAGCTGGAATAGTAGGAACTTCTTTAAGATTATTAATTCGAGCT"
        "GAATTAGGAAATCCTGGTTCTTTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCACATGCTT"
        "TTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGTAATTGACTTGTTCCTTTAAT"
        "ATTAGGAGCACCAGATATAGCTTTCCCTCGAATAAATAATATAAGTTTTTGACTTCTTCCCCCCTCATTA"
        "ACTCTTTTAATTTCAAGAAGAATTGTAGAAAATGGAGCAGGTACAGGATGAACAGTTTATCCCCCACTTT"
        "CATCTAATATTGCTCATGGAGGAAGATCCGTTGATTTAGCTATTTTCTCCTTACATTTAGCTGGTATTTC"
        "TTCTATTTTAGGAGCTATTAACTTTATTACCACAATTATTAATATACGATTAAATAGATTATCTTTCGAT"
        "CAAATACCATTATTTATTTGAGCAGTAGGAATTACAGCATTTTTATTATTATTATTATCTTTACCTGTAT"
        "TAGCTGGAGCTATTACTATACTTTTAACAGATCGAAATTTAAATACATCTTTTTTTGACCCTGCAGGAGG"
        "AGGAGACCCTATTTTATATCAACATTTATTT",
    ),
}

# Forward tag+primer FASTA (Step 4 cheat sheet).
PRIMERS_F = {
    "F1_WaspExF_Tab1": "AACAAGCCCCTTTATCWTSWRRWWTTGS",
    "F2_WaspExF_Tab2": "GGAATGAGTCCTTTATCWTSWRRWWTTGS",
    "F3_WaspExF_Tab3": "AATTGCCGGTCCTTTATCWTSWRRWWTTGS",
}

# Reverse tag+primer FASTA (Step 4 cheat sheet).
PRIMERS_R = {
    "R1_LuthienR_Tab29": "GAGTAACCACTTCWGGRTGWCCAAARAAYCA",
    "R2_LuthienR_Tab54": "CGATGAGTTACTTCWGGRTGWCCAAARAAYCA",
}

IUPAC = {
    "W": "AT", "S": "CG", "R": "AG", "Y": "CT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def resolve_iupac(seq, rng):
    """Collapse ambiguity codes to one concrete base, as a real synthesized primer molecule would carry."""
    return "".join(rng.choice(IUPAC[b]) if b in IUPAC else b for b in seq)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_dna(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def apply_ont_errors(seq, rng, sub_rate=0.05, ins_rate=0.015, del_rate=0.015):
    """Simulate a typical ONT per-base error profile (substitutions + indels)."""
    out = []
    for base in seq:
        r = rng.random()
        if r < del_rate:
            continue  # deletion: skip base
        elif r < del_rate + ins_rate:
            out.append(rng.choice("ACGT"))  # insertion before this base
            out.append(base)
        elif r < del_rate + ins_rate + sub_rate:
            out.append(rng.choice([b for b in "ACGT" if b != base]))  # substitution
        else:
            out.append(base)
    return "".join(out)


def quality_string(length, rng, mean_q=14, low_q=6):
    """Rough ONT-like quality string: mostly middling quality with noisy dips."""
    return "".join(chr(33 + max(2, min(40, int(rng.gauss(mean_q, 4))))) for _ in range(length))


def make_read(fwd_seq, rev_seq, insert, rng, error_rate_scale=1.0):
    # Calibrated against the pipeline's actual amplicon_sorter (default
    # similar_species=85%/similar_consensus=96% identity thresholds): two
    # independently-noisy copies of the same true sequence have expected
    # pairwise identity of roughly 1 - 2*error_rate*(1-error_rate). The
    # original 8% total error rate put that at ~85%, right on the
    # similar_species edge, so with only a handful of reads per well
    # amplicon_sorter assigned 0/N reads to any cluster. At ~4% total error
    # pairwise identity is ~92%, and with 100 reads/well amplicon_sorter
    # reliably clusters and rebuilds an exact-length, 100%-identity
    # consensus (verified locally against this exact container/version).
    raw = fwd_seq + insert + revcomp(rev_seq)
    noisy = apply_ont_errors(
        raw, rng,
        sub_rate=0.03 * error_rate_scale,
        ins_rate=0.005 * error_rate_scale,
        del_rate=0.005 * error_rate_scale,
    )
    qual = quality_string(len(noisy), rng)
    return noisy, qual


def write_fastq_record(fh, read_id, seq, qual):
    fh.write(f"@{read_id}\n{seq}\n+\n{qual}\n".encode())


# ---------------------------------------------------------------------------
# Well design: one entry per (forward tag, reverse tag) combination. Each
# well's `composition` is a list of (species_or_None, read_count) pairs --
# most wells are a single species, but a well can mix two, matching what
# real gut-content wells sometimes show (see the real WaspEx run behind
# this dataset: e.g. one real well split between Providencia rettgeri and
# a second Belonogaster species). The wasp/POS_CON wells total 150 reads so
# amplicon_sorter has enough redundancy for a stable consensus per species
# even after a mixed well is split (calibrated down to ~40 reads/species
# minimum -- see make_read()). EXT_NEG and BLANK both stay at 2 tiny junk
# reads with no real insert -- a clean extraction blank or PCR no-template
# control shouldn't produce a consensus at all (no tissue went in, so
# there's nothing for even a non-exclusion primer to amplify), and these
# fall well below amplicon_sorter's -min 300bp filter so it correctly
# drops them rather than fabricating a species-level hit.
#
# Woodland (plate01) is deliberately the more diverse site: 5 real prey
# species (Operophtera_brumata, Tortrix_viridana, Erannis_defoliaria,
# Thecla_betulae, Colotois_pennaria) mixed across its wasp wells, including
# a 3-way split (calibrated down to 40 reads/species minimum -- see
# make_read()). Meadow (plate02) is deliberately less diverse: one dominant
# species (Pieris_brassicae) across all three wasp wells, with two others
# (Autographa_gamma, Noctua_pronuba) only as minor traces in one well each
# -- giving the two plates' community matrices genuinely different
# richness/evenness instead of just different sample names.
# ---------------------------------------------------------------------------

WELLS_BY_PLATE = {
    "plate01": [
        # (fwd_header, rev_header, sample_name, composition, error_scale)
        ("F1_WaspExF_Tab1", "R1_LuthienR_Tab29", "EXT_NEG_F1_R1",
         [(None, 2)], 1.0),
        ("F2_WaspExF_Tab2", "R1_LuthienR_Tab29", "WD_wasp01",
         [("Operophtera_brumata", 60), ("Tortrix_viridana", 50), ("Thecla_betulae", 40)], 1.0),
        ("F3_WaspExF_Tab3", "R1_LuthienR_Tab29", "WD_wasp02",
         [("Tortrix_viridana", 90), ("Colotois_pennaria", 60)], 1.0),
        ("F1_WaspExF_Tab1", "R2_LuthienR_Tab54", "WD_wasp03",
         [("Operophtera_brumata", 90), ("Erannis_defoliaria", 60)], 1.0),
        ("F2_WaspExF_Tab2", "R2_LuthienR_Tab54", "POS_CON_F2_R2",
         [("Drosophila_melanogaster", 150)], 1.0),
        ("F3_WaspExF_Tab3", "R2_LuthienR_Tab54", "BLANK_F3_R2",
         [(None, 2)], 1.0),
    ],
    "plate02": [
        ("F1_WaspExF_Tab1", "R1_LuthienR_Tab29", "EXT_NEG_F1_R1",
         [(None, 2)], 1.0),
        ("F2_WaspExF_Tab2", "R1_LuthienR_Tab29", "MW_wasp01",
         [("Pieris_brassicae", 110), ("Noctua_pronuba", 40)], 1.0),
        ("F3_WaspExF_Tab3", "R1_LuthienR_Tab29", "MW_wasp02",
         [("Pieris_brassicae", 150)], 1.0),
        ("F1_WaspExF_Tab1", "R2_LuthienR_Tab54", "MW_wasp03",
         [("Pieris_brassicae", 110), ("Autographa_gamma", 40)], 1.0),
        ("F2_WaspExF_Tab2", "R2_LuthienR_Tab54", "POS_CON_F2_R2",
         [("Drosophila_melanogaster", 150)], 1.0),
        ("F3_WaspExF_Tab3", "R2_LuthienR_Tab54", "BLANK_F3_R2",
         [(None, 2)], 1.0),
    ],
}

PLATES = [
    # (plate id, ONT barcode dir, fastq filename)
    ("plate01", "barcode01", "plate1_combined.fastq.gz"),
    ("plate02", "barcode02", "plate2_combined.fastq.gz"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=".", help="Output directory (default: current dir)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument(
        "--fastq-prefix", default="wasp_test_data",
        help="Prefix written into samplesheet.csv's fastq column. Nextflow/nf-schema resolves "
             "relative samplesheet paths against the pipeline's launch/project directory, not "
             "the samplesheet's own folder, so this must match this dataset's path from repo root.",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    # --- primer-tag FASTAs ---
    with open(outdir / "primers_f.fasta", "w") as f:
        for header, seq in PRIMERS_F.items():
            f.write(f">{header}\n{seq}\n")
    with open(outdir / "primers_r.fasta", "w") as f:
        for header, seq in PRIMERS_R.items():
            f.write(f">{header}\n{seq}\n")

    # --- reference DB for --custom_db ---
    # Accession must come first: BLAST's default -outfmt 6 sseqid is just the
    # first whitespace-delimited token of the subject header, and
    # assign_taxonomy.R looks up that sseqid via taxonomizr's accessionToTaxa()
    # -- if the species name came first, BLAST would report the name as
    # sseqid instead of the accession, and every taxid lookup would silently
    # fail (every ASV coming back "Unassigned" despite perfect BLAST hits).
    with open(outdir / "custom_db.fasta", "w") as f:
        for name, (accession, seq) in REFERENCES.items():
            f.write(f">{accession} {name}\n{seq}\n")

    # --- samplesheet + metadata + FASTQs ---
    metadata_rows = ["id,primer_comb,sample"]
    samplesheet_rows = ["id,fastq"]

    for plate_id, barcode_dir, fastq_name in PLATES:
        plate_dir = outdir / barcode_dir
        plate_dir.mkdir(exist_ok=True)
        samplesheet_rows.append(f"{plate_id},{args.fastq_prefix}/{barcode_dir}/{fastq_name}")

        wells = WELLS_BY_PLATE[plate_id]
        read_counter = 0
        species_counts = {}

        with gzip.open(plate_dir / fastq_name, "wb") as fq:
            for fwd_header, rev_header, sample, composition, err_scale in wells:
                primer_comb = f"{fwd_header}_{rev_header}"
                metadata_rows.append(f"{plate_id},{primer_comb},{sample}")

                fwd_resolved = resolve_iupac(PRIMERS_F[fwd_header], rng)
                rev_resolved = resolve_iupac(PRIMERS_R[rev_header], rng)

                for insert_source, n_reads in composition:
                    if insert_source is not None:
                        species_counts[insert_source] = species_counts.get(insert_source, 0) + n_reads
                    for _ in range(n_reads):
                        if insert_source is None:
                            # PCR no-template control: primer-dimer-like artifact, no real insert.
                            insert = random_dna(rng.randint(0, 6), rng)
                        else:
                            insert = REFERENCES[insert_source][1]

                        seq, qual = make_read(fwd_resolved, rev_resolved, insert, rng, err_scale)
                        read_counter += 1
                        # Opaque read ID, matching real ONT output: demultiplexing must be
                        # worked out from the primer-tag sequence in the read, not the header.
                        read_id = f"{plate_id}-{read_counter:08x}-sim runid=wasp_test"
                        write_fastq_record(fq, read_id, seq, qual)

        species_summary = ", ".join(f"{sp}={n}" for sp, n in species_counts.items())
        print(f"{plate_id}: wrote {read_counter} reads ({species_summary}) -> {plate_dir / fastq_name}")

    with open(outdir / "metadata.csv", "w") as f:
        f.write("\n".join(metadata_rows) + "\n")
    with open(outdir / "samplesheet.csv", "w") as f:
        f.write("\n".join(samplesheet_rows) + "\n")

    print(f"\nDone. Test dataset written to {outdir.resolve()}")


if __name__ == "__main__":
    main()
