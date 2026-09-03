#!/usr/bin/env python3
"""
Shrink a full taxonomizr sqlite database down to just what's needed for a
small --custom_db BLAST reference FASTA.

bin/assign_taxonomy.R resolves taxonomy in two steps against --sql_db:
    accessionToTaxa(sseqid, sql_db)   -> reads table accessionTaxa(base, accession, taxa),
                                          matching sseqid (a FULL versioned accession, e.g.
                                          "OR614726.1") against the `accession` column -- NOT
                                          `base`, which holds the version-stripped form. This
                                          is taxonomizr's default (accessionToTaxa()'s `version`
                                          arg defaults to "version", which -- confusingly --
                                          means "match the full accession.version string against
                                          the column named `accession`"; verified directly against
                                          taxonomizr 0.11.1's source and this project's actual
                                          accessionTaxa.sql schema, since neither is obvious
                                          from the column names alone.
    getTaxonomy(taxid, sql_db)        -> walks nodes(id, parent, rank) to root,
                                          reading names(id, name, scientific) along the way

Since every possible BLAST hit against custom_db.fasta is necessarily one of
its own sequences, accessionTaxa only ever needs rows for those accessions.
nodes/names can't be pruned down to just the leaf taxids though -- getTaxonomy
needs the full ancestor chain (species -> genus -> family -> ... -> root) to
build a lineage, so this script resolves each accession's taxid, walks
nodes.parent up to the root, and keeps every node/name on that path.

This does NOT build the source database -- it filters an existing one, so you
need a local copy of the full taxonomizr db first (e.g. via R's
taxonomizr::prepareDatabase(), which downloads NCBI's accession2taxid +
taxdump and builds it).

Your source db's accessionTaxa table is a dated snapshot of NCBI's
accession2taxid files, so recently-submitted accessions (newer than that
snapshot) won't be in it -- accessionToTaxa() will find nothing for them
even though the accession is real and the nodes/names tables already have
the (usually long-established) taxid. Rather than re-downloading a fresh
multi-GB snapshot just for that, pass --extra-taxid ACCESSION:TAXID (repeatable)
to manually supply the mapping for specific accessions; look the taxid up at
https://www.ncbi.nlm.nih.gov/taxonomy once, e.g. via the accession's GenBank
record's /db_xref="taxon:NNNNN" field.

Usage:
    python3 shrink_taxonomizr_db.py \\
        --custom-db custom_db.fasta \\
        --source-db /path/to/full/accessionTaxa.sql \\
        --output-db nameNode.small.sqlite \\
        --extra-taxid OZ481787.1:7227 \\
        --extra-taxid PQ729218.1:32391
"""

import argparse
import re
import sqlite3
from pathlib import Path

ACCESSION_RE = re.compile(r"([A-Za-z0-9_]+\.\d+)")


def extract_accessions(fasta_path):
    accessions = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            m = ACCESSION_RE.search(line[1:].strip())
            if m:
                accessions.append(m.group(1))
    return accessions


def parse_extra_taxids(pairs):
    overrides = {}
    for pair in pairs or []:
        acc, _, taxid = pair.partition(":")
        if not taxid:
            raise SystemExit(f"--extra-taxid must be ACCESSION:TAXID, got {pair!r}")
        overrides[acc] = int(taxid)
    return overrides


def find_taxids(cur, accessions, extra_taxids=None):
    extra_taxids = extra_taxids or {}
    cols = {row[1] for row in cur.execute("PRAGMA table_info(accessionTaxa)")}

    taxids = set()
    missing = []
    manual = {}  # accession -> taxid, for rows not present in the source db at all
    for acc in accessions:
        if acc in extra_taxids:
            taxids.add(extra_taxids[acc])
            manual[acc] = extra_taxids[acc]
            continue
        # `accession` holds the full versioned string -- match that first.
        row = cur.execute("SELECT taxa FROM accessionTaxa WHERE accession = ?", (acc,)).fetchone()
        if row is None and "base" in cols:
            # fall back to the version-stripped `base` column, in case the source db
            # has a different version of the same accession (NCBI bumped it since).
            row = cur.execute(
                "SELECT taxa FROM accessionTaxa WHERE base = ?", (acc.split(".")[0],)
            ).fetchone()
        if row is None:
            missing.append(acc)
        else:
            taxids.add(row[0])
    return taxids, missing, manual


def collect_ancestors(cur, taxids):
    ancestors = set()
    frontier = set(taxids)
    while frontier:
        ancestors |= frontier
        placeholders = ",".join("?" * len(frontier))
        rows = cur.execute(
            f"SELECT id, parent FROM nodes WHERE id IN ({placeholders})", tuple(frontier)
        ).fetchall()
        frontier = {parent for _id, parent in rows if parent not in ancestors and parent != _id}
    return ancestors


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--custom-db", required=True, help="custom_db.fasta to read accessions from")
    parser.add_argument("--source-db", required=True, help="full taxonomizr sqlite db to filter")
    parser.add_argument("--output-db", required=True, help="path to write the shrunk sqlite db")
    parser.add_argument(
        "--extra-taxid", action="append", metavar="ACCESSION:TAXID",
        help="manually supply the taxid for an accession missing from --source-db "
             "(repeatable); see the module docstring for how to look one up",
    )
    args = parser.parse_args()
    extra_taxids = parse_extra_taxids(args.extra_taxid)

    accessions = extract_accessions(args.custom_db)
    print(f"{len(accessions)} accession(s) in {args.custom_db}: {accessions}")
    if not accessions:
        raise SystemExit("No accessions found -- expected FASTA headers like '>name ACCESSION.version'.")

    src = sqlite3.connect(args.source_db)
    cur = src.cursor()

    taxids, missing, manual = find_taxids(cur, accessions, extra_taxids)
    if missing:
        raise SystemExit(
            f"{len(missing)} accession(s) not found in accessionTaxa and no --extra-taxid given "
            f"for them: {missing}"
        )
    if manual:
        print(f"Manually patched {len(manual)} accession(s) via --extra-taxid: {manual}")
    if not taxids:
        raise SystemExit("No taxids resolved -- check --source-db and accession formats.")
    print(f"Resolved {len(taxids)} taxid(s): {sorted(taxids)}")

    ancestors = collect_ancestors(cur, taxids)
    print(f"{len(ancestors)} taxid(s) after walking the ancestor chain to root")

    out_path = Path(args.output_db)
    if out_path.exists():
        out_path.unlink()
    src.execute("ATTACH DATABASE ? AS out", (str(out_path),))

    src.execute("CREATE TABLE out.accessionTaxa AS SELECT * FROM accessionTaxa WHERE 0")
    db_accessions = [acc for acc in accessions if acc not in manual]
    if db_accessions:
        acc_placeholders = ",".join("?" * len(db_accessions))
        src.execute(
            f"INSERT INTO out.accessionTaxa SELECT * FROM accessionTaxa WHERE accession IN ({acc_placeholders})",
            db_accessions,
        )
    if manual:
        # these accessions don't exist in --source-db at all, so there's no row to
        # copy -- build one directly from the --extra-taxid mapping instead. `accession`
        # gets the full versioned string (what accessionToTaxa() actually matches on by
        # default), `base` gets the version-stripped form.
        cols = [row[1] for row in cur.execute("PRAGMA table_info(accessionTaxa)")]
        for acc, taxid in manual.items():
            values = {"base": acc.split(".")[0], "accession": acc, "taxa": taxid}
            row = [values.get(col) for col in cols]
            placeholders = ",".join("?" * len(cols))
            src.execute(f"INSERT INTO out.accessionTaxa VALUES ({placeholders})", row)

    src.execute("CREATE TABLE out.nodes AS SELECT * FROM nodes WHERE 0")
    src.execute("CREATE TABLE out.names AS SELECT * FROM names WHERE 0")
    anc_placeholders = ",".join("?" * len(ancestors))
    src.execute(
        f"INSERT INTO out.nodes SELECT * FROM nodes WHERE id IN ({anc_placeholders})", tuple(ancestors)
    )
    src.execute(
        f"INSERT INTO out.names SELECT * FROM names WHERE id IN ({anc_placeholders})", tuple(ancestors)
    )

    src.execute("CREATE INDEX out.idx_accessionTaxa_accession ON accessionTaxa(accession)")
    src.execute("CREATE INDEX out.idx_nodes_id ON nodes(id)")
    src.execute("CREATE INDEX out.idx_names_id ON names(id)")

    src.commit()
    src.execute("DETACH DATABASE out")
    src.close()

    size_kb = out_path.stat().st_size / 1024
    print(f"\nWrote {out_path} ({size_kb:.1f} KiB): "
          f"{len(accessions)} accession row(s), {len(ancestors)} node/name row(s)")


if __name__ == "__main__":
    main()
