"""Package domain-screen pilot TSVs as sorted, reproducible gzip files.

Usage::

    python -I scripts/domain_screen_pilot_package.py SRC.tsv OUT.tsv.gz all|nontrivial

Rows are sorted by cohort, class, gene pair, breakpoints, Pfam family and
sample, so identical junctions sit together and compress well; the gzip
header carries no timestamp. ``nontrivial`` keeps only clean, ambiguous,
tandem, gap, out-of-frame, junction-stop and splice-site-unresolved rows.
"""

import csv
import gzip
import io
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)
src, out, mode = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
NT = {
    "clean_reconstitution",
    "ambiguous_junction_register",
    "tandem_duplication_or_insertion",
    "gap_or_deletion",
    "out_of_frame",
    "premature_stop_at_junction",
    "splice_site_unresolved",
}
with src.open(newline="") as h:
    r = csv.reader(h, delimiter="\t")
    header = next(r)
    rows = list(r)
ix = {c: i for i, c in enumerate(header)}
if mode == "nontrivial":
    rows = [x for x in rows if x[ix["reconstitution_class"]] in NT]


def key(x):
    return (
        x[ix["cohort"]],
        x[ix["reconstitution_class"]],
        x[ix["five_prime_gene"]],
        x[ix["three_prime_gene"]],
        x[ix["five_prime_chromosome"]],
        x[ix["five_prime_breakpoint"]],
        x[ix["three_prime_chromosome"]],
        x[ix["three_prime_breakpoint"]],
        x[ix["pfam_family"]],
        x[ix["sample_id"]],
    )


rows.sort(key=key)
with (
    gzip.GzipFile(out, "wb", compresslevel=9, mtime=0) as raw,
    io.TextIOWrapper(raw, newline="", encoding="utf-8") as text,
):
    writer = csv.writer(text, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
print(out.name, len(rows))
