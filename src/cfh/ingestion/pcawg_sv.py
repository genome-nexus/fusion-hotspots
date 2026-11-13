"""Strict, data-only ingestion of the public PCAWG TCGA consensus BEDPE archive.

Coordinates remain BED half-open intervals until conversion to the existing
MSK SV model (end = 1-based anchor). A BEDPE '+' retains the left genomic
side, '-' the right: a gene contributes 5' iff its strand equals that sign.
Every compatible ordered pair is emitted; multiple overlapping genes are
explicitly ambiguous. Intergenic and incompatible calls stay in the audit.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import tarfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, TextIO

HEADER = (
    "chrom1",
    "start1",
    "end1",
    "chrom2",
    "start2",
    "end2",
    "sv_id",
    "pe_support",
    "strand1",
    "strand2",
    "svclass",
    "svmethod",
)
CALLERS = frozenset({"SNOWMAN", "BRASS", "dRANGER", "DELLY"})
CHROMS = frozenset([str(i) for i in range(1, 23)] + ["X", "Y", "MT", "hs37d5"])
MAX_UNPACKED = 3_000_000_000
BIN = 100_000


@dataclass(frozen=True)
class Gene:
    identifier: str
    symbol: str
    chromosome: str
    start: int  # zero-based inclusive
    end: int  # zero-based exclusive
    strand: str
    biotype: str


class GeneIndex:
    def __init__(self, genes: Iterable[Gene]):
        self.bins: dict[tuple[str, int], list[Gene]] = defaultdict(list)
        for gene in genes:
            if gene.start < 0 or gene.end <= gene.start or gene.strand not in {"+", "-"}:
                raise ValueError(f"invalid gene interval: {gene}")
            for bucket in range(gene.start // BIN, (gene.end - 1) // BIN + 1):
                self.bins[gene.chromosome, bucket].append(gene)

    def overlap(self, chromosome: str, start: int, end: int) -> list[Gene]:
        genes = {
            gene
            for bucket in range(start // BIN, (end - 1) // BIN + 1)
            for gene in self.bins.get((chromosome, bucket), [])
            if gene.start < end and start < gene.end
        }
        return sorted(genes, key=lambda g: (g.symbol, g.identifier))


def read_gencode(comp: Path, attrs: Path) -> GeneIndex:
    """Union transcript spans per GENCODE v19 gene from UCSC hg19 tables.

    All biotypes are retained in the overlap audit, including antisense genes.
    Noncoding genes will not acquire an invented protein/domain model.
    """
    annotations = {}
    with gzip.open(attrs, "rt") as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) < 5:
                raise ValueError("malformed GENCODE attrs row")
            annotations[row[4]] = (row[0], row[1], row[2])
    spans: dict[tuple[str, str, str, str, str], tuple[int, int]] = {}
    with gzip.open(comp, "rt") as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) != 16 or row[1] not in annotations:
                raise ValueError("malformed/unmatched GENCODE genePred row")
            identifier, symbol, biotype = annotations[row[1]]
            key = (identifier, symbol, row[2].removeprefix("chr"), row[3], biotype)
            start, end = int(row[4]), int(row[5])
            old = spans.get(key, (start, end))
            spans[key] = min(start, old[0]), max(end, old[1])
    return GeneIndex(
        Gene(gid, sym, chrom, s, e, strand, bio)
        for (gid, sym, chrom, strand, bio), (s, e) in spans.items()
    )


def sample_index(path: Path) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    required = {
        "aliquot_id",
        "donor_unique_id",
        "donor_wgs_exclusion_white_gray",
        "dcc_project_code",
        "icgc_specimen_id",
        "submitter_specimen_id",
        "icgc_donor_id",
        "submitter_donor_id",
        "library_strategy",
        "dcc_specimen_type",
    }
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not required <= set(reader.fieldnames or []):
            raise ValueError("sample sheet missing required columns")
        for line, row in enumerate(reader, 2):
            if None in row or any(not row.get(k) for k in required):
                raise ValueError(f"sample sheet malformed row {line}")
            if row["donor_wgs_exclusion_white_gray"] not in {"Whitelist", "Graylist", "Excluded"}:
                raise ValueError(f"sample sheet unknown status at {line}")
            key = row["aliquot_id"]
            if key in index and index[key] != row:
                raise ValueError(f"conflicting sample sheet aliquot {key}")
            index[key] = row
    return index


def parse_bedpe(handle: TextIO, source: str) -> list[dict[str, Any]]:
    reader = csv.DictReader(handle, delimiter="\t")
    if tuple(reader.fieldnames or ()) != HEADER:
        raise ValueError(f"{source}: unexpected BEDPE header {reader.fieldnames}")
    result = []
    identifiers: set[str] = set()
    for line, row in enumerate(reader, 2):
        try:
            if None in row or any(row[k] in (None, "") for k in HEADER):
                raise ValueError("missing/surplus field")
            record: dict[str, Any] = dict(row)
            for key in ("start1", "end1", "start2", "end2", "pe_support"):
                record[key] = int(row[key])
            for side in (1, 2):
                if (
                    record[f"chrom{side}"] not in CHROMS
                    or record[f"start{side}"] < 0
                    or record[f"end{side}"] != record[f"start{side}"] + 1
                    or record[f"strand{side}"] not in {"+", "-"}
                ):
                    raise ValueError(f"invalid breakend {side} (requires single-base interval)")
            methods = row["svmethod"].split("_")
            if (
                not set(methods) <= CALLERS
                or len(set(methods)) != len(methods)
                or record["pe_support"] < 0
            ):
                raise ValueError("invalid caller/support field")
            expected = {
                "DEL": ("+", "-"),
                "DUP": ("-", "+"),
                "h2hINV": ("+", "+"),
                "t2tINV": ("-", "-"),
            }
            if row["svclass"] not in {*expected, "TRA"}:
                raise ValueError("unknown SV class")
            if row["svclass"] in expected and (
                (row["strand1"], row["strand2"]) != expected[row["svclass"]]
                or row["chrom1"] != row["chrom2"]
            ):
                raise ValueError("SV class inconsistent with strands/chromosomes")
            if row["sv_id"] in identifiers:
                raise ValueError("duplicate SV id")
            identifiers.add(row["sv_id"])
            record.update(caller_support=len(methods), source_file=source, source_line=line)
            result.append(record)
        except (ValueError, TypeError) as error:
            raise ValueError(f"{source}:{line}: {error}") from error
    return result


def read_archive(
    path: Path, samples: dict[str, dict[str, str]]
) -> tuple[list[dict[str, Any]], dict]:
    """Never extract paths or execute archive contents; cap nested inflation too."""
    result = []
    donors = set()
    seen = set()
    inflated = 0
    file_row_counts: dict[str, int] = {}
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        member_bytes = sum(m.size for m in members)
        if member_bytes > MAX_UNPACKED:
            raise ValueError("STOP: archive exceeds 3 GB unpacked")
        for member in members:
            if not member.isfile() and not member.isdir():
                raise ValueError(f"unsafe archive member {member.name}")
            if not member.name.endswith(".bedpe.gz"):
                continue
            aliquot = Path(member.name).name.split(".")[0]
            if aliquot in seen or aliquot not in samples:
                raise ValueError(f"duplicate or unmapped aliquot {aliquot}")
            seen.add(aliquot)
            sample = samples[aliquot]
            if (
                not sample["dcc_project_code"].endswith("-US")
                or sample["library_strategy"] != "WGS"
                or "tumour" not in sample["dcc_specimen_type"].lower()
            ):
                raise ValueError(f"not TCGA tumour WGS: {aliquot}")
            stream = archive.extractfile(member)
            assert stream is not None
            with gzip.GzipFile(fileobj=stream) as zipped:
                payload = zipped.read(MAX_UNPACKED - inflated + 1)
            inflated += len(payload)
            if inflated > MAX_UNPACKED:
                raise ValueError("STOP: nested BEDPE exceeds 3 GB unpacked")
            rows = parse_bedpe(io.StringIO(payload.decode("utf-8")), member.name)
            file_row_counts[aliquot] = len(rows)
            for row in rows:
                row.update(sample)
            result.extend(rows)
            donors.add(sample["donor_unique_id"])
    return result, {
        "archive_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "tar_member_bytes": member_bytes,
        "nested_bedpe_bytes": inflated,
        "bedpe_files": len(seen),
        "donors": len(donors),
        "sv_rows": len(result),
        "file_row_counts": file_row_counts,
        "empty_bedpe_files": sum(count == 0 for count in file_row_counts.values()),
    }


def pair_genes(
    row: dict[str, Any], genes: GeneIndex
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    overlaps = [genes.overlap(row[f"chrom{s}"], row[f"start{s}"], row[f"end{s}"]) for s in (1, 2)]
    ambiguous = any(len({g.identifier for g in hits}) > 1 for hits in overlaps)
    pairs = []
    for first in overlaps[0]:
        for second in overlaps[1]:
            if first.symbol == second.symbol:
                continue
            first_five = first.strand == row["strand1"]
            second_five = second.strand == row["strand2"]
            if first_five == second_five:
                continue
            five, three = (first, second) if first_five else (second, first)
            five_site, three_site = (1, 2) if first_five else (2, 1)
            pairs.append(
                {
                    **row,
                    "Sample_Id": row["aliquot_id"],
                    "Patient_Id": row["icgc_donor_id"],
                    "Site1_Hugo_Symbol": five.symbol,
                    "Site2_Hugo_Symbol": three.symbol,
                    "Site1_Chromosome": five.chromosome,
                    "Site2_Chromosome": three.chromosome,
                    "Site1_Position": row[f"end{five_site}"],
                    "Site2_Position": row[f"end{three_site}"],
                    "Connection_Type": "5to3",
                    "NCBI_Build": "GRCh37",
                    "SV_Status": "SOMATIC",
                    "Event_Info": f"protein fusion {{{five.symbol}:{three.symbol}}}",
                    "Annotation": None,
                    "Breakpoint_Type": None,
                    "Site2_Effect_On_Frame": None,
                    "Tumor_Split_Read_Count": None,
                    "Tumor_Paired_End_Read_Count": None,
                    "gene_pairing_status": "ambiguous" if ambiguous else "unique",
                    "five_gene_id": five.identifier,
                    "three_gene_id": three.identifier,
                    "five_gene_strand": five.strand,
                    "three_gene_strand": three.strand,
                }
            )
    audit = {
        **row,
        "genes1": ";".join(f"{g.symbol}:{g.strand}:{g.biotype}" for g in overlaps[0]),
        "genes2": ";".join(f"{g.symbol}:{g.strand}:{g.biotype}" for g in overlaps[1]),
        "pair_count": len(pairs),
        "pairing_status": (
            "intergenic"
            if not all(overlaps)
            else "ambiguous"
            if ambiguous
            else "sense_pairs"
            if pairs
            else "no_sense_pair"
        ),
    }
    return pairs, audit
