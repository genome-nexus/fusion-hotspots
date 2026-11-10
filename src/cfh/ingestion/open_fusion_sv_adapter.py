"""Present normalized open fusion tables as structural-variant rows.

:func:`cfh.ingestion.open_fusion_tables.normalize_table` writes one schema for
all four open tables. The domain screen reads cBioPortal structural-variant
rows (see ``sv_parser``). This adapter maps one onto the other without
filtering, deduplicating or guessing:

* ``gene_5p``/``gene_3p`` become Site1/Site2, and ``Event_Info`` records the
  explicit 5'->3' order (``{GENE5:GENE3} fusion``), so orientation never
  depends on a connection type.
* Positions, chromosomes and the build are copied verbatim. STAR-Fusion
  reports the last 5' base and the first 3' base of the junction, which is the
  screen's default ``both_bases_retained`` convention.
* Junction (split) reads become ``Tumor_Split_Read_Count`` and spanning
  fragments ``Tumor_Paired_End_Read_Count``. There is no precise/imprecise
  label in these tables, so ``Breakpoint_Type`` stays empty.
* Rows with a missing gene, or the same gene on both sides, are kept; the
  screen itself decides they are not screenable.

The ``Open_*`` columns carry the table's own fields through the screen.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from cfh.ingestion.open_fusion_tables import SCHEMA

OPEN_PASSTHROUGH_COLUMNS = (
    "Open_source_table",
    "Open_tissue_class",
    "Open_tumor_type",
    "Open_strand_5p",
    "Open_strand_3p",
    "Open_ffpm",
    "Open_splice_type",
    "Open_caller",
    "Open_fi_validated",
)
"""Open-table fields carried onto screen rows (``cfh domain-screen`` writes them)."""


def open_table_sv_rows(path: str | Path) -> pd.DataFrame:
    """Read a normalized open-table TSV as structural-variant rows.

    The header must be exactly :data:`~cfh.ingestion.open_fusion_tables.SCHEMA`.
    Output has one row per input row, in input order.
    """
    with Path(path).open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = tuple(next(reader, ()))
        if header != SCHEMA:
            raise ValueError(f"{path} is not a normalized open fusion table (unexpected header)")
        records = []
        for number, values in enumerate(reader, start=1):
            if len(values) != len(SCHEMA):
                raise ValueError(f"{path} row {number} has {len(values)} fields")
            row = dict(zip(SCHEMA, values, strict=True))
            records.append(_sv_row(row, number))
    return pd.DataFrame.from_records(records)


def _blank(value: str) -> str | None:
    return value if value != "" else None


def _count(value: str) -> int | None:
    """A nonnegative integer (integral scientific notation such as ``1e+05`` allowed)."""
    if value == "":
        return None
    number = float(value)
    if not number.is_integer() or number < 0:
        raise ValueError(f"count {value!r} is not a nonnegative integer")
    return int(number)


def _sv_row(row: dict[str, str], number: int) -> dict[str, object]:
    five, three = _blank(row["gene_5p"]), _blank(row["gene_3p"])
    event_info = f"{{{five}:{three}}} fusion" if five and three else "fusion (genes unavailable)"
    return {
        "Sample_Id": _blank(row["sample_id"]),
        "Patient_Id": _blank(row["patient_id"]),
        "Site1_Hugo_Symbol": five,
        "Site1_Chromosome": _blank(row["chrom_5p"]),
        "Site1_Position": _count(row["pos_5p"]),
        "Site2_Hugo_Symbol": three,
        "Site2_Chromosome": _blank(row["chrom_3p"]),
        "Site2_Position": _count(row["pos_3p"]),
        "Event_Info": event_info,
        "Annotation": None,
        "NCBI_Build": _blank(row["genome_build"]),
        "Tumor_Split_Read_Count": _count(row["junction_reads"]),
        "Tumor_Paired_End_Read_Count": _count(row["spanning_reads"]),
        "Breakpoint_Type": None,
        "Source_row_number": number,
        "Open_source_table": row["source_table"],
        "Open_tissue_class": _blank(row["tissue_class"]),
        "Open_tumor_type": _blank(row["tumor_type"]),
        "Open_strand_5p": _blank(row["strand_5p"]),
        "Open_strand_3p": _blank(row["strand_3p"]),
        "Open_ffpm": _blank(row["ffpm"]),
        "Open_splice_type": _blank(row["splice_type"]),
        "Open_caller": _blank(row["caller"]),
        "Open_fi_validated": _blank(row["fi_validated"]),
    }
