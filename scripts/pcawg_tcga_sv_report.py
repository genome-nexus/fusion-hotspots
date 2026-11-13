"""Summarize the completed PCAWG TCGA run and emit focused ETS/kinase tables."""

from __future__ import annotations

import json
from collections import Counter

from pcawg_tcga_sv_screen import NONTRIVIAL, OUT, read_rows, write_rows

ETS = ("ERG", "FLI1", "ETV1", "ETV4", "ETV5", "GABPA", "ETS1", "ETS2", "FEV")
KINASE = {"PF00069", "PF07714"}


def table(counts: dict) -> str:
    return "| Outcome | Count |\n|---|---:|\n" + "".join(
        f"| {key} | {value:,} |\n" for key, value in sorted(counts.items())
    )


def main() -> None:
    summary = json.loads((OUT / "summary.json").read_text())
    screen = read_rows(OUT / "domain_screen.tsv.gz")
    audit = read_rows(OUT / "breakend_gene_audit.tsv.gz")
    candidates = read_rows(OUT / "nontrivial_candidates.tsv")
    ets_screen = [r for r in screen if {r["five_prime_gene"], r["three_prime_gene"]} & set(ETS)]
    kinases = [r for r in screen if r["pfam_family"] in KINASE]
    write_rows(OUT / "ets_sense_pairs.tsv", ets_screen)
    write_rows(OUT / "kinase_shared_family.tsv", kinases)
    ets_audit = []
    for row in audit:
        hits = {
            part.split(":")[0]
            for side in ("genes1", "genes2")
            for part in row[side].split(";")
            if part
        }
        if hits & set(ETS):
            ets_audit.append(row)
    write_rows(OUT / "ets_breakend_audit.tsv", ets_audit)
    ets_table = (
        "| Gene | SVs overlapping gene | Sense pairs | Non-trivial family rows |\n"
        "|---|---:|---:|---:|\n"
    )
    for gene in ETS:
        overlap = sum(
            any(
                gene == part.split(":")[0]
                for side in ("genes1", "genes2")
                for part in row[side].split(";")
            )
            for row in audit
        )
        own = [r for r in ets_screen if gene in (r["five_prime_gene"], r["three_prime_gene"])]
        nontrivial = sum(r["reconstitution_class"] in NONTRIVIAL for r in own)
        ets_table += f"| {gene} | {overlap} | {len(own)} | {nontrivial} |\n"
    priority = {
        "clean_reconstitution": 0,
        "gap_or_deletion": 1,
        "tandem_duplication_or_insertion": 2,
        "ambiguous_junction_register": 3,
    }
    ranked = sorted(
        candidates,
        key=lambda r: (
            {"likely_intact": 0, "uncertain": 1, "likely_broken": 2}.get(r["sequence_tier"], 3),
            priority.get(r["reconstitution_class"], 4),
            r["gene_pairing_status"] != "unique",
            r["paralog_sequence_flag"] == "flagged",
            -int(r["caller_support"]),
            r["fusion_name"],
            r["sample_id"],
            r["pfam_family"],
        ),
    )
    strongest = []
    seen = set()
    for row in ranked:
        key = row["fusion_name"], row["pfam_family"]
        if key in seen:
            continue
        seen.add(key)
        strongest.append(row)
        if len(strongest) == 10:
            break
    write_rows(OUT / "strongest_candidates.tsv", strongest)
    bullets = (
        "\n".join(
            f"- **{r['fusion_name']} / {r['pfam_family']}**, {r['icgc_donor_id']} "
            f"({r['dcc_project_code']}, {r['caller_support']}/4 callers): "
            f"{r['reconstitution_class']}; sequence-only tier {r['sequence_tier']}; "
            f"paralog check {r['paralog_sequence_flag']}, gene pairing {r['gene_pairing_status']}; "
            f"overall **{r['tier']}** (DNA/read/repeat/RNA validation incomplete)."
            for r in strongest
        )
        or "No non-trivial candidate rows."
    )
    flag_text = "\n".join(
        f"### {name}\n\n{table(counts)}" for name, counts in summary["flag_counts"].items()
    )
    kinase_nt = [r for r in kinases if r["reconstitution_class"] in NONTRIVIAL]
    metrics = read_rows(OUT / "domain_function_metrics.tsv.gz")
    sequence_table = (
        "| Fusion | Indel aa | Domain bits | GA margin | Best/union HMM coverage | "
        "Lost key features | Sequence tier |\n"
        "|---|---:|---:|---:|---:|---:|---|\n"
    )
    for metric in metrics:
        sequence_table += (
            f"| {metric['fusion']} | {metric['pilot_indel_residues']} | "
            f"{metric['chimera_domain_bits']} | {metric['chimera_ga_margin']} | "
            f"{metric['chimera_match_coverage']}/{metric['chimera_match_coverage_union']} | "
            f"{metric['key_lost']} | {metric['tier']} |\n"
        )
    report = f"""# PCAWG TCGA consensus DNA SV domain-completion screen — 2026-11-13

This first genome-wide DNA screen examines **{summary["ingest"]["sv_rows"]:,} PASS consensus SVs**,
from **822 TCGA donor BEDPE files** (796 whitelisted, 26 greylisted). **789 donors**
have at least one SV; **33 files contain only a header**. The sample sheet contains
828 eligible TCGA donors (801 whitelist + 27 greylist); six have no archive file.
Empty files and absent files are different observations, not negative biological
results. See `sample_inventory.tsv` and `donors_without_bedpe.tsv`.

There are **zero clean reconstitutions**, **{len(candidates)} non-trivial event-family rows**
from **12 distinct SVs**, and
**{summary["sequence_records"]} distinct sequence records**. Overall candidate tiers
remain fail-closed when evidence is absent. The #155 sequence-only tier is a
separate heuristic about the predicted protein; it does not validate a DNA call.

## Data, build and provenance

The [PCAWG open archive README](https://object.genomeinformatics.org/icgc25k-open/PCAWG/consensus_sv/README.md)
publishes archive MD5 `650981d737522350503bba74e6d2c65d`, which matches this download.
Archive size is 3,128,049 bytes; tar members total 3,053,586 bytes, and nested
BEDPE text totals 11,353,243 bytes (below the 3 GB stopping threshold).
All downloads reside in the new gitignored `data/raw/pcawg_sv/` directory.
Tar paths are never extracted and no downloaded content is executed.
`input_manifest.json` records source URLs, byte sizes and SHA256 hashes.

**Build: GRCh37 / hs37d5.** The [PCAWG primary paper](https://www.nature.com/articles/s41586-020-1969-6)
identifies hs37d5 for alignment in Extended Data Fig. 1. The saved Ensembl
assembly response identifies GRCh37.p13; GENCODE v19 uses GRCh37 (UCSC hg19).
The public BEDPE itself has no assembly header, so this is dataset-level build
provenance, not a per-row declaration. No liftOver or GRCh38 cache is used.
Transcript/Pfam models use the existing MSK SV screen and MSKCC isoform override
on `https://www.genomenexus.org`; CDS and reference windows use
`https://grch37.rest.ensembl.org`. Build-stamped caches live under
`runs/.pcawg_tcga_sv_work/transcripts_GRCh37/` and `flags_GRCh37/`.

## Ingestion and gene pairing

The strict adapter checks the exact header, field counts, integer coordinates,
single-base half-open intervals, legal chromosomes/strands/SV types, type/strand
consistency, caller vocabulary and duplicates. Malformed data stops the run.
Aliquot IDs must map to TCGA tumour WGS samples. Donor and specimen IDs,
whitelist/greylist, cancer project, SV type, strands, four-caller support and raw
`pe_support` are retained. The source has no split-read or confidence-interval
columns; single-base intervals do not establish breakpoint precision. Read and
precision confidence remain **not assessed**. `pe_support` is source metadata,
never substituted for split reads or used to pass an evidence gate.

GENCODE v19 transcript spans are unioned per gene and strand, retaining all
biotypes in the overlap audit. BEDPE `+` retains the left genomic side; `-`
retains the right side (confirmed by the DEL/DUP and inversion sign patterns).
A gene contributes 5′ when gene strand equals retained-side sign; otherwise it
contributes 3′. Only distinct genes with complementary roles become ordered
sense-sense fusions. Reverse order is handled symmetrically. Every compatible
combination is retained; multiple overlapping gene IDs mark the assignment
ambiguous. Intergenic ends and genic calls with no compatible pair remain in
`breakend_gene_audit.tsv.gz`, without an invented partner. Intragenic calls are
not two-partner fusions. Noncoding/unmapped genes receive no invented model.

{table(summary["gene_pairing_counts"])}

These are SV-level categories; a call can yield multiple gene pairs.
**{summary["sense_gene_pairs"]:,} ordered pairs** enter the existing any-Pfam screen.
The BEDPE half-open `end` is the 1-based anchor passed into the MSK input model.
The existing junction predictor evaluates its standard retained-base convention
and its alternative 3′-base-dropped convention; both frame results are saved.
DNA calls retain mid-exon cuts and do not impose an exon-boundary filter.
Predicted products are reference-only hypotheses: public BEDPE has no inserted
junction sequence and does not demonstrate transcription or the splicing outcome.

## Domain classes

Counts are **event-family rows**, including one no-shared/model row per pair
when no family can be classified; they are not counts of unique SVs or donors.
An absent canonical lookup is explicitly `no_transcript_model` / `not_assessed`,
not a successful no-shared-family result. Missing model/CDS also cannot pass.
Canonical gene IDs and exon strands must agree with GENCODE at each breakend;
otherwise the row is not assessed. Successful batch omissions are cached as absent
mappings; failed batches are
listed in `summary.json` and are never persisted as biological negatives.
There are **20,928 outcome rows**, with **1,073 annotated shared-family rows**
from **598 gene-pair events**; only **149 shared-family rows** have usable models
and reach artefact assessment. The 9,907 model-unassessed outcomes limit sensitivity.

{table(summary["screen"]["class_counts"])}

## Sequence checks and final tiers

The unchanged #155 machinery rebuilds clean/repeat/gap/ambiguous-register
sequences, scores Pfam HMMs, evaluates UniProt key residues/motifs, indel position,
and parent structure context. Out-of-frame, premature-stop and unresolved-splice
rows remain explicit but have no fabricated sequence metrics or tier. Required
missing HMM/annotation/structure inputs fail closed. `domain_function_metrics.tsv.gz`
contains every sequence metric and missing-input reason; `chimera_records.json.gz`
links sequence records to sample observations. Candidate rows carry the sequence
record ID, sequence-only tier, final tier, junction peptide, both frames, donor,
cancer project, coordinates, caller support and all assessed flags.

Final candidate tiers (DNA evidence plus sequence):

{table(summary["candidate_tiers"])}

Sequence-only tiers on candidate rows (unscorable rows = not assessed):

{table(summary["candidate_sequence_tiers"])}

## Artefact checks

The #154 paralog/single-transcript explanation and genomic-window checks run on
all translated/shared-family classes, including `not_split_on_both`, for context.
Reference similarities at ±150 and ±500 bp retain identity and aligned length.
No segdup/RepeatMasker tracks were supplied: repeat/duplication flags are
`not_assessed`, never a pass. Missing read support also remains `not_assessed`.
Consequently an unflagged reference-window comparison cannot by itself pass
the combined genomic/evidence flag. Normal RNA recurrence is not assessed.
Shared-family call burden includes annotated shared families whose transcript
model is unavailable, and is counted per donor by distinct aliquot + SV ID,
without multiplying Pfam families or gene pairings. `donor_paralog_calls` counts
distinct donor SVs whose junction can be explained by a single transcript or
paralog sequence under #154; it is not the total DNA SV burden. These tests
identify suspicion and cannot establish the origin of a DNA rearrangement.
`donor_burden.tsv` preserves the complete donor inventory and observed burden
lower bounds. A donor with unannotated pairs cannot pass a low-burden flag;
positive high-burden flags can still be raised from observed calls.

{flag_text}

## Strongest candidates for follow-up

Sorted by sequence-only functional tier, junction class, unambiguous pairing,
absence of a positive paralog flag, then caller support; no candidate is promoted through missing
validation data. Full rows are in `strongest_candidates.tsv` and every non-trivial
row, including failed/uncertain frames, is in `nontrivial_candidates.tsv`.

{bullets}

## Explicit ETS-family examination

`ets_breakend_audit.tsv` includes all SVs overlapping any requested ETS gene,
including intergenic/antisense/intragenic cases.
No requested ETS gene yields a non-trivial domain-completion candidate; in
particular, there is no ERG::FLI1 or FLI1::ERG consensus call in this dataset.
`ets_sense_pairs.tsv` includes
all compatible sense gene pairs and their screen outcome. This DNA approach
can examine mid-exon ETS breaks that RNA exon-boundary fusion tables omit.

{ets_table}

## Kinase pairs

Shared kinase-family rows are defined explicitly as Pfam **PF00069** (protein
kinase) or **PF07714** (tyrosine kinase), not by gene-name heuristics.
There are **{len(kinases)} shared kinase-family rows**, **{len(kinase_nt)} non-trivial**.
`kinase_shared_family.tsv` contains every such pair with donor and breakpoints;
non-trivial rows also appear in the complete candidate table. This definition
does not assert exhaustive coverage of atypical kinase families or pairs assigned
to different Pfam families.

{table(dict(Counter(r["reconstitution_class"] for r in kinases)))}

## Limits and validation

This is a PASS-only consensus screen. Homologous regions are likely under-called;
absence of a consensus SV is not evidence that a fusion cannot occur. Canonical
isoforms and available Pfam annotations limit sensitivity; noncoding genes,
missing mappings, and alternate transcripts remain annotation limitations. No
RNA confirmation, inserted-junction sequence, split reads, controlled VCF
confidence intervals, repeat/segdup track assessment or matched normal RNA
validation is available here. Sequence-only plausibility does not imply an
expressed, stable or functional fusion protein, nor establish oncogenicity.

The ingest/strand suite contains 59 collected cases across 18 test functions
(parametrized strand and
malformed-input coverage), including missing-lookup and unavailable-evidence
regressions. Final repository checks and exact suite counts are recorded in
`validation.json` and the PR description. Commands run from this worktree:

```sh
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev mypy src
uv run --extra dev pytest -q
```

Reproduce after downloading the manifest inputs to the ignored raw directory:

```sh
uv run python scripts/pcawg_tcga_sv_screen.py --stage screen
uv run python scripts/pcawg_tcga_sv_screen.py --stage flags
uv run python scripts/pcawg_tcga_sv_screen.py --stage sequence
uv run python scripts/pcawg_tcga_sv_report.py
```
"""
    (OUT / "REPORT.md").write_text(report)
    print(f"report: {len(candidates)} non-trivial rows; {len(kinases)} kinase rows")


if __name__ == "__main__":
    main()
