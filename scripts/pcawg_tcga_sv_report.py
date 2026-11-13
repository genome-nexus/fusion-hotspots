"""Summarize the completed PCAWG TCGA run and emit focused ETS/kinase tables."""

from __future__ import annotations

import json
from collections import Counter

from pcawg_tcga_sv_screen import NONTRIVIAL, OUT, read_rows, save_summary, write_rows

ETS = ("ERG", "FLI1", "ETV1", "ETV4", "ETV5", "GABPA", "ETS1", "ETS2", "FEV")
KINASE = {"PF00069", "PF07714"}


def table(counts: dict) -> str:
    return "| Outcome | Count |\n|---|---:|\n" + "".join(
        f"| {key} | {value:,} |\n" for key, value in sorted(counts.items())
    )


def candidate_key(row: dict) -> tuple:
    return tuple(
        row[k]
        for k in ("aliquot_id", "sv_id", "fusion_name", "pfam_family", "reconstitution_class")
    )


def comparison(
    summary: dict, screen: list[dict], candidates: list[dict], kinases: list[dict]
) -> str:
    baseline = json.loads((OUT / "baseline_539b788.json").read_text())
    old_keys = {candidate_key(r) for r in baseline["nontrivial_candidates"]}
    new = [r for r in candidates if candidate_key(r) not in old_keys]
    write_rows(OUT / "new_candidates_after_lookup_fix.tsv", new)
    old_classes = baseline["screen"]["class_counts"]
    new_classes = summary["screen"]["class_counts"]
    values = {
        "no_transcript_model class rows (all reasons)": (
            old_classes["no_transcript_model"],
            new_classes["no_transcript_model"],
        ),
        "Missing canonical mappings": (
            baseline["model_outcome_counts"]["missing_canonical_mapping"],
            summary["model_outcome_counts"].get("missing_canonical_mapping", 0),
        ),
        "Annotated shared-family rows": (
            baseline["screen"]["family_rows"],
            summary["screen"]["family_rows"],
        ),
        "Rows reaching artefact assessment": (baseline["flag_rows"], summary["flag_rows"]),
        "Non-trivial family rows": (len(baseline["nontrivial_candidates"]), len(candidates)),
        "Sequence records": (baseline["sequence_records"], summary["sequence_records"]),
        "Clean reconstitutions": (
            old_classes["clean_reconstitution"],
            new_classes["clean_reconstitution"],
        ),
        "ETS non-trivial rows": (
            0,
            sum(bool({r["five_prime_gene"], r["three_prime_gene"]} & set(ETS)) for r in candidates),
        ),
        "Kinase non-trivial rows": (
            1,
            sum(r["reconstitution_class"] in NONTRIVIAL for r in kinases),
        ),
    }
    save_summary(
        before_after={
            "baseline_commit": baseline["source_commit"],
            "counts": {k: {"before": a, "after": b} for k, (a, b) in values.items()},
            "new_nontrivial_rows": len(new),
            "new_fusions": sorted({r["fusion_name"] for r in new}),
            "final_tiers_before": baseline["candidate_tiers"],
            "final_tiers_after": summary["candidate_tiers"],
        }
    )
    text = "| Measure | Before 539b788 | Corrected |\n|---|---:|---:|\n"
    text += "".join(f"| {k} | {a:,} | {b:,} |\n" for k, (a, b) in values.items())
    text += "\nNew non-trivial candidates after restoring single-gene fallback:\n\n"
    text += (
        "\n".join(
            f"- {r['fusion_name']} / {r['pfam_family']}, {r['icgc_donor_id']}: "
            f"{r['reconstitution_class']}; sequence {r['sequence_tier']}; final DNA {r['tier']}."
            for r in new
        )
        or "None."
    )
    text += (
        "\n\nBaseline final DNA tiers: "
        + str(baseline["candidate_tiers"])
        + "; corrected: "
        + str(summary["candidate_tiers"])
        + ".\n"
    )
    return text


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
            f"overall **{r['tier']}** (DNA read/repeat validation incomplete)."
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
    comparison_text = comparison(summary, screen, candidates, kinases)
    distinct_svs = len({(r["aliquot_id"], r["sv_id"]) for r in candidates})
    categories = summary["model_outcome_counts"]
    overridden_ambiguous = summary["ambiguous_breakpoint_shared_rows"] - categories.get(
        "base_screen_ambiguous_breakpoints", 0
    )
    report = f"""# PCAWG TCGA consensus DNA SV domain-completion screen — 2026-11-13

This first genome-wide DNA screen examines **{summary["ingest"]["sv_rows"]:,} PASS consensus SVs**,
from **822 TCGA donor BEDPE files** (796 whitelisted, 26 greylisted). **789 donors**
have at least one SV; **33 files contain only a header**. The sample sheet contains
828 eligible TCGA donors (801 whitelist + 27 greylist); six have no archive file.
Empty files and absent files are different observations, not negative biological
results. See `sample_inventory.tsv` and `donors_without_bedpe.tsv`.

There are **{summary["screen"]["class_counts"]["clean_reconstitution"]} clean reconstitutions**,
**{len(candidates)} non-trivial event-family rows** from **{distinct_svs} distinct SVs**, and
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
otherwise the row is not assessed. Batch omissions do **not** establish absence:
the batch helper matches returned `hugoSymbols`, which can be different or missing.
Every omitted/unavailable batch result falls back to the base single-gene GET.
Only a definitive single-gene not-found response is negatively cached; network
and server failures never are. Removed **3,657 poisoned omission markers** from
the prior run. `pre_annotation_class` and `pre_annotation_reason` preserve the
base screen result when the identity/strand check overrides it.
There are **{len(screen):,} outcome rows**, with
**{summary["screen"]["family_rows"]:,} annotated shared-family rows** from
**{summary["screen"]["events_with_a_shared_family"]:,} gene-pair events**;
**{summary["flag_rows"]:,} rows** reach artefact assessment.

The inherited `no_transcript_model` class includes both missing data and base-screen
product/locus verdicts. It must not be read as a count of missing models alone:

{table(categories)}

`missing_canonical_mapping` is a definitive single-gene lookup absence;
`missing_or_invalid_transcript_model` covers missing CDS/exon/valid model data.
Annotation ID/strand mismatches remain separate. Base-screen coding-product
verdicts include past-stop, outside-transcript and before-CDS breaks; outside-locus
and ambiguous-breakpoint decisions also remain separate from missing models.

## Lookup correction: before and after

{comparison_text}

Both runs examine the same 129,273 SVs and 20,462 ordered gene pairs. Recovery
of a partner model can introduce several Pfam family rows in place of one
unassessed placeholder, so family-row counts need not change one-for-one.

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

Final candidate tiers (DNA-applicable evidence plus sequence):

Normal RNA recurrence is **not applicable** to this DNA screen and is excluded
from the final DNA gate. The remaining absent **split-read support** and
**repeat/segmental-duplication tracks** keep the final DNA tiers `not_assessed`.
Positive artefact flags and ambiguous gene assignments can add further limitations.

{table(summary["candidate_tiers"])}

Sequence-only tiers on candidate rows (unscorable rows = not assessed):

{table(summary["candidate_sequence_tiers"])}

{sequence_table}

## Artefact checks

The #154 paralog/single-transcript explanation and genomic-window checks run on
all translated/shared-family classes, including `not_split_on_both`, for context.
Reference similarities at ±150 and ±500 bp retain identity and aligned length.
No segdup/RepeatMasker tracks were supplied: repeat/duplication flags are
`not_assessed`, never a pass. Missing read support also remains `not_assessed`.
Consequently an unflagged reference-window comparison cannot by itself pass
the combined genomic/evidence flag. Normal RNA recurrence is reported separately
as `not_applicable`; it does not participate in the DNA final tier.
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
The table below reports all requested ETS outcomes after lookup recovery.
The complete breakend audit contains no ERG/FLI1 consensus SV in either direction.
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
isoforms and available Pfam annotations limit sensitivity; genuine missing mappings,
missing/invalid models, identity mismatches and alternate transcripts are separate
limitations (see the reason-count table above). The prior run lost **722 shared-family
rows across 52 SVs**, mostly PCDHA/UGT1A, to **ambiguous genomic breakpoints** in
the base screen. The corrected run has **{summary["ambiguous_breakpoint_shared_rows"]:,} such rows
across {summary["ambiguous_breakpoint_svs"]} SVs**. These are sensitivity losses from
overlapping target loci, not missing transcript lookups; base-branch breakpoint
selection logic is unchanged here. Of those corrected rows,
{categories.get("base_screen_ambiguous_breakpoints", 0)} retain ambiguity as their
primary reason; {overridden_ambiguous}
additional rows also fail the annotation identity/strand check (the original
ambiguity is retained in `pre_annotation_reason`). No
RNA confirmation, inserted-junction sequence, split reads, controlled VCF
confidence intervals, repeat/segdup track assessment or matched normal RNA
validation is available here. Sequence-only plausibility does not imply an
expressed, stable or functional fusion protein, nor establish oncogenicity.

The focused suite covers malformed inputs and all strand combinations, plus
batch-omission fallback, wrong-symbol responses, definitive single-gene absence,
network/server errors, DNA-only tier gating and missing-data reason categories.
Final repository checks and exact suite counts are recorded in
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
