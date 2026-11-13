# PCAWG TCGA consensus DNA SV domain-completion screen — 2026-11-13

This first genome-wide DNA screen examines **129,273 PASS consensus SVs**,
from **822 TCGA donor BEDPE files** (796 whitelisted, 26 greylisted). **789 donors**
have at least one SV; **33 files contain only a header**. The sample sheet contains
828 eligible TCGA donors (801 whitelist + 27 greylist); six have no archive file.
Empty files and absent files are different observations, not negative biological
results. See `sample_inventory.tsv` and `donors_without_bedpe.tsv`.

There are **0 clean reconstitutions**,
**16 non-trivial event-family rows** from **15 distinct SVs**, and
**5 distinct sequence records**. Overall candidate tiers
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

| Outcome | Count |
|---|---:|
| ambiguous | 8,209 |
| intergenic | 76,908 |
| no_sense_pair | 32,343 |
| sense_pairs | 11,813 |


These are SV-level categories; a call can yield multiple gene pairs.
**20,462 ordered pairs** enter the existing any-Pfam screen.
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
There are **21,004 outcome rows**, with
**1,247 annotated shared-family rows** from
**694 gene-pair events**;
**161 rows** reach artefact assessment.

The inherited `no_transcript_model` class includes both missing data and base-screen
product/locus verdicts. It must not be read as a count of missing models alone:

| Outcome | Count |
|---|---:|
| annotation_identity_or_strand_mismatch | 482 |
| base_screen_ambiguous_breakpoints | 807 |
| base_screen_no_coding_product | 84 |
| base_screen_outside_canonical_locus | 101 |
| missing_canonical_mapping | 6,833 |
| missing_or_invalid_transcript_model | 5 |


`missing_canonical_mapping` is a definitive single-gene lookup absence;
`missing_or_invalid_transcript_model` covers missing CDS/exon/valid model data.
Annotation ID/strand mismatches remain separate. Base-screen coding-product
verdicts include past-stop, outside-transcript and before-CDS breaks; outside-locus
and ambiguous-breakpoint decisions also remain separate from missing models.

## Lookup correction: before and after

| Measure | Before 539b788 | Corrected |
|---|---:|---:|
| no_transcript_model class rows (all reasons) | 9,907 | 8,312 |
| Missing canonical mappings | 8,792 | 6,833 |
| Annotated shared-family rows | 1,073 | 1,247 |
| Rows reaching artefact assessment | 149 | 161 |
| Non-trivial family rows | 13 | 16 |
| Sequence records | 4 | 5 |
| Clean reconstitutions | 0 | 0 |
| ETS non-trivial rows | 0 | 0 |
| Kinase non-trivial rows | 1 | 1 |

New non-trivial candidates after restoring single-gene fallback:

- KIAA1919::SLC16A10 / PF07690, DO2783: gap_or_deletion; sequence likely_broken; final DNA not_assessed.
- ARSE::ARSD / PF00884, DO44462: out_of_frame; sequence not_assessed; final DNA not_assessed.
- SERPINB11::SERPINB12 / PF00079, DO11854: out_of_frame; sequence not_assessed; final DNA not_assessed.

Baseline final DNA tiers: {'not_assessed': 13}; corrected: {'not_assessed': 16}.


Both runs examine the same 129,273 SVs and 20,462 ordered gene pairs. Recovery
of a partner model can introduce several Pfam family rows in place of one
unassessed placeholder, so family-row counts need not change one-for-one.

| Outcome | Count |
|---|---:|
| ambiguous_junction_register | 0 |
| clean_reconstitution | 0 |
| gap_or_deletion | 4 |
| no_shared_family | 12,531 |
| no_transcript_model | 8,312 |
| not_split_on_both | 145 |
| out_of_frame | 7 |
| premature_stop_at_junction | 0 |
| splice_site_unresolved | 4 |
| tandem_duplication_or_insertion | 1 |


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

| Outcome | Count |
|---|---:|
| not_assessed | 16 |


Sequence-only tiers on candidate rows (unscorable rows = not assessed):

| Outcome | Count |
|---|---:|
| likely_broken | 4 |
| not_assessed | 11 |
| uncertain | 1 |


| Fusion | Indel aa | Domain bits | GA margin | Best/union HMM coverage | Lost key features | Sequence tier |
|---|---:|---:|---:|---:|---:|---|
| ABCA6::ABCA10 | -184 | -2.03 | -33.63 | 0.012/0.012 | 0 | likely_broken |
| AKR1C2::AKR1C3 | -198 | 50.26 | 26.26 | 0.388/0.388 | 32 | likely_broken |
| KIAA1919::SLC16A10 | -229 | 28.13 | -5.27 | 0.233/0.233 | not_assessed | likely_broken |
| GNA14::GNAQ | -138 | 137.37 | 88.57 | 0.324/0.577 | 12 | likely_broken |
| MAP3K3::STRADA | 69 | 186.61 | 154.91 | 0.695/0.996 | 0 | uncertain |


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

### breakpoint_repeat_flag

| Outcome | Count |
|---|---:|
| not_assessed | 161 |

### genomic_similarity_flag

| Outcome | Count |
|---|---:|
| flagged | 85 |
| not_assessed | 76 |

### normal_recurrence_flag

| Outcome | Count |
|---|---:|
| not_applicable | 161 |

### paralog_sequence_flag

| Outcome | Count |
|---|---:|
| flagged | 12 |
| not_assessed | 4 |
| pass | 145 |

### read_evidence_flag

| Outcome | Count |
|---|---:|
| not_assessed | 161 |

### sample_burden_flag

| Outcome | Count |
|---|---:|
| flagged | 77 |
| not_assessed | 79 |
| pass | 5 |


## Strongest candidates for follow-up

Sorted by sequence-only functional tier, junction class, unambiguous pairing,
absence of a positive paralog flag, then caller support; no candidate is promoted through missing
validation data. Full rows are in `strongest_candidates.tsv` and every non-trivial
row, including failed/uncertain frames, is in `nontrivial_candidates.tsv`.

- **MAP3K3::STRADA / PF00069**, DO5626 (BRCA-US, 4/4 callers): tandem_duplication_or_insertion; sequence-only tier uncertain; paralog check pass, gene pairing ambiguous; overall **not_assessed** (DNA read/repeat validation incomplete).
- **GNA14::GNAQ / PF00503**, DO2341 (BRCA-US, 4/4 callers): gap_or_deletion; sequence-only tier likely_broken; paralog check pass, gene pairing unique; overall **not_assessed** (DNA read/repeat validation incomplete).
- **KIAA1919::SLC16A10 / PF07690**, DO2783 (BRCA-US, 4/4 callers): gap_or_deletion; sequence-only tier likely_broken; paralog check pass, gene pairing unique; overall **not_assessed** (DNA read/repeat validation incomplete).
- **ABCA6::ABCA10 / PF12698**, DO52171 (LIHC-US, 3/4 callers): gap_or_deletion; sequence-only tier likely_broken; paralog check pass, gene pairing unique; overall **not_assessed** (DNA read/repeat validation incomplete).
- **AKR1C2::AKR1C3 / PF00248**, DO20928 (KIRP-US, 3/4 callers): gap_or_deletion; sequence-only tier likely_broken; paralog check pass, gene pairing unique; overall **not_assessed** (DNA read/repeat validation incomplete).
- **GPC6::GPC5 / PF01153**, DO41157 (UCEC-US, 4/4 callers): out_of_frame; sequence-only tier not_assessed; paralog check pass, gene pairing unique; overall **not_assessed** (DNA read/repeat validation incomplete).
- **NFASC::NTRK1 / PF13895**, DO12484 (GBM-US, 4/4 callers): splice_site_unresolved; sequence-only tier not_assessed; paralog check not_assessed, gene pairing unique; overall **not_assessed** (DNA read/repeat validation incomplete).
- **SERPINB11::SERPINB12 / PF00079**, DO11854 (GBM-US, 4/4 callers): out_of_frame; sequence-only tier not_assessed; paralog check pass, gene pairing unique; overall **not_assessed** (DNA read/repeat validation incomplete).
- **ZNF69::ZNF844 / PF01352**, DO220140 (SARC-US, 4/4 callers): splice_site_unresolved; sequence-only tier not_assessed; paralog check not_assessed, gene pairing unique; overall **not_assessed** (DNA read/repeat validation incomplete).
- **LMNB2::KRT85 / PF00038**, DO219964 (SARC-US, 3/4 callers): out_of_frame; sequence-only tier not_assessed; paralog check pass, gene pairing unique; overall **not_assessed** (DNA read/repeat validation incomplete).

## Explicit ETS-family examination

`ets_breakend_audit.tsv` includes all SVs overlapping any requested ETS gene,
including intergenic/antisense/intragenic cases.
The table below reports all requested ETS outcomes after lookup recovery.
The complete breakend audit contains no ERG/FLI1 consensus SV in either direction.
`ets_sense_pairs.tsv` includes
all compatible sense gene pairs and their screen outcome. This DNA approach
can examine mid-exon ETS breaks that RNA exon-boundary fusion tables omit.

| Gene | SVs overlapping gene | Sense pairs | Non-trivial family rows |
|---|---:|---:|---:|
| ERG | 72 | 25 | 0 |
| FLI1 | 7 | 0 | 0 |
| ETV1 | 7 | 2 | 0 |
| ETV4 | 7 | 2 | 0 |
| ETV5 | 4 | 1 | 0 |
| GABPA | 5 | 0 | 0 |
| ETS1 | 13 | 5 | 0 |
| ETS2 | 1 | 1 | 0 |
| FEV | 0 | 0 | 0 |


## Kinase pairs

Shared kinase-family rows are defined explicitly as Pfam **PF00069** (protein
kinase) or **PF07714** (tyrosine kinase), not by gene-name heuristics.
There are **29 shared kinase-family rows**, **1 non-trivial**.
`kinase_shared_family.tsv` contains every such pair with donor and breakpoints;
non-trivial rows also appear in the complete candidate table. This definition
does not assert exhaustive coverage of atypical kinase families or pairs assigned
to different Pfam families.

| Outcome | Count |
|---|---:|
| no_transcript_model | 21 |
| not_split_on_both | 7 |
| tandem_duplication_or_insertion | 1 |


## Limits and validation

This is a PASS-only consensus screen. Homologous regions are likely under-called;
absence of a consensus SV is not evidence that a fusion cannot occur. Canonical
isoforms and available Pfam annotations limit sensitivity; genuine missing mappings,
missing/invalid models, identity mismatches and alternate transcripts are separate
limitations (see the reason-count table above). The prior run lost **722 shared-family
rows across 52 SVs**, mostly PCDHA/UGT1A, to **ambiguous genomic breakpoints** in
the base screen. The corrected run has **859 such rows
across 74 SVs**. These are sensitivity losses from
overlapping target loci, not missing transcript lookups; base-branch breakpoint
selection logic is unchanged here. Of those corrected rows,
807 retain ambiguity as their
primary reason; 52
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
