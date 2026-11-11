# Sequence-level functional-domain checks for chimeric domains

`cfh.domain_function` (driver: `scripts/domain_function_seq.py`) asks, for every
in-frame chimeric domain the junction screen produced, whether there is
*sequence-level evidence* that the rebuilt domain could still work. It is the
first of three stacked steps: structure prediction of chimeras (PR 2) and
stability / protein-language-model scores (PR 3: PROVEAN, ESM, INDELi, FoldX)
consume the records written here. None of those are computed here.

**Nothing here is a calibrated predictor of activity.** The tier is a
transparent, documented heuristic over reported metrics; calibrating an
intact/broken cutoff would need assayed chimeras.

## Inputs and records

Input: the pilot's `combined_candidates.tsv.gz` (all five cohorts). Rows in the
four translated shared-domain classes are used: `clean_reconstitution`,
`tandem_duplication_or_insertion`, `gap_or_deletion`,
`ambiguous_junction_register`. `out_of_frame`, `premature_stop_at_junction` and
`splice_site_unresolved` have no in-frame chimeric domain.

Each junction is rebuilt with the screen's own `predict_junction` +
`classify_family` on the same cached transcript models, and the rebuilt class
must equal the screen's class (a mismatch is counted in the manifest, never
accepted). Rows whose chimeric domain is identical (same Pfam family, parent
transcripts and domain sequence with ±30 flanking residues) collapse into one
**record** (`chimera_records.json.gz`, `schema_version =
"domain_function_seq/1"`). A record carries: cohorts and per-cohort row counts,
fusion name(s), Pfam accession, 5'/3' parent gene, transcript, full protein,
parent domain span and sequence, junction residues, hybrid-codon residue, the
full translated fusion protein, the chimeric domain span (+ flanked sequence),
the screen's junction class/reason and indel (size, residues), read support
(max/total split reads, max paired-end reads, sample count) and the per-row
counts of every pilot artefact flag (`artefact_clean_rows` = rows where all six
flags pass). `cfh.domain_function.record.read_records_json` loads it and checks
the schema version.

## Metrics (`domain_function_metrics.tsv.gz`, one row per record)

### Profile-HMM integrity (HMMER via pyhmmer)

* The family HMM is fetched per family from the InterPro API
  (`/api/entry/pfam/<acc>/?annotation=hmm`) and cached; the Pfam/InterPro
  release is in `run_manifest.json`. No Pfam-A database is downloaded.
* The **translated fusion sequence** is scanned: the chimeric domain plus 30
  flanking residues of the actual fusion protein (not a pre-cut domain slice).
  Each parent is scanned the same way (parent domain ±30 residues).
* HMMER runs with filters off (`--max`), every hit reported, `Z = domZ = 1`
  (E-values are per single comparison; thresholds use bit scores).
* Reported per sequence (chimera / parent5 / parent3): best domain hit
  overlapping the domain (`*_domain_bits`, `*_i_evalue`, `*_ga_margin` = bits −
  domain GA, `*_passes_ga`, `*_match_coverage` = HMM match states spanned /
  model length, `*_hmm_span`, `*_ali_span` in protein coordinates), full
  sequence bits, `delta_bits_vs_5p/3p`, number of domain hits in the chimera,
  whether the chimera's hit spans the junction, and per sequence the number of hits
  overlapping the domain (`*_hits_in_domain`) with the union of their match-state
  coverage (`*_match_coverage_union`), so a domain HMMER splits around an indel is
  visible rather than judged on one fragment.
* **Junction residues in insert states**: junction residues are the 3 chimera
  residues each side of the junction plus every junction-inserted residue
  (below); `junction_residues_in_inserts` counts those HMMER aligned to insert
  states and `junction_residues_outside_alignment` those outside the hit.
* **HMMvar-style score** (Liu et al. 2014, BMC Bioinformatics 15:5,
  doi:10.1186/1471-2105-15-5): `hmmvar_delta_vs_5p/3p` = chimera minus parent
  HMMER full-sequence log-odds score (bits) of the scanned windows. It is the
  unnormalised log-likelihood-ratio difference against the profile; the windows
  of an indel chimera differ in length from the parent's.

Caveats. Passing GA supports *family membership*, not a correct fold or
activity. HMMER absorbs extra residues in insert states, so a junction
insertion can sit in insert states and still score well — hence the insert
count. A low score in a family whose parents also score low (`U1`) is
uninformative.

### Alignment model of the junction

The two parent windows are aligned (Biopython `PairwiseAligner`, BLOSUM62, gap
−11/−1, free end gaps). Chimera residues take their parent's column (the hybrid
residue takes the 3' column of residue F−1). A residue whose column is not
after the previous one is **junction-inserted** (a tandem copy of aligned
columns); columns skipped at the junction are **junction-deleted**
(`aligned_inserted_residues`, `aligned_inserted_sequence`,
`aligned_deleted_columns`, `aligned_deleted_5p/3p`). The screen's own indel
call is kept alongside (`pilot_indel_residues`, `pilot_indel_sequence`).

### Key residues (UniProtKB/Swiss-Prot)

Each parent gene maps to the reviewed human entry with that primary gene name.
Key features: Active site, Binding site (incl. metal), Site (excluding
translocation breakpoints), DNA binding, Motif, Region (excluding Disordered).
UniProt positions map to the parent transcript protein by global alignment;
below 90 % identity over the parent window → `not_assessed`. Each feature
overlapping the parent window gets its worst residue status: **retained**
(chimera keeps the parent residue), **swapped** (chimera has the other
parent's aligned residue — the expected outcome of a paralogue chimera),
**altered** (differs from both), **lost** (column junction-deleted).
`junction_to_key_5p/3p` = residues from the junction to the nearest key *site*
(point features and spans ≤ 10 residues; a whole-domain DNA-binding span would
make every junction "near"); `junction_near_key` = within 5.

### Secondary structure at the junction (AlphaFold DB)

For each parent's UniProt accession the AFDB model (`/api/prediction/<acc>`,
current version, fragment F1) is downloaded once; only a compact derived JSON is
cached and the mmCIF is deleted. Secondary structure is the file's own
`_struct_conf` category (AFDB writes DSSP there): H/G/I → helix, E → strand,
anything else (turn, bend, PPII, coil) → loop. A model without `_struct_conf`,
a missing model, or a model whose sequence differs from UniProt gives
`not_assessed` — never an assumed loop. Reported: DSSP context of 7 residues
around each parent's junction, `junction_ss_5p/3p`, the classes of the
junction-indel residues (`indel_ss`, `indel_in_ss`; `no_indel` for none) and the
minimum pLDDT of those residues.

### SCHEMA-like broken contacts

After Voigt et al. 2002 (Nat Struct Biol 9:553; see also Meyer et al. 2003,
Protein Sci 12:1686, doi:10.1110/ps.0306603). Contacts in each parent model:
any heavy-atom pair within 4.5 Å between residues ≥ 3 apart in sequence,
restricted to that parent's Pfam domain and to residues with pLDDT ≥ 70
(`schema_plddt_coverage_*` = fraction of domain residues kept). A contact is
broken (counts toward **E**) when the chimera's residue pair at those alignment
columns matches neither parent's pair. `schema_E_core_*` restricts to pairs of
buried residues (FreeSASA relative SASA < 0.20 in the full-length model).
Standard SCHEMA assumes aligned blocks, so junction indels are reported
separately: `schema_lost_to_deletion_*` (contacts with a residue in a deleted
column) and `schema_touching_inserted_*` (contacts of the source copy of an
inserted residue). E is reported against each parent structure; **there is no
universal E threshold** and E is not used for tiering.

## HEURISTIC tier

Rules (also in `run_manifest.json` → `tier_rules`, and per record in
`tier_rules` / `tier_missing_inputs`):

| rule | tier | condition |
|---|---|---|
| B1_below_GA | likely_broken | chimera best hit < domain GA (or no hit) while both parents ≥ GA |
| B2_low_coverage | likely_broken | chimera hits in the domain together span < 50 % of match states while both parents ≥ GA |
| B3_point_feature_lost | likely_broken | an Active site / Binding site / Site residue is in a deleted column |
| U1_parent_below_GA | uncertain | a parent scores below GA |
| U2_bit_drop | uncertain | chimera bits > 10 below the weaker parent |
| U3_partial_coverage | uncertain | chimera coverage < 80 % |
| U4_point_feature_altered | uncertain | a point key residue differs from both parents |
| U5_indel_in_ss | uncertain | a junction-indel residue is helix/strand in a parent model |
| U6_indel_near_key | uncertain | indel present and junction ≤ 5 residues from a key site |
| U7_junction_outside_hit | uncertain | the chimera's best hit does not span the junction |
| U8_domain_split | uncertain | HMMER splits the chimeric domain into more hits than either parent |
| N1_missing | not_assessed | not broken, no U rule fired, but an input needed for intact is missing |
| I1_all_pass | likely_intact | everything assessed and nothing above fired |

Broken rules fire on whatever is available; `likely_intact` requires every
input (HMM, key residues for both parents, indel secondary structure when there
is an indel). Thresholds (50 %, 80 %, 10 bits, 5 residues) are documented
heuristics, not calibrated values. `metric_disagreements` lists metric pairs
that point different ways (e.g. passes GA but indel in a helix; passes GA with
junction residues in insert states; below GA but key sites retained; intact but
broken core contacts) instead of hiding them.

## Shortlist for structure prediction

`shortlist_for_structure.tsv`: best record per fusion and family, records with
≥ 1 artefact-clean row only, ranked by tier (intact, uncertain, not_assessed,
broken), then lowest `schema_E_max` (max over the two parent structures;
missing last), then highest GA margin, then max split reads and sample count;
top 25, followed by the named cases (FLI1::ERG, STAT5B::STAT5A, GPC4::GPC3,
GPC6::GPC5, GPC3::GPC4, LPO::MPO, PIM3::ULK1) when they did not make the cut,
marked in `shortlist_reason`.
