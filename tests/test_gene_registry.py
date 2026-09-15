import pytest

from cfh.genes.registry import GeneConfig, KeyDomain, derive_gene_config_defaults, load_gene_config
from cfh.mapping.genome_nexus_source import CanonicalTranscript, ExonRecord, PfamDomain


def test_loads_braf_yaml_into_validated_gene_config():
    config = load_gene_config("braf")
    assert isinstance(config, GeneConfig)
    assert config.canonical_transcript_id == "NM_004333"
    assert config.protein_id == "P15056"
    assert "kinase" in config.key_domains[0].name.lower()


def test_loads_ret_yaml_with_live_genome_nexus_identifiers():
    config = load_gene_config("ret")

    assert isinstance(config, GeneConfig)
    assert config.gene_symbol == "RET"
    assert config.canonical_transcript_id == "NM_020975"
    assert config.protein_id == "P07949"
    assert config.entrez_gene_id == 5979
    assert config.key_domains[0].accession == "PF07714"
    assert config.benchmark_reference.fusion_count == 7


@pytest.mark.parametrize(
    ("gene", "transcript", "protein", "entrez"),
    [
        ("alk", "NM_004304", "Q9UM73", 238),
        ("ntrk1", "NM_002529", "P04629", 4914),
    ],
)
def test_loads_alk_and_ntrk1_yaml_with_genome_nexus_canonical_identifiers(
    gene, transcript, protein, entrez
):
    """Curated IDs are the live Genome Nexus canonical-transcript values.

    Both targets use the catalytic protein-kinase Pfam family returned by
    that endpoint; its live residue bounds are deliberately resolved by the
    production Genome Nexus mapping path rather than copied into YAML.
    """
    config = load_gene_config(gene)

    assert config.canonical_transcript_id == transcript
    assert config.protein_id == protein
    assert config.entrez_gene_id == entrez
    assert config.key_domains[0].accession == "PF07714"
    assert config.key_domains[0].source == "genome_nexus"


def test_unknown_gene_raises():
    with pytest.raises(FileNotFoundError):
        load_gene_config("not_a_real_gene")


def test_braf_disruption_required_domains_are_configured_with_real_pfam_accessions():
    """BRAF opts into the domain-disruption test with its real N-terminal
    autoinhibitory module (RAS-binding + cysteine-rich domains), sourced
    from the same live Genome Nexus data used for the kinase domain."""
    config = load_gene_config("braf")
    accessions = {domain.accession for domain in config.disruption_required_domains}
    assert accessions == {"PF02196", "PF00130"}
    assert all(domain.source == "genome_nexus" for domain in config.disruption_required_domains)


def test_ret_disruption_required_domains_are_configured_with_real_pfam_accessions():
    """RET opts into the domain-disruption test with its real N-terminal
    Cadherin domain, sourced from the same live Genome Nexus data used for
    the kinase domain.

    Before this fix, ``ret.yaml`` left ``disruption_required_domains``
    unset, so it was silently re-derived from a live Genome Nexus call on
    every run instead of being curated (see
    ``tests/test_cohort_auto_config.py`` and the domain_disruption
    benchmark tests for how that made RET's committed run artifacts
    nondeterministic run-to-run). Curating it here, like BRAF, removes that
    live dependency entirely.
    """
    config = load_gene_config("ret")
    accessions = {domain.accession for domain in config.disruption_required_domains}
    assert accessions == {"PF00028"}
    assert all(domain.source == "genome_nexus" for domain in config.disruption_required_domains)


def test_fgfr2_disruption_required_domains_uses_a_coordinate_override():
    """Regression: fgfr2.yaml's disruption_required_domains used to name
    "Immunoglobulin D1 domain" -- a name that never matches UniProt's own
    feature for that region ("Ig-like C2-type 1"), so domain_disruption
    failed outright on every real FGFR2 run
    (see cfh.mapping.feature_mapper's KeyDomain.start_aa docstring). It now
    targets the better-evidenced C-terminal exon-18 tail via an explicit
    coordinate override, which needs no live-source name match at all."""
    config = load_gene_config("fgfr2")
    assert len(config.disruption_required_domains) == 1
    domain = config.disruption_required_domains[0]
    assert domain.source == "literature"
    assert domain.start_aa == 768
    assert domain.end_aa == 821
    assert domain.accession is None


def test_alk_and_ntrk1_explicitly_opt_out_of_disruption_required_domains():
    """ALK and NTRK1 have real N-terminal extracellular domains preceding
    their kinase domain (verified live), but neither is curated as a
    ``disruption_required_domains`` entry -- an explicit empty list opts
    both genes out of domain_disruption for good, rather than leaving the
    field unset and letting it silently auto-derive (and potentially flip
    between runs) from whatever a live Genome Nexus call happens to return.
    """
    for gene in ("alk", "ntrk1"):
        config = load_gene_config(gene)
        assert config.disruption_required_domains == []
        assert "disruption_required_domains" in config.model_fields_set


def test_disruption_required_domains_defaults_to_empty_list():
    """Opt-in field: a config that never mentions it must gracefully default
    to no-op, not error, the same as the existing key_domains pattern."""
    config = GeneConfig(gene_symbol="FAKE", canonical_transcript_id="NM_1", protein_id="P1")
    assert config.disruption_required_domains == []


def test_key_domain_coordinate_override_defaults_to_unset():
    domain = KeyDomain(name="Some domain", source="test")
    assert domain.start_aa is None
    assert domain.end_aa is None


def test_key_domain_coordinate_override_accepts_a_matched_pair():
    domain = KeyDomain(name="C-terminal tail", source="literature", start_aa=768, end_aa=821)
    assert domain.start_aa == 768
    assert domain.end_aa == 821


@pytest.mark.parametrize(("start_aa", "end_aa"), [(768, None), (None, 821)])
def test_key_domain_coordinate_override_rejects_only_one_bound_set(start_aa, end_aa):
    with pytest.raises(ValueError, match="start_aa and end_aa must be set together"):
        KeyDomain(name="Bad", source="test", start_aa=start_aa, end_aa=end_aa)


def test_key_domain_coordinate_override_rejects_end_before_start():
    with pytest.raises(ValueError, match="end_aa must be >= start_aa"):
        KeyDomain(name="Bad", source="test", start_aa=821, end_aa=768)


def test_eml4_alk_and_a_synthetic_tmprss2_erg_pair_leave_disruption_domains_unset():
    """Negative controls: EML4-ALK (joint-partner mechanism) and TMPRSS2-ERG
    (promoter-swap/expression-driven mechanism) must not configure this
    field, so domain_disruption gracefully no-ops for both."""
    eml4_alk = load_gene_config("eml4-alk")
    assert eml4_alk.disruption_required_domains == []

    tmprss2_erg = GeneConfig(gene_pair=("TMPRSS2", "ERG"))
    assert tmprss2_erg.disruption_required_domains == []


def test_loads_tmprss2_erg_yaml_as_a_curated_gene_pair_config():
    """TMPRSS2-ERG is a second worked ``gene_pair`` example alongside
    EML4-ALK, but a genuinely different oncogenic mechanism: a promoter-swap/
    expression-driven fusion rather than a domain-retention story. The
    curated config must configure only ``joint_partner_dependency`` (no
    ``disruption_required_domains``/``key_domains``, which don't apply)."""
    config = load_gene_config("tmprss2-erg")

    assert isinstance(config, GeneConfig)
    assert config.gene_symbol is None
    assert config.gene_pair == ("TMPRSS2", "ERG")
    assert config.disruption_required_domains == []
    assert config.key_domains == []
    assert "joint_partner_dependency" in config.analysis_modes
    assert config.mechanism_note is not None
    assert "promoter-swap" in config.mechanism_note.lower()
    assert "not domain-retention" in config.mechanism_note.lower()


def test_loads_eml4_alk_yaml_with_curated_mechanism_note():
    """EML4-ALK's real config curates a domain-retention/ligand-independent-
    dimerization mechanism_note -- the same family of mechanism as the
    single-gene alk.yaml config, unlike TMPRSS2-ERG's promoter-swap note."""
    config = load_gene_config("eml4-alk")
    assert config.mechanism_note is not None
    assert "dimerization" in config.mechanism_note.lower()


def test_loads_erg_yaml_with_live_genome_nexus_identifiers():
    """ERG is TMPRSS2-ERG's curated single-gene component -- mirroring how
    EML4-ALK curates only ALK (its 3' partner) -- so the pair benchmark can
    live-fetch real structural-variant data for the pair (see
    ``_partner_component_configs``). Curating ERG's real Ets-domain here
    does not imply domain retention drives the TMPRSS2-ERG *pair*'s
    oncogenicity -- that config deliberately leaves key_domains unset."""
    config = load_gene_config("erg")

    assert isinstance(config, GeneConfig)
    assert config.gene_symbol == "ERG"
    assert config.canonical_transcript_id == "NM_001136154"
    assert config.protein_id == "P11308"
    assert config.entrez_gene_id == 2078
    assert config.key_domains[0].accession == "PF00178"
    assert "ets" in config.key_domains[0].name.lower()


def test_mechanism_note_defaults_to_none_and_is_opt_in():
    """Opt-in field: a config that never sets it must default to ``None``
    rather than requiring every config to supply one."""
    config = GeneConfig(gene_pair=("FAKE5", "FAKE6"), analysis_modes=["joint_partner_dependency"])
    assert config.mechanism_note is None


def test_mechanism_note_verified_defaults_to_true():
    """Every currently-curated gene's mechanism_note was human-cross-checked
    against primary literature -- defaulting to True means none of them
    need to set this explicitly to keep their "Curated mechanism"
    report/viewer treatment."""
    config = GeneConfig(
        gene_symbol="FAKE",
        canonical_transcript_id="NM_1",
        protein_id="P1",
        mechanism_note="Some curated text.",
    )
    assert config.mechanism_note_verified is True


def test_mechanism_note_verified_can_be_set_false_for_an_ai_drafted_note():
    config = GeneConfig(
        gene_symbol="FAKE",
        canonical_transcript_id="NM_1",
        protein_id="P1",
        mechanism_note="Some AI-drafted text, not yet independently verified.",
        mechanism_note_verified=False,
    )
    assert config.mechanism_note_verified is False


def test_derive_defaults_uses_most_n_terminal_key_domain_and_complete_pfam_list():
    config = GeneConfig(
        gene_symbol="FAKE",
        canonical_transcript_id="NM_1",
        protein_id="P1",
        key_domains=[
            KeyDomain(name="later", source="genome_nexus", accession="PF00003"),
            KeyDomain(name="earlier", source="genome_nexus", accession="PF00002"),
        ],
    )
    canonical = CanonicalTranscript(
        transcript_id="ENST1",
        refseq_mrna_id="NM_1",
        protein_id="P1",
        protein_length=150,
        uniprot_id=None,
        pfam_domains=[
            PfamDomain("PF00001", 10, 40),
            PfamDomain("PF00002", 60, 90),
            PfamDomain("PF00003", 110, 140),
        ],
        exons=[
            ExonRecord("e1", 1, 150, 1, 1),
            ExonRecord("e2", 151, 300, 2, 1),
            ExonRecord("e3", 301, 450, 3, 1),
        ],
        utrs=[],
    )

    derived = derive_gene_config_defaults(
        config,
        canonical,
        domain_name_resolver=lambda accession: {"PF00001": "Regulator"}[accession],
    )

    assert derived.expected_retained_exon_hint == "2"
    assert [domain.accession for domain in derived.disruption_required_domains] == ["PF00001"]
    assert derived.disruption_required_domains[0].name == "Regulator"


def test_derive_defaults_never_replaces_explicit_values():
    explicit_domain = KeyDomain(name="Curated", source="curator", accession="PF99999")
    config = GeneConfig(
        gene_symbol="FAKE",
        canonical_transcript_id="NM_1",
        protein_id="P1",
        key_domains=[KeyDomain(name="key", source="genome_nexus", accession="PF00002")],
        disruption_required_domains=[explicit_domain],
        expected_retained_exon_hint="exon 99",
    )
    canonical = CanonicalTranscript(
        transcript_id="ENST1",
        refseq_mrna_id="NM_1",
        protein_id="P1",
        protein_length=100,
        uniprot_id=None,
        pfam_domains=[PfamDomain("PF00001", 10, 20), PfamDomain("PF00002", 40, 80)],
        exons=[ExonRecord("e1", 1, 300, 1, 1)],
        utrs=[],
    )

    derived = derive_gene_config_defaults(config, canonical)

    assert derived.expected_retained_exon_hint == "exon 99"
    assert derived.disruption_required_domains == [explicit_domain]


def test_derive_exon_default_uses_closest_preceding_boundary_when_start_is_in_a_gap():
    config = GeneConfig(
        gene_symbol="FAKE",
        canonical_transcript_id="NM_1",
        protein_id="P1",
        key_domains=[KeyDomain(name="key", source="genome_nexus", accession="PF00001")],
    )
    canonical = CanonicalTranscript(
        transcript_id="ENST1",
        refseq_mrna_id="NM_1",
        protein_id="P1",
        protein_length=100,
        uniprot_id=None,
        pfam_domains=[PfamDomain("PF00001", 60, 80)],
        exons=[ExonRecord("e1", 1, 150, 7, 1)],
        utrs=[],
    )

    derived = derive_gene_config_defaults(config, canonical)

    assert derived.expected_retained_exon_hint == "7"


def test_derive_defaults_respects_explicit_empty_disruption_domain_list():
    config = GeneConfig.model_validate(
        {
            "gene_symbol": "FAKE",
            "canonical_transcript_id": "NM_1",
            "protein_id": "P1",
            "key_domains": [{"name": "key", "source": "genome_nexus", "accession": "PF00002"}],
            "disruption_required_domains": [],
        }
    )
    canonical = CanonicalTranscript(
        transcript_id="ENST1",
        refseq_mrna_id="NM_1",
        protein_id="P1",
        protein_length=100,
        uniprot_id=None,
        pfam_domains=[PfamDomain("PF00001", 10, 20), PfamDomain("PF00002", 40, 80)],
        exons=[],
        utrs=[],
    )

    assert derive_gene_config_defaults(config, canonical).disruption_required_domains == []
