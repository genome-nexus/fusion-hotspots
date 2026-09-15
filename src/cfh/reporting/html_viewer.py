"""Backend-less static-site HTML viewer for an already-written ``cfh`` run
directory (Milestone 4, "Option A"): plain HTML/CSS/vanilla JS that reads
the already-generated ``results.json``/``summary.json`` and SVG artifacts
and re-renders them into one or more self-contained ``.html`` files under
``<run_directory>/viewer/``.

This module is purely a *presentation* layer over artifacts some other part
of the pipeline already wrote to disk -- it re-reads ``results.json``/
``summary.json``/``*.svg`` exactly as committed and never recomputes,
re-derives, or overrides a single statistic. Nothing here can change a
run's numbers. The one exception is the curated-vs-auto-configured badge
(see :func:`_config_source_for_gene`): it checks whether a curated YAML
exists for the gene, which is metadata about provenance, not a statistic.

No server and no JS build toolchain: every page is one flat ``.html`` file
with its CSS/JS inlined and every SVG figure inlined directly into the page
markup (not linked via ``<img src="...svg">`` or fetched with
``fetch()``) so it opens correctly straight off disk via a ``file://`` URL
-- ``fetch()`` of a sibling file is blocked by CORS under ``file://`` in
some browsers, but an inlined ``<svg>`` element (and inlined
``<script type="application/json">`` data) is just part of the page's own
DOM and always works.

Two run shapes are recognized (see :func:`build_run_viewer`):

* A single-gene run directory (``cfh real-benchmark``/``cfh analyze``
  output): one gene page, written to ``<run_directory>/viewer/index.html``.
* A cohort-scan run directory (``cfh cohort-scan`` output, identified by
  the presence of ``cohort_scan/summary.json``): a landing page at
  ``<run_directory>/viewer/index.html`` embedding the genome-wide
  Manhattan plot with a click-through handler, plus one gene page per
  ``cohort_scan/gene_reports/<gene>/`` entry under
  ``<run_directory>/viewer/genes/<gene>.html``.

A single-gene run may itself be one of two page shapes: an ordinary gene
page (domain-retention/fusion-schematic figures, a per-event table) or, for
a ``gene_pair`` config (e.g. EML4-ALK, TMPRSS2-ERG), a joint-partner page --
observed/expected/enrichment-p-value, not domain retention -- rendered by a
dedicated template (:func:`_render_gene_pair_page`) rather than force-fit
into the single-gene layout.

Several sections are optional and gated on the underlying data actually
being present in a given run's ``results.json`` (or, for cross-cohort
concordance, a sibling ``cross_cohort_concordance.json`` next to it):
the ``composite_score`` evidence-ranking table, the genomic-coordinate
clustering view (``genomic_position_recurrence``), and the cross-cohort CMH
panel (``cross_cohort_concordance``). Many genes/older runs will legitimately
lack one or more of these -- the corresponding section is simply omitted,
never rendered empty or broken.
"""

from __future__ import annotations

import json
from html import escape as _esc
from pathlib import Path
from typing import Any

from cfh.genes.registry import load_gene_config

# Mirrors cfh.reporting.fusion_schematic._status_word's display language so
# the viewer's toggle uses the same words as the schematic it controls.
_DOMAIN_STATUS_LABELS = {
    "retained": "Retained",
    "disrupted": "Truncated",
    "lost": "Lost",
}

_EVENT_FIELDS = (
    "event_id",
    "sample_id",
    "tumor_type",
    "oncotree_code",
    "partner_gene",
    "target_role",
    "breakpoint_protein_position",
    "domain_status",
)

_COMPOSITE_SCORE_COLUMNS = (
    ("Rank", "Rank"),
    ("Partner_gene", "Partner gene"),
    ("Event_count", "Events"),
    ("Composite_score", "Composite score"),
    ("Recurrence_score", "Recurrence"),
    ("Domain_retention_score", "Domain retention"),
    ("Domain_disruption_score", "Domain disruption"),
    ("Cutpoint_proximity_score", "Cutpoint proximity"),
    ("Confidence_certainty_score", "Confidence"),
)


def build_run_viewer(run_directory: Path | str) -> Path | None:
    """Build the static-site viewer bundle for ``run_directory`` (a
    top-level ``runs/<run_id>/`` directory), writing it to
    ``<run_directory>/viewer/``.

    Returns the path to the written ``viewer/index.html``, or ``None`` when
    ``run_directory`` has neither a single-gene ``results.json`` nor a
    cohort-scan ``cohort_scan/summary.json`` to build from (nothing is
    written in that case).
    """
    run_directory = Path(run_directory)
    if (run_directory / "cohort_scan" / "summary.json").exists():
        return _build_cohort_scan_viewer(run_directory)
    if (run_directory / "results.json").exists():
        return _build_single_gene_viewer(
            run_directory, run_directory / "viewer" / "index.html", back_link=None
        )
    return None


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def _read_json_or_none(path: Path) -> Any | None:
    """Like :func:`_read_json`, but tolerant of a missing or unparsable
    sibling artifact -- an optional data file (e.g.
    ``cross_cohort_concordance.json``) must never take down viewer
    generation for the run it sits next to."""
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _read_text_or_none(path: Path) -> str | None:
    return path.read_text() if path.exists() else None


def _config_source_for_gene(gene_symbol: str) -> str | None:
    """Whether ``gene_symbol`` resolves to a hand-curated YAML under
    ``genes/configs/`` ("curated") or would be auto-derived live from
    Genome Nexus ("auto") -- mirrors ``cfh.cohort.scan._resolve_configs``'s
    own curated-first resolution order exactly (a successful
    :func:`~cfh.genes.registry.load_gene_config` load is the definitive
    "curated" signal there too). Never raises: an unexpected error here
    must not take down the rest of the page."""
    if not gene_symbol:
        return None
    try:
        load_gene_config(gene_symbol)
        return "curated"
    except FileNotFoundError:
        return "auto"
    except Exception:  # noqa: BLE001 - a badge must never crash viewer generation
        return None


def _config_source_by_gene(summary: dict) -> dict[str, str]:
    """Per-gene ``config_source`` ("curated"/"auto"/"unresolved") already
    computed and recorded by a cohort scan, keyed by upper-cased gene
    symbol -- reused as-is rather than re-derived, since the scan already
    knows definitively (including "unresolved" genes that never even got a
    ``GeneConfig``)."""
    return {
        str(row["gene_symbol"]).upper(): row["config_source"]
        for row in summary.get("genes") or []
        if isinstance(row, dict) and row.get("gene_symbol") and row.get("config_source")
    }


def _config_source_badge(config_source: str | None) -> str:
    if config_source == "curated":
        return (
            '<span class="badge badge-curated" '
            'title="Loaded from a hand-curated YAML under genes/configs/">'
            "Curated config</span>"
        )
    if config_source == "auto":
        return (
            '<span class="badge badge-auto" '
            'title="Derived live from Genome Nexus, no curated YAML on file">'
            "Auto-configured</span>"
        )
    return ""


def _algorithm_result(payload: dict, algorithm_name: str) -> dict[str, Any] | None:
    """Find one named entry in ``payload["algorithm_results"]`` (a flat
    list of ``AlgorithmResult``-shaped dicts, one per registered algorithm
    that ran for this gene) -- or ``None`` if that algorithm didn't run for
    this gene/run at all, which is the common, expected case for most of
    the optional sections below."""
    for item in payload.get("algorithm_results") or []:
        if isinstance(item, dict) and item.get("Algorithm") == algorithm_name:
            return item
    return None


def _algorithm_failed(result: dict[str, Any]) -> bool:
    """Mirrors ``cfh.algorithms.composite_score._failed``'s orchestrator-
    wrapped-exception convention, read back from the plain JSON dict
    instead of a reconstructed ``AlgorithmResult``."""
    return any(
        str(warning).startswith("Algorithm failed") for warning in (result.get("Warnings") or [])
    )


def _cross_cohort_data(payload: dict, gene_run_dir: Path | None) -> dict[str, Any] | None:
    """A gene's cross-cohort CMH concordance report, if one exists for this
    run: either a ``cross_cohort_concordance`` entry in
    ``algorithm_results`` (forward-compatible, should the orchestrator ever
    embed it there directly) or, as produced today by
    ``cfh.algorithms.cross_cohort_concordance.compare_cohort_runs``/
    ``cfh compare-cohorts``, a sibling ``cross_cohort_concordance.json``
    committed next to this gene's ``results.json``. Neither is required --
    most genes/runs legitimately have neither, and this returns ``None``
    for both."""
    result = _algorithm_result(payload, "cross_cohort_concordance")
    if result is not None and not _algorithm_failed(result):
        data = result.get("Summary") or result.get("Tables")
        if isinstance(data, dict) and data.get("cmh"):
            return data
    if gene_run_dir is not None:
        sibling = _read_json_or_none(gene_run_dir / "cross_cohort_concordance.json")
        if isinstance(sibling, dict) and sibling.get("cmh"):
            return sibling
    return None


def _event_records(payload: dict) -> list[dict[str, Any]]:
    return [
        {field: event.get(field) for field in _EVENT_FIELDS}
        for event in payload.get("events") or []
    ]


def _scalar_summary_items(summary: dict) -> list[tuple[str, Any]]:
    """The subset of a run's ``summary`` dict worth showing as a simple
    stat grid: scalar values only (partner-count tables, contingency
    tables, and other structured fields already have a dedicated home in
    report.md/report.pdf). Generic over every run type -- no gene-specific
    field names are hardcoded here."""
    return [
        (key, value)
        for key, value in summary.items()
        if value is not None and not isinstance(value, (list, dict))
    ]


def _format_summary_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def _humanize(key: str) -> str:
    return key.replace("_", " ").title()


_PAGE_STYLE = """
:root { color-scheme: light; }
* { box-sizing: border-box; }
body {
  margin: 0; padding: 24px 20px 60px; background: #f7f7f8; color: #1a1a1a;
  font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
}
.wrap { max-width: 1080px; margin: 0 auto; }
h1 { font-size: 22px; margin: 0 0 4px; }
h1 .badge { position: relative; top: -2px; }
h2 {
  font-size: 16px; margin: 32px 0 10px;
  border-bottom: 1px solid #ddd; padding-bottom: 6px;
}
.subtitle { color: #555; font-size: 13px; margin: 0 0 20px; }
.back-link {
  display: inline-block; margin-bottom: 14px; font-size: 13px;
  color: #2878b5; text-decoration: none;
}
.back-link:hover { text-decoration: underline; }
.card {
  background: #fff; border: 1px solid #e2e2e2;
  border-radius: 8px; padding: 16px 18px;
}
.stat-grid {
  display: grid; grid-template-columns: repeat(auto-fill, minmax(170px, 1fr)); gap: 10px;
}
.stat {
  background: #fff; border: 1px solid #e2e2e2;
  border-radius: 8px; padding: 10px 12px;
}
.stat .label {
  font-size: 11px; color: #666; text-transform: uppercase; letter-spacing: 0.02em;
}
.stat .value { font-size: 18px; font-weight: 600; margin-top: 2px; }
.svg-frame {
  width: 100%; overflow-x: auto; background: #fff;
  border: 1px solid #e2e2e2; border-radius: 8px; padding: 8px;
}
.svg-frame svg { display: block; max-width: 100%; height: auto; }
.controls {
  display: flex; flex-wrap: wrap; gap: 20px; align-items: center;
  margin: 12px 0 16px; padding: 12px 14px; background: #fff;
  border: 1px solid #e2e2e2; border-radius: 8px;
}
.control-group { display: flex; align-items: center; gap: 8px; flex-wrap: wrap; }
.control-group label.title { font-size: 12px; font-weight: 600; color: #333; }
.control-group select { font-size: 13px; padding: 4px 6px; }
.status-toggle { display: inline-flex; align-items: center; gap: 4px; font-size: 13px; }
#events-summary { font-size: 12px; color: #555; margin-left: auto; }
table.events { width: 100%; border-collapse: collapse; font-size: 12.5px; }
table.events th, table.events td {
  text-align: left; padding: 5px 8px;
  border-bottom: 1px solid #eee; white-space: nowrap;
}
table.events th { position: sticky; top: 0; background: #fafafa; }
table.events tr.row-highlight { background: #fff3cd; transition: background 1.2s ease; }
.lollipop-tooltip {
  position: fixed; z-index: 1000; pointer-events: none;
  background: #1a1a1a; color: #fff; font-size: 12px; line-height: 1.3;
  padding: 4px 8px; border-radius: 4px; white-space: nowrap;
}
.events-table-frame {
  max-height: 420px; overflow: auto;
  border: 1px solid #e2e2e2; border-radius: 8px; background: #fff;
}
.gene-links { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 14px; }
.gene-links a {
  font-size: 12.5px; padding: 4px 9px; border-radius: 999px; text-decoration: none;
  border: 1px solid #ccc; color: #1a1a1a; background: #fff;
}
.gene-links a:hover { background: #f0f0f0; }
.gene-links a.significant { border-color: #d62728; color: #d62728; }
.note { font-size: 12px; color: #777; margin-top: 8px; }
.panel-note { font-size: 12px; color: #555; margin: 8px 0 0; padding-left: 18px; }
.badge {
  display: inline-block; font-size: 10.5px; font-weight: 600; padding: 2px 9px;
  border-radius: 999px; margin-left: 10px; vertical-align: middle;
  text-transform: uppercase; letter-spacing: 0.03em;
}
.badge-curated { background: #e6f4ea; color: #1a7f37; border: 1px solid #b7e0c2; }
.badge-auto { background: #fff4e5; color: #9a6700; border: 1px solid #f0dca6; }
table.data-table { width: 100%; border-collapse: collapse; font-size: 12.5px; }
table.data-table th, table.data-table td {
  text-align: left; padding: 5px 8px;
  border-bottom: 1px solid #eee; white-space: nowrap;
}
table.data-table th { position: sticky; top: 0; background: #fafafa; }
table.sortable th[data-sort-key] { cursor: pointer; user-select: none; }
table.sortable th[data-sort-key]:hover { background: #f0f0f0; }
table.sortable th[data-sort-key]::after { content: ""; margin-left: 4px; color: #999; }
table.sortable th[data-sort-dir="asc"]::after { content: "\\25B2"; margin-left: 4px; color: #555; }
table.sortable th[data-sort-dir="desc"]::after { content: "\\25BC"; margin-left: 4px; color: #555; }
.data-table-frame {
  max-height: 380px; overflow: auto;
  border: 1px solid #e2e2e2; border-radius: 8px; background: #fff;
}
.two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
@media (max-width: 720px) { .two-col { grid-template-columns: 1fr; } }
"""


def _landing_stat_grid(summary: dict) -> str:
    items = [
        ("Study", summary.get("study_id")),
        (
            "Genes after gating",
            f"{summary.get('genes_after_gating')} / {summary.get('total_genes_before_gating')}",
        ),
        ("FDR-significant genes", len(summary.get("significant_genes") or [])),
        ("Significance level (q)", summary.get("significance_level")),
        ("Curated gene configs", summary.get("curated_gene_count")),
        ("Auto-configured genes", summary.get("auto_config_gene_count")),
    ]
    cells = "".join(
        f'<div class="stat"><div class="label">{_esc(str(label))}</div>'
        f'<div class="value">{_esc(_format_summary_value(value))}</div></div>'
        for label, value in items
        if value is not None
    )
    return f'<div class="stat-grid">{cells}</div>'


def _gene_link_list(summary: dict, available_genes: set[str]) -> str:
    significant = set(summary.get("significant_genes") or [])
    links = []
    for gene_symbol in sorted(available_genes, key=str.upper):
        css_class = " significant" if gene_symbol.upper() in significant else ""
        links.append(
            f'<a class="gene-link{css_class}" href="genes/{_esc(gene_symbol)}.html">'
            f"{_esc(gene_symbol.upper())}</a>"
        )
    return f'<div class="gene-links">{"".join(links)}</div>'


_LANDING_SCRIPT = """
document.addEventListener("DOMContentLoaded", function () {
  var available = new Set(JSON.parse(document.getElementById("available-genes").textContent));
  document.querySelectorAll("svg circle[data-gene]").forEach(function (circle) {
    var gene = (circle.getAttribute("data-gene") || "").toLowerCase();
    if (!available.has(gene)) { return; }
    circle.style.cursor = "pointer";
    circle.addEventListener("click", function () {
      window.location.href = "genes/" + gene + ".html";
    });
  });
});
"""


def _render_landing_page(summary: dict, manhattan_svg: str, *, available_genes: set[str]) -> str:
    study_id = _esc(str(summary.get("study_id", "")))
    available_genes_json = _safe_json_script(sorted(available_genes))
    subtitle = (
        f"{study_id} &mdash; click a significant gene's point below (or use the list) "
        "to open its full per-gene view."
    )
    clickable_note = (
        "Points are clickable only for genes with a full per-gene report generated "
        "(FDR-significant, curated, or an honorable mention); every other scanned gene "
        "still carries a hover tooltip on its point, but is left as plain (non-clickable) "
        "text/point rather than a dead link."
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{study_id} fusion-hotspot cohort scan</title>
<style>{_PAGE_STYLE}</style>
</head>
<body>
<div class="wrap">
  <h1>Genome-wide fusion-hotspot cohort scan</h1>
  <p class="subtitle">{subtitle}</p>
  {_landing_stat_grid(summary)}
  <h2>Genome-wide summary</h2>
  <div class="svg-frame">{manhattan_svg}</div>
  <p class="note">{clickable_note}</p>
  <h2>Genes with a full report ({len(available_genes)})</h2>
  {_gene_link_list(summary, available_genes)}
</div>
<script type="application/json" id="available-genes">{available_genes_json}</script>
<script>{_LANDING_SCRIPT}</script>
</body>
</html>
"""


def _safe_json_script(data: Any) -> str:
    """JSON-serialize ``data`` for safe embedding inside an inline
    ``<script type="application/json">`` block -- escaping ``</`` so a
    value containing a literal ``</script>`` cannot prematurely close the
    tag (``"\\/"`` is a valid JSON escape for ``/``, so this round-trips
    through ``JSON.parse`` unchanged)."""
    return json.dumps(data, allow_nan=False).replace("</", "<\\/")


def _build_cohort_scan_viewer(run_directory: Path) -> Path:
    cohort_dir = run_directory / "cohort_scan"
    summary = _read_json(cohort_dir / "summary.json")
    manhattan_svg = _read_text_or_none(cohort_dir / "manhattan.svg") or ""
    config_source_map = _config_source_by_gene(summary)

    gene_reports_dir = cohort_dir / "gene_reports"
    available_genes = sorted(
        entry.name
        for entry in (gene_reports_dir.iterdir() if gene_reports_dir.exists() else [])
        if entry.is_dir() and (entry / "results.json").exists()
    )

    viewer_dir = run_directory / "viewer"
    genes_dir = viewer_dir / "genes"
    genes_dir.mkdir(parents=True, exist_ok=True)
    for gene_lower in available_genes:
        _build_single_gene_viewer(
            gene_reports_dir / gene_lower,
            genes_dir / f"{gene_lower}.html",
            back_link="../index.html",
            config_source=config_source_map.get(gene_lower.upper()),
        )

    index_path = viewer_dir / "index.html"
    index_path.write_text(
        _render_landing_page(summary, manhattan_svg, available_genes=set(available_genes))
    )
    return index_path


def _is_gene_pair_payload(payload: dict) -> bool:
    return bool((payload.get("summary") or {}).get("gene_pair"))


def _build_single_gene_viewer(
    gene_run_dir: Path,
    output_path: Path,
    *,
    back_link: str | None,
    config_source: str | None = None,
) -> Path:
    payload = _read_json(gene_run_dir / "results.json")
    domain_svg = _read_text_or_none(
        gene_run_dir / "visualizations" / "domain_retention_outliers.svg"
    )
    fusion_svg = _read_text_or_none(gene_run_dir / "visualizations" / "fusion_schematic.svg")
    resolved_config_source = config_source or _config_source_for_gene(
        str(payload.get("gene_symbol", ""))
    )
    cross_cohort = _cross_cohort_data(payload, gene_run_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if _is_gene_pair_payload(payload):
        html = _render_gene_pair_page(
            payload,
            back_link=back_link,
            config_source=resolved_config_source,
            cross_cohort=cross_cohort,
        )
    else:
        html = _render_gene_page(
            payload,
            domain_svg=domain_svg,
            fusion_svg=fusion_svg,
            back_link=back_link,
            config_source=resolved_config_source,
            cross_cohort=cross_cohort,
        )
    output_path.write_text(html)
    return output_path


_GENE_PAGE_SCRIPT = r"""
document.addEventListener("DOMContentLoaded", function () {
  var events = JSON.parse(document.getElementById("events-data").textContent);
  var statusLabels = JSON.parse(document.getElementById("status-labels").textContent);
  var eventsById = {};
  events.forEach(function (e) {
    if (e.event_id) { eventsById[e.event_id] = e; }
  });

  function distinctValues(field) {
    var values = events.map(function (e) { return e[field]; }).filter(Boolean);
    return Array.from(new Set(values)).sort();
  }
  var tumorTypes = distinctValues("tumor_type");
  var domainStatuses = distinctValues("domain_status");

  var tumorSelect = document.getElementById("tumor-type-filter");
  tumorSelect.innerHTML = '<option value="">All tumor types</option>' + tumorTypes.map(
    function (t) {
      return '<option value="' + t.replace(/"/g, "&quot;") + '">' + t + "</option>";
    }
  ).join("");

  var statusContainer = document.getElementById("status-filter");
  statusContainer.innerHTML = domainStatuses.map(function (s) {
    var label = statusLabels[s] || (s.charAt(0).toUpperCase() + s.slice(1));
    var checkbox = '<input type="checkbox" value="' + s + '" checked> ';
    return '<label class="status-toggle">' + checkbox + label + "</label>";
  }).join("");

  // Group fusion-schematic rows by the (partner, breakpoint-aa) key baked
  // into every row shape's <title> text by cfh.reporting.fusion_schematic
  // (e.g. "...BRAF–KIAA1549 fusion, breakpoint BRAF aa 381; ...") --
  // read back here, never recomputed or duplicated from the algorithm side.
  var schematicGroups = new Map();
  events.forEach(function (e) {
    if (!e.partner_gene || e.breakpoint_protein_position == null) { return; }
    var key = e.partner_gene + "|" + e.breakpoint_protein_position;
    var group = schematicGroups.get(key);
    if (!group) {
      group = { tumorTypes: new Set(), domainStatuses: new Set() };
      schematicGroups.set(key, group);
    }
    if (e.tumor_type) { group.tumorTypes.add(e.tumor_type); }
    if (e.domain_status) { group.domainStatuses.add(e.domain_status); }
  });
  var TITLE_RE = /–([^;]+?) fusion, breakpoint \S+ aa (\d+);/;
  var schematicSvg = document.getElementById("fusion-schematic-svg");

  var summaryEl = document.getElementById("events-summary");
  var rows = Array.from(document.querySelectorAll("#events-table tbody tr"));

  function applyFilters() {
    var tumor = tumorSelect.value;
    var checked = statusContainer.querySelectorAll("input:checked");
    var activeStatuses = new Set(Array.from(checked).map(function (i) { return i.value; }));

    var shown = 0;
    rows.forEach(function (tr) {
      var matchesTumor = !tumor || tr.dataset.tumorType === tumor;
      var matchesStatus = !tr.dataset.domainStatus || activeStatuses.has(tr.dataset.domainStatus);
      var visible = matchesTumor && matchesStatus;
      tr.hidden = !visible;
      if (visible) { shown += 1; }
    });
    if (summaryEl) {
      summaryEl.textContent = "Showing " + shown + " of " + rows.length + " events";
    }

    if (schematicSvg) {
      schematicSvg.querySelectorAll("title").forEach(function (titleEl) {
        var match = TITLE_RE.exec(titleEl.textContent || "");
        var shape = titleEl.parentElement;
        if (!match) { return; }
        var group = schematicGroups.get(match[1] + "|" + match[2]);
        if (!group) { shape.style.opacity = ""; return; }
        var groupMatchesTumor = !tumor || group.tumorTypes.has(tumor);
        var groupMatchesStatus = Array.from(group.domainStatuses).some(function (s) {
          return activeStatuses.has(s);
        });
        shape.style.opacity = (groupMatchesTumor && groupMatchesStatus) ? "1" : "0.12";
      });
    }
  }

  tumorSelect.addEventListener("change", applyFilters);
  statusContainer.addEventListener("change", applyFilters);
  applyFilters();

  // Hovering or clicking a lollipop-track point (data-event-id, set by
  // cfh.real_benchmark._domain_track_svg) shows a sample-id tooltip and/or
  // jumps to that event's row in the table below -- clearing any active
  // filters first so the row can't be hidden, since the point itself is
  // drawn independent of the table's tumor-type/domain-status filters.
  var domainTrackFrame = document.getElementById("domain-track-svg");
  var tooltip = document.getElementById("lollipop-tooltip");
  if (domainTrackFrame) {
    domainTrackFrame.querySelectorAll("circle[data-event-id]").forEach(function (circle) {
      var eventId = circle.getAttribute("data-event-id");
      if (!eventId) { return; }
      var matchedEvent = eventsById[eventId];
      if (tooltip && matchedEvent && matchedEvent.sample_id) {
        circle.addEventListener("mouseenter", function () {
          tooltip.textContent = "Sample " + matchedEvent.sample_id;
          tooltip.hidden = false;
        });
        circle.addEventListener("mousemove", function (mouseEvent) {
          tooltip.style.left = mouseEvent.clientX + 14 + "px";
          tooltip.style.top = mouseEvent.clientY + 14 + "px";
        });
        circle.addEventListener("mouseleave", function () {
          tooltip.hidden = true;
        });
      }
      var row = document.querySelector(
        '#events-table tbody tr[data-event-id="' + CSS.escape(eventId) + '"]'
      );
      if (!row) { return; }
      circle.style.cursor = "pointer";
      circle.addEventListener("click", function () {
        tumorSelect.value = "";
        statusContainer.querySelectorAll("input").forEach(function (i) { i.checked = true; });
        applyFilters();
        row.scrollIntoView({ behavior: "smooth", block: "center" });
        row.classList.add("row-highlight");
        window.setTimeout(function () { row.classList.remove("row-highlight"); }, 1500);
      });
    });
  }
});
"""

# Attaches click-to-sort behavior to every `table.sortable` on the page
# (currently just the composite-score evidence-ranking table). Sorting is
# purely a DOM re-order of the existing <tr> elements -- it never touches
# the embedded JSON data, and each numeric cell carries its raw value in a
# `data-sort` attribute so sorting is correct even when the visible text is
# formatted (e.g. "0.3337" vs. the 4-significant-figure display string).
_SORTABLE_TABLE_SCRIPT = r"""
document.addEventListener("DOMContentLoaded", function () {
  document.querySelectorAll("table.sortable").forEach(function (table) {
    var tbody = table.tBodies[0];
    if (!tbody) { return; }
    var headers = Array.from(table.querySelectorAll("thead th[data-sort-key]"));
    var currentKey = null;
    var currentDir = 1;
    headers.forEach(function (th, index) {
      th.addEventListener("click", function () {
        var key = th.getAttribute("data-sort-key");
        currentDir = (currentKey === key) ? -currentDir : 1;
        currentKey = key;
        headers.forEach(function (h) { h.removeAttribute("data-sort-dir"); });
        th.setAttribute("data-sort-dir", currentDir === 1 ? "asc" : "desc");

        var rows = Array.from(tbody.querySelectorAll("tr"));
        rows.sort(function (a, b) {
          var aCell = a.children[index];
          var bCell = b.children[index];
          var aRaw = aCell ? aCell.getAttribute("data-sort") : null;
          var bRaw = bCell ? bCell.getAttribute("data-sort") : null;
          var aNum = aRaw === null || aRaw === "" ? NaN : Number(aRaw);
          var bNum = bRaw === null || bRaw === "" ? NaN : Number(bRaw);
          var cmp;
          if (!isNaN(aNum) && !isNaN(bNum)) {
            cmp = aNum - bNum;
          } else {
            cmp = String(aRaw).localeCompare(String(bRaw));
          }
          return cmp * currentDir;
        });
        rows.forEach(function (row) { tbody.appendChild(row); });
      });
    });
  });
});
"""


def _events_controls_html() -> str:
    return """
<div class="controls">
  <div class="control-group">
    <label class="title" for="tumor-type-filter">Tumor type</label>
    <select id="tumor-type-filter"><option value="">All tumor types</option></select>
  </div>
  <div class="control-group">
    <span class="title">Domain status</span>
    <span id="status-filter"></span>
  </div>
  <span id="events-summary"></span>
</div>
"""


def _render_events_table(events: list[dict[str, Any]]) -> str:
    if not events:
        return '<p class="note">No per-event data available for this run.</p>'
    header = (
        "<tr><th>Event</th><th>Sample</th><th>Tumor type</th><th>Oncotree</th>"
        "<th>Partner gene</th><th>Role</th><th>Breakpoint (aa)</th><th>Domain status</th></tr>"
    )
    body_rows = []
    for event in events:
        tumor_type = event.get("tumor_type") or ""
        domain_status = event.get("domain_status") or ""
        event_id_value = str(event.get("event_id") or "")
        body_rows.append(
            '<tr data-tumor-type="{tumor_attr}" data-domain-status="{status_attr}" '
            'data-event-id="{event_id_attr}">'
            "<td>{event_id}</td><td>{sample_id}</td><td>{tumor_type}</td><td>{oncotree}</td>"
            "<td>{partner}</td><td>{role}</td><td>{breakpoint}</td><td>{status}</td></tr>".format(
                tumor_attr=_esc(tumor_type),
                status_attr=_esc(domain_status),
                event_id_attr=_esc(event_id_value),
                event_id=_esc(event_id_value),
                sample_id=_esc(str(event.get("sample_id") or "")),
                tumor_type=_esc(tumor_type),
                oncotree=_esc(str(event.get("oncotree_code") or "")),
                partner=_esc(str(event.get("partner_gene") or "")),
                role=_esc(str(event.get("target_role") or "")),
                breakpoint=_esc(str(event.get("breakpoint_protein_position") or "")),
                status=_esc(_DOMAIN_STATUS_LABELS.get(domain_status, domain_status)),
            )
        )
    return (
        '<div class="events-table-frame"><table class="events" id="events-table">'
        f"<thead>{header}</thead><tbody>{''.join(body_rows)}</tbody></table></div>"
    )


def _composite_score_section(payload: dict) -> str:
    """A sortable table of ``composite_score``'s per-partner ranked results
    (``algorithm_results[*].Tables.composite_evidence_ranking``), omitted
    entirely when that algorithm didn't run or produced no ranked partner
    for this gene/run."""
    result = _algorithm_result(payload, "composite_score")
    if result is None or _algorithm_failed(result):
        return ""
    rows = ((result.get("Tables") or {}).get("composite_evidence_ranking")) or []
    if not rows:
        return ""

    header_cells = "".join(
        f'<th data-sort-key="{_esc(key)}">{_esc(label)}</th>'
        for key, label in _COMPOSITE_SCORE_COLUMNS
    )
    body_rows = []
    for row in rows:
        cells = []
        for key, _label in _COMPOSITE_SCORE_COLUMNS:
            value = row.get(key)
            sort_value = "" if value is None else str(value)
            display = _format_summary_value(value) if value is not None else ""
            cells.append(f'<td data-sort="{_esc(sort_value)}">{_esc(display)}</td>')
        body_rows.append(f"<tr>{''.join(cells)}</tr>")

    return (
        "<h2>Composite evidence ranking (per partner gene)</h2>"
        f'<p class="panel-note">Click a column header to sort ascending/descending. '
        f"{len(rows)} partner gene(s) ranked by <code>composite_score</code>.</p>"
        '<div class="data-table-frame">'
        '<table class="data-table sortable" id="composite-score-table">'
        f"<thead><tr>{header_cells}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody></table></div>"
    )


def _joint_partner_section(payload: dict) -> str:
    """The ``joint_partner`` algorithm's own pair-enrichment table
    (observed/expected/p-value/odds-ratio) -- distinct from a single-gene
    domain-retention result, so it gets its own small table rather than
    being coerced into the composite-score/events layout."""
    result = _algorithm_result(payload, "joint_partner")
    if result is None or _algorithm_failed(result):
        return ""
    rows = (result.get("Tables") or {}).get("pair_results") or []
    if not rows:
        return ""
    header = (
        "<tr><th>5&prime; gene</th><th>3&prime; gene</th><th>Eligible events</th>"
        "<th>Observed</th><th>Expected</th><th>p-value</th><th>Odds ratio</th></tr>"
    )
    body = "".join(
        "<tr><td>{g5}</td><td>{g3}</td><td>{eligible}</td><td>{observed}</td>"
        "<td>{expected}</td><td>{p_value}</td><td>{odds_ratio}</td></tr>".format(
            g5=_esc(str(row.get("gene5") or "")),
            g3=_esc(str(row.get("gene3") or "")),
            eligible=_esc(str(row.get("eligible_event_count", ""))),
            observed=_esc(str(row.get("observed_count", ""))),
            expected=_esc(_format_summary_value(row.get("expected_count"))),
            p_value=_esc(_format_summary_value(row.get("p_value"))),
            odds_ratio=_esc(
                "—" if row.get("odds_ratio") is None else _format_summary_value(row["odds_ratio"])
            ),
        )
        for row in rows
    )
    return (
        "<h2>Gene-pair enrichment test</h2>"
        '<div class="data-table-frame"><table class="data-table">'
        f"<thead>{header}</thead><tbody>{body}</tbody></table></div>"
    )


def _genomic_clustering_section(payload: dict) -> str:
    """The genomic-coordinate clustering view (``genomic_position_recurrence``),
    linked to the protein-position picture via its own
    ``protein_position_genomic_spread`` cross-reference table: for every
    protein junction position shared by >=2 events, whether their genomic
    breakpoints are one exact DNA coordinate or scattered across a wider
    span. Gated on ``Summary.determinable`` -- most genes/older runs have
    no genomic breakpoints supplied at all, and this must degrade to
    nothing rendered, never an empty/broken section."""
    result = _algorithm_result(payload, "genomic_position_recurrence")
    if result is None or _algorithm_failed(result):
        return ""
    summary = result.get("Summary") or {}
    if not summary.get("determinable"):
        return ""
    tables = result.get("Tables") or {}
    spread_rows = tables.get("protein_position_genomic_spread") or []
    bin_rows = (tables.get("genomic_bin_recurrence") or [])[:10]
    note = summary.get("genomic_vs_protein_clustering_note") or ""
    reference_build = summary.get("reference_build")

    spread_html = ""
    if spread_rows:
        header = (
            "<tr><th>Protein position (aa)</th><th>Events</th>"
            "<th>Distinct genomic positions</th><th>Genomic span (bp)</th>"
            "<th>Spans &gt;1 chromosome</th></tr>"
        )
        body = "".join(
            "<tr><td>{aa}</td><td>{n_events}</td><td>{distinct}</td><td>{span}</td>"
            "<td>{multi}</td></tr>".format(
                aa=_esc(str(row.get("protein_position_aa"))),
                n_events=_esc(str(row.get("n_events"))),
                distinct=_esc(str(row.get("n_distinct_genomic_positions"))),
                span=_esc(
                    "—" if row.get("genomic_span_bp") is None else str(row["genomic_span_bp"])
                ),
                multi=_esc("Yes" if row.get("spans_multiple_chromosomes") else "No"),
            )
            for row in spread_rows
        )
        spread_html = (
            "<p><strong>Protein position &harr; genomic breakpoint spread</strong> "
            "(linked view: each row is a protein junction position also shown in the "
            "fusion schematic/events table above).</p>"
            '<div class="data-table-frame"><table class="data-table">'
            f"<thead>{header}</thead><tbody>{body}</tbody></table></div>"
        )

    bin_html = ""
    if bin_rows:
        header = "<tr><th>Chromosome</th><th>Bin start</th><th>Bin end</th><th>Events</th></tr>"
        body = "".join(
            "<tr><td>{chrom}</td><td>{start}</td><td>{end}</td><td>{n_events}</td></tr>".format(
                chrom=_esc(str(row.get("chromosome"))),
                start=_esc(str(row.get("bin_start"))),
                end=_esc(str(row.get("bin_end"))),
                n_events=_esc(str(row.get("n_events"))),
            )
            for row in bin_rows
        )
        bin_html = (
            "<p><strong>Top genomic bins</strong> "
            f"(fixed {_esc(str(summary.get('bin_size_bp', '')))} bp width):</p>"
            '<div class="data-table-frame"><table class="data-table">'
            f"<thead>{header}</thead><tbody>{body}</tbody></table></div>"
        )

    build_note = f" (build {_esc(str(reference_build))})" if reference_build else ""
    tables_html = (
        f'<div class="two-col">{spread_html}{bin_html}</div>'
        if spread_html and bin_html
        else spread_html + bin_html
    )
    return (
        f"<h2>Genomic-coordinate clustering{build_note}</h2>"
        f'<p class="note">{_esc(str(note))}</p>'
        f"{tables_html}"
    )


def _cross_cohort_panel(cross_cohort: dict[str, Any] | None) -> str:
    """A small panel showing the pooled Cochran-Mantel-Haenszel
    statistic/p-value/common-odds-ratio from a gene's cross-cohort
    concordance report, when one exists for this run (see
    :func:`_cross_cohort_data`). Omitted entirely otherwise -- this is not
    computed for most genes."""
    if not cross_cohort:
        return ""
    cmh = cross_cohort.get("cmh") or {}
    strata = cross_cohort.get("strata") or []
    significant = cross_cohort.get("nominal_pooled_significant")
    stat_cells = "".join(
        f'<div class="stat"><div class="label">{_esc(label)}</div>'
        f'<div class="value">{_esc(_format_summary_value(value))}</div></div>'
        for label, value in (
            ("Pooled CMH statistic", cmh.get("statistic")),
            ("Pooled p-value", cmh.get("p_value")),
            ("Common odds ratio", cmh.get("common_odds_ratio")),
            ("Informative strata", f"{cmh.get('informative_strata')} / {cmh.get('total_strata')}"),
        )
        if value is not None
    )
    strata_rows = "".join(
        "<tr><td>{study}</td><td>{direction}</td><td>{informative}</td></tr>".format(
            study=_esc(str(stratum.get("study_id") or "")),
            direction=_esc(str(stratum.get("association_direction") or "")),
            informative=_esc("Yes" if stratum.get("informative") else "No"),
        )
        for stratum in strata
    )
    strata_html = (
        (
            '<div class="data-table-frame" style="margin-top:12px;"><table class="data-table">'
            "<thead><tr><th>Cohort</th><th>Association direction</th><th>Informative</th></tr>"
            f"</thead><tbody>{strata_rows}</tbody></table></div>"
        )
        if strata
        else ""
    )
    warnings_items = "".join(
        f"<li>{_esc(str(warning))}</li>" for warning in (cross_cohort.get("warnings") or [])
    )
    warnings_html = f'<ul class="panel-note">{warnings_items}</ul>' if warnings_items else ""
    significance_note = (
        "Nominally significant" if significant else "Not nominally significant"
    ) + " at &alpha;=0.05."
    return (
        "<h2>Cross-cohort concordance (Cochran&ndash;Mantel&ndash;Haenszel)</h2>"
        f'<div class="card"><div class="stat-grid">{stat_cells}</div>'
        f'<p class="note">{significance_note}</p>'
        f"{strata_html}"
        f"{warnings_html}"
        "</div>"
    )


def _page_head(title: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{_PAGE_STYLE}</style>
</head>
<body>
<div class="wrap">
"""


def _render_gene_page(
    payload: dict,
    *,
    domain_svg: str | None,
    fusion_svg: str | None,
    back_link: str | None,
    config_source: str | None = None,
    cross_cohort: dict[str, Any] | None = None,
) -> str:
    gene_symbol = _esc(str(payload.get("gene_symbol", "")))
    study_id = _esc(str(payload.get("study_id", "")))
    summary = payload.get("summary") or {}
    events = _event_records(payload)

    stat_cells = "".join(
        f'<div class="stat"><div class="label">{_esc(_humanize(key))}</div>'
        f'<div class="value">{_esc(_format_summary_value(value))}</div></div>'
        for key, value in _scalar_summary_items(summary)
    )

    sections = []
    back_html = (
        f'<a class="back-link" href="{_esc(back_link)}">&larr; Back to cohort-scan summary</a>'
        if back_link
        else ""
    )

    if domain_svg:
        sections.append(
            "<h2>Domain retention / outlier lollipop track</h2>"
            '<p class="panel-note">Click a point to jump to its event in the table below.</p>'
            f'<div class="svg-frame" id="domain-track-svg">{domain_svg}</div>'
        )

    sections.append(_composite_score_section(payload))

    controls_html = ""
    if fusion_svg:
        controls_html = _events_controls_html()
        sections.append(
            "<h2>Fusion-transcript schematic</h2>"
            + controls_html
            + f'<div class="svg-frame" id="fusion-schematic-svg">{fusion_svg}</div>'
        )

    sections.append(_genomic_clustering_section(payload))
    sections.append(_cross_cohort_panel(cross_cohort))

    sections.append("<h2>Events" + ("" if fusion_svg else " (filters below)") + "</h2>")
    if not fusion_svg:
        sections.append(_events_controls_html())
    sections.append(_render_events_table(events))
    events_json = _safe_json_script(events)
    status_labels_json = _safe_json_script(_DOMAIN_STATUS_LABELS)

    badge_html = _config_source_badge(config_source)
    return (
        _page_head(f"{gene_symbol} fusion-hotspot report")
        + f"""  {back_html}
  <h1>{gene_symbol}{badge_html}</h1>
  <p class="subtitle">{study_id}</p>
  <div class="stat-grid">{stat_cells}</div>
  {"".join(sections)}
</div>
<div class="lollipop-tooltip" id="lollipop-tooltip" hidden></div>
<script type="application/json" id="events-data">{events_json}</script>
<script type="application/json" id="status-labels">{status_labels_json}</script>
<script>{_GENE_PAGE_SCRIPT}</script>
<script>{_SORTABLE_TABLE_SCRIPT}</script>
</body>
</html>
"""
    )


def _render_gene_pair_page(
    payload: dict,
    *,
    back_link: str | None,
    config_source: str | None = None,
    cross_cohort: dict[str, Any] | None = None,
) -> str:
    """Dedicated template for a ``gene_pair`` config's run (e.g. EML4-ALK,
    TMPRSS2-ERG): the result shape is observed/expected/enrichment-p-value
    from ``joint_partner``, not a single-gene domain-retention result, so
    this does not reuse :func:`_render_gene_page`'s domain/fusion-schematic
    layout -- only the sections that make sense for both (composite score,
    genomic clustering, cross-cohort CMH, the events table/tumor filter)
    are shared."""
    gene_symbol = _esc(str(payload.get("gene_symbol", "")))
    study_id = _esc(str(payload.get("study_id", "")))
    summary = payload.get("summary") or {}
    events = _event_records(payload)

    stat_cells = "".join(
        f'<div class="stat"><div class="label">{_esc(_humanize(key))}</div>'
        f'<div class="value">{_esc(_format_summary_value(value))}</div></div>'
        for key, value in _scalar_summary_items(summary)
    )

    back_html = (
        f'<a class="back-link" href="{_esc(back_link)}">&larr; Back to cohort-scan summary</a>'
        if back_link
        else ""
    )

    sections = [
        _joint_partner_section(payload),
        _composite_score_section(payload),
        _genomic_clustering_section(payload),
        _cross_cohort_panel(cross_cohort),
        "<h2>Events (filters below)</h2>",
        _events_controls_html(),
        _render_events_table(events),
    ]
    events_json = _safe_json_script(events)
    status_labels_json = _safe_json_script(_DOMAIN_STATUS_LABELS)

    badge_html = _config_source_badge(config_source)
    return (
        _page_head(f"{gene_symbol} gene-pair fusion-hotspot report")
        + f"""  {back_html}
  <h1>{gene_symbol}{badge_html}</h1>
  <p class="subtitle">{study_id} &mdash; gene-pair (joint-partner) fusion enrichment</p>
  <div class="stat-grid">{stat_cells}</div>
  {"".join(sections)}
</div>
<script type="application/json" id="events-data">{events_json}</script>
<script type="application/json" id="status-labels">{status_labels_json}</script>
<script>{_GENE_PAGE_SCRIPT}</script>
<script>{_SORTABLE_TABLE_SCRIPT}</script>
</body>
</html>
"""
    )


__all__ = ["build_run_viewer"]
