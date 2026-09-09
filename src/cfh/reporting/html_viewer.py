"""Backend-less static-site HTML viewer for an already-written ``cfh`` run
directory (Milestone 4, "Option A"): plain HTML/CSS/vanilla JS that reads
the already-generated ``results.json``/``summary.json`` and SVG artifacts
and re-renders them into one or more self-contained ``.html`` files under
``<run_directory>/viewer/``.

This module is purely a *presentation* layer over artifacts some other part
of the pipeline already wrote to disk -- it re-reads ``results.json``/
``summary.json``/``*.svg`` exactly as committed and never recomputes,
re-derives, or overrides a single statistic. Nothing here can change a
run's numbers.

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

Any other/incomplete run directory (e.g. a ``gene_pair`` run, which has no
``gene_track`` and so no domain/fusion SVGs) still gets a best-effort gene
page -- whatever sections have real data to show (summary stats, the
events table) are rendered; sections with no underlying artifact are
simply omitted, never fabricated.
"""

from __future__ import annotations

import json
from html import escape as _esc
from pathlib import Path
from typing import Any

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


def _read_text_or_none(path: Path) -> str | None:
    return path.read_text() if path.exists() else None


def _build_cohort_scan_viewer(run_directory: Path) -> Path:
    cohort_dir = run_directory / "cohort_scan"
    summary = _read_json(cohort_dir / "summary.json")
    manhattan_svg = _read_text_or_none(cohort_dir / "manhattan.svg") or ""

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
        )

    index_path = viewer_dir / "index.html"
    index_path.write_text(
        _render_landing_page(summary, manhattan_svg, available_genes=set(available_genes))
    )
    return index_path


def _build_single_gene_viewer(
    gene_run_dir: Path, output_path: Path, *, back_link: str | None
) -> Path:
    payload = _read_json(gene_run_dir / "results.json")
    domain_svg = _read_text_or_none(
        gene_run_dir / "visualizations" / "domain_retention_outliers.svg"
    )
    fusion_svg = _read_text_or_none(gene_run_dir / "visualizations" / "fusion_schematic.svg")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        _render_gene_page(
            payload, domain_svg=domain_svg, fusion_svg=fusion_svg, back_link=back_link
        )
    )
    return output_path


def _safe_json_script(data: Any) -> str:
    """JSON-serialize ``data`` for safe embedding inside an inline
    ``<script type="application/json">`` block -- escaping ``</`` so a
    value containing a literal ``</script>`` cannot prematurely close the
    tag (``"\\/"`` is a valid JSON escape for ``/``, so this round-trips
    through ``JSON.parse`` unchanged)."""
    return json.dumps(data, allow_nan=False).replace("</", "<\\/")


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
        "still carries a hover tooltip on its point."
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


_GENE_PAGE_SCRIPT = r"""
document.addEventListener("DOMContentLoaded", function () {
  var events = JSON.parse(document.getElementById("events-data").textContent);
  var statusLabels = JSON.parse(document.getElementById("status-labels").textContent);

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
});
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
        body_rows.append(
            '<tr data-tumor-type="{tumor_attr}" data-domain-status="{status_attr}">'
            "<td>{event_id}</td><td>{sample_id}</td><td>{tumor_type}</td><td>{oncotree}</td>"
            "<td>{partner}</td><td>{role}</td><td>{breakpoint}</td><td>{status}</td></tr>".format(
                tumor_attr=_esc(tumor_type),
                status_attr=_esc(domain_status),
                event_id=_esc(str(event.get("event_id") or "")),
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


def _render_gene_page(
    payload: dict, *, domain_svg: str | None, fusion_svg: str | None, back_link: str | None
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
            f'<div class="svg-frame">{domain_svg}</div>'
        )

    controls_html = ""
    if fusion_svg:
        controls_html = """
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
        sections.append(
            "<h2>Fusion-transcript schematic</h2>"
            + controls_html
            + f'<div class="svg-frame" id="fusion-schematic-svg">{fusion_svg}</div>'
        )

    sections.append("<h2>Events" + ("" if fusion_svg else " (filters below)") + "</h2>")
    if not fusion_svg:
        sections.append(
            """
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
        )
    sections.append(_render_events_table(events))
    events_json = _safe_json_script(events)
    status_labels_json = _safe_json_script(_DOMAIN_STATUS_LABELS)

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{gene_symbol} fusion-hotspot report</title>
<style>{_PAGE_STYLE}</style>
</head>
<body>
<div class="wrap">
  {back_html}
  <h1>{gene_symbol}</h1>
  <p class="subtitle">{study_id}</p>
  <div class="stat-grid">{stat_cells}</div>
  {"".join(sections)}
</div>
<script type="application/json" id="events-data">{events_json}</script>
<script type="application/json" id="status-labels">{status_labels_json}</script>
<script>{_GENE_PAGE_SCRIPT}</script>
</body>
</html>
"""


__all__ = ["build_run_viewer"]
