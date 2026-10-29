"""Landing page for the published GitHub Pages site.

Indexes every gene that any published run covers, so a reader can search
for a gene and land on the page with the most to say about it, with the
other runs that include it listed alongside. Built offline from committed
run artifacts only (``cohort_scan/summary.json``, ``results.json``), like
the per-run viewers it links to.

Which page is "most meaningful" for a gene, best first:

1. its full report in the newest cohort scan (current code, cohort-wide
   FDR, tumor-type adjustment, fragility, exon/intron view);
2. a dedicated single-gene run on that same study;
3. the newest cohort scan's summary, for a gene analyzed without a full
   report;
4. single-gene runs on other studies (replication cohorts);
5. gene-pair runs.

Ties go to the newer run.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from html import escape as _esc
from pathlib import Path
from typing import Any

_RUN_DIR_PATTERN = re.compile(r"^(?P<run_type>.+)_(?P<timestamp>\d{8}T\d{6}Z)$")

_COHORT_REPORT, _SAME_STUDY_RUN, _COHORT_SUMMARY, _OTHER_STUDY_RUN, _PAIR_RUN = range(5)


@dataclass(frozen=True)
class _Run:
    name: str
    kind: str  # "cohort" | "gene" | "pair"
    study_id: str
    date: str  # YYYY-MM-DD
    timestamp: str
    url: str


def study_label(study_id: str) -> str:
    """A readable name for a cBioPortal study id."""
    tcga = re.fullmatch(r"([a-z]+)_tcga_pan_can_atlas_(\d{4})", study_id)
    if tcga:
        return f"TCGA {tcga.group(1).upper()} PanCancer Atlas"
    impact = re.fullmatch(r"msk_impact(?:_(\d+k))?_(\d{4})", study_id)
    if impact:
        size = f" {impact.group(1).upper()}" if impact.group(1) else ""
        return f"MSK-IMPACT{size} ({impact.group(2)})"
    return study_id


def _read_json(path: Path) -> Any | None:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def _run_meta(run_dir: Path, kind: str, study_id: str) -> _Run:
    match = _RUN_DIR_PATTERN.match(run_dir.name)
    timestamp = match.group("timestamp") if match else ""
    date = f"{timestamp[:4]}-{timestamp[4:6]}-{timestamp[6:8]}" if timestamp else ""
    return _Run(run_dir.name, kind, study_id, date, timestamp, f"runs/{run_dir.name}/index.html")


def _positive_controls(cohort_dir: Path) -> tuple[int, int] | None:
    path = cohort_dir / "positive_controls.tsv"
    if not path.exists():
        return None
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return sum(row.get("status") == "pass" for row in rows), len(rows)


def collect_site_index(runs_dir: Path, published: list[str]) -> dict:
    """Gene and run index for the landing page, from published runs only."""
    runs: list[dict] = []
    entries: dict[str, list[tuple[tuple, dict]]] = {}
    cohorts: list[tuple[_Run, dict, Path]] = []
    singles: list[tuple[_Run, dict]] = []

    for name in published:
        run_dir = runs_dir / name
        summary = _read_json(run_dir / "cohort_scan" / "summary.json")
        if summary is not None:
            meta = _run_meta(run_dir, "cohort", summary.get("study_id", ""))
            cohorts.append((meta, summary, run_dir / "cohort_scan"))
            continue
        payload = _read_json(run_dir / "results.json")
        if payload is None:
            continue
        pair = (payload.get("summary") or {}).get("gene_pair")
        meta = _run_meta(run_dir, "pair" if pair else "gene", payload.get("study_id", ""))
        singles.append((meta, payload))

    newest_cohort = max(cohorts, key=lambda item: item[0].timestamp, default=None)
    primary_study = newest_cohort[0].study_id if newest_cohort else None

    def add(gene: str, priority: int, meta: _Run, entry: dict) -> None:
        # Newer runs first within a priority: negate the timestamp's digits.
        newness = -int(meta.timestamp[:8] + meta.timestamp[9:15]) if meta.timestamp else 0
        entries.setdefault(gene.upper(), []).append(
            (
                (priority, newness),
                {
                    "run": meta.name,
                    "kind": meta.kind,
                    "study": study_label(meta.study_id),
                    "date": meta.date,
                    **entry,
                },
            )
        )

    for meta, summary, cohort_dir in cohorts:
        genes = summary.get("genes") or []
        analyzed = [row for row in genes if row.get("status") == "ok"]
        controls = _positive_controls(cohort_dir)
        reconstitution = summary.get("domain_reconstitution") or {}
        runs.append(
            {
                "kind": "cohort",
                "title": f"Cohort scan: {study_label(meta.study_id)}",
                "study": study_label(meta.study_id),
                "date": meta.date,
                "url": meta.url,
                "newest": newest_cohort is not None and meta == newest_cohort[0],
                "genes_analyzed": len(analyzed),
                "significant_genes": len(summary.get("significant_genes") or []),
                "positive_controls": list(controls) if controls else None,
                "domain_completing_fusions": reconstitution.get("complementary"),
            }
        )
        reports = {
            entry.name
            for entry in (cohort_dir / "gene_reports").glob("*")
            if (entry / "results.json").exists()
        }
        honorable = {row.get("gene_symbol") for row in summary.get("honorable_mentions") or []}
        is_newest = newest_cohort is not None and meta == newest_cohort[0]
        for row in genes:
            gene = str(row.get("gene_symbol") or "")
            if not gene:
                continue
            has_report = gene.lower() in reports
            priority = (
                _COHORT_REPORT
                if has_report and is_newest
                else _COHORT_SUMMARY
                if is_newest
                else _OTHER_STUDY_RUN
            )
            add(
                gene,
                priority,
                meta,
                {
                    "label": "Full gene report" if has_report else "Cohort summary",
                    "url": (
                        f"runs/{meta.name}/genes/{gene.lower()}.html" if has_report else meta.url
                    ),
                    "events": row.get("n_events_analyzed"),
                    "patients": row.get("distinct_patient_count"),
                    "q": row.get("min_fdr_adjusted_q_value"),
                    "significant": row.get("fdr_significant") is True,
                    "honorable": gene in honorable,
                    "analyzed": row.get("status") == "ok",
                },
            )

    for meta, payload in singles:
        summary = payload.get("summary") or {}
        pair = summary.get("gene_pair")
        genes = list(pair) if pair else [str(payload.get("gene_symbol") or "")]
        title = "–".join(genes) if pair else genes[0]
        runs.append(
            {
                "kind": meta.kind,
                "title": f"{title} {'gene-pair' if pair else 'single-gene'} run",
                "study": study_label(meta.study_id),
                "date": meta.date,
                "url": meta.url,
                "events": summary.get("total_fusions"),
            }
        )
        priority = (
            _PAIR_RUN
            if pair
            else _SAME_STUDY_RUN
            if meta.study_id == primary_study
            else _OTHER_STUDY_RUN
        )
        for gene in genes:
            if not gene:
                continue
            add(
                gene,
                priority,
                meta,
                {
                    "label": f"{title} pair run" if pair else "Single-gene run",
                    "url": meta.url,
                    "events": summary.get("total_fusions"),
                    "patients": summary.get("known_patient_count"),
                    "q": None,
                    "significant": False,
                    "honorable": False,
                    "analyzed": True,
                },
            )

    gene_index = []
    for gene, ranked in sorted(entries.items()):
        ordered = [entry for _, entry in sorted(ranked, key=lambda item: item[0])]
        best = ordered[0]
        cohort_entry = next((entry for entry in ordered if entry["kind"] == "cohort"), None)
        flags = cohort_entry or best
        gene_index.append(
            {
                "gene": gene,
                "best": best,
                "also": ordered[1:],
                "events": flags.get("events"),
                "q": flags.get("q"),
                "significant": flags.get("significant", False),
                "honorable": flags.get("honorable", False),
                "analyzed": flags.get("analyzed", True),
            }
        )
    order = {"cohort": 0, "gene": 1, "pair": 2}
    runs.sort(key=lambda run: (order[run["kind"]], not run.get("newest"), run["title"]))
    return {"runs": runs, "genes": gene_index}


def render_site_index(index: dict) -> str:
    """The landing page HTML."""
    newest = next((run for run in index["runs"] if run.get("newest")), None)
    analyzed = sum(1 for gene in index["genes"] if gene["analyzed"])
    significant = sorted(gene["gene"] for gene in index["genes"] if gene["significant"])
    stats = []
    if newest:
        stats.append(("Genes analysed", f"{newest['genes_analyzed']:,}"))
        stats.append(("FDR-significant genes", str(newest["significant_genes"])))
        if newest.get("positive_controls"):
            passed, total = newest["positive_controls"]
            stats.append(("Textbook fusions recovered", f"{passed} of {total}"))
        if newest.get("domain_completing_fusions") is not None:
            stats.append(
                ("Truncated domains completed by partner", str(newest["domain_completing_fusions"]))
            )
    stats_html = "".join(
        f'<div class="stat"><div class="stat-value">{_esc(value)}</div>'
        f'<div class="stat-label">{_esc(label)}</div></div>'
        for label, value in stats
    )
    chips = "".join(
        f'<button type="button" class="chip" data-gene="{_esc(gene)}">{_esc(gene)}</button>'
        for gene in significant
    )
    run_cards = "".join(_run_card(run) for run in index["runs"])
    subtitle = (
        f"Latest cohort scan: {_esc(newest['study'])}, {_esc(newest['date'])}" if newest else ""
    )
    data = json.dumps(index["genes"], separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Fusion Hotspots</title>
<style>{_STYLE}</style>
</head>
<body>
<header class="hero">
  <div class="wrap">
    <p class="eyebrow">genome-nexus / fusion-hotspots</p>
    <h1>Fusion hotspots</h1>
    <p class="lede">Where recurrent gene fusions break, which protein domains they keep or
    lose, and whether that pattern is selected for. {subtitle}.</p>
    <div class="stats">{stats_html}</div>
  </div>
</header>
<main class="wrap">
  <section aria-labelledby="search-heading">
    <h2 id="search-heading">Find a gene</h2>
    <p class="hint">{analyzed:,} genes analysed. Each opens the page with the most detail for
    that gene; other runs that include it are listed beside it.</p>
    <div class="search">
      <input id="gene-search" type="search" autocomplete="off" spellcheck="false"
        placeholder="Search a gene, e.g. ALK, FGFR2, ERG" aria-label="Search genes">
      <span class="count" id="result-count" aria-live="polite"></span>
    </div>
    <div class="chips" aria-label="FDR-significant genes">
      <span class="chips-label">FDR-significant:</span>{chips}
    </div>
    <div class="table-frame">
      <table id="gene-table">
        <thead><tr>
          <th scope="col">Gene</th><th scope="col">Best page</th>
          <th scope="col" class="num">Events</th><th scope="col" class="num">Best q</th>
          <th scope="col">Also in</th>
        </tr></thead>
        <tbody></tbody>
      </table>
    </div>
    <p class="hint" id="more-hint" hidden></p>
  </section>
  <section aria-labelledby="runs-heading">
    <h2 id="runs-heading">Runs</h2>
    <div class="cards">{run_cards}</div>
  </section>
</main>
<footer class="wrap">Built from the committed <code>runs/</code> artifacts of
<a href="https://github.com/genome-nexus/fusion-hotspots">genome-nexus/fusion-hotspots</a>.
q-values are Benjamini–Hochberg across the cohort scan's tests.</footer>
<script type="application/json" id="gene-data">{data}</script>
<script>{_SCRIPT}</script>
</body>
</html>
"""


def _run_card(run: dict) -> str:
    details = []
    if run["kind"] == "cohort":
        details.append(f"{run['genes_analyzed']:,} genes analysed")
        details.append(f"{run['significant_genes']} FDR-significant")
        if run.get("positive_controls"):
            passed, total = run["positive_controls"]
            details.append(f"{passed}/{total} textbook fusions recovered")
    elif run.get("events") is not None:
        details.append(f"{run['events']:,} fusions")
    badge = '<span class="badge">Latest</span>' if run.get("newest") else ""
    kind_label = {"cohort": "Cohort scan", "gene": "Single gene", "pair": "Gene pair"}[run["kind"]]
    return (
        f'<a class="card card-{run["kind"]}" href="{_esc(run["url"])}">'
        f'<span class="card-kind">{kind_label}</span>{badge}'
        f'<span class="card-title">{_esc(run["title"])}</span>'
        f'<span class="card-meta">{_esc(run["study"])} · {_esc(run["date"])}</span>'
        f'<span class="card-detail">{_esc(" · ".join(details))}</span></a>'
    )


_STYLE = """
:root {
  --bg: #f7f8fa; --surface: #ffffff; --text: #17202a; --muted: #5b6770;
  --line: #e1e5ea; --accent: #0e7c7b; --accent-soft: #e3f2f1; --sig: #b4441d;
  --sig-soft: #fbe9e2; --hero: #0f2b3a; --hero-text: #f1f6f8;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #11161b; --surface: #182029; --text: #e6ecf0; --muted: #9aa7b1;
    --line: #2a3540; --accent: #4fc0bd; --accent-soft: #173533; --sig: #f08a64;
    --sig-soft: #3a2219; --hero: #0b1c26; --hero-text: #e6ecf0;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: var(--text);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif;
}
a { color: var(--accent); }
.wrap { max-width: 68rem; margin: 0 auto; padding: 0 16px; }
.hero { background: var(--hero); color: var(--hero-text); padding: 40px 0 32px; }
.eyebrow { margin: 0; font-size: 12px; letter-spacing: .08em; text-transform: uppercase;
  opacity: .7; }
h1 { margin: 4px 0 8px; font-size: 34px; letter-spacing: -.01em; }
.lede { margin: 0; max-width: 46rem; opacity: .85; }
.stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: 12px; margin-top: 24px; }
.stat { background: rgba(255,255,255,.07); border: 1px solid rgba(255,255,255,.12);
  border-radius: 10px; padding: 12px 14px; }
.stat-value { font-size: 24px; font-weight: 650; }
.stat-label { font-size: 12.5px; opacity: .75; }
main section { margin: 32px 0; }
h2 { font-size: 20px; margin: 0 0 4px; }
.hint { color: var(--muted); margin: 0 0 12px; font-size: 13.5px; }
.search { position: relative; }
#gene-search { width: 100%; font: inherit; font-size: 17px; padding: 12px 110px 12px 14px;
  border: 1px solid var(--line); border-radius: 10px; background: var(--surface);
  color: var(--text); }
#gene-search:focus { outline: 2px solid var(--accent); outline-offset: 1px; }
.count { position: absolute; right: 14px; top: 50%; transform: translateY(-50%);
  color: var(--muted); font-size: 13px; }
.chips { display: flex; flex-wrap: wrap; gap: 6px; align-items: center; margin: 12px 0; }
.chips-label { color: var(--muted); font-size: 13px; margin-right: 2px; }
.chip { font: inherit; font-size: 13px; border: 1px solid var(--line); background: var(--surface);
  color: var(--text); border-radius: 999px; padding: 2px 10px; cursor: pointer; }
.chip:hover { border-color: var(--accent); color: var(--accent); }
.table-frame { overflow-x: auto; background: var(--surface); border: 1px solid var(--line);
  border-radius: 10px; }
table { width: 100%; border-collapse: collapse; font-size: 14px; }
th, td { text-align: left; padding: 8px 12px; border-bottom: 1px solid var(--line);
  vertical-align: top; }
th { font-size: 12px; text-transform: uppercase; letter-spacing: .04em; color: var(--muted);
  font-weight: 600; }
tbody tr:last-child td { border-bottom: 0; }
tbody tr:hover { background: var(--accent-soft); }
.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
.gene { font-weight: 650; white-space: nowrap; }
.tag { display: inline-block; font-size: 11px; font-weight: 600; border-radius: 4px;
  padding: 0 6px; margin-left: 6px; vertical-align: 1px; }
.tag-sig { background: var(--sig-soft); color: var(--sig); }
.tag-hm { background: var(--accent-soft); color: var(--accent); }
.tag-na { background: var(--line); color: var(--muted); }
.best { min-width: 13rem; }
.best a { font-weight: 600; text-decoration: none; white-space: nowrap; }
.best a:hover { text-decoration: underline; }
.sub { display: block; color: var(--muted); font-size: 12.5px; }
.also a { display: block; font-size: 13px; white-space: nowrap; text-decoration: none; }
.also a:hover { text-decoration: underline; }
.cards { display: grid; grid-template-columns: repeat(auto-fill, minmax(240px, 1fr)); gap: 12px; }
.card { display: flex; flex-direction: column; gap: 2px; padding: 14px 16px;
  background: var(--surface); border: 1px solid var(--line); border-radius: 10px;
  color: var(--text); text-decoration: none; position: relative; }
.card:hover { border-color: var(--accent); }
.card-cohort { border-left: 4px solid var(--accent); }
.card-kind { font-size: 11.5px; text-transform: uppercase; letter-spacing: .05em;
  color: var(--muted); }
.badge { position: absolute; top: 12px; right: 12px; font-size: 11px; font-weight: 600;
  background: var(--accent); color: var(--surface); border-radius: 4px; padding: 1px 6px; }
.card-title { font-weight: 650; }
.card-meta, .card-detail { font-size: 13px; color: var(--muted); }
footer { color: var(--muted); font-size: 13px; padding: 16px 16px 40px; }
@media (max-width: 600px) {
  h1 { font-size: 28px; }
  th:nth-child(5), td:nth-child(5) { display: none; }
  th, td { padding: 8px; }
  .best { min-width: 0; }
  .tag { display: table; margin: 2px 0 0; }
  .best a { white-space: normal; }
}
"""

_SCRIPT = r"""
(function () {
  var genes = JSON.parse(document.getElementById("gene-data").textContent);
  var input = document.getElementById("gene-search");
  var tbody = document.querySelector("#gene-table tbody");
  var countEl = document.getElementById("result-count");
  var moreHint = document.getElementById("more-hint");
  var LIMIT = 60;

  function esc(text) {
    return String(text == null ? "" : text).replace(/[&<>"']/g, function (c) {
      return {"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c];
    });
  }
  function fmtQ(q) {
    if (q == null) { return "–"; }
    return q < 0.001 ? q.toExponential(1) : q.toPrecision(2);
  }
  function rank(gene, query) {
    if (!query) { return gene.significant ? 0 : gene.honorable ? 1 : 2; }
    if (gene.gene === query) { return 0; }
    if (gene.gene.indexOf(query) === 0) { return 1; }
    return gene.gene.indexOf(query) >= 0 ? 2 : -1;
  }
  function row(gene) {
    var tags = "";
    if (gene.significant) { tags += '<span class="tag tag-sig">FDR-significant</span>'; }
    else if (gene.honorable) { tags += '<span class="tag tag-hm">Honorable mention</span>'; }
    if (!gene.analyzed) { tags += '<span class="tag tag-na">Not analysed</span>'; }
    var best = gene.best;
    var also = gene.also.map(function (entry) {
      var kind = {cohort: "cohort scan", gene: "single gene", pair: "gene pair"}[entry.kind];
      return '<a href="' + esc(entry.url) + '" title="' + esc(entry.label + ", " + entry.date) +
        '">' + esc(entry.study) + " · " + kind + "</a>";
    }).join("");
    return "<tr><td class=\"gene\">" + esc(gene.gene) + tags + "</td>" +
      '<td class="best"><a href="' + esc(best.url) + '">' + esc(best.label) + " →</a>" +
      '<span class="sub">' + esc(best.study) + " · " + esc(best.date) + "</span></td>" +
      '<td class="num">' + (gene.events == null ? "–" : esc(gene.events)) + "</td>" +
      '<td class="num">' + fmtQ(gene.q) + "</td>" +
      '<td class="also">' + (also || '<span class="sub">–</span>') + "</td></tr>";
  }
  function render() {
    var query = input.value.trim().toUpperCase();
    var matches = genes
      .map(function (gene) { return [rank(gene, query), gene]; })
      .filter(function (pair) { return pair[0] >= 0; })
      .sort(function (a, b) {
        if (a[0] !== b[0]) { return a[0] - b[0]; }
        var qa = a[1].q == null ? 2 : a[1].q, qb = b[1].q == null ? 2 : b[1].q;
        return qa !== qb ? qa - qb : a[1].gene.localeCompare(b[1].gene);
      })
      .map(function (pair) { return pair[1]; });
    tbody.innerHTML = matches.slice(0, LIMIT).map(row).join("") ||
      '<tr><td colspan="5" class="sub">No studied gene matches “' + esc(input.value) +
      "”.</td></tr>";
    countEl.textContent = query
      ? matches.length + (matches.length === 1 ? " gene" : " genes")
      : "";
    moreHint.hidden = matches.length <= LIMIT;
    moreHint.textContent = "Showing the first " + LIMIT + " of " + matches.length +
      " — keep typing to narrow it down.";
    return matches;
  }
  input.addEventListener("input", render);
  input.addEventListener("keydown", function (event) {
    if (event.key !== "Enter") { return; }
    var matches = render();
    if (matches.length) { window.location.href = matches[0].best.url; }
  });
  document.querySelectorAll(".chip[data-gene]").forEach(function (chip) {
    chip.addEventListener("click", function () {
      input.value = chip.getAttribute("data-gene");
      render();
      input.focus();
    });
  });
  var initial = new URLSearchParams(window.location.search).get("gene");
  if (initial) { input.value = initial; }
  render();
})();
"""
