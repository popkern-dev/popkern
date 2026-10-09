// The feature tree: one row per feature of the build report with its
// operation, outcome and time. A blend row states the continuity it
// measured against the continuity asked for, and how it knows.

const $ = (id) => document.getElementById(id);

/** Render `envelope.report` into the report table, scoped to `part`
 *  (an opened assembly part) when it is not null. */
export function renderReport(envelope, part) {
  const table = $("report");
  table.innerHTML = "";
  const scoped = envelope.assembly != null;
  const features = (envelope.report?.features ?? []).filter((f) => part == null || f.part === part);
  for (const f of features) {
    const row = document.createElement("tr");
    const label = f.feature?.label ?? f.label ?? "?";
    row.dataset.feature = label;
    const operation = f.feature?.operation ?? f.operation ?? "?";
    const outcome = typeof f.outcome === "string" ? f.outcome : Object.keys(f.outcome ?? {})[0] ?? "?";
    const elapsed = typeof f.elapsed_ms === "number" ? `${f.elapsed_ms.toFixed(1)} ms` : "";
    const downgraded = f.blend?.downgraded === true;
    const continuity = f.blend?.continuity ?? null;
    let verdict = outcome;
    if (continuity && continuity.requested) {
      // A blend row always states what was measured against what
      // was asked — a downgrade is never a bare success mark — and
      // how: certified over every seam, or from samples (pending
      // certification, or where a certificate could not decide).
      const how = continuity.deferred_seams > 0 ? "sampled, certifying…" : continuity.certified_achieved ? "certified" : "measured";
      verdict = downgraded
        ? `downgraded: ${continuity.achieved} ${how}, ${continuity.requested} requested`
        : `${outcome} · ${continuity.achieved} ${how}`;
    }
    const shown = scoped && part == null && f.part ? `${f.part}: ${label}` : label;
    // Labels and messages come from the script: text, never markup.
    const verdictClass = downgraded ? "cell-warn" : outcome === "Served" ? "cell-ok" : "cell-bad";
    for (const [text, className] of [[shown, ""], [operation, ""], [verdict, verdictClass], [elapsed, ""]]) {
      const cell = document.createElement("td");
      cell.textContent = text;
      if (className) cell.className = className;
      row.append(cell);
    }
    if (Array.isArray(f.warnings) && f.warnings.length) {
      const warn = document.createElement("tr");
      warn.className = "row-warn";
      // A structured warning reads as its message; the rest of its
      // payload stays in the report (export, console), not the table.
      const cell = document.createElement("td");
      cell.colSpan = 4;
      cell.textContent = f.warnings.map((w) => (typeof w === "string" ? w : (w?.message ?? JSON.stringify(w)))).join("; ");
      warn.append(cell);
      table.append(row, warn);
    } else {
      table.append(row);
    }
  }
  const summary = envelope.report
    ? `${envelope.report.failures} failure(s), ${envelope.report.downgraded_blends} downgraded blend(s)${envelope.report.is_clean ? ", clean" : ", with warnings"}`
    : "no report";
  $("report-summary").textContent = summary;
}

/** Deferred seams across a report's blend rows (the certification
 *  work still to do after the build served). */
export function deferredSeams(report) {
  return (report?.features ?? []).reduce((n, f) => n + (f.blend?.continuity?.deferred_seams ?? 0), 0);
}
