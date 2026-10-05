// The viewport overlay for a build that did not serve: the outcome
// class (script error, kernel decline, gate failure), its rendered
// message with the remediation for a budget stop, and a jump to the
// script line it names. The last good model stays on screen beneath.
import { describeLimit } from "./time-limit.js";
import { setStatus } from "./status-line.js";

const $ = (id) => document.getElementById(id);

/** Show `outcome` (a non-served envelope or a failed reply's outcome)
 *  over the viewport, mark its line in `editor`, and state it on the
 *  status line. `timeLimit` words the budget remediation;
 *  `keptModel` says the last good model is still on screen. */
export function showOutcome(outcome, { editor, timeLimit, keptModel, openScript = null }) {
  const overlay = $("overlay");
  overlay.hidden = false;
  const lines = [];
  if (outcome.outcome === "gate_failure") {
    lines.push(outcome.step === "cancelled" ? outcome.detail : `${outcome.step} failed: ${outcome.detail}`);
  } else {
    lines.push(outcome.rendered ?? `${outcome.kind}`);
  }
  // A project migration that failed: name it, and the upstream
  // migrations changed since it last built (where a vanished label
  // usually comes from).
  const failed = outcome.project?.error;
  if (failed?.migration) {
    lines.push("", `in ${failed.migration}${failed.line ? `, line ${failed.line}` : ""}`);
    if (failed.upstream_changed?.length) lines.push(`changed upstream since it last built: ${failed.upstream_changed.join(", ")}`);
  }
  if (outcome.kind === "BudgetExceeded" || outcome.step === "cancelled") {
    lines.push("", timeLimit == null ? "Press run to try again." : `The time limit is ${describeLimit(timeLimit)}: raise it in the top bar and run again.`);
  }
  $("overlay-text").textContent = lines.join("\n");
  $("overlay-class").textContent = outcome.outcome.replace("_", " ");
  const script = /at ([A-Za-z0-9_.-]+\.rhai):\d+:\d+:/.exec(outcome.rendered ?? "");
  if (script && openScript) openScript(script[1]);
  const position = /:(\d+):(\d+):/.exec(outcome.rendered ?? "");
  if (position) {
    editor.mark(Number(position[1]));
    $("overlay-goto").hidden = false;
    $("overlay-goto").onclick = () => editor.goTo(Number(position[1]), Number(position[2]));
  } else {
    editor.mark(null);
    $("overlay-goto").hidden = true;
  }
  setStatus(`${outcome.outcome.replace("_", " ")}: ${outcome.kind ?? outcome.step}${keptModel ? " (showing the last good model)" : ""}`, true);
}

/** Hide the overlay and clear the editor's error mark (a build served). */
export function clearOutcome(editor) {
  editor.mark(null);
  $("overlay").hidden = true;
}
