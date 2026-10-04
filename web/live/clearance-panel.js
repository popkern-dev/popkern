// The interference and clearance list (PK-630): one row per occurrence
// pair the check found — interfering with its overlap volume, touching,
// or clear with its certified distance — worst first, as the worker
// orders them. A row selects its two occurrences.

const $ = (id) => document.getElementById(id);

/** A clearance row's verdict. */
export function clearanceText(pair) {
  const n = (x) => (Math.abs(x) >= 1e-3 || x === 0 ? x.toPrecision(4) : x.toExponential(3));
  if (pair.status === "interfering") return pair.volume != null ? `interfering: ${n(pair.volume)} overlap` : `interfering (volume unknown: ${pair.decline ?? "declined"})`;
  if (pair.status === "touching") return "touching";
  return `clear: ${n(pair.distance)}${pair.refined ? "" : " (mesh)"}, ≥ ${n(pair.lower_bound)}`;
}

/** Render an `interference` reply for the check run at `clearance`;
 *  clicking a row calls `select([a, b])`. Returns the summary text, or
 *  null when the check failed (the summary then says why). */
export function renderClearance(reply, clearance, select) {
  const list = $("clearance-list");
  list.innerHTML = "";
  list.hidden = false;
  $("clearance-summary").hidden = false;
  if (!reply.ok) {
    $("clearance-summary").textContent = reply.outcome?.detail ?? reply.outcome?.rendered ?? "the check failed";
    return null;
  }
  const { pairs, occurrences, candidate_pairs } = reply.report;
  const clashes = pairs.filter((p) => p.status === "interfering").length;
  $("clearance-summary").textContent = `${occurrences} occurrence(s), ${candidate_pairs} pair(s) near: ${clashes} interfering, ${pairs.length - clashes} within ${clearance}`;
  for (const pair of pairs) {
    const row = document.createElement("li");
    row.className = `tree-row ${pair.status === "interfering" ? "clash" : pair.status === "touching" ? "contact" : ""}`;
    const name = document.createElement("span");
    name.className = "tree-name";
    name.textContent = `${pair.a} ↔ ${pair.b}`;
    const what = document.createElement("span");
    what.className = "tree-kind";
    what.textContent = clearanceText(pair);
    row.title = (pair.faces ?? []).join("  ↔  ") || what.textContent;
    row.append(name, what);
    row.addEventListener("click", () => select([pair.a, pair.b]));
    list.append(row);
  }
  return $("clearance-summary").textContent;
}

/** The clearance button: build the live model if needed, run the check
 *  at the clearance the viewer bar names, and list the pairs. Resolves
 *  to the worker's reply (null when no live model could be built). */
export async function checkClearance(page) {
  if (!(await page.ensureLive())) return null;
  const clearance = Math.max(0, Number($("clearance").value) || 0);
  page.status(`checking clearance within ${clearance}…`);
  const reply = await page.host.interference(clearance);
  const summary = renderClearance(reply, clearance, (paths) => page.viewer.select(paths));
  if (summary != null) page.status(`clearance checked: ${summary}`);
  return reply;
}
