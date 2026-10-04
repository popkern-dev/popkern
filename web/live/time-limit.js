// The build time limit the user picks, in seconds (null = none). The
// kernel's wall budget stops the script at the limit less a grace, and
// the page terminates the worker at the limit itself (the backstop for
// a long kernel call the script cannot interrupt). The operation
// budget, a guard against runaway script loops, scales with the time.
import { readTimeLimit, writeTimeLimit } from "./local-store.js";

/** The presets the top bar offers, shortest first; null is "no limit". */
export const TIME_LIMITS = [10, 24, 60, 120, 300, 900, null];
/** The limit a first visit builds under. */
export const DEFAULT_TIME_LIMIT = 24;
const OPERATIONS_PER_SECOND = 1_000_000;

/** The worker budgets (PROTOCOL.md `init` options) for `limit`. */
export function budgetFor(limit) {
  if (limit == null) return { max_operations: null, timeout_ms: null, continuity_certification: "Deferred" };
  const kernel = limit - Math.min(4, limit / 6);
  return { max_operations: Math.round(kernel * OPERATIONS_PER_SECOND), timeout_ms: kernel * 1000, continuity_certification: "Deferred" };
}

/** "24 s", "5 min" or "no limit". */
export function describeLimit(limit) {
  if (limit == null) return "no limit";
  return limit < 60 ? `${limit} s` : `${limit / 60} min`;
}

/** The remembered limit, or the default when none (or an unknown
 *  value) is stored. */
export function storedLimit() {
  const stored = readTimeLimit();
  if (stored === "none") return null;
  const value = Number(stored);
  return TIME_LIMITS.includes(value) ? value : DEFAULT_TIME_LIMIT;
}

/** Change `page`'s build time limit for later runs (remembered here). */
export async function applyTimeLimit(page, limit) {
  page.timeLimit = limit;
  writeTimeLimit(limit);
  const budget = budgetFor(limit);
  // A restarted worker (stop, crash) starts under the new limit too.
  page.host.options = budget;
  if (page.host.info) await page.host.setBudget(budget);
}

/** Fill the time-limit `select` with the presets, `current` selected. */
export function renderTimeLimitOptions(select, current) {
  select.innerHTML = "";
  for (const limit of TIME_LIMITS) {
    const option = document.createElement("option");
    option.value = limit == null ? "none" : String(limit);
    option.textContent = describeLimit(limit) + (limit === DEFAULT_TIME_LIMIT ? " (default)" : "");
    option.selected = limit === current;
    select.append(option);
  }
}
