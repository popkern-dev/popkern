// The history dock's drawing (scrubber.js owns the timeline and what
// is on screen): one card per step — its number, operation icon and
// label, flagged when its feature failed or warned — with the steps
// before the playhead built, the current one lit and the later ones
// dimmed; migration chapters; the label naming the step on screen; on
// a phone the step card under the strip (the step, its operation and a
// jump to its script line); and the empty state when there is no
// history to step through.
import { iconSvg, operationIcon } from "./icons.js";

const $ = (id) => document.getElementById(id);
// Operations drawn as minor (narrow) cards.
const MINOR = new Set(["sketch", "reference"]);
// Auto labels (`extrude-1`, `boolean-3`) are not named marks.
const AUTO_LABEL = /^[a-z_]+-\d+$/;

/** "Extrude", "Tapped hole", "Mate" for an operation name. */
export function operationTitle(operation) {
  const text = String(operation ?? "step").replace(/_/g, " ");
  return text.charAt(0).toUpperCase() + text.slice(1);
}

export class HistoryStrip {
  /** `go(index)` shows step `index` (1-based). */
  constructor(go) {
    this.list = $("scrubber-ticks");
    this.list.addEventListener("click", (e) => {
      const card = e.target.closest("[data-index]");
      if (card) go(Number(card.dataset.index));
    });
  }

  /** No timeline: say why, offering the live build when `offerLive`. */
  empty(text, offerLive = false) {
    $("scrubber").hidden = true;
    $("step-card").hidden = true;
    $("history-empty").hidden = false;
    $("history-empty-text").textContent = text;
    $("history").hidden = !offerLive;
    $("scrubber-label").textContent = "";
    this.list.replaceChildren();
  }

  /** Draw one card per step. `flags` maps a feature label to "failed"
   *  or "warn" (from the build report); `assembly` draws construction
   *  steps (mates as named marks). `chapters` are the migrations. */
  render(steps, { flags = new Map(), assembly = false, chapters = [], goChapter }) {
    $("history-empty").hidden = true;
    $("scrubber").hidden = false;
    const range = $("scrubber-range");
    range.min = "1";
    range.max = String(steps.length);
    range.value = String(steps.length);
    const width = String(steps.length).length;
    this.list.replaceChildren();
    for (const s of steps) {
      const li = document.createElement("li");
      li.className = "step tick";
      li.dataset.label = s.label;
      if (MINOR.has(s.operation)) li.classList.add("minor");
      if (s.suppressed) li.classList.add("suppressed");
      const flag = s.failed ? "failed" : flags.get(s.label) ?? null;
      if (flag) li.classList.add(flag);
      if (assembly ? s.kind === "mate" : !AUTO_LABEL.test(s.label) && !MINOR.has(s.operation)) li.classList.add("named");
      const card = document.createElement("button");
      card.type = "button";
      card.dataset.index = String(s.index);
      card.title = `${s.index}. ${s.label} (${s.operation})${flag === "failed" ? " — failed" : flag === "warn" ? " — with warnings" : ""}`;
      const index = document.createElement("span");
      index.className = "step-index";
      index.textContent = String(s.index).padStart(Math.max(2, width), "0");
      const icon = document.createElement("span");
      icon.className = "step-icon";
      // Constant icon markup; the label below is text, never markup.
      icon.innerHTML = operationIcon(s.operation);
      const name = document.createElement("span");
      name.className = "step-name";
      name.textContent = s.label;
      card.append(index, icon, name);
      if (flag) {
        const mark = document.createElement("span");
        mark.className = "step-flag";
        mark.textContent = flag === "failed" ? "✕" : "▲";
        mark.setAttribute("aria-hidden", "true");
        card.append(mark);
      }
      li.append(card);
      this.list.append(li);
    }
    const bar = $("scrubber-chapters");
    bar.replaceChildren();
    bar.hidden = chapters.length < 2;
    for (const [k, c] of chapters.entries()) {
      const chip = document.createElement("button");
      chip.type = "button";
      chip.className = "chapter";
      chip.dataset.first = String(c.first);
      chip.dataset.last = String(c.last);
      chip.textContent = c.name;
      chip.title = `${c.name}: steps ${c.first}–${c.last}`;
      chip.addEventListener("click", () => goChapter(c, k));
      bar.append(chip);
    }
  }

  /** Light step `at` of `steps`: the cards, the label, the range, the
   *  chapter, and the phone's step card. `line` is the step's line in
   *  the open script (null: not in it); `note` trails the label. */
  paint(steps, at, { line = null, note = "" } = {}) {
    const s = steps[at - 1];
    if (!s) return;
    const cards = this.list.children;
    for (let i = 0; i < cards.length; i++) {
      cards[i].classList.toggle("current", i === at - 1);
      cards[i].classList.toggle("after", i > at - 1);
    }
    this.reveal(cards[at - 1]);
    $("scrubber-range").value = String(at);
    for (const chip of $("scrubber-chapters").children) {
      chip.classList.toggle("current", at >= Number(chip.dataset.first) && at <= Number(chip.dataset.last));
    }
    const where = s.migration ? ` · ${s.migration}` : "";
    const lineText = line != null ? ` · line ${line}` : "";
    $("scrubber-label").textContent = `${s.label} · step ${at} of ${steps.length}${where}${lineText}${note}`;
    const card = $("step-card");
    card.hidden = false;
    card.querySelector(".step-card-icon").innerHTML = operationIcon(s.operation);
    card.querySelector(".step-card-name").textContent = s.label;
    card.querySelector(".step-card-what").textContent = `${operationTitle(s.operation)} · ${at}/${steps.length}${where}`;
    const jump = $("step-card-line");
    jump.hidden = line == null;
    if (line == null) jump.replaceChildren();
    else {
      jump.textContent = `line ${line}`;
      jump.insertAdjacentHTML("beforeend", iconSvg("next"));
    }
    jump.dataset.line = line == null ? "" : String(line);
  }

  /** Scroll the strip (only the strip) so `card` is in view. */
  reveal(card) {
    if (!card) return;
    const strip = this.list;
    const left = card.offsetLeft;
    const right = left + card.offsetWidth;
    if (left < strip.scrollLeft || right > strip.scrollLeft + strip.clientWidth) {
      strip.scrollLeft = Math.max(0, left - (strip.clientWidth - card.offsetWidth) / 2);
    }
    // A dock dragged to another height wraps the cards into rows that
    // scroll down (panel-resize.js).
    const top = card.offsetTop;
    const bottom = top + card.offsetHeight;
    if (top < strip.scrollTop || bottom > strip.scrollTop + strip.clientHeight) {
      strip.scrollTop = Math.max(0, top - (strip.clientHeight - card.offsetHeight) / 2);
    }
  }
}
