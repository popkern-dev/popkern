// The film stage (docs/architecture/presentation-films.md §6): while a
// film is open the model view shows it, beside the page's usual panels
// (the history dock, the script and the inspector stay as they were).
// The stage has one way out (Close film, top right; Esc), one control
// bar along its bottom — play / pause, the time, a timeline of the
// film's steps as chapters, the current step and the export menu — the
// caption over the film, a guide showing the film's frame above the
// bar, and while the film is being prepared, what it waits on.
import { PRESETS, codecSupport, exportFilm } from "./film-export.js";
import { saveBytes } from "./exports.js";

const $ = (id) => document.getElementById(id);

/** `seconds` as m:ss. */
export function clockText(seconds) {
  const s = Math.max(0, Math.floor(seconds + 1e-6));
  return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, "0")}`;
}

/** The step of `film` playing at time `t` (the last one at the end). */
function stepAt(film, t) {
  return film.steps.find((s) => t < s.start + s.duration) ?? film.steps.at(-1) ?? null;
}

export function wireFilmBar(page) {
  const player = page.film;
  const preset = $("film-preset");
  for (const [key, p] of Object.entries(PRESETS)) {
    const option = document.createElement("option");
    option.value = key;
    option.textContent = p.label;
    preset.append(option);
  }
  codecSupport().then((support) => {
    for (const option of preset.options) {
      const format = PRESETS[option.value].format;
      if ((format === "mp4" && !support.mp4) || (format === "webm" && !support.webm)) {
        option.disabled = true;
        option.textContent += " — not in this browser";
      }
    }
    if (preset.selectedOptions[0]?.disabled) preset.value = "gif";
  });

  const toggle = () => (player.playing ? player.pause() : player.play());
  $("film-play").addEventListener("click", toggle);
  const closeFilm = () => {
    // Back to the geometry: open the last migration.
    const migrations = Object.keys(page.projects.current?.files ?? {}).filter((p) => p.startsWith("migrations/")).sort();
    if (page.projects.filmOpen && migrations.length) page.projects.select(migrations.at(-1));
    // However the film was opened (or while it is being prepared), it closes.
    player.close();
  };
  $("film-exit").addEventListener("click", closeFilm);
  const typing = (t) => t && (t.isContentEditable || ["INPUT", "TEXTAREA", "SELECT"].includes(t.tagName) || t.closest?.("#editor"));
  document.addEventListener("keydown", (e) => {
    if (!(player.active || player.preparing) || typing(e.target)) return;
    if (e.key === "Escape") {
      e.preventDefault();
      closeFilm();
    } else if (e.key === " " && player.active && e.target.closest?.(".pane-viewer")) {
      e.preventDefault();
      toggle();
    }
  });

  wireTimeline(player);
  wireExport(page, preset);

  // The preparing card counts the seconds it has waited.
  let ticker = 0;
  const paint = () => {
    const open = player.active;
    const preparing = player.preparing;
    const staged = open || Boolean(preparing);
    document.body.classList.toggle("film-stage", staged);
    $("film-exit").hidden = !staged;
    $("film-bar").hidden = !open;
    $("film-caption").hidden = !open;
    $("film-frame").hidden = !open;
    const loading = $("film-loading");
    loading.hidden = !preparing;
    if (preparing) {
      loading.classList.toggle("film-loading-quiet", preparing.quiet);
      $("film-loading-title").textContent = preparing.quiet ? "Updating the film" : "Preparing the film";
      const waited = Math.floor((performance.now() - preparing.started) / 1000);
      $("film-loading-stage").textContent = waited >= 1 ? `${preparing.stage}… ${waited} s` : `${preparing.stage}…`;
      if (!ticker) ticker = setInterval(paint, 500);
    } else if (ticker) {
      clearInterval(ticker);
      ticker = 0;
    }
    page.rail?.renderFilms?.();
    if (!open) return;
    const film = player.film;
    document.body.classList.toggle("film-playing", player.playing);
    $("film-play").setAttribute("aria-label", player.playing ? "Pause the film" : "Play the film");
    $("film-time").textContent = `${clockText(player.time)} / ${clockText(film.duration)}`;
    const step = stepAt(film, player.time);
    $("film-step").textContent = step?.title || film.name;
    paintTimeline(player);
    const caption = player.caption;
    const box = $("film-caption");
    box.style.opacity = caption ? String(caption.opacity) : "0";
    box.classList.toggle("film-caption-title", caption?.kind === "title");
    box.querySelector("strong").textContent = caption?.title ?? "";
    box.querySelector("span").textContent = caption?.subtitle ?? "";
    box.classList.toggle("film-caption-dark", luminance(film.background) > 0.5);
    placeFrame(film);
  };
  player.onChange(paint);
  // The frame follows the viewport: a phone's sheet tabs, a rotation or
  // a window resize change the pane's shape under an open film.
  // The canvas resizing (the viewer then renders at the new aspect)
  // re-frames the film too.
  new ResizeObserver(() => player.active && (placeFrame(player.film), player.apply(player.time))).observe($("viewer"));
  paint();
}

/** The timeline: one chapter per step, sized by its duration, filled up
 *  to the time; a playhead; drag or click to seek (a chapter's start
 *  jumps the script to its line too), arrow keys step a second. */
function wireTimeline(player) {
  const line = $("film-timeline");
  const seekAt = (clientX) => {
    if (!player.active) return;
    const t = timeAtX(player.film, clientX);
    if (t != null) player.seek(t);
  };
  let dragging = false;
  line.addEventListener("pointerdown", (e) => {
    dragging = true;
    // Capture keeps a drag that leaves the bar seeking; a pointer the
    // browser no longer tracks (a pen lifted mid-event) cannot be.
    try {
      line.setPointerCapture(e.pointerId);
    } catch {
      // seek on this press alone
    }
    line.classList.add("film-timeline-drag");
    seekAt(e.clientX);
  });
  line.addEventListener("pointermove", (e) => dragging && seekAt(e.clientX));
  const release = () => {
    dragging = false;
    line.classList.remove("film-timeline-drag");
  };
  line.addEventListener("pointerup", release);
  line.addEventListener("pointercancel", release);
  line.addEventListener("keydown", (e) => {
    if (!player.active) return;
    const step = { ArrowLeft: -1, ArrowRight: 1, PageDown: -5, PageUp: 5 }[e.key];
    if (e.key === "Home" || e.key === "End") player.seek(e.key === "Home" ? 0 : player.film.duration);
    else if (step != null) player.seek(player.time + step);
    else return;
    e.preventDefault();
  });
}

/** The film time under client x on the timeline: the chapter under it
 *  and the fraction across that chapter (gaps snap to the nearer). */
function timeAtX(film, clientX) {
  const chapters = [...$("film-chapters").children];
  if (!chapters.length) return null;
  let best = null;
  for (const [i, el] of chapters.entries()) {
    const r = el.getBoundingClientRect();
    const step = film.steps[i];
    const d = clientX < r.left ? r.left - clientX : clientX > r.right ? clientX - r.right : 0;
    if (!best || d < best.d) best = { d, t: step.start + step.duration * Math.min(1, Math.max(0, (clientX - r.left) / Math.max(1, r.width))) };
  }
  return Math.min(film.duration, best.t);
}

function paintTimeline(player) {
  const film = player.film;
  const box = $("film-chapters");
  const key = `${film.name}:${film.duration}:${film.steps.map((s) => `${s.start}/${s.duration}/${s.title}`).join(",")}`;
  if (box.dataset.film !== key) {
    box.dataset.film = key;
    box.replaceChildren(
      ...film.steps.map((s) => {
        const chapter = document.createElement("span");
        chapter.className = `film-chapter film-chapter-${s.kind}`;
        chapter.style.flexGrow = String(Math.max(s.duration, 1e-3));
        chapter.title = s.title || (s.kind === "hold" ? "pause" : "");
        chapter.append(document.createElement("i"));
        return chapter;
      }),
    );
  }
  const line = $("film-timeline");
  line.setAttribute("aria-valuemax", film.duration.toFixed(1));
  line.setAttribute("aria-valuenow", player.time.toFixed(1));
  line.setAttribute("aria-valuetext", `${clockText(player.time)} of ${clockText(film.duration)}`);
  let head = 0;
  const lineBox = line.getBoundingClientRect();
  [...box.children].forEach((chapter, i) => {
    const s = film.steps[i];
    const f = Math.min(1, Math.max(0, (player.time - s.start) / Math.max(s.duration, 1e-9)));
    chapter.firstChild.style.width = `${f * 100}%`;
    chapter.classList.toggle("film-chapter-current", f > 0 && f < 1);
    if (player.time >= s.start - 1e-9) {
      const r = chapter.getBoundingClientRect();
      head = r.left - lineBox.left + f * r.width;
    }
  });
  line.querySelector(".film-playhead").style.transform = `translateX(${head}px)`;
}

/** The export menu: a format and Export (stop while it renders). */
function wireExport(page, preset) {
  const player = page.film;
  let controller = null;
  $("film-export").addEventListener("click", async () => {
    if (controller) {
      controller.abort();
      return;
    }
    if (!player.active) return;
    controller = new AbortController();
    const button = $("film-export");
    button.textContent = "Stop";
    const state = $("film-export-state");
    const started = performance.now();
    try {
      const out = await exportFilm(page, preset.value, {
        signal: controller.signal,
        progress: (done, total) => (state.textContent = `frame ${done} of ${total}`),
      });
      const name = `${player.film.name}.${out.extension}`;
      saveBytes(out.bytes, name, out.mime);
      const mb = out.bytes.length / 1e6;
      state.textContent = `${name}: ${mb.toFixed(1)} MB, ${out.frames} frame(s) in ${((performance.now() - started) / 1000).toFixed(1)} s`;
    } catch (e) {
      state.textContent = e.name === "AbortError" ? "export stopped" : `export failed: ${e.message}`;
    } finally {
      controller = null;
      button.textContent = "Export";
    }
  });
}

/** The film's frame inside the canvas (the pane can be any shape; on a
 *  film's stage the canvas stops above the control bar's strip): the
 *  frame element covers the canvas and clips its hole's shade, so the
 *  dimming outside the film's frame stays inside the model view. */
function placeFrame(film) {
  const canvas = $("viewer");
  const guide = $("film-frame").querySelector(".film-hole");
  const caption = $("film-caption");
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (!w || !h) return;
  const aspect = film.size[0] / film.size[1];
  let fw = w;
  let fh = w / aspect;
  if (fh > h) {
    fh = h;
    fw = h * aspect;
  }
  const left = (w - fw) / 2;
  const top = (h - fh) / 2;
  Object.assign(guide.style, { left: `${left}px`, top: `${top}px`, width: `${fw}px`, height: `${fh}px` });
  const unit = Math.min(fw, fh) / 1080;
  Object.assign(caption.style, { left: `${left + 56 * unit}px`, top: `${top + 56 * unit}px`, fontSize: `${Math.max(9, 30 * unit)}px` });
}

function luminance(hex) {
  const n = parseInt(String(hex).slice(1), 16);
  if (!Number.isFinite(n)) return 1;
  return (0.2126 * ((n >> 16) & 255) + 0.7152 * ((n >> 8) & 255) + 0.0722 * (n & 255)) / 255;
}
