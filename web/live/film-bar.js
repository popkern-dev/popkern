// The film bar (docs/architecture/presentation-films.md §6): while a film
// is open it takes the history dock's place — play / pause, the time, a
// scrub range with a marker per step, the caption over the model, a
// guide showing the film's frame, and the export menu.
import { PRESETS, codecSupport, exportFilm } from "./film-export.js";
import { saveBytes } from "./exports.js";

const $ = (id) => document.getElementById(id);
const RANGE_STEPS = 1000;

export function wireFilmBar(page) {
  const player = page.film;
  const range = $("film-range");
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

  $("film-play").addEventListener("click", () => (player.playing ? player.pause() : player.play()));
  range.addEventListener("input", () => {
    if (!player.active) return;
    player.seek((Number(range.value) / RANGE_STEPS) * player.film.duration);
  });
  $("film-close").addEventListener("click", () => {
    // Back to the geometry: open the last migration.
    const migrations = Object.keys(page.projects.current?.files ?? {}).filter((p) => p.startsWith("migrations/")).sort();
    if (migrations.length) page.projects.select(migrations.at(-1));
    else player.close();
  });
  let controller = null;
  $("film-export").addEventListener("click", async () => {
    if (controller) {
      controller.abort();
      return;
    }
    if (!player.active) return;
    controller = new AbortController();
    const button = $("film-export");
    button.textContent = "stop";
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
  // Space toggles play while the bar has focus.
  $("film-bar").addEventListener("keydown", (e) => {
    if (e.key === " " && e.target.tagName !== "SELECT") {
      e.preventDefault();
      if (player.playing) player.pause();
      else player.play();
    }
  });

  const paint = () => {
    const open = player.active;
    $("film-bar").hidden = !open;
    $("film-caption").hidden = !open;
    $("film-frame").hidden = !open;
    if (!open) return;
    const film = player.film;
    document.body.classList.toggle("film-playing", player.playing);
    $("film-play").setAttribute("aria-label", player.playing ? "Pause the film" : "Play the film");
    $("film-name").textContent = film.name;
    $("film-time").textContent = `${player.time.toFixed(1)} / ${film.duration.toFixed(1)} s`;
    if (document.activeElement !== range) range.value = String(Math.round((player.time / film.duration) * RANGE_STEPS));
    const marks = $("film-marks");
    if (marks.dataset.film !== `${film.name}:${film.duration}:${film.steps.length}`) {
      marks.dataset.film = `${film.name}:${film.duration}:${film.steps.length}`;
      marks.replaceChildren(
        ...film.steps.map((s) => {
          const mark = document.createElement("button");
          mark.type = "button";
          mark.className = `film-mark film-mark-${s.kind}`;
          mark.style.left = `${(s.start / film.duration) * 100}%`;
          mark.title = s.title || "pause";
          mark.addEventListener("click", () => {
            player.seek(s.start + 1e-6);
            if (s.line) page.editor.goTo(s.line, 1);
          });
          return mark;
        }),
      );
    }
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
  paint();
}

/** The film's frame inside the viewport (the pane can be any shape). */
function placeFrame(film) {
  const pane = document.querySelector(".pane-viewer");
  const guide = $("film-frame");
  const caption = $("film-caption");
  if (!pane) return;
  const w = pane.clientWidth;
  const h = pane.clientHeight;
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
