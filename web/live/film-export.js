// Film export (docs/architecture/presentation-films.md §7): frame by frame
// at a fixed timestep on a renderer of the export size — the viewer's
// scene and a film camera — the caption composited on a 2D canvas, then
// encoded: MP4 (H.264) and WebM (VP9) through WebCodecs and the vendored
// muxers, GIF through gifenc, or one PNG frame. Deterministic: frame k
// is the film at k / fps, whatever the machine's speed.
import * as THREE from "./vendor/three.module.min.js";
import { ArrayBufferTarget as Mp4Target, Muxer as Mp4Muxer } from "./vendor/mp4-muxer.js";
import { ArrayBufferTarget as WebmTarget, Muxer as WebmMuxer } from "./vendor/webm-muxer.js";
import { GIFEncoder, applyPalette, quantize } from "./vendor/gifenc.esm.js";
import { FilmDetail } from "./film-lod.js";
import { frameTimes } from "./film-state.js";

/** The export menu (the design's preset table). */
export const PRESETS = {
  "mp4-1080p": { label: "MP4 1080p (Slides, Discord, social)", format: "mp4", size: [1920, 1080], bitrate: 8e6 },
  "mp4-vertical": { label: "MP4 vertical 9:16 (Shorts, Reels)", format: "mp4", size: [1080, 1920], bitrate: 8e6 },
  "mp4-720p": { label: "MP4 720p (smaller)", format: "mp4", size: [1280, 720], bitrate: 4e6 },
  webm: { label: "WebM 1080p (web)", format: "webm", size: [1920, 1080], bitrate: 6e6 },
  gif: { label: "GIF 640 wide (READMEs)", format: "gif", width: 640, fps: 15 },
  png: { label: "PNG of this frame", format: "png" },
};

/** The frame size of `preset` for `film` (a GIF keeps the film's
 *  aspect at its width). */
export function presetSize(preset, film) {
  if (preset.size) return preset.size;
  if (preset.width) {
    const h = Math.round((preset.width * film.size[1]) / film.size[0] / 2) * 2;
    return [preset.width, h];
  }
  return film.size;
}

/** What this browser can export: `{ mp4, webm }` codec support. */
export async function codecSupport() {
  if (typeof VideoEncoder === "undefined") return { mp4: false, webm: false };
  const test = async (codec) => {
    try {
      return (await VideoEncoder.isConfigSupported({ codec, width: 1920, height: 1080, bitrate: 8e6, framerate: 30 })).supported === true;
    } catch {
      return false;
    }
  };
  return { mp4: await test("avc1.640028"), webm: await test("vp09.00.40.08") };
}

/** The caption drawn into a frame, scaled to its width. */
export function drawCaption(ctx, caption, width, height, background) {
  if (!caption || caption.opacity <= 0) return;
  const unit = Math.min(width, height) / 1080;
  const dark = luminance(background) > 0.5;
  ctx.save();
  ctx.globalAlpha = caption.opacity;
  ctx.fillStyle = dark ? "#111318" : "#f5f7fa";
  ctx.textBaseline = "top";
  const pad = 56 * unit;
  const titleSize = (caption.kind === "title" ? 64 : 48) * unit;
  ctx.font = `600 ${titleSize}px Inter, "Segoe UI", system-ui, sans-serif`;
  ctx.fillText(caption.title, pad, pad, width - 2 * pad);
  if (caption.subtitle) {
    ctx.globalAlpha = caption.opacity * 0.8;
    ctx.font = `400 ${30 * unit}px Inter, "Segoe UI", system-ui, sans-serif`;
    ctx.fillText(caption.subtitle, pad, pad + titleSize * 1.25, width - 2 * pad);
  }
  ctx.restore();
}

function luminance(hex) {
  const n = parseInt(String(hex).slice(1), 16);
  if (!Number.isFinite(n)) return 1;
  return (0.2126 * ((n >> 16) & 255) + 0.7152 * ((n >> 8) & 255) + 0.0722 * (n & 255)) / 255;
}

/** Export the open film with `presetKey`: resolves to `{ bytes, mime,
 *  extension, frames }`. `progress(done, total)` is called per frame;
 *  `signal` (an AbortSignal) stops it. */
export async function exportFilm(page, presetKey, { progress = () => {}, signal = null } = {}) {
  const player = page.film;
  const film = player.film;
  if (!film) throw new Error("no film is open");
  const preset = PRESETS[presetKey];
  if (!preset) throw new Error(`no export preset ${presetKey}`);
  const [width, height] = presetSize(preset, film);
  const viewer = page.viewer;
  const canvas = document.createElement("canvas");
  canvas.width = width;
  canvas.height = height;
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, stencil: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(1);
  renderer.setSize(width, height, false);
  renderer.localClippingEnabled = true;
  player.look.prepare(renderer);
  const camera = new THREE.PerspectiveCamera(30, width / height, 0.01, 1000);
  const frame = document.createElement("canvas");
  frame.width = width;
  frame.height = height;
  const ctx = frame.getContext("2d", { willReadFrequently: preset.format === "gif" });
  // Screen-space detail: each frame draws every part at the triangle
  // density its size on screen needs (film-lod.js).
  const detail = new FilmDetail(viewer);
  const startTime = player.time;
  const wasPlaying = player.playing;
  player.pause();
  const draw = async (t) => {
    await player.apply(t, { camera, aspect: width / height, render: false });
    // The cut outline follows a moving plane once per animation frame on
    // screen; an export frame needs it now.
    if (viewer.section.active) viewer.section.outlineNow();
    // Fat lines (section outlines) measure their width in this
    // renderer's pixels; the gizmo stays off.
    for (const line of viewer.section.overlays()) renderer.getDrawingBufferSize(line.material.uniforms.viewport.value);
    detail.apply(camera, height);
    renderer.render(viewer.scene, camera);
    ctx.drawImage(canvas, 0, 0);
    drawCaption(ctx, player.caption, width, height, film.background);
  };
  try {
    if (preset.format === "png") {
      await draw(startTime);
      const blob = await new Promise((resolve) => frame.toBlob(resolve, "image/png"));
      return { bytes: new Uint8Array(await blob.arrayBuffer()), mime: "image/png", extension: "png", frames: 1 };
    }
    const fps = preset.fps ?? film.fps;
    const times = frameTimes(film, fps);
    if (preset.format === "gif") return await encodeGif(times, fps, draw, ctx, width, height, progress, signal);
    return await encodeVideo(preset, times, fps, draw, frame, width, height, progress, signal);
  } finally {
    detail.restore();
    renderer.dispose();
    for (const line of viewer.section.overlays()) viewer.renderer.getDrawingBufferSize(line.material.uniforms.viewport.value);
    await player.apply(startTime);
    if (wasPlaying) player.play();
  }
}

async function encodeVideo(preset, times, fps, draw, frame, width, height, progress, signal) {
  const mp4 = preset.format === "mp4";
  const target = mp4 ? new Mp4Target() : new WebmTarget();
  const muxer = mp4
    ? new Mp4Muxer({ target, video: { codec: "avc", width, height, frameRate: fps }, fastStart: "in-memory" })
    : new WebmMuxer({ target, video: { codec: "V_VP9", width, height, frameRate: fps } });
  let failure = null;
  const encoder = new VideoEncoder({
    output: (chunk, meta) => muxer.addVideoChunk(chunk, meta),
    error: (e) => (failure = e),
  });
  // High profile 4.0 holds 1080p and 1080x1920 at 30 fps.
  const codec = mp4 ? (width * height > 1280 * 720 ? "avc1.640028" : "avc1.64001f") : "vp09.00.40.08";
  const config = { codec, width, height, bitrate: preset.bitrate, framerate: fps, latencyMode: "quality" };
  if (mp4) config.avc = { format: "avc" };
  encoder.configure(config);
  const usec = 1e6 / fps;
  for (let k = 0; k < times.length; k++) {
    if (signal?.aborted) throw new DOMException("export stopped", "AbortError");
    if (failure) throw failure;
    await draw(times[k]);
    const videoFrame = new VideoFrame(frame, { timestamp: Math.round(k * usec), duration: Math.round(usec) });
    encoder.encode(videoFrame, { keyFrame: k % (2 * fps) === 0 });
    videoFrame.close();
    while (encoder.encodeQueueSize > 4) await new Promise((r) => setTimeout(r, 0));
    progress(k + 1, times.length);
  }
  await encoder.flush();
  if (failure) throw failure;
  encoder.close();
  muxer.finalize();
  return {
    bytes: new Uint8Array(target.buffer),
    mime: mp4 ? "video/mp4" : "video/webm",
    extension: mp4 ? "mp4" : "webm",
    frames: times.length,
  };
}

async function encodeGif(times, fps, draw, ctx, width, height, progress, signal) {
  const gif = GIFEncoder();
  const delay = Math.round(1000 / fps);
  for (let k = 0; k < times.length; k++) {
    if (signal?.aborted) throw new DOMException("export stopped", "AbortError");
    await draw(times[k]);
    const { data } = ctx.getImageData(0, 0, width, height);
    const palette = quantize(data, 256);
    gif.writeFrame(applyPalette(data, palette), width, height, { palette, delay });
    progress(k + 1, times.length);
    // Yield so the page stays responsive between frames.
    await new Promise((r) => setTimeout(r, 0));
  }
  gif.finish();
  return { bytes: gif.bytes(), mime: "image/gif", extension: "gif", frames: times.length };
}
