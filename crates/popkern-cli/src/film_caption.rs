//! Film captions on a frame (`film-export.js` `drawCaption`): the title
//! at the top left (64 px for a title step, 48 for a step, weight 600)
//! and the subtitle under it (30 px, 80 % opacity), scaled by the
//! frame's shorter side / 1080, in near-black on a light background
//! and near-white on a dark one, faded by the caption's opacity.
//!
//! The page draws with the browser's font (Inter, else the system UI
//! font); the command line uses `--font` or the first system font it
//! finds (Inter, DejaVu Sans, Liberation Sans, Noto Sans, Helvetica,
//! Arial, Segoe UI). Text wider than the frame's padded width shrinks
//! to fit (canvas `fillText`'s `maxWidth`).

use std::path::{Path, PathBuf};

use fontdue::{Font, FontSettings};
use popkern_rhai::film::state::Caption;

/// The caption fonts: a regular and a bold face.
pub struct Fonts {
    regular: Font,
    bold: Font,
    /// Where they came from (for the report).
    pub source: PathBuf,
}

/// `(regular, bold)` file names to look for, in preference order.
const CANDIDATES: &[(&str, &str)] = &[
    ("Inter-Regular.ttf", "Inter-SemiBold.ttf"),
    ("Inter-Regular.otf", "Inter-SemiBold.otf"),
    ("DejaVuSans.ttf", "DejaVuSans-Bold.ttf"),
    ("LiberationSans-Regular.ttf", "LiberationSans-Bold.ttf"),
    ("NotoSans-Regular.ttf", "NotoSans-Bold.ttf"),
    ("Helvetica.ttc", "Helvetica.ttc"),
    ("Arial.ttf", "Arial Bold.ttf"),
    ("arial.ttf", "arialbd.ttf"),
    ("segoeui.ttf", "seguisb.ttf"),
];

fn font_dirs() -> Vec<PathBuf> {
    let mut dirs: Vec<PathBuf> = [
        "/usr/share/fonts",
        "/usr/local/share/fonts",
        "/System/Library/Fonts",
        "/Library/Fonts",
        "C:\\Windows\\Fonts",
    ]
    .iter()
    .map(PathBuf::from)
    .collect();
    if let Some(home) = std::env::var_os("HOME") {
        let home = PathBuf::from(home);
        dirs.push(home.join(".local/share/fonts"));
        dirs.push(home.join(".fonts"));
        dirs.push(home.join("Library/Fonts"));
    }
    dirs
}

/// Every file under `dir` (four levels deep), name → path.
fn files_under(dir: &Path, depth: usize, out: &mut Vec<(String, PathBuf)>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    let mut entries: Vec<_> = entries.flatten().collect();
    entries.sort_by_key(|e| e.file_name());
    for entry in entries {
        let path = entry.path();
        if path.is_dir() {
            if depth > 0 {
                files_under(&path, depth - 1, out);
            }
        } else if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            out.push((name.to_string(), path.clone()));
        }
    }
}

fn load(path: &Path, index: u32) -> Result<Font, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Font::from_bytes(
        bytes,
        FontSettings {
            collection_index: index,
            ..FontSettings::default()
        },
    )
    .map_err(|e| format!("{}: {e}", path.display()))
}

impl Fonts {
    /// `--font` (one file for both weights), else the first system pair.
    pub fn find(explicit: Option<&Path>) -> Result<Fonts, String> {
        if let Some(path) = explicit {
            return Ok(Fonts {
                regular: load(path, 0)?,
                bold: load(path, 0)?,
                source: path.to_path_buf(),
            });
        }
        let mut files = Vec::new();
        for dir in font_dirs() {
            files_under(&dir, 4, &mut files);
        }
        let find = |name: &str| {
            files
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, p)| p.clone())
        };
        for (regular, bold) in CANDIDATES {
            let (Some(r), Some(b)) = (find(regular), find(bold)) else {
                continue;
            };
            // A collection (Helvetica.ttc) holds the bold face second.
            let bold_index = u32::from(r == b);
            if let (Ok(regular), Ok(bold)) = (load(&r, 0), load(&b, bold_index)) {
                return Ok(Fonts {
                    regular,
                    bold,
                    source: r,
                });
            }
        }
        Err("no caption font found: pass --font <file.ttf>".into())
    }
}

/// Relative luminance of `#rrggbb` (1 when unreadable, as the page).
fn luminance(background: &str) -> f32 {
    let Ok(n) = u32::from_str_radix(background.trim_start_matches('#'), 16) else {
        return 1.0;
    };
    (0.2126 * ((n >> 16) & 255) as f32
        + 0.7152 * ((n >> 8) & 255) as f32
        + 0.0722 * (n & 255) as f32)
        / 255.0
}

/// Draw `text` with its em box's top at `top`, from `left`, at `size`
/// px (shrunk to `max_width`), in `color` at `alpha`.
#[allow(clippy::too_many_arguments)]
fn draw_text(
    frame: &mut [u8],
    width: usize,
    height: usize,
    font: &Font,
    text: &str,
    left: f32,
    top: f32,
    size: f32,
    max_width: f32,
    color: [u8; 3],
    alpha: f32,
) {
    let advance = |size: f32| -> f32 {
        text.chars()
            .map(|c| font.metrics(c, size).advance_width)
            .sum()
    };
    let full = advance(size);
    let size = if full > max_width && full > 0.0 {
        size * max_width / full
    } else {
        size
    };
    let ascent = font
        .horizontal_line_metrics(size)
        .map_or(size * 0.8, |m| m.ascent);
    // Canvas `textBaseline = "top"` puts the em box's top at `top`; the
    // em box sits about a tenth of the size above the ascender line.
    let baseline = top + ascent;
    let mut pen = left;
    for c in text.chars() {
        let (metrics, coverage) = font.rasterize(c, size);
        let x0 = (pen + metrics.xmin as f32).round() as isize;
        let y0 = (baseline - metrics.height as f32 - metrics.ymin as f32).round() as isize;
        for gy in 0..metrics.height {
            for gx in 0..metrics.width {
                let (x, y) = (x0 + gx as isize, y0 + gy as isize);
                if x < 0 || y < 0 || x >= width as isize || y >= height as isize {
                    continue;
                }
                let a = alpha * coverage[gy * metrics.width + gx] as f32 / 255.0;
                if a <= 0.0 {
                    continue;
                }
                let o = (y as usize * width + x as usize) * 4;
                for k in 0..3 {
                    let dst = frame[o + k] as f32;
                    frame[o + k] = (dst + (color[k] as f32 - dst) * a).round() as u8;
                }
            }
        }
        pen += metrics.advance_width;
    }
}

/// Draw `caption` over `frame` (RGBA, `width × height`).
pub fn draw_caption(
    frame: &mut [u8],
    width: usize,
    height: usize,
    caption: &Caption,
    background: &str,
    fonts: &Fonts,
) {
    if caption.opacity <= 0.0 {
        return;
    }
    let unit = width.min(height) as f32 / 1080.0;
    let color = if luminance(background) > 0.5 {
        [0x11, 0x13, 0x18]
    } else {
        [0xf5, 0xf7, 0xfa]
    };
    let pad = 56.0 * unit;
    let title_size = if caption.kind == "title" { 64.0 } else { 48.0 } * unit;
    let max_width = width as f32 - 2.0 * pad;
    let opacity = caption.opacity as f32;
    draw_text(
        frame,
        width,
        height,
        &fonts.bold,
        &caption.title,
        pad,
        pad,
        title_size,
        max_width,
        color,
        opacity,
    );
    if !caption.subtitle.is_empty() {
        draw_text(
            frame,
            width,
            height,
            &fonts.regular,
            &caption.subtitle,
            pad,
            pad + title_size * 1.25,
            30.0 * unit,
            max_width,
            color,
            opacity * 0.8,
        );
    }
}
