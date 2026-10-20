//! Gaussian and mean curvature colour maps with a stated range
//! (PK-504). Each pixel is evaluated on its own surface at its own
//! `(u, v)` (the zebra frame's per-pixel surface location), not
//! interpolated from vertices, and drawn unlit: lighting would change
//! what a colour means. Mean curvature is signed against the OUTWARD
//! normal (the pixel's facet normal), Gaussian curvature needs no sign.
//!
//! The range is part of the reading, so it is chosen explicitly and
//! reported: [`CurvatureRange::Fitted`] fits the sampled extreme |value|
//! (never below a floor), [`CurvatureRange::Fixed`] holds one range
//! across models and counts the pixels it clips. The extrema reported
//! are SAMPLES at pixel centres, not certified bounds of the surface. A
//! legend (the colour bar and its ± range, printed) can be drawn on the
//! image itself.

use popkern_core::topology::BRep;
use popkern_measure::{CurvatureKind, curvature_at};
use popkern_mesh::TriMesh;
use serde::Serialize;

use crate::camera::Camera;
use crate::surface_analysis::diverging_color;
use crate::zebra::{SurfaceHit, rasterize_zebra_frame};

const BACKGROUND: [u8; 4] = [38, 38, 38, 255];

/// How a curvature map's colour range is chosen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CurvatureRange {
    /// The sampled extreme |value| over the rendered pixels, never
    /// below `floor` (a flat part maps to white instead of noise).
    Fitted {
        /// The smallest half-range, in 1/length (curvature) units.
        floor: f64,
    },
    /// One half-range for every model compared under it; values
    /// beyond it clip and are counted.
    Fixed(f64),
}

impl Default for CurvatureRange {
    /// Fitted, floor `1e-12` per length unit: far under any inspectable
    /// curvature, only separating flat from curved.
    fn default() -> Self {
        CurvatureRange::Fitted { floor: 1e-12 }
    }
}

/// What a curvature map shows, and how.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvatureMapOptions {
    /// The scalar mapped.
    pub kind: CurvatureKind,
    /// The colour range policy.
    pub range: CurvatureRange,
    /// Draw the colour bar and its range on the image.
    pub legend: bool,
}

/// What a rendered curvature map used and found.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CurvatureMapReport {
    /// `"gaussian"`, `"mean"` or `"max_principal_abs"`.
    pub kind: &'static str,
    /// The half-range mapped to full colour (blue -h .. white 0 .. red +h).
    pub half_range: f64,
    /// Whether the range was fixed by the caller (else fitted).
    pub fixed: bool,
    /// Pixels the surface covers.
    pub pixels: usize,
    /// Of those, pixels whose value lies beyond the half-range (fixed
    /// ranges only; a fitted range clips none).
    pub clipped: usize,
    /// Of those, pixels no surface value was evaluated for (no surface
    /// hit, a degenerate chart): drawn grey.
    pub unevaluated: usize,
    /// The smallest and largest SAMPLED value (pixel centres; not a
    /// certified bound of the surface).
    pub sampled_min: f64,
    /// See `sampled_min`.
    pub sampled_max: f64,
}

/// Render a curvature map of `brep` (meshed as `mesh`) seen through
/// `camera`: RGBA pixels and the report of the range used.
pub fn render_curvature_map(
    brep: &BRep,
    mesh: &TriMesh,
    camera: &Camera,
    width: u32,
    height: u32,
    options: CurvatureMapOptions,
) -> (Vec<u8>, CurvatureMapReport) {
    let frame = rasterize_zebra_frame(brep, mesh, camera, width, height);
    let mut values: Vec<Option<f64>> = Vec::with_capacity(frame.samples.len());
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut pixels, mut unevaluated) = (0usize, 0usize);
    for sample in &frame.samples {
        let Some(s) = sample else {
            values.push(None);
            continue;
        };
        pixels += 1;
        let surface = brep
            .faces
            .get(s.surface_face)
            .and_then(|f| f.surface.as_ref());
        let value = match (s.hit, surface) {
            (SurfaceHit::NoSurface, _) | (_, None) => None,
            (_, Some(surface)) => {
                curvature_at(&surface.geometry, s.uv.0, s.uv.1, options.kind).map(|v| {
                    // `curvature_at` signs H against the CHART normal
                    // du x dv; flip it where that points inward of the
                    // pixel's outward facet normal.
                    use popkern_core::geometry::surface::SurfaceEval;
                    let (su, sv) = (
                        surface.geometry.partial_u(s.uv.0, s.uv.1),
                        surface.geometry.partial_v(s.uv.0, s.uv.1),
                    );
                    let chart = su.cross(&sv);
                    let inward =
                        crate::zebra::dot([chart.x, chart.y, chart.z], s.mesh_normal) < 0.0;
                    if options.kind == CurvatureKind::Mean && inward {
                        -v
                    } else {
                        v
                    }
                })
            }
        };
        match value.filter(|v| v.is_finite()) {
            Some(v) => {
                lo = lo.min(v);
                hi = hi.max(v);
                values.push(Some(v));
            }
            None => {
                unevaluated += 1;
                values.push(Some(f64::NAN));
            }
        }
    }
    if !lo.is_finite() {
        (lo, hi) = (0.0, 0.0);
    }
    let (half_range, fixed) = match options.range {
        CurvatureRange::Fitted { floor } => (lo.abs().max(hi.abs()).max(floor.max(0.0)), false),
        CurvatureRange::Fixed(h) => (h.abs(), true),
    };
    let mut clipped = 0usize;
    let mut rgba = Vec::with_capacity(values.len() * 4);
    for value in &values {
        let px = match value {
            None => BACKGROUND,
            Some(v) if v.is_nan() => [128, 128, 128, 255],
            Some(v) => {
                if v.abs() > half_range {
                    clipped += 1;
                }
                diverging_color(*v, half_range)
            }
        };
        rgba.extend_from_slice(&px);
    }
    if options.legend {
        draw_legend(&mut rgba, width, height, half_range, options.kind);
    }
    let report = CurvatureMapReport {
        kind: match options.kind {
            CurvatureKind::Gaussian => "gaussian",
            CurvatureKind::Mean => "mean",
            CurvatureKind::MaxPrincipalAbs => "max_principal_abs",
        },
        half_range,
        fixed,
        pixels,
        clipped,
        unevaluated,
        sampled_min: lo,
        sampled_max: hi,
    };
    (rgba, report)
}

/// The colour bar down the right edge, `+h` printed at its top, `0` at
/// its middle and `-h` at its bottom, the scalar's letter above.
fn draw_legend(rgba: &mut [u8], width: u32, height: u32, half: f64, kind: CurvatureKind) {
    let (w, h) = (width as i64, height as i64);
    if w < 60 || h < 60 {
        return;
    }
    let scale = (h / 160).max(1);
    let (bar_w, top, bottom) = (8 * scale, 12 * scale, h - 12 * scale);
    let x0 = w - bar_w - 4 * scale;
    let mut put = |x: i64, y: i64, c: [u8; 4]| {
        if (0..w).contains(&x) && (0..h).contains(&y) {
            let i = ((y * w + x) * 4) as usize;
            rgba[i..i + 4].copy_from_slice(&c);
        }
    };
    for y in top..bottom {
        let t = 1.0 - 2.0 * (y - top) as f64 / (bottom - top - 1).max(1) as f64;
        let c = diverging_color(t * half.max(f64::MIN_POSITIVE), half.max(f64::MIN_POSITIVE));
        for x in x0..x0 + bar_w {
            put(x, y, c);
        }
    }
    let letter = match kind {
        CurvatureKind::Gaussian => "K",
        CurvatureKind::Mean => "H",
        CurvatureKind::MaxPrincipalAbs => "k",
    };
    let white = [236, 236, 236, 255];
    let labels = [
        (letter.to_string(), top - 8 * scale),
        (format!("+{half:.2e}"), top),
        ("0".to_string(), (top + bottom) / 2 - 3 * scale),
        (format!("-{half:.2e}"), bottom - 6 * scale),
    ];
    for (text, y) in labels {
        let width_px = text.chars().count() as i64 * 4 * scale;
        let x = x0 - width_px - 3 * scale;
        // A dark backing: the label reads over the model too.
        for by in y - scale..y + 6 * scale {
            for bx in x - scale..x + width_px {
                put(bx, by, BACKGROUND);
            }
        }
        for (k, ch) in text.chars().enumerate() {
            let glyph = glyph(ch);
            for (row, bits) in glyph.iter().enumerate() {
                for col in 0..3 {
                    if bits & (0b100 >> col) != 0 {
                        for dy in 0..scale {
                            for dx in 0..scale {
                                put(
                                    x + (k as i64 * 4 + col) * scale + dx,
                                    y + row as i64 * scale + dy,
                                    white,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

/// A 3x5 bitmap glyph (rows top to bottom, bit 2 the left column).
fn glyph(c: char) -> [u8; 5] {
    match c {
        '0' => [0b111, 0b101, 0b101, 0b101, 0b111],
        '1' => [0b010, 0b110, 0b010, 0b010, 0b111],
        '2' => [0b111, 0b001, 0b111, 0b100, 0b111],
        '3' => [0b111, 0b001, 0b111, 0b001, 0b111],
        '4' => [0b101, 0b101, 0b111, 0b001, 0b001],
        '5' => [0b111, 0b100, 0b111, 0b001, 0b111],
        '6' => [0b111, 0b100, 0b111, 0b101, 0b111],
        '7' => [0b111, 0b001, 0b010, 0b010, 0b010],
        '8' => [0b111, 0b101, 0b111, 0b101, 0b111],
        '9' => [0b111, 0b101, 0b111, 0b001, 0b111],
        '.' => [0b000, 0b000, 0b000, 0b000, 0b010],
        '-' => [0b000, 0b000, 0b111, 0b000, 0b000],
        '+' => [0b000, 0b010, 0b111, 0b010, 0b000],
        'e' => [0b000, 0b111, 0b111, 0b100, 0b111],
        'K' => [0b101, 0b110, 0b100, 0b110, 0b101],
        'H' => [0b101, 0b101, 0b111, 0b101, 0b101],
        'k' => [0b100, 0b101, 0b110, 0b110, 0b101],
        _ => [0; 5],
    }
}

#[cfg(test)]
#[path = "curvature_map_tests.rs"]
mod tests;
