//! A PNG of one per-vertex field on the analysis mesh, through the
//! pinned `popkern-viz` software rasteriser (no GPU, so the bytes are
//! the same on every run): a blue–green–red ramp over the field's range
//! and a legend bar with the minimum and maximum.

use std::path::Path;

use glam::Vec3;
use popkern_analysis::model::AnalysisModel;
use popkern_viz::camera::Camera;
use popkern_viz::cpu_render::render_prepared;
use popkern_viz::headless::write_png;
use popkern_viz::mesh_data::MeshVertex;

/// Image size, pixels.
pub const WIDTH: u32 = 800;
/// Image height, pixels.
pub const HEIGHT: u32 = 600;

/// Ramp stops: blue, cyan, green, yellow, red.
const RAMP: [[f32; 3]; 5] = [
    [0.10, 0.20, 0.85],
    [0.10, 0.75, 0.90],
    [0.20, 0.80, 0.25],
    [0.95, 0.85, 0.15],
    [0.90, 0.15, 0.10],
];

/// The ramp colour of `t` in [0, 1].
fn ramp(t: f64) -> [f32; 3] {
    let t = if t.is_finite() {
        t.clamp(0.0, 1.0)
    } else {
        0.0
    };
    let x = t as f32 * (RAMP.len() - 1) as f32;
    let i = (x.floor() as usize).min(RAMP.len() - 2);
    let f = x - i as f32;
    std::array::from_fn(|k| RAMP[i][k] * (1.0 - f) + RAMP[i + 1][k] * f)
}

/// 3 × 5 glyphs (rows top to bottom, bits left to right) for the
/// legend's numbers.
fn glyph(c: char) -> Option<[u8; 5]> {
    Some(match c {
        '0' => [7, 5, 5, 5, 7],
        '1' => [2, 6, 2, 2, 7],
        '2' => [7, 1, 7, 4, 7],
        '3' => [7, 1, 7, 1, 7],
        '4' => [5, 5, 7, 1, 1],
        '5' => [7, 4, 7, 1, 7],
        '6' => [7, 4, 7, 5, 7],
        '7' => [7, 1, 1, 1, 1],
        '8' => [7, 5, 7, 5, 7],
        '9' => [7, 5, 7, 1, 7],
        '.' => [0, 0, 0, 0, 2],
        '-' => [0, 0, 7, 0, 0],
        '+' => [0, 2, 7, 2, 0],
        'e' => [0, 7, 7, 4, 7],
        _ => return None,
    })
}

/// Draws `text` with its top-left corner at (x, y), each font pixel
/// `scale` image pixels.
fn draw_text(rgba: &mut [u8], x: u32, y: u32, text: &str, scale: u32) {
    let mut cx = x;
    for c in text.chars() {
        if let Some(rows) = glyph(c) {
            for (r, bits) in rows.iter().enumerate() {
                for b in 0..3u32 {
                    if bits & (4 >> b) != 0 {
                        for dy in 0..scale {
                            for dx in 0..scale {
                                put(
                                    rgba,
                                    cx + b * scale + dx,
                                    y + r as u32 * scale + dy,
                                    [235; 3],
                                );
                            }
                        }
                    }
                }
            }
        }
        cx += 4 * scale;
    }
}

fn put(rgba: &mut [u8], x: u32, y: u32, c: [u8; 3]) {
    if x >= WIDTH || y >= HEIGHT {
        return;
    }
    let i = 4 * (y as usize * WIDTH as usize + x as usize);
    if let Some(p) = rgba.get_mut(i..i + 4) {
        p.copy_from_slice(&[c[0], c[1], c[2], 255]);
    }
}

/// Unit outward normal at each vertex of every solid (area-weighted).
fn vertex_normals(model: &AnalysisModel) -> Vec<Vec<[f32; 3]>> {
    model
        .solids
        .iter()
        .map(|s| {
            let mut n = vec![[0.0f64; 3]; s.vertices.len()];
            for t in &s.triangles {
                let [a, b, c] = t.map(|i| s.vertices[i as usize]);
                let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
                let x = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                for &i in t {
                    for k in 0..3 {
                        n[i as usize][k] += x[k];
                    }
                }
            }
            n.iter()
                .map(|x| {
                    let l = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
                    if l > 0.0 {
                        x.map(|c| (c / l) as f32)
                    } else {
                        [0.0, 0.0, 1.0]
                    }
                })
                .collect()
        })
        .collect()
}

/// Renders `field` (per solid, per vertex) on `model` and writes the
/// PNG to `path`; returns the field's (min, max).
pub fn write_field_png(
    model: &AnalysisModel,
    field: &[Vec<f64>],
    path: &Path,
) -> Result<(f64, f64), String> {
    let (lo, hi) = field
        .iter()
        .flatten()
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &v| {
            (l.min(v), h.max(v))
        });
    if !(lo.is_finite() && hi.is_finite()) {
        return Err("the field has no finite values".into());
    }
    let span = if hi > lo { hi - lo } else { 1.0 };
    let normals = vertex_normals(model);
    let mut vertices = Vec::new();
    let mut indices = Vec::new();
    let (mut bmin, mut bmax) = (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY));
    for ((s, values), n) in model.solids.iter().zip(field).zip(&normals) {
        let base = vertices.len() as u32;
        for ((p, v), nv) in s.vertices.iter().zip(values).zip(n) {
            let q = Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32);
            bmin = bmin.min(q);
            bmax = bmax.max(q);
            let c = ramp((v - lo) / span);
            vertices.push(MeshVertex {
                position: q.to_array(),
                normal: *nv,
                color: [c[0], c[1], c[2], 1.0],
            });
        }
        // The rasteriser's screen y points down, which mirrors the
        // winding: outward counter-clockwise triangles are given reversed.
        indices.extend(
            s.triangles
                .iter()
                .flat_map(|t| [t[0], t[2], t[1]])
                .map(|i| base + i),
        );
    }
    let mut camera = Camera::default_iso();
    camera.aspect = WIDTH as f32 / HEIGHT as f32;
    camera.fit_to_bounds(bmin, bmax);
    let mut rgba = render_prepared(&vertices, &indices, &[], &camera, WIDTH, HEIGHT, false);
    // Legend: the ramp from the minimum (bottom) to the maximum (top).
    let (x0, top, bottom) = (WIDTH - 70, 60u32, HEIGHT - 60);
    for y in top..=bottom {
        let c = ramp(f64::from(bottom - y) / f64::from(bottom - top));
        for x in x0..x0 + 20 {
            put(&mut rgba, x, y, c.map(|v| (v * 255.0).round() as u8));
        }
    }
    draw_text(&mut rgba, x0 - 60, top - 24, &format!("{hi:.3e}"), 2);
    draw_text(&mut rgba, x0 - 60, bottom + 10, &format!("{lo:.3e}"), 2);
    write_png(path, WIDTH, HEIGHT, &rgba).map_err(|e| e.to_string())?;
    Ok((lo, hi))
}
