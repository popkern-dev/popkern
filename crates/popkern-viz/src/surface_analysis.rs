//! Surface-quality analysis rendering (PK-520): zebra stripes and
//! curvature false-color over the CPU rasterizer.
//!
//! Zebra stripes are reflection lines of a family of parallel light
//! strips: the classic class-A continuity inspection. A G1 break
//! kinks the stripes across a seam; a G2 break kinks their
//! DERIVATIVE (the stripe family bends abruptly while staying
//! connected). Stripes are computed PER PIXEL from the interpolated
//! surface normal — per-vertex banding would alias the exact
//! discontinuities this view exists to show.
//!
//! Curvature false-color maps [`popkern_measure::curvature_field`]
//! (analytic per-vertex Gaussian/mean curvature) through a diverging
//! palette, unlit — lighting would distort the value reading. The
//! range is robust-symmetric (98th percentile of |value|) so a
//! single degenerate spike cannot wash out the map.

use std::path::Path;

use popkern_core::topology::BRep;
use popkern_measure::{CurvatureKind, curvature_field};
use popkern_mesh::TriMesh;

use crate::camera::Camera;
use crate::cpu_render::rasterize::{Framebuffer, ScreenVertex, rasterize_triangle};
use crate::error::VizError;

/// Which analysis view to render.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurfaceAnalysisMode {
    /// Reflection-line zebra stripes; `stripe_count` families across
    /// the reflection hemisphere (12 is a reasonable default).
    Zebra {
        /// Number of stripe periods across the reflection range.
        stripe_count: u32,
    },
    /// False-color curvature map.
    Curvature {
        /// Which curvature scalar to map.
        kind: CurvatureKind,
    },
}

/// Render an analysis view of `brep` to an RGBA buffer
/// (`width * height * 4`).
pub fn render_surface_analysis(
    brep: &BRep,
    mesh: &TriMesh,
    camera: &Camera,
    width: u32,
    height: u32,
    mode: SurfaceAnalysisMode,
) -> Result<Vec<u8>, VizError> {
    match mode {
        SurfaceAnalysisMode::Zebra { stripe_count } => {
            Ok(render_zebra(mesh, camera, width, height, stripe_count))
        }
        SurfaceAnalysisMode::Curvature { kind } => {
            render_curvature(brep, mesh, camera, width, height, kind)
        }
    }
}

/// Tessellate and render an analysis view straight to a PNG.
pub fn render_surface_analysis_to_png(
    brep: &BRep,
    camera: &Camera,
    width: u32,
    height: u32,
    mode: SurfaceAnalysisMode,
    output_path: &Path,
) -> Result<(), VizError> {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .map_err(|e| VizError::SoftwareRender(format!("tessellation failed: {e}")))?;
    let pixels = render_surface_analysis(brep, &mesh, camera, width, height, mode)?;
    crate::headless::write_png(output_path, width, height, &pixels)
}

/// Map a curvature value to a diverging blue–white–red color over a
/// symmetric range: negative blue, zero white, positive red.
pub fn diverging_color(value: f64, half_range: f64) -> [u8; 4] {
    let t = if half_range > 0.0 {
        (value / half_range).clamp(-1.0, 1.0)
    } else {
        0.0
    };
    let (cold, warm) = ([59.0, 76.0, 192.0], [180.0, 4.0, 38.0]);
    let white = [242.0, 242.0, 242.0];
    let lerp = |a: [f64; 3], b: [f64; 3], f: f64| -> [u8; 4] {
        [
            (a[0] + (b[0] - a[0]) * f) as u8,
            (a[1] + (b[1] - a[1]) * f) as u8,
            (a[2] + (b[2] - a[2]) * f) as u8,
            255,
        ]
    };
    if t < 0.0 {
        lerp(white, cold, -t)
    } else {
        lerp(white, warm, t)
    }
}

/// Robust symmetric half-range: 1.25× the median of the NONZERO
/// magnitudes. Flat regions (exactly zero curvature) dominate most
/// parts and would drag a plain percentile to zero; degenerate-chart
/// spikes would inflate a max-based range until the working band
/// reads as faint tint — the nonzero median saturates the curved
/// features (a constant-radius fillet paints at full strength) while
/// staying immune to both populations.
pub fn robust_half_range(values: &[f64]) -> f64 {
    // 1e-12 (1/model-length) is a DISPLAY floor separating "flat"
    // from "curved" — orders of magnitude below any inspectable
    // curvature at sane model scale, not a geometric tolerance.
    let mut mags: Vec<f64> = values
        .iter()
        .map(|v| v.abs())
        .filter(|m| m.is_finite() && *m > 1e-12)
        .collect();
    if mags.is_empty() {
        return 0.0;
    }
    mags.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    (1.25 * mags[mags.len() / 2]).max(f64::MIN_POSITIVE)
}

fn render_curvature(
    brep: &BRep,
    mesh: &TriMesh,
    camera: &Camera,
    width: u32,
    height: u32,
    kind: CurvatureKind,
) -> Result<Vec<u8>, VizError> {
    let field = curvature_field(brep, mesh, kind)
        .map_err(|e| VizError::SoftwareRender(format!("curvature field: {e}")))?;
    let half_range = robust_half_range(&field.per_vertex);
    let vp = camera.view_projection();
    let (w, h) = (width as f32, height as f32);
    let mut fb = Framebuffer::new(width, height);
    fb.clear([38, 38, 38, 255]);
    let near_w = camera.near.max(1.0e-4);
    for tri in &mesh.triangles {
        let mut screen: [Option<ScreenVertex>; 3] = [None, None, None];
        let mut behind = false;
        for (slot, &vi) in screen.iter_mut().zip(tri.iter()) {
            let mv = &mesh.vertices[vi];
            let clip = vp
                * glam::Vec4::new(
                    mv.position.x as f32,
                    mv.position.y as f32,
                    mv.position.z as f32,
                    1.0,
                );
            if clip.w <= near_w {
                behind = true;
                break;
            }
            let ndc = clip / clip.w;
            *slot = Some(ScreenVertex {
                x: (ndc.x * 0.5 + 0.5) * w,
                y: (1.0 - (ndc.y * 0.5 + 0.5)) * h,
                z: ndc.z,
                color: diverging_color(field.per_vertex[vi], half_range),
            });
        }
        if behind {
            continue;
        }
        if let (Some(v0), Some(v1), Some(v2)) =
            (screen[0].as_ref(), screen[1].as_ref(), screen[2].as_ref())
        {
            rasterize_triangle(&mut fb, v0, v1, v2);
        }
    }
    Ok(fb.to_rgba())
}

fn render_zebra(
    mesh: &TriMesh,
    camera: &Camera,
    width: u32,
    height: u32,
    stripe_count: u32,
) -> Vec<u8> {
    let vp = camera.view_projection();
    let (w, h) = (width as f32, height as f32);
    let eye = camera.eye;
    // Stripe family axis: the camera's right vector — stripes read as
    // vertical bands of a horizon of parallel lights.
    let forward = (camera.target - camera.eye).normalize_or_zero();
    let right = forward.cross(camera.up).normalize_or_zero();
    let stripes = stripe_count.max(1) as f32;
    let near_w = camera.near.max(1.0e-4);

    let mut fb = Framebuffer::new(width, height);
    fb.clear([38, 38, 38, 255]);
    for tri in &mesh.triangles {
        let mut sv = [[0.0f32; 3]; 3]; // screen x, y, depth
        let mut world = [glam::Vec3::ZERO; 3];
        let mut normals = [glam::Vec3::ZERO; 3];
        let mut behind = false;
        for (k, &vi) in tri.iter().enumerate() {
            let mv = &mesh.vertices[vi];
            let p = glam::Vec3::new(
                mv.position.x as f32,
                mv.position.y as f32,
                mv.position.z as f32,
            );
            let clip = vp * glam::Vec4::new(p.x, p.y, p.z, 1.0);
            if clip.w <= near_w {
                behind = true;
                break;
            }
            let ndc = clip / clip.w;
            sv[k] = [
                (ndc.x * 0.5 + 0.5) * w,
                (1.0 - (ndc.y * 0.5 + 0.5)) * h,
                ndc.z,
            ];
            world[k] = p;
            normals[k] =
                glam::Vec3::new(mv.normal.x as f32, mv.normal.y as f32, mv.normal.z as f32);
        }
        if behind {
            continue;
        }
        rasterize_zebra_triangle(&mut fb, &sv, &world, &normals, eye, right, stripes);
    }
    fb.to_rgba()
}

/// Rasterize one triangle with per-pixel normal interpolation and
/// reflection-line banding.
fn rasterize_zebra_triangle(
    fb: &mut Framebuffer,
    sv: &[[f32; 3]; 3],
    world: &[glam::Vec3; 3],
    normals: &[glam::Vec3; 3],
    eye: glam::Vec3,
    stripe_axis: glam::Vec3,
    stripes: f32,
) {
    let (w, h) = (fb.width as i32, fb.height as i32);
    let min_x = sv
        .iter()
        .map(|v| v[0])
        .fold(f32::INFINITY, f32::min)
        .floor() as i32;
    let max_x = sv
        .iter()
        .map(|v| v[0])
        .fold(f32::NEG_INFINITY, f32::max)
        .ceil() as i32;
    let min_y = sv
        .iter()
        .map(|v| v[1])
        .fold(f32::INFINITY, f32::min)
        .floor() as i32;
    let max_y = sv
        .iter()
        .map(|v| v[1])
        .fold(f32::NEG_INFINITY, f32::max)
        .ceil() as i32;
    let (min_x, max_x) = (min_x.max(0), max_x.min(w - 1));
    let (min_y, max_y) = (min_y.max(0), max_y.min(h - 1));
    if min_x > max_x || min_y > max_y {
        return;
    }
    let area = (sv[1][0] - sv[0][0]) * (sv[2][1] - sv[0][1])
        - (sv[2][0] - sv[0][0]) * (sv[1][1] - sv[0][1]);
    if area.abs() < 1e-12 {
        return;
    }
    for py in min_y..=max_y {
        for px in min_x..=max_x {
            let (fx, fy) = (px as f32 + 0.5, py as f32 + 0.5);
            let b0 = ((sv[1][0] - fx) * (sv[2][1] - fy) - (sv[2][0] - fx) * (sv[1][1] - fy)) / area;
            let b1 = ((sv[2][0] - fx) * (sv[0][1] - fy) - (sv[0][0] - fx) * (sv[2][1] - fy)) / area;
            let b2 = 1.0 - b0 - b1;
            if b0 < 0.0 || b1 < 0.0 || b2 < 0.0 {
                continue;
            }
            let depth = b0 * sv[0][2] + b1 * sv[1][2] + b2 * sv[2][2];
            let idx = (py as usize) * (fb.width as usize) + (px as usize);
            if depth >= fb.depth[idx] {
                continue;
            }
            let n = (normals[0] * b0 + normals[1] * b1 + normals[2] * b2).normalize_or_zero();
            let p = world[0] * b0 + world[1] * b1 + world[2] * b2;
            let view = (p - eye).normalize_or_zero();
            let r = view - n * (2.0 * view.dot(n));
            // Band on the reflected direction's component along the
            // stripe axis: s ∈ [-1, 1] → stripe phase.
            let s = r.dot(stripe_axis);
            let phase = (s * 0.5 + 0.5) * stripes;
            let band = phase.fract() < 0.5;
            let shade = if band { 24u8 } else { 236u8 };
            fb.depth[idx] = depth;
            fb.color[idx] = [shade, shade, shade, 255];
        }
    }
}

#[cfg(test)]
#[path = "surface_analysis_tests.rs"]
mod tests;
