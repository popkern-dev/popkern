//! Surface-quality analysis rendering (PK-520): zebra stripes and
//! curvature false-color over the CPU rasterizer.
//!
//! Zebra stripes are reflection lines of a family of parallel light
//! strips: the classic class-A continuity inspection. A G1 break
//! kinks the stripes across a seam; a G2 break kinks their
//! DERIVATIVE (the stripe family bends abruptly while staying
//! connected). Stripes are computed PER PIXEL — per-vertex banding
//! would alias the exact discontinuities this view exists to show —
//! from either the interpolated mesh normal or the surface's own
//! normal (PK-503, [`crate::zebra`]), so a stripe defect can be
//! attributed to the geometry or to the tessellation.
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
use crate::zebra::{ZebraModeDisagreement, ZebraNormals, rasterize_zebra_frame};

/// Which analysis view to render.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurfaceAnalysisMode {
    /// Reflection-line zebra stripes; `stripe_count` families across
    /// the reflection hemisphere (12 is a reasonable default).
    Zebra {
        /// Number of stripe periods across the reflection range.
        stripe_count: u32,
        /// Which normal the stripes reflect in (see [`crate::zebra`]).
        normals: ZebraNormals,
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
        SurfaceAnalysisMode::Zebra {
            stripe_count,
            normals,
        } => {
            Ok(rasterize_zebra_frame(brep, mesh, camera, width, height)
                .shade(normals, stripe_count))
        }
        SurfaceAnalysisMode::Curvature { kind } => {
            render_curvature(brep, mesh, camera, width, height, kind)
        }
    }
}

/// Tessellate (display quality — see
/// [`popkern_mesh::TessellationOptions::display`]) and render an
/// analysis view straight to a PNG.
pub fn render_surface_analysis_to_png(
    brep: &BRep,
    camera: &Camera,
    width: u32,
    height: u32,
    mode: SurfaceAnalysisMode,
    output_path: &Path,
) -> Result<(), VizError> {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::display())
        .map_err(|e| VizError::SoftwareRender(format!("tessellation failed: {e}")))?;
    let pixels = render_surface_analysis(brep, &mesh, camera, width, height, mode)?;
    crate::headless::write_png(output_path, width, height, &pixels)
}

/// Tessellate with `tessellation` (viewers use
/// [`popkern_mesh::TessellationOptions::display`]), render a zebra
/// view in `normals` to a PNG, and return the mesh-vs-surface normal
/// disagreement of the same frame — the number that says how much of
/// the image's stripe structure the mesh owns.
#[allow(clippy::too_many_arguments)] // one render request: view, frame, stripes, source, mesh, file
pub fn render_zebra_to_png(
    brep: &BRep,
    camera: &Camera,
    width: u32,
    height: u32,
    stripe_count: u32,
    normals: ZebraNormals,
    tessellation: &popkern_mesh::TessellationOptions,
    output_path: &Path,
) -> Result<ZebraModeDisagreement, VizError> {
    let mesh = popkern_mesh::tessellate(brep, tessellation)
        .map_err(|e| VizError::SoftwareRender(format!("tessellation failed: {e}")))?;
    let frame = rasterize_zebra_frame(brep, &mesh, camera, width, height);
    crate::headless::write_png(
        output_path,
        width,
        height,
        &frame.shade(normals, stripe_count),
    )?;
    Ok(frame.mode_disagreement(stripe_count))
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

#[cfg(test)]
#[path = "surface_analysis_tests.rs"]
mod tests;
