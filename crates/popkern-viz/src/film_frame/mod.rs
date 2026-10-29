//! Presentation frames for films and still renders
//! (docs/architecture/cli-render-export.md §4): a CPU rasterizer that
//! draws a scene the way popkern live's film look does — per-pixel
//! shading in each part's colour (`studio`/`cad`, or the `toon`
//! three-tone ramp), black face-boundary edges, fades, the hover-blue
//! highlight, a section plane with closed-looking cuts, and a soft
//! contact shadow on a ground plane under the model — on a plain
//! background, 2× supersampled.
//!
//! Separate from [`crate::cpu_render`], whose vertex shading the zebra
//! and debug renders keep. Deterministic: the frame is split into row
//! bands rendered on threads, but each pixel's computation is fixed, so
//! the same scene gives the same bytes.

mod raster;
mod shade;
mod shadow;

use glam::{DMat4, DVec3, Mat4, Vec3, Vec4};
use rayon::prelude::*;

pub use shade::Look;

/// One drawn body: its display mesh in its own frame, placed.
#[derive(Debug, Clone)]
pub struct FrameItem<'a> {
    /// `x, y, z` per vertex.
    pub positions: &'a [f32],
    /// Unit normal per vertex.
    pub normals: &'a [f32],
    /// Three vertex indices per triangle.
    pub indices: &'a [u32],
    /// Edge segments (two points each), body frame.
    pub edges: &'a [[f32; 3]],
    /// Body frame to world (column-major, glam).
    pub world: DMat4,
    /// sRGB colour, 0..1.
    pub color: [f32; 3],
    /// 0..1 (below 1 the item blends over what is behind it).
    pub opacity: f32,
    /// 0..1 hover-blue emphasis.
    pub highlight: f32,
}

/// A perspective camera.
#[derive(Debug, Clone, Copy)]
pub struct FrameCamera {
    /// Eye position.
    pub eye: DVec3,
    /// The point looked at.
    pub target: DVec3,
    /// Up direction.
    pub up: DVec3,
    /// Vertical field of view, degrees.
    pub vfov_deg: f64,
    /// Near plane distance.
    pub near: f64,
    /// Far plane distance.
    pub far: f64,
}

/// A whole frame's scene.
#[derive(Debug, Clone)]
pub struct FrameScene<'a> {
    /// The bodies.
    pub items: Vec<FrameItem<'a>>,
    /// The view.
    pub camera: FrameCamera,
    /// Shading.
    pub look: Look,
    /// Background, sRGB bytes.
    pub background: [u8; 3],
    /// A section plane `normal · x = offset`, removing the side the
    /// normal points to.
    pub section: Option<(DVec3, f64)>,
    /// The model's static bounds: the ground and its shadow fit them
    /// (`None`: no ground).
    pub ground: Option<[DVec3; 2]>,
}

/// Supersampling factor per axis.
pub const SUPERSAMPLE: usize = 2;

/// A transformed triangle ready to rasterize.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Tri {
    /// Clip-space corners.
    pub clip: [Vec4; 3],
    /// World positions.
    pub world: [Vec3; 3],
    /// World normals.
    pub normal: [Vec3; 3],
    /// Item index.
    pub item: u32,
}

/// The view-projection matrix of `camera` at `aspect`.
fn view_projection(camera: &FrameCamera, aspect: f64) -> Mat4 {
    let view = DMat4::look_at_rh(camera.eye, camera.target, camera.up);
    let proj = DMat4::perspective_rh_gl(
        camera.vfov_deg.to_radians(),
        aspect,
        camera.near,
        camera.far,
    );
    (proj * view).as_mat4()
}

/// Every visible item's triangles in clip space, split into opaque
/// and translucent (the latter far to near by item centre).
fn triangles(scene: &FrameScene<'_>, vp: Mat4) -> (Vec<Tri>, Vec<Tri>) {
    let mut opaque = Vec::new();
    let mut translucent: Vec<(f64, Vec<Tri>)> = Vec::new();
    for (index, item) in scene.items.iter().enumerate() {
        if item.opacity <= 0.001 {
            continue;
        }
        let world = item.world.as_mat4();
        let normal_matrix = glam::Mat3::from_mat4(world).inverse().transpose();
        let count = item.positions.len() / 3;
        let wp: Vec<Vec3> = (0..count)
            .map(|v| {
                world.transform_point3(Vec3::new(
                    item.positions[3 * v],
                    item.positions[3 * v + 1],
                    item.positions[3 * v + 2],
                ))
            })
            .collect();
        let wn: Vec<Vec3> = (0..count)
            .map(|v| {
                (normal_matrix
                    * Vec3::new(
                        item.normals.get(3 * v).copied().unwrap_or(0.0),
                        item.normals.get(3 * v + 1).copied().unwrap_or(0.0),
                        item.normals.get(3 * v + 2).copied().unwrap_or(1.0),
                    ))
                .normalize_or_zero()
            })
            .collect();
        let clip: Vec<Vec4> = wp.iter().map(|p| vp * p.extend(1.0)).collect();
        let mut tris = Vec::with_capacity(item.indices.len() / 3);
        for t in item.indices.as_chunks::<3>().0.iter() {
            let [a, b, c] = [t[0] as usize, t[1] as usize, t[2] as usize];
            if a >= count || b >= count || c >= count {
                continue;
            }
            tris.push(Tri {
                clip: [clip[a], clip[b], clip[c]],
                world: [wp[a], wp[b], wp[c]],
                normal: [wn[a], wn[b], wn[c]],
                item: index as u32,
            });
        }
        if item.opacity >= 0.999 {
            opaque.extend(tris);
        } else {
            let centre = wp.iter().fold(Vec3::ZERO, |s, p| s + *p) / (count.max(1) as f32);
            let depth = centre.as_dvec3().distance(scene.camera.eye);
            translucent.push((depth, tris));
        }
    }
    translucent.sort_by(|a, b| b.0.total_cmp(&a.0));
    (
        opaque,
        translucent.into_iter().flat_map(|(_, t)| t).collect(),
    )
}

/// Render `scene` at `width × height`: sRGB RGBA bytes, row-major from
/// the top.
pub fn render_frame(scene: &FrameScene<'_>, width: u32, height: u32) -> Vec<u8> {
    let (w, h) = (width.max(1) as usize, height.max(1) as usize);
    let (sw, sh) = (w * SUPERSAMPLE, h * SUPERSAMPLE);
    let vp = view_projection(&scene.camera, w as f64 / h as f64);
    let (opaque, translucent) = triangles(scene, vp);
    let edges = raster::edge_segments(scene, vp);
    let shadow = scene
        .ground
        .map(|bounds| shadow::ShadowMap::build(scene, &opaque, &translucent, bounds));
    let context = raster::Context {
        scene,
        vp_inverse: vp.inverse(),
        width: sw,
        height: sh,
        shadow: shadow.as_ref(),
    };
    // Row bands, each rendered whole by one thread.
    let band = 16;
    let bands: Vec<Vec<[f32; 3]>> = (0..sh.div_ceil(band))
        .into_par_iter()
        .map(|b| {
            let rows = b * band..((b + 1) * band).min(sh);
            raster::render_band(&context, &opaque, &translucent, &edges, rows)
        })
        .collect();
    let linear: Vec<[f32; 3]> = bands.into_iter().flatten().collect();
    // Box filter down to the output size, then sRGB bytes.
    let mut out = vec![0u8; w * h * 4];
    let n = (SUPERSAMPLE * SUPERSAMPLE) as f32;
    for y in 0..h {
        for x in 0..w {
            let mut sum = [0.0f32; 3];
            for dy in 0..SUPERSAMPLE {
                for dx in 0..SUPERSAMPLE {
                    let p = linear[(y * SUPERSAMPLE + dy) * sw + x * SUPERSAMPLE + dx];
                    for c in 0..3 {
                        sum[c] += p[c];
                    }
                }
            }
            let o = (y * w + x) * 4;
            for c in 0..3 {
                out[o + c] =
                    (shade::linear_to_srgb(sum[c] / n) * 255.0 + 0.5).clamp(0.0, 255.0) as u8;
            }
            out[o + 3] = 255;
        }
    }
    out
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
