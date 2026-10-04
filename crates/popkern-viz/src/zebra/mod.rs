//! Reflection-line zebra rendering with attributable normals (PK-503).
//!
//! A zebra stripe defect has two possible owners: the SURFACE (a real
//! G1/G2 break, a pinched patch) or the TESSELLATION (vertex normals
//! interpolated linearly across facets that turn by up to the facet
//! bound). One image cannot tell them apart. This module rasterizes
//! the mesh ONCE into a deferred frame — per pixel: the visible
//! triangle's face, perspective-correct position, interpolated vertex
//! normal and interpolated native `(u, v)` — and shades that frame in
//! either normal source:
//!
//! - [`ZebraNormals::Mesh`]: the interpolated vertex normal (what the
//!   popkern live's GPU shader shows);
//! - [`ZebraNormals::Surface`]: the face surface's own normal
//!   evaluated at the interpolated `(u, v)`, in the face's facing
//!   sense, with the one-sided limit at chart singularities
//!   ([`popkern_mesh::surface_limit_normal`]).
//!
//! Both views share triangles, camera, visibility, view vectors and
//! stripe set; only the normal differs. Attribution rule: a stripe
//! defect present in the surface view is a surface defect; one present
//! only in the mesh view is tessellation. [`ZebraFrame::mode_disagreement`]
//! measures the per-pixel angle between the two normals so the
//! tessellation's share of any image is a number, not a judgement.

use popkern_core::geometry::surface::SurfaceEval as _;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::TriMesh;

use crate::camera::Camera;

mod facing;
mod footprint;
mod locate;
mod reassign;
mod report;

pub use report::{
    ZebraFaceDisagreement, ZebraModeDisagreement, ZebraSeamCrossing, ZebraWorstPixel,
};

/// Which normal a zebra view reflects in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ZebraNormals {
    /// Interpolated mesh vertex normals (the display pipeline's view).
    #[default]
    Mesh,
    /// The face surface's normal evaluated at the interpolated native
    /// `(u, v)` (the geometry's own view, tessellation-independent up
    /// to the chord error of the sample point).
    Surface,
}

/// Background shade of an uncovered pixel.
const BACKGROUND: [u8; 4] = [38, 38, 38, 255];
/// Dark and light stripe shades.
const DARK: u8 = 24;
const LIGHT: u8 = 236;

/// One visible pixel of a [`ZebraFrame`].
#[derive(Debug, Clone, Copy)]
pub struct ZebraSample {
    /// Face owning the visible triangle (the mesh view's face).
    pub face: FaceHandle,
    /// Index of the visible mesh triangle (into the rendered mesh's
    /// `triangles`), for attributing a pixel's reading to its facet.
    pub triangle: usize,
    /// Face whose surface the pixel's ray actually meets (the surface
    /// view's face). Differs from `face` within about a pixel of a
    /// seam: the display mesh's seam runs through shared edge samples
    /// while the true seam bulges between them, so a triangle of one
    /// face can cover pixels that show its neighbour's surface.
    pub surface_face: FaceHandle,
    /// Native `(u, v)` of the located surface point on `surface_face`.
    pub uv: (f64, f64),
    /// Perspective-correct interpolated world position on the facet.
    pub position: [f64; 3],
    /// Unit interpolated mesh vertex normal.
    pub mesh_normal: [f64; 3],
    /// How `uv` was located on the surface.
    pub hit: SurfaceHit,
    /// The located surface point lies more than one facet-width outside
    /// the visible triangle's own 3D footprint (a silhouette
    /// sliver, PK-578). Near a silhouette the mesh's chord outline and
    /// the true one differ, so a grazing ray through the facet meets
    /// the surface elsewhere. The surface view is right there, but it
    /// does not show the point the facet approximates. The report
    /// counts these pixels apart from the angle statistics.
    pub outside_facet: bool,
    /// Unit surface normal at `uv` in the face's facing sense; `None`
    /// when the face carries no surface or the chart has no regular
    /// point nearby (the surface view then falls back to the mesh
    /// normal and the report counts the pixel as unevaluated).
    pub surface_normal: Option<[f64; 3]>,
    /// Surface normal at the surface point closest to the pixel's
    /// facet point: the point the facet APPROXIMATES, as opposed to the
    /// point the eye sees along the ray. Comparing it with
    /// [`Self::mesh_normal`] measures tessellation accuracy (PK-578);
    /// near grazing the two surface points differ, legitimately. The
    /// facet's interpolated parameter point is not a substitute: it is
    /// undefined at chart poles (ball_bearing read 19° there), and it
    /// hides a chord that crosses iso-lines. coil_spring's helix seam
    /// chords lie on the tube 0.016 in `u` off their vertices'
    /// iso-line, where the true normal is 5.2° from the interpolated
    /// one (PK-587).
    pub facet_surface_normal: Option<[f64; 3]>,
}

/// How a pixel's surface parameter was located.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceHit {
    /// The pixel's view ray intersects the face surface at `uv`
    /// (Newton on `S(u, v) = eye + t·d`, seeded at the facet's
    /// interpolated parameter): the normal is the one the eye sees.
    Ray,
    /// The ray passes just outside the face surface's parameter domain
    /// (a pixel beside a seam whose ray meets the surface's analytic
    /// continuation): `uv` is the domain point nearest the ray, on the
    /// domain boundary.
    DomainBoundary,
    /// The ray misses the claiming face's domain but meets a
    /// neighbouring face's surface inside ITS domain (a pixel between
    /// the mesh seam and the true seam): `uv` is on `surface_face`.
    NeighbourRay,
    /// The ray found no intersection near the seed (grazing, or a
    /// concave facet whose ray passes outside the surface): `uv` is the
    /// closest surface point to the facet point.
    ClosestPoint,
    /// Neither refinement converged: `uv` is the facet's interpolated
    /// parameter.
    Interpolated,
    /// The face carries no surface.
    NoSurface,
}

/// A deferred zebra frame: the visible sample of every pixel, plus
/// the camera-derived eye and stripe axis both views share.
#[derive(Debug, Clone)]
pub struct ZebraFrame {
    /// Frame width in pixels.
    pub width: u32,
    /// Frame height in pixels.
    pub height: u32,
    /// Row-major visible samples (`None` = background).
    pub samples: Vec<Option<ZebraSample>>,
    /// Eye position.
    pub eye: [f64; 3],
    /// Stripe family axis (the camera's right vector): stripes read as
    /// vertical bands of a horizon of parallel lights.
    pub stripe_axis: [f64; 3],
    /// Faces whose mesh vertex normals disagree in sign with their own
    /// surface chart (a mixed facing vote): the mesh's facing on these
    /// faces is internally inconsistent, so the surface view's sense
    /// is a majority decision.
    pub mixed_facing_faces: Vec<FaceHandle>,
}

/// How far outside its facet's 3D footprint a located point may
/// sit, in barycentric units, before the pixel counts as a silhouette
/// sliver (one full facet width: surface curvature moves genuine hits
/// only slightly past an edge).
const FACET_FOOTPRINT_SLACK: f64 = 1.0;

pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let l = dot(a, a).sqrt();
    (l > 1e-300 && l.is_finite()).then(|| [a[0] / l, a[1] / l, a[2] / l])
}

pub(crate) fn angle_between(a: [f64; 3], b: [f64; 3]) -> f64 {
    // atan2 of |a×b| against a·b stays accurate at tiny angles, where
    // acos of a near-unit dot loses half the digits.
    let c = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    dot(c, c).sqrt().atan2(dot(a, b))
}

/// Rasterize `mesh` into a deferred [`ZebraFrame`] seen by `camera`.
///
/// Visibility is resolved first (nearest triangle per pixel, depth in
/// NDC); only visible pixels are then resolved into samples, so each
/// surface evaluation is paid once per final pixel regardless of
/// overdraw.
pub fn rasterize_zebra_frame(
    brep: &BRep,
    mesh: &TriMesh,
    camera: &Camera,
    width: u32,
    height: u32,
) -> ZebraFrame {
    let vp = camera.view_projection();
    let (w, h) = (width as f32, height as f32);
    let near_w = camera.near.max(1.0e-4);
    let n_px = (width as usize) * (height as usize);
    let mut depth = vec![f32::INFINITY; n_px];
    // Visible triangle and its perspective-correct barycentrics.
    let mut visible: Vec<Option<(usize, [f64; 3])>> = vec![None; n_px];

    for (ti, tri) in mesh.triangles.iter().enumerate() {
        let mut sv = [[0.0f32; 3]; 3];
        let mut inv_w = [0.0f64; 3];
        let mut behind = false;
        for (k, &vi) in tri.iter().enumerate() {
            let p = mesh.vertices[vi].position;
            let clip = vp * glam::Vec4::new(p.x as f32, p.y as f32, p.z as f32, 1.0);
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
            inv_w[k] = 1.0 / clip.w as f64;
        }
        if behind {
            continue;
        }
        let area = (sv[1][0] - sv[0][0]) * (sv[2][1] - sv[0][1])
            - (sv[2][0] - sv[0][0]) * (sv[1][1] - sv[0][1]);
        if area.abs() < 1e-12 {
            continue;
        }
        let min_x = (sv
            .iter()
            .map(|v| v[0])
            .fold(f32::INFINITY, f32::min)
            .floor() as i32)
            .max(0);
        let max_x = (sv
            .iter()
            .map(|v| v[0])
            .fold(f32::NEG_INFINITY, f32::max)
            .ceil() as i32)
            .min(width as i32 - 1);
        let min_y = (sv
            .iter()
            .map(|v| v[1])
            .fold(f32::INFINITY, f32::min)
            .floor() as i32)
            .max(0);
        let max_y = (sv
            .iter()
            .map(|v| v[1])
            .fold(f32::NEG_INFINITY, f32::max)
            .ceil() as i32)
            .min(height as i32 - 1);
        for py in min_y..=max_y {
            for px in min_x..=max_x {
                let (fx, fy) = (px as f32 + 0.5, py as f32 + 0.5);
                let b0 =
                    ((sv[1][0] - fx) * (sv[2][1] - fy) - (sv[2][0] - fx) * (sv[1][1] - fy)) / area;
                let b1 =
                    ((sv[2][0] - fx) * (sv[0][1] - fy) - (sv[0][0] - fx) * (sv[2][1] - fy)) / area;
                let b2 = 1.0 - b0 - b1;
                if b0 < 0.0 || b1 < 0.0 || b2 < 0.0 {
                    continue;
                }
                let z = b0 * sv[0][2] + b1 * sv[1][2] + b2 * sv[2][2];
                let idx = (py as usize) * (width as usize) + (px as usize);
                if z >= depth[idx] {
                    continue;
                }
                let q = [
                    b0 as f64 * inv_w[0],
                    b1 as f64 * inv_w[1],
                    b2 as f64 * inv_w[2],
                ];
                let qs = q[0] + q[1] + q[2];
                if qs.abs() < 1e-30 {
                    continue;
                }
                depth[idx] = z;
                visible[idx] = Some((ti, [q[0] / qs, q[1] / qs, q[2] / qs]));
            }
        }
    }

    let eye = [
        camera.eye.x as f64,
        camera.eye.y as f64,
        camera.eye.z as f64,
    ];
    let (signs, mixed_facing_faces) = facing::facing_signs(brep, mesh);
    let sign_of = |face: FaceHandle| -> f64 {
        signs
            .iter()
            .find(|(f, _)| *f == face)
            .map_or(1.0, |(_, s)| *s)
    };
    let resolved: Vec<Option<(ZebraSample, f64)>> = visible
        .into_iter()
        .map(|hit| {
            let (ti, c) = hit?;
            let tri = mesh.triangles[ti];
            let face = mesh.face_map[ti];
            let vs = [
                &mesh.vertices[tri[0]],
                &mesh.vertices[tri[1]],
                &mesh.vertices[tri[2]],
            ];
            let lerp3 = |f: &dyn Fn(usize) -> [f64; 3]| -> [f64; 3] {
                let (a, b, d) = (f(0), f(1), f(2));
                [
                    c[0] * a[0] + c[1] * b[0] + c[2] * d[0],
                    c[0] * a[1] + c[1] * b[1] + c[2] * d[1],
                    c[0] * a[2] + c[1] * b[2] + c[2] * d[2],
                ]
            };
            let position = lerp3(&|k| {
                let p = vs[k].position;
                [p.x, p.y, p.z]
            });
            let mesh_normal = unit(lerp3(&|k| {
                let n = vs[k].normal;
                [n.x, n.y, n.z]
            }))?;
            let surface = brep
                .faces
                .get(face)
                .and_then(|f| f.surface.as_ref())
                .map(|sof| &sof.geometry);
            let uvs = locate::unwrapped_vertex_uvs(surface, &vs);
            let seed = locate::interpolate_uv(&uvs, c);
            let (uv, hit, facet_uv) = match surface {
                Some(s) => {
                    let query = locate::PixelQuery {
                        eye,
                        facet: position,
                        reach: locate::triangle_reach(&vs),
                        sense: sign_of(face),
                    };
                    let (uv, hit) =
                        locate::locate_on_surface(s, &[seed, uvs[0], uvs[1], uvs[2]], &query);
                    (uv, hit, query.closest_point(s, seed))
                }
                None => (seed, SurfaceHit::NoSurface, None),
            };
            let corners = [0, 1, 2].map(|k| {
                let p = vs[k].position;
                [p.x, p.y, p.z]
            });
            let outside_facet = matches!(hit, SurfaceHit::Ray | SurfaceHit::DomainBoundary)
                && surface.is_some_and(|s| {
                    let q = s.point_at(uv.0, uv.1);
                    footprint::footprint_excess(&corners, [q.x, q.y, q.z]) > FACET_FOOTPRINT_SLACK
                });
            let normal_at = |at: (f64, f64)| {
                surface.and_then(|s| {
                    let n = popkern_mesh::surface_limit_normal(s, at.0, at.1)?;
                    let sg = sign_of(face);
                    Some([sg * n.x, sg * n.y, sg * n.z])
                })
            };
            let surface_normal = normal_at(uv);
            let facet_surface_normal = facet_uv.and_then(normal_at);
            Some((
                ZebraSample {
                    face,
                    triangle: ti,
                    surface_face: face,
                    uv,
                    position,
                    mesh_normal,
                    hit,
                    outside_facet,
                    surface_normal,
                    facet_surface_normal,
                },
                locate::triangle_reach(&vs),
            ))
        })
        .collect::<Vec<_>>();
    let samples = reassign::reassign_seam_pixels(brep, width, &resolved, eye, &sign_of);

    let forward = (camera.target - camera.eye).normalize_or_zero();
    let right = forward.cross(camera.up).normalize_or_zero();
    ZebraFrame {
        width,
        height,
        samples,
        eye,
        stripe_axis: [right.x as f64, right.y as f64, right.z as f64],
        mixed_facing_faces,
    }
}

impl ZebraSample {
    /// The face the given view shows at this pixel.
    pub fn face_in(&self, normals: ZebraNormals) -> FaceHandle {
        match normals {
            ZebraNormals::Mesh => self.face,
            ZebraNormals::Surface => self.surface_face,
        }
    }

    /// The normal the given view reflects in (the surface view falls
    /// back to the mesh normal on unevaluated pixels).
    pub fn normal(&self, normals: ZebraNormals) -> [f64; 3] {
        match normals {
            ZebraNormals::Mesh => self.mesh_normal,
            ZebraNormals::Surface => self.surface_normal.unwrap_or(self.mesh_normal),
        }
    }
}

impl ZebraFrame {
    /// Reflection coordinate `s ∈ [-1, 1]` of pixel `(x, y)` in the
    /// given view: the reflected view direction's component along the
    /// stripe axis. The stripe phase is `(s/2 + 1/2) · stripe_count`;
    /// a G1 break is a jump in `s` along an image path, a G2 break a
    /// jump in its derivative.
    pub fn reflection_coordinate(&self, x: u32, y: u32, normals: ZebraNormals) -> Option<f64> {
        if x >= self.width || y >= self.height {
            return None;
        }
        let sample = self.samples[(y as usize) * (self.width as usize) + (x as usize)]?;
        Some(self.reflection_of(&sample, normals))
    }

    pub(crate) fn reflection_of(&self, sample: &ZebraSample, normals: ZebraNormals) -> f64 {
        let n = sample.normal(normals);
        let p = sample.position;
        let view = unit([p[0] - self.eye[0], p[1] - self.eye[1], p[2] - self.eye[2]])
            .unwrap_or([0.0, 0.0, 0.0]);
        let vn = 2.0 * dot(view, n);
        let r = [
            view[0] - n[0] * vn,
            view[1] - n[1] * vn,
            view[2] - n[2] * vn,
        ];
        dot(r, self.stripe_axis)
    }

    pub(crate) fn band(
        &self,
        sample: &ZebraSample,
        normals: ZebraNormals,
        stripe_count: u32,
    ) -> bool {
        let phase = (self.reflection_of(sample, normals) * 0.5 + 0.5) * stripe_count.max(1) as f64;
        phase.fract() < 0.5
    }

    /// Shade the frame as an RGBA zebra image (`width * height * 4`).
    pub fn shade(&self, normals: ZebraNormals, stripe_count: u32) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.samples.len() * 4);
        for sample in &self.samples {
            let px = match sample {
                None => BACKGROUND,
                Some(s) => {
                    let shade = if self.band(s, normals, stripe_count) {
                        DARK
                    } else {
                        LIGHT
                    };
                    [shade, shade, shade, 255]
                }
            };
            out.extend_from_slice(&px);
        }
        out
    }
}

#[cfg(test)]
mod tests;
