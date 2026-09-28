//! Seam-class grouping of zebra seam crossings (PK-503): the
//! attribution harness and its tests read image-space continuity per
//! kind of seam — cap–cap, cap rim (cap–strip), strip–support,
//! strip–strip — rather than one number over the whole part.

use std::collections::HashSet;

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};
use popkern_viz::camera::Camera;
use popkern_viz::zebra::{
    ZebraFrame, ZebraModeDisagreement, ZebraNormals, ZebraSeamCrossing, rasterize_zebra_frame,
};

/// The role of a face on either side of a seam.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SeamClass(u8, u8);

impl SeamClass {
    /// Vertex-blend cap face role.
    pub const CAP: u8 = 0;
    /// Edge-blend strip face role.
    pub const STRIP: u8 = 1;
    /// Support (planar) face role.
    pub const SUPPORT: u8 = 2;

    /// The unordered class of a seam between faces of roles `a`, `b`.
    pub fn of(a: u8, b: u8) -> SeamClass {
        SeamClass(a.min(b), a.max(b))
    }

    /// Human-readable class name.
    pub fn label(&self) -> &'static str {
        match (self.0, self.1) {
            (0, 0) => "cap-cap",
            (0, 1) => "cap-rim",
            (0, 2) => "cap-support",
            (1, 1) => "strip-strip",
            (1, 2) => "strip-support",
            _ => "support-support",
        }
    }
}

/// Group crossings by seam class (`role` maps a face to
/// [`SeamClass::CAP`] / [`SeamClass::STRIP`] / [`SeamClass::SUPPORT`]),
/// classes in order.
pub fn classify_crossings(
    crossings: &[ZebraSeamCrossing],
    role: impl Fn(FaceHandle) -> u8,
) -> Vec<(SeamClass, Vec<ZebraSeamCrossing>)> {
    let mut groups: Vec<(SeamClass, Vec<ZebraSeamCrossing>)> = Vec::new();
    for c in crossings {
        let class = SeamClass::of(role(c.left), role(c.right));
        match groups.iter_mut().find(|(k, _)| *k == class) {
            Some((_, v)) => v.push(*c),
            None => groups.push((class, vec![*c])),
        }
    }
    groups.sort_by_key(|(k, _)| *k);
    groups
}

/// Summary of a group of seam crossings. Crossings of DIFFERENT seams
/// form different populations (a round's two rims see different
/// one-sided slopes), so a pooled median sits on a mode boundary and
/// jumps when a handful of rows change; each statistic is therefore
/// taken per seam (unordered face pair) first, and the class reads its
/// worst seam.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SeamSummary {
    /// Number of row crossings.
    pub rows: usize,
    /// Number of distinct seams (face pairs) with at least
    /// [`MIN_SEAM_ROWS`] crossings.
    pub seams: usize,
    /// Worst per-seam median value jump.
    pub value_p50: f64,
    /// Worst per-seam 95th-percentile value jump.
    pub value_p95: f64,
    /// Worst per-seam median normal jump (radians).
    pub normal_p50: f64,
    /// Worst per-seam 95th-percentile normal jump (radians).
    pub normal_p95: f64,
    /// Worst per-seam median slope jump (per pixel).
    pub slope_p50: f64,
    /// Worst per-seam 95th-percentile slope jump (per pixel).
    pub slope_p95: f64,
}

/// Seams crossed by fewer image rows than this carry too few samples
/// for a percentile and are left out of the summary.
pub const MIN_SEAM_ROWS: usize = 8;

fn percentile(mut v: Vec<f64>, q: f64) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(f64::total_cmp);
    v[((v.len() - 1) as f64 * q).round() as usize]
}

/// Summarize a group of crossings (see [`SeamSummary`]).
pub fn summarize(rows: &[ZebraSeamCrossing]) -> SeamSummary {
    let mut seams: Vec<((FaceHandle, FaceHandle), Vec<&ZebraSeamCrossing>)> = Vec::new();
    for c in rows {
        let key = if c.left < c.right {
            (c.left, c.right)
        } else {
            (c.right, c.left)
        };
        match seams.iter_mut().find(|(k, _)| *k == key) {
            Some((_, v)) => v.push(c),
            None => seams.push((key, vec![c])),
        }
    }
    let mut out = SeamSummary {
        rows: rows.len(),
        ..SeamSummary::default()
    };
    for (_, crossings) in seams.iter().filter(|(_, v)| v.len() >= MIN_SEAM_ROWS) {
        let values: Vec<f64> = crossings.iter().map(|c| c.value_jump).collect();
        let slopes: Vec<f64> = crossings.iter().map(|c| c.slope_jump).collect();
        let normals: Vec<f64> = crossings.iter().map(|c| c.normal_jump_rad).collect();
        out.seams += 1;
        out.normal_p50 = out.normal_p50.max(percentile(normals.clone(), 0.5));
        out.normal_p95 = out.normal_p95.max(percentile(normals, 0.95));
        out.value_p50 = out.value_p50.max(percentile(values.clone(), 0.5));
        out.value_p95 = out.value_p95.max(percentile(values, 0.95));
        out.slope_p50 = out.slope_p50.max(percentile(slopes.clone(), 0.5));
        out.slope_p95 = out.slope_p95.max(percentile(slopes, 0.95));
    }
    out
}

/// A vertex-blend fixture framed on its cap.
pub struct AttributionFixture {
    /// Fixture name.
    pub name: &'static str,
    /// The blended body.
    pub brep: BRep,
    /// The vertex-blend cap faces.
    pub caps: HashSet<FaceHandle>,
    /// Close-up camera on the cap (square frame).
    pub camera: Camera,
}

fn fillet_vertex(body: &BRep, vertex: VertexHandle, r: f64) -> (BRep, HashSet<FaceHandle>) {
    fillet_vertex_radii(body, vertex, &[r])
}

/// [`fillet_vertex`] with per-edge radii, cycled over the vertex's
/// edges in handle order (the setback family's convention).
fn fillet_vertex_radii(
    body: &BRep,
    vertex: VertexHandle,
    radii: &[f64],
) -> (BRep, HashSet<FaceHandle>) {
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .into_iter()
        .zip(radii.iter().cycle())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let result =
        fillet_edges(body, &specs, &options, Tolerance::default(), None).expect("vertex fillet");
    let caps = result.corner_patches.iter().map(|p| p.face).collect();
    (result.brep, caps)
}

fn apex(body: &BRep) -> VertexHandle {
    body.vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex")
}

fn look(eye: [f32; 3], target: [f32; 3], up: glam::Vec3) -> Camera {
    Camera {
        eye: glam::Vec3::from(eye),
        target: glam::Vec3::from(target),
        up,
        fov_y: 0.5,
        aspect: 1.0,
        near: 0.005,
        far: 100.0,
    }
}

/// The attribution fixtures: the pentagon apex flagship (the
/// pyramid_apex gallery geometry), pentagon5_blunt at the G2 request,
/// and the cube-corner KP cap.
pub fn attribution_fixtures() -> Vec<AttributionFixture> {
    let tol = Tolerance::default();
    let mut out = Vec::new();
    // The pyramid_apex gallery geometry / PK-519 flagship.
    let pyr = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
    let (brep, caps) = fillet_vertex(&pyr, apex(&pyr), 0.16);
    out.push(AttributionFixture {
        name: "pentagon_flagship",
        brep,
        caps,
        camera: look([0.55, -0.45, 1.75], [0.0, 0.0, 1.0], glam::Vec3::Z),
    });
    // pentagon5_blunt at the G2 request.
    let blunt = BRep::n_gonal_pyramid(5, 1.0, 3.0, tol).expect("blunt");
    let (brep, caps) = fillet_vertex(&blunt, apex(&blunt), 0.05);
    out.push(AttributionFixture {
        name: "pentagon5_blunt_g2",
        brep,
        caps,
        camera: look([0.18, -0.15, 3.1], [0.0, 0.0, 2.85], glam::Vec3::Z),
    });
    // pentagon5_variable at the G2 request (1:3 radii): its rim's end
    // layer at the spoke junctions (PK-581) is the largest rim
    // residual in the setback family.
    let variable = BRep::n_gonal_pyramid(5, 1.0, 0.8, tol).expect("variable");
    let (brep, caps) =
        fillet_vertex_radii(&variable, apex(&variable), &[0.02, 0.04, 0.06, 0.04, 0.02]);
    out.push(AttributionFixture {
        name: "pentagon5_variable_g2",
        brep,
        caps,
        camera: look([0.16, -0.13, 0.98], [0.0, 0.0, 0.74], glam::Vec3::Z),
    });
    // Cube corner KP cap.
    let cube = BRep::unit_cube(tol).expect("cube");
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| {
            v.position
                .distance_to(&popkern_core::geometry::Point3::new(1.0, 1.0, 1.0))
                < 1e-9
        })
        .map(|(h, _)| h)
        .expect("corner");
    let (brep, caps) = fillet_vertex(&cube, corner, 0.2);
    out.push(AttributionFixture {
        name: "cube_corner",
        brep,
        caps,
        camera: look([1.9, 1.7, 1.8], [0.9, 0.9, 0.9], glam::Vec3::Z),
    });
    out
}

/// One fixture's attribution: the mesh-vs-surface normal comparison
/// and, per seam class, the seam readings in both views.
pub struct Attribution {
    /// The deferred frame both views shade.
    pub frame: ZebraFrame,
    /// Mesh-vs-surface normal disagreement.
    pub disagreement: ZebraModeDisagreement,
    /// Per seam class: (class, surface-view summary, mesh-view summary).
    pub classes: Vec<(SeamClass, SeamSummary, SeamSummary)>,
}

impl Attribution {
    /// The summaries of one seam class, if the frame crosses it.
    pub fn class(&self, label: &str) -> Option<(SeamSummary, SeamSummary)> {
        self.classes
            .iter()
            .find(|(c, _, _)| c.label() == label)
            .map(|(_, s, m)| (*s, *m))
    }
}

/// Render `fx` into a `size` × `size` deferred frame on the display
/// mesh and measure its attribution for `stripes` stripes.
pub fn attribute(fx: &AttributionFixture, size: u32, stripes: u32) -> Attribution {
    let mesh = popkern_mesh::tessellate(&fx.brep, &popkern_mesh::TessellationOptions::display())
        .expect("display mesh");
    let frame = rasterize_zebra_frame(&fx.brep, &mesh, &fx.camera, size, size);
    let disagreement = frame.mode_disagreement(stripes);
    let is_plane = |f: FaceHandle| {
        matches!(
            fx.brep
                .faces
                .get(f)
                .and_then(|x| x.surface.as_ref())
                .map(|s| &s.geometry),
            Some(SurfaceGeometry::Plane { .. })
        )
    };
    let role = |f: FaceHandle| {
        if fx.caps.contains(&f) {
            SeamClass::CAP
        } else if is_plane(f) {
            SeamClass::SUPPORT
        } else {
            SeamClass::STRIP
        }
    };
    let surf = classify_crossings(&frame.seam_crossings(ZebraNormals::Surface), role);
    let mesh_x = classify_crossings(&frame.seam_crossings(ZebraNormals::Mesh), role);
    let classes = surf
        .iter()
        .map(|(class, rows)| {
            let m_rows = mesh_x
                .iter()
                .find(|(c, _)| c == class)
                .map(|(_, r)| r.as_slice())
                .unwrap_or(&[]);
            (*class, summarize(rows), summarize(m_rows))
        })
        .collect();
    Attribution {
        frame,
        disagreement,
        classes,
    }
}
