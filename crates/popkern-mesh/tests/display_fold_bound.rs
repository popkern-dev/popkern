//! In-face crease bound of the display mesh (PK-626).
//!
//! The angular facet bound (`TessellationOptions::max_facet_angle`) is
//! enforced per triangle on VERTEX-normal spread, but what a viewer's
//! crease filter, an STL consumer or flat shading sees is the dihedral
//! between neighbouring facets of one face. The two diverged on the
//! boss plate — the reference part (plate ∪ boss, 0.08 G1 rounds on the
//! seam, the top perimeter and the vertical edges):
//!
//! - The torus seam round's tensor strip was fixed up by the angular
//!   refine pass, whose diagonal bisection staggered half-step rows;
//!   the thin triangles tilted by the long chord's sagitta and creased
//!   their neighbours 25–29° (96 edges), every vertex spread in bound.
//! - The sphere corner patches carried a pole vertex whose outgoing
//!   meridian segment was labelled with the incoming edge's span, so
//!   split rounds halved the wrong edge toward the pole: 4x-graded
//!   samples, 1e-13-area slivers, 76° folds.
//!
//! The gate: welded watertight, every in-face dihedral within the bound
//! (plus a tolerance for the ruled-lattice rows the bound does not
//! govern) away from chart-singular vertices, and no corner-fan sliver
//! below 1e-8 of its face's area.

use std::collections::HashMap;

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::{DISPLAY_MAX_FACET_ANGLE, TessellationOptions, TriMesh, welded_open_edges};
use popkern_sweep::{SweepOptions, SweepProfile, sweep};

/// Closed polyline wire through `points` (degree 1, closed by repeating
/// the first point).
fn closed_polyline(points: &[Point3]) -> NurbsCurve {
    let mut control_points = points.to_vec();
    control_points.push(points[0]);
    let n = control_points.len();
    let mut knots = vec![0.0];
    knots.extend((0..n).map(|i| i as f64 / (n - 1) as f64));
    knots.push(1.0);
    NurbsCurve {
        degree: 1,
        knots,
        control_points,
        weights: vec![1.0; n],
    }
}

/// Rational quadratic full circle of `radius` about the z-axis at `z`.
fn circle(radius: f64, z: f64) -> NurbsCurve {
    let r = radius;
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let control_points = [
        (r, 0.0),
        (r, r),
        (0.0, r),
        (-r, r),
        (-r, 0.0),
        (-r, -r),
        (0.0, -r),
        (r, -r),
        (r, 0.0),
    ]
    .map(|(x, y)| Point3::new(x, y, z))
    .to_vec();
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points,
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

fn extrude(profile: NurbsCurve, z0: f64, z1: f64) -> BRep {
    let path = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, z0), Point3::new(0.0, 0.0, z1)],
        weights: vec![1.0, 1.0],
    };
    sweep(
        &SweepProfile::Wire(profile),
        &path,
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("extrusion")
    .brep
}

/// The boss plate: a 2.4 x 1.6 x 0.4 plate united with an r 0.35 boss
/// (z 0.2..1.0), then 0.08 G1 rounds on the boss seam, the top
/// perimeter and the four vertical edges.
fn boss_plate() -> BRep {
    let tol = Tolerance::default();
    let (hw, hd) = (1.2, 0.8);
    let plate = extrude(
        closed_polyline(&[
            Point3::new(-hw, -hd, 0.0),
            Point3::new(hw, -hd, 0.0),
            Point3::new(hw, hd, 0.0),
            Point3::new(-hw, hd, 0.0),
        ]),
        0.0,
        0.4,
    );
    let boss = extrude(circle(0.35, 0.0), 0.2, 1.0);
    let body = popkern_bool::union(plate, boss, tol, None).expect("union");
    let mut edges = Vec::new();
    for (eh, edge) in body.edges.iter() {
        let he = &body.half_edges[edge.half_edge];
        let a = body.vertices[he.origin].position;
        let b = body.vertices[body.half_edges[he.twin].origin].position;
        let mid = edge.curve.as_ref().map(|c| {
            let curve: &dyn CurveEval = &c.geometry;
            curve.point_at(0.5 * (c.t_start + c.t_end))
        });
        let on_top = |p: Point3| (p.z - 0.4).abs() < 1e-6;
        let seam = mid.is_some_and(|m| on_top(m) && (m.x.hypot(m.y) - 0.35).abs() < 1e-6);
        let top_perimeter = on_top(a) && on_top(b) && a.distance_to(&b) > 1.0;
        let vertical = (a.z - b.z).abs() > 0.3
            && (a.x - b.x).abs() < 1e-9
            && (a.y - b.y).abs() < 1e-9
            && (a.x.abs() - hw).abs() < 1e-9;
        if seam || top_perimeter || vertical {
            edges.push(eh);
        }
    }
    assert_eq!(edges.len(), 9, "seam + 4 top perimeter + 4 vertical edges");
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.08),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G1,
        ..Default::default()
    };
    fillet_edges(&body, &specs, &options, tol, None)
        .expect("rounds serve")
        .brep
}

fn facet_normal(p: [Point3; 3]) -> Vec3 {
    let a = Vec3::new(p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z);
    let b = Vec3::new(p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z);
    a.cross(&b)
}

/// Per-face crease census: the worst in-face dihedral (radians, edges
/// incident to a chart-singular vertex excluded) and the number of
/// triangles below `1e-8` of their face's area, with the face kind.
struct FaceCreases {
    face: FaceHandle,
    kind: &'static str,
    worst: f64,
    worst_at: Point3,
    slivers: usize,
}

fn crease_census(brep: &BRep, mesh: &TriMesh) -> Vec<FaceCreases> {
    let mut census = Vec::new();
    for fm in &mesh.face_meshes {
        let surface = &brep.faces[fm.face]
            .surface
            .as_ref()
            .expect("surface")
            .geometry;
        // Weld the face's own vertices: the pole row of a sphere chart
        // repeats the pole once per station.
        let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
        let ids: Vec<usize> = fm
            .vertices
            .iter()
            .map(|v| {
                let p = v.position;
                let k = (
                    (p.x * 1e9).round() as i64,
                    (p.y * 1e9).round() as i64,
                    (p.z * 1e9).round() as i64,
                );
                let next = key.len();
                *key.entry(k).or_insert(next)
            })
            .collect();
        let singular: Vec<bool> = fm
            .vertices
            .iter()
            .map(|v| pole_at(surface, v.position))
            .collect();
        let positions = |t: [usize; 3]| t.map(|i| fm.vertices[i].position);
        let areas: Vec<f64> = fm
            .triangles
            .iter()
            .map(|&t| 0.5 * facet_normal(positions(t)).length())
            .collect();
        let face_area: f64 = areas.iter().sum();
        let mut adjacency: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (ti, t) in fm.triangles.iter().enumerate() {
            let w = t.map(|i| ids[i]);
            if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
                continue; // a pole-row triangle: zero area by construction
            }
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if singular[a] || singular[b] {
                    continue;
                }
                let (x, y) = (ids[a], ids[b]);
                adjacency.entry((x.min(y), x.max(y))).or_default().push(ti);
            }
        }
        let (mut worst, mut worst_at) = (0.0_f64, Point3::new(0.0, 0.0, 0.0));
        for triangles in adjacency.values() {
            let &[t0, t1] = triangles.as_slice() else {
                continue;
            };
            let (n0, n1) = (
                facet_normal(positions(fm.triangles[t0])),
                facet_normal(positions(fm.triangles[t1])),
            );
            let angle = (n0.dot(&n1) / (n0.length() * n1.length()))
                .clamp(-1.0, 1.0)
                .acos();
            if angle > worst {
                worst = angle;
                worst_at = fm.vertices[fm.triangles[t0][0]].position;
            }
        }
        census.push(FaceCreases {
            face: fm.face,
            kind: surface.kind_name(),
            worst,
            worst_at,
            slivers: areas.iter().filter(|&&a| a < 1e-8 * face_area).count(),
        });
    }
    census
}

/// Whether `p` is a sphere's chart pole (where the chart's normal is
/// undefined and a crease is a chart artefact, not a facet defect).
fn pole_at(surface: &SurfaceGeometry, p: Point3) -> bool {
    let SurfaceGeometry::Sphere { radius, .. } = surface else {
        return false;
    };
    let (_, (v0, v1)) = surface.parameter_range();
    [v0, v1]
        .iter()
        .any(|&v| surface.point_at(0.0, v).distance_to(&p) < 1e-9 * radius.max(1.0))
}

#[test]
fn boss_plate_display_mesh_creases_within_the_facet_bound() {
    let brep = boss_plate();
    let mesh = popkern_mesh::tessellate(&brep, &TessellationOptions::display()).expect("mesh");
    assert_eq!(welded_open_edges(&mesh), 0, "watertight");
    let census = crease_census(&brep, &mesh);
    // Tensor rows and boundary fans may crease slightly past one
    // cell's turn; a staggered refine row or a pole sliver fan reads
    // 2x the bound or more.
    let bound = DISPLAY_MAX_FACET_ANGLE + 0.5_f64.to_radians();
    let over: Vec<String> = census
        .iter()
        .filter(|f| f.worst > bound)
        .map(|f| {
            format!(
                "{:?} {} {:.1}° at ({:.4},{:.4},{:.4})",
                f.face,
                f.kind,
                f.worst.to_degrees(),
                f.worst_at.x,
                f.worst_at.y,
                f.worst_at.z
            )
        })
        .collect();
    assert!(over.is_empty(), "in-face creases over the bound: {over:?}");
    let slivers: Vec<String> = census
        .iter()
        .filter(|f| f.slivers > 0)
        .map(|f| format!("{:?} {} x{}", f.face, f.kind, f.slivers))
        .collect();
    assert!(slivers.is_empty(), "sliver triangles: {slivers:?}");
    // The corner patches and the seam round are present: the census
    // covers the classes the gate is about.
    for kind in ["Torus", "Sphere"] {
        assert!(census.iter().any(|f| f.kind == kind), "no {kind} face");
    }
}
