//! Marched mitre corners (PK-528).
//!
//! PK-522 mitred only equal-radius circular-cylinder corners, where the
//! strip axes are concurrent and every seam is an exact bisector
//! ellipse. Mixed radii declined (`MismatchedRadii`), and so did
//! curvature-continuous (G2) or variable-radius strips
//! (`NonCylindricalStrip`), although a point on every strip exists
//! there too. Those strips now get a mitre point solved on the surfaces
//! themselves and seams marched onto both strips.

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, VertexMode};
use popkern_core::geometry::{CurveGeometry, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

#[path = "support/mesh_gauges.rs"]
mod mesh_gauges;

fn mitred(continuity: ContinuityLevel) -> FilletOptions {
    FilletOptions {
        propagate_tangent: false,
        vertex_mode: VertexMode::Mitred,
        continuity,
        ..FilletOptions::default()
    }
}

/// The unit cube and its three edges at the origin corner, in x, y,
/// z order.
fn cube_corner_edges() -> (BRep, [EdgeHandle; 3]) {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let edge_to = |far: Point3| -> EdgeHandle {
        cube.edges
            .iter()
            .find(|(_, e)| {
                let he = &cube.half_edges[e.half_edge];
                let ends = [
                    cube.vertices[he.origin].position,
                    cube.vertices[cube.half_edges[he.twin].origin].position,
                ];
                let at = |p: Point3| ends.iter().any(|q| (*q - p).length() < 1e-12);
                at(Point3::new(0.0, 0.0, 0.0)) && at(far)
            })
            .map(|(h, _)| h)
            .expect("cube edge")
    };
    let edges = [
        edge_to(Point3::new(1.0, 0.0, 0.0)),
        edge_to(Point3::new(0.0, 1.0, 0.0)),
        edge_to(Point3::new(0.0, 0.0, 1.0)),
    ];
    (cube, edges)
}

/// Volume of the mitred cube: the intersection of the three rounded
/// solids, one per fillet. Beyond the corner box `[0, R]³`
/// (R = largest radius) each strip removes `(1 − π/4)·r²` per unit
/// length. Inside the box, the removed union is integrated by
/// midpoint quadrature.
fn mitred_cube_volume(radii: [f64; 3]) -> f64 {
    let big = radii.iter().copied().fold(0.0, f64::max);
    let strips: f64 = radii
        .iter()
        .map(|r| (1.0 - std::f64::consts::FRAC_PI_4) * r * r * (1.0 - big))
        .sum();
    let [rx, ry, rz] = radii;
    let removed =
        |a: f64, b: f64, r: f64| a < r && b < r && (a - r).powi(2) + (b - r).powi(2) > r * r;
    let n = 240;
    let h = big / n as f64;
    let mut box_removed = 0.0;
    for i in 0..n {
        let x = (i as f64 + 0.5) * h;
        for j in 0..n {
            let y = (j as f64 + 0.5) * h;
            for k in 0..n {
                let z = (k as f64 + 0.5) * h;
                if removed(y, z, rx) || removed(x, z, ry) || removed(x, y, rz) {
                    box_removed += h * h * h;
                }
            }
        }
    }
    1.0 - strips - box_removed
}

/// Enclosed volume from per-face fine tessellations (the combined-body
/// path does not honour the chord height on narrow strips, PK-505).
fn flux_volume(brep: &BRep) -> f64 {
    let options = popkern_mesh::TessellationOptions {
        chord_height: 1e-6,
        min_subdivisions: 8,
        max_subdivisions: 1024,
        ..Default::default()
    };
    let mut volume = 0.0;
    for face in brep.faces.keys() {
        let mesh = popkern_mesh::tessellate_face(brep, face, &options).expect("face meshes");
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]].position,
                mesh.vertices[t[1]].position,
                mesh.vertices[t[2]].position,
            );
            volume += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x);
        }
    }
    volume.abs() / 6.0
}

/// The single vertex no input vertex explains: the mitre point, with
/// its incident edge curves.
fn mitre_vertex_curves(input: &BRep, out: &BRep) -> Vec<CurveGeometry> {
    let (vertex, _) = out
        .vertices
        .iter()
        .filter(|(_, v)| v.position.x.max(v.position.y).max(v.position.z) < 0.2)
        .find(|(_, v)| {
            input
                .vertices
                .values()
                .all(|w| (w.position - v.position).length() > 1e-6)
                && [v.position.x, v.position.y, v.position.z]
                    .iter()
                    .all(|c| *c > 1e-9)
        })
        .expect("mitre vertex strictly inside the corner box");
    out.edges
        .values()
        .filter(|e| {
            let he = &out.half_edges[e.half_edge];
            he.origin == vertex || out.half_edges[he.twin].origin == vertex
        })
        .map(|e| e.curve.as_ref().expect("curve").geometry.clone())
        .collect()
}

fn serve(radii: [FilletRadius; 3], continuity: ContinuityLevel) -> (BRep, BRep) {
    let (cube, edges) = cube_corner_edges();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii)
        .map(|(&edge, radius)| FilletSpec { edge, radius })
        .collect();
    let result = popkern_blend::fillet_edges(
        &cube,
        &specs,
        &mitred(continuity),
        Tolerance::default(),
        None,
    )
    .unwrap_or_else(|e| panic!("marched mitre must serve: {e}"));
    let report = &result.repair_report.mitre_corner;
    assert_eq!(report.corners_built, 1, "{report:?}");
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    assert_eq!(
        mesh_gauges::welded_open_edges(&result.brep),
        0,
        "welded mesh closes"
    );
    (cube, result.brep)
}

/// Mixed radii: the strip axes miss one another, so no bisector
/// ellipse exists, yet the three cylinders share a point. The served
/// solid must be the intersection of the three rounded solids.
#[test]
fn mixed_radius_mitres_are_the_rounded_solid_intersection() {
    for radii in [[0.05, 0.08, 0.05], [0.05, 0.08, 0.12]] {
        let (cube, out) = serve(radii.map(FilletRadius::Constant), ContinuityLevel::G1);
        let curves = mitre_vertex_curves(&cube, &out);
        assert_eq!(curves.len(), 3, "{radii:?}: valence-3 mitre vertex");
        assert!(
            curves.iter().all(|c| matches!(c, CurveGeometry::Nurbs(_))),
            "{radii:?}: marched seams"
        );
        let big = radii.iter().copied().fold(0.0, f64::max);
        let volume = flux_volume(&out);
        let expected = mitred_cube_volume(radii);
        assert!(
            (volume - expected).abs() < 0.01 * big * big * big + 2e-6,
            "{radii:?}: volume {volume:.9} vs rounded-solid intersection {expected:.9}"
        );
    }
}

/// Curvature-continuous (G2) strips are not circular cylinders.
#[test]
fn g2_strip_corner_mitres() {
    let (cube, out) = serve([0.08; 3].map(FilletRadius::Constant), ContinuityLevel::G2);
    assert_eq!(mitre_vertex_curves(&cube, &out).len(), 3);
}

/// Variable-radius strips are swept rational surfaces.
#[test]
fn variable_radius_corner_mitres() {
    let law = |end: f64| FilletRadius::Variable { start: 0.06, end };
    let (cube, out) = serve([law(0.12), law(0.1), law(0.08)], ContinuityLevel::G1);
    assert_eq!(mitre_vertex_curves(&cube, &out).len(), 3);
}

/// Every edge filleted with radii cycling 0.04 / 0.06 / 0.08: every
/// corner mixes radii, and the seams cross the cylinder charts' period
/// seam (a projection landing across it once jumped the pcurve by 2π:
/// 3.7e-3 off), and meet at crossings where both strips are tangent to
/// the support plane (alternating projection stalled 1.2e-5..3.5e-5
/// off there; the slice Newton does not).
#[test]
fn whole_bodies_mitre_every_mixed_radius_corner() {
    use popkern_core::geometry::Vec3;
    let tol = Tolerance::default();
    let oblique = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.3, 0.7, 0.0),
            Point3::new(-0.2, 0.6, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        tol,
    )
    .expect("prism");
    for (label, body, corners) in [
        ("cube", BRep::unit_cube(tol).expect("cube"), 8),
        ("oblique prism", oblique, 8),
        (
            "triangular pyramid",
            BRep::n_gonal_pyramid(3, 1.0, 1.2, tol).expect("pyramid"),
            4,
        ),
    ] {
        let specs: Vec<FilletSpec> = body
            .edges
            .keys()
            .enumerate()
            .map(|(i, edge)| FilletSpec {
                edge,
                radius: FilletRadius::Constant([0.04, 0.06, 0.08][i % 3]),
            })
            .collect();
        let result =
            popkern_blend::fillet_edges(&body, &specs, &mitred(ContinuityLevel::G1), tol, None)
                .unwrap_or_else(|e| {
                    panic!("{label}: mixed-radius whole-body mitre must serve: {e}")
                });
        assert_eq!(
            result.repair_report.mitre_corner.corners_built, corners,
            "{label}"
        );
        let errors = result.brep.validate_all();
        assert!(errors.is_empty(), "{label}: validate_all: {errors:?}");
        assert_eq!(mesh_gauges::welded_open_edges(&result.brep), 0, "{label}");
    }
}
