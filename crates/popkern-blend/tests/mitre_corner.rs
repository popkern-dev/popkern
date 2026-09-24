//! Mitred corners of three or more fillets (PK-522).
//!
//! `VertexMode::Mitred` once built no corner geometry: the corner kept
//! the face-trim residue (support-plane flaps with edges off their
//! strips, 3.5e-2 at the unit-cube corner, r = 0.05), and the
//! corner-resolution gate (PK-516) declined it. The mitre-corner
//! settlement now trims every adjacent strip pair to its exact
//! bisector-ellipse intersection, with all seams meeting at the
//! strips' common mitre point. Where no such point exists (mismatched
//! radii, an irregular corner of four or more strips) the corner
//! declines with a typed reason.

use popkern_blend::{
    BlendError, BlendFailure, FilletOptions, FilletRadius, FilletSpec, VertexMode,
};
use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveGeometry, CurveOnEdge, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

#[path = "support/mesh_gauges.rs"]
mod mesh_gauges;

fn mitred() -> FilletOptions {
    FilletOptions {
        propagate_tangent: false,
        vertex_mode: VertexMode::Mitred,
        ..FilletOptions::default()
    }
}

fn incident_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    brep.edges
        .iter()
        .filter(|(_, e)| {
            let he = &brep.half_edges[e.half_edge];
            he.origin == vertex || brep.half_edges[he.twin].origin == vertex
        })
        .map(|(h, _)| h)
        .collect()
}

fn specs(edges: &[EdgeHandle], radii: &[f64]) -> Vec<FilletSpec> {
    edges
        .iter()
        .zip(radii)
        .map(|(&edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect()
}

fn highest_vertex(brep: &BRep) -> VertexHandle {
    brep.vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("vertex")
}

/// The served vertex at `at`, with the curves of its edges.
fn vertex_star(brep: &BRep, at: Point3) -> (VertexHandle, Vec<CurveGeometry>) {
    let (vertex, _) = brep
        .vertices
        .iter()
        .find(|(_, v)| (v.position - at).length() < 1e-9)
        .unwrap_or_else(|| panic!("no vertex at the mitre point {at:?}"));
    let curves = incident_edges(brep, vertex)
        .into_iter()
        .map(|e| {
            brep.edges[e]
                .curve
                .as_ref()
                .expect("curve")
                .geometry
                .clone()
        })
        .collect();
    (vertex, curves)
}

/// Enclosed volume as the sum of per-face divergence fluxes over each
/// face's own fine tessellation. Each face is meshed independently to
/// the chord height (measured: strip areas within 1.5e-7 of exact).
/// The combined-body tessellation does not honour the chord height on
/// narrow strips (recorded on PK-505).
fn flux_volume(brep: &BRep) -> f64 {
    let options = popkern_mesh::TessellationOptions {
        chord_height: 1e-6,
        min_subdivisions: 8,
        max_subdivisions: 1024,
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

/// Unit cube, the three fillets at the origin corner: each strip
/// removes `(1 − π/4)·r²` per unit length beyond the corner box, and
/// inside the box the solid keeps one octant of the Steinmetz
/// tricylinder, `(2 − √2)·r³`.
fn mitred_cube_volume(r: f64) -> f64 {
    let strips = 3.0 * (1.0 - std::f64::consts::FRAC_PI_4) * r * r * (1.0 - r);
    let corner_box = (std::f64::consts::SQRT_2 - 1.0) * r * r * r;
    1.0 - strips - corner_box
}

fn cube_at_origin_corner() -> (BRep, Vec<EdgeHandle>) {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| (v.position - Point3::new(0.0, 0.0, 0.0)).length() < 1e-12)
        .map(|(h, _)| h)
        .expect("origin corner");
    let edges = incident_edges(&cube, corner);
    assert_eq!(edges.len(), 3);
    (cube, edges)
}

#[test]
fn cube_corner_mitre_is_exact() {
    for r in [0.05, 0.2] {
        let (cube, edges) = cube_at_origin_corner();
        let result = popkern_blend::fillet_edges(
            &cube,
            &specs(&edges, &[r; 3]),
            &mitred(),
            Tolerance::default(),
            None,
        )
        .unwrap_or_else(|e| panic!("r={r}: mitred cube corner must serve: {e}"));
        let report = &result.repair_report.mitre_corner;
        assert_eq!(report.corners_built, 1, "r={r}: {report:?}");
        assert!(report.declines.is_empty(), "r={r}: {report:?}");
        let errors = result.brep.validate_all();
        assert!(errors.is_empty(), "r={r}: validate_all: {errors:?}");
        assert_eq!(result.brep.faces.len(), 9, "r={r}: 6 planes + 3 strips");

        // The three seams meet at the tricylinder vertex.
        let m = r - r / std::f64::consts::SQRT_2;
        let (_, curves) = vertex_star(&result.brep, Point3::new(m, m, m));
        assert_eq!(curves.len(), 3, "r={r}: valence-3 mitre vertex");
        assert!(
            curves
                .iter()
                .all(|c| matches!(c, CurveGeometry::Ellipse { .. })),
            "r={r}: every seam is an exact bisector ellipse"
        );

        let volume = flux_volume(&result.brep);
        let expected = mitred_cube_volume(r);
        // Gate: 1% of the corner-box term (a wrong mitre is a fraction
        // of it) over the oracle's floor (planar faces lose chord area
        // along their arcs: measured 3.3e-6 at r = 0.2).
        let corner_term = (std::f64::consts::SQRT_2 - 1.0) * r * r * r;
        assert!(
            (volume - expected).abs() < 0.01 * corner_term + 1e-6,
            "r={r}: volume {volume:.9} vs analytic {expected:.9}"
        );
    }
}

#[test]
fn regular_pyramid_apexes_mitre() {
    for n in 3..=6 {
        let pyramid = BRep::n_gonal_pyramid(n, 1.0, 1.2, Tolerance::default()).expect("pyramid");
        let edges = incident_edges(&pyramid, highest_vertex(&pyramid));
        let result = popkern_blend::fillet_edges(
            &pyramid,
            &specs(&edges, &vec![0.16; n]),
            &mitred(),
            Tolerance::default(),
            None,
        )
        .unwrap_or_else(|e| panic!("n={n}: regular apex mitre must serve: {e}"));
        assert_eq!(result.repair_report.mitre_corner.corners_built, 1, "n={n}");
        let errors = result.brep.validate_all();
        assert!(errors.is_empty(), "n={n}: validate_all: {errors:?}");
        // By symmetry the mitre point is on the axis.
        let on_axis = result
            .brep
            .vertices
            .values()
            .find(|v| v.position.x.hypot(v.position.y) < 1e-9 && v.position.z > 0.5)
            .map(|v| v.position)
            .expect("mitre vertex on the pyramid axis");
        let (_, curves) = vertex_star(&result.brep, on_axis);
        assert_eq!(curves.len(), n, "n={n}: one seam per strip pair");
        assert!(
            curves
                .iter()
                .all(|c| matches!(c, CurveGeometry::Ellipse { .. })),
            "n={n}: every seam is an exact bisector ellipse"
        );
    }
}

/// The square pyramid with its apex moved off the axis: the four
/// side planes no longer share a point at distance r from all of them,
/// so the fillet axes are not concurrent and no mitre point exists.
fn skewed_square_pyramid(dx: f64, dy: f64) -> BRep {
    let mut b = BRep::n_gonal_pyramid(4, 1.0, 1.2, Tolerance::default()).expect("pyramid");
    let apex = highest_vertex(&b);
    b.vertices[apex].position.x += dx;
    b.vertices[apex].position.y += dy;
    let edges: Vec<EdgeHandle> = b.edges.keys().collect();
    for e in edges {
        let he = b.edges[e].half_edge;
        let s = b.vertices[b.half_edges[he].origin].position;
        let t = b.vertices[b.half_edges[b.half_edges[he].twin].origin].position;
        b.edges[e].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: s,
                direction: t - s,
            },
            t_start: 0.0,
            t_end: 1.0,
        });
    }
    let faces: Vec<_> = b.faces.keys().collect();
    for f in faces {
        let h0 = b.loops[b.faces[f].outer_loop].half_edge;
        let h1 = b.half_edges[h0].next;
        let h2 = b.half_edges[h1].next;
        let [p0, p1, p2] = [h0, h1, h2].map(|h| b.vertices[b.half_edges[h].origin].position);
        let normal = (p1 - p0).cross(&(p2 - p0)).normalized().expect("plane");
        let u_axis = (p1 - p0).normalized().expect("plane");
        let v_axis = normal.cross(&u_axis);
        b.faces[f].surface.as_mut().expect("surface").geometry = SurfaceGeometry::Plane {
            origin: p0,
            normal,
            u_axis,
            v_axis,
        };
    }
    b.attach_planar_pcurves().expect("pcurves");
    b.validate().expect("skewed pyramid validates");
    b
}

fn declined_corner_reason(
    input: &BRep,
    result: Result<popkern_blend::BlendResult, BlendError>,
) -> String {
    let Err(BlendError::PartialFailure {
        partial_result,
        failures,
        ..
    }) = result
    else {
        panic!("expected a transactional corner decline");
    };
    assert_eq!(
        partial_result.faces.len(),
        input.faces.len(),
        "input handed back"
    );
    failures
        .iter()
        .find_map(|(_, f)| match f {
            BlendFailure::VertexBlendFailed { reason, .. } => Some(reason.clone()),
            _ => None,
        })
        .unwrap_or_else(|| panic!("expected VertexBlendFailed: {failures:?}"))
}

#[test]
fn irregular_valence_four_apex_declines_without_a_mitre_point() {
    let pyramid = skewed_square_pyramid(0.2, 0.15);
    let edges = incident_edges(&pyramid, highest_vertex(&pyramid));
    let result = popkern_blend::fillet_edges(
        &pyramid,
        &specs(&edges, &[0.16; 4]),
        &mitred(),
        Tolerance::default(),
        None,
    );
    let reason = declined_corner_reason(&pyramid, result);
    assert!(
        reason.contains("mitre declined: no point lies on every strip"),
        "typed no-common-point decline: {reason}"
    );
}

#[test]
fn mismatched_radii_mitre_declines_typed() {
    let (cube, edges) = cube_at_origin_corner();
    let result = popkern_blend::fillet_edges(
        &cube,
        &specs(&edges, &[0.05, 0.08, 0.05]),
        &mitred(),
        Tolerance::default(),
        None,
    );
    let reason = declined_corner_reason(&cube, result);
    assert!(
        reason.contains("mitre declined: strip radii differ"),
        "typed mismatched-radii decline: {reason}"
    );
}

fn triangular_prism() -> BRep {
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 0.8, 0.0),
    ];
    BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .expect("prism")
}

/// A quad prism with 72°/108° profile corners and one skew side.
fn oblique_prism() -> BRep {
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.3, 0.7, 0.0),
        Point3::new(-0.2, 0.6, 0.0),
    ];
    BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .expect("prism")
}

/// Oblique trihedral corners (no symmetry): the residue carries stub
/// chains instead of quad flaps, and the reduction to the
/// rail-crossing polygon handles both. Before the reduction was
/// generic every one of these corners declined.
#[test]
fn oblique_trihedral_corners_mitre() {
    for (label, body) in [
        ("triangular", triangular_prism()),
        ("oblique", oblique_prism()),
    ] {
        let tops: Vec<VertexHandle> = body
            .vertices
            .iter()
            .filter(|(_, v)| v.position.z > 0.5)
            .map(|(h, _)| h)
            .collect();
        for vertex in tops {
            let edges = incident_edges(&body, vertex);
            let result = popkern_blend::fillet_edges(
                &body,
                &specs(&edges, &[0.08; 3]),
                &mitred(),
                Tolerance::default(),
                None,
            )
            .unwrap_or_else(|e| panic!("{label} {vertex:?}: oblique mitre must serve: {e}"));
            assert_eq!(
                result.repair_report.mitre_corner.corners_built, 1,
                "{label}"
            );
            let errors = result.brep.validate_all();
            assert!(
                errors.is_empty(),
                "{label} {vertex:?}: validate_all: {errors:?}"
            );
        }
    }
}

/// Every edge of a body filleted mitred: every corner in one
/// operation, acute (pyramid base), obtuse (prism) and valence 4
/// (pyramid apex) alike.
#[test]
fn whole_bodies_mitre_every_corner() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let pyramid = BRep::n_gonal_pyramid(4, 1.0, 1.2, Tolerance::default()).expect("pyramid");
    for (label, body, corners) in [
        ("cube", cube, 8),
        ("triangular prism", triangular_prism(), 6),
        ("oblique prism", oblique_prism(), 8),
        ("square pyramid", pyramid, 5),
    ] {
        let edges: Vec<EdgeHandle> = body.edges.keys().collect();
        let result = popkern_blend::fillet_edges(
            &body,
            &specs(&edges, &vec![0.05; edges.len()]),
            &mitred(),
            Tolerance::default(),
            None,
        )
        .unwrap_or_else(|e| panic!("{label}: whole-body mitre must serve: {e}"));
        let report = &result.repair_report.mitre_corner;
        assert_eq!(report.corners_built, corners, "{label}: {report:?}");
        let errors = result.brep.validate_all();
        assert!(errors.is_empty(), "{label}: validate_all: {errors:?}");
        assert_eq!(
            mesh_gauges::welded_open_edges(&result.brep),
            0,
            "{label}: welded mesh closes"
        );
    }
}
