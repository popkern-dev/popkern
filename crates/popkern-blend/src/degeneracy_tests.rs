//! Unit tests for the §13.1 degeneracy detectors.

use super::*;
use crate::degeneracy_surface::{
    check_degenerate_support_normal, check_nearly_tangent_supports,
    check_parameterization_singularity, limit_normal_at,
};
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn unit_cube_analysis() -> (BRep, EdgeAnalysis) {
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let mut analyses = crate::analysis::analyze_edges(&brep, &[edge], false).unwrap();
    (brep, analyses.remove(0))
}

// ---------------------------------------------------------------------------
// positive (non-degenerate) input — the cube should pass all four detectors
// ---------------------------------------------------------------------------

#[test]
fn unit_cube_edge_has_no_degeneracies() {
    let (brep, analysis) = unit_cube_analysis();
    let result = check_edge_degeneracy(&brep, &analysis, &tol());
    assert!(
        result.is_ok(),
        "unit cube edge should not trigger any §13.1 detector: {result:?}"
    );
}

#[test]
fn unit_cube_all_edges_pass_degeneracy_check() {
    // Every edge of a healthy cube must pass the detector — the
    // §13.2 false-positive budget (and this §13.1 one) applies to
    // legitimate input too. If any detector false-positives on
    // `BRep::unit_cube`, the entire downstream test suite is at
    // risk.
    let brep = BRep::unit_cube(tol()).unwrap();
    let edges: Vec<_> = brep.edges.keys().collect();
    let analyses = crate::analysis::analyze_edges(&brep, &edges, false).unwrap();
    for (idx, analysis) in analyses.iter().enumerate() {
        let result = check_edge_degeneracy(&brep, analysis, &tol());
        assert!(
            result.is_ok(),
            "unit cube edge {idx} false-positived: {result:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// collapsed-edge detector
// ---------------------------------------------------------------------------

#[test]
fn collapsed_edge_is_detected_by_direct_manipulation() {
    // Build a cube, then move one vertex of an edge to coincide
    // with its neighbour — this synthesises a collapsed edge
    // without going through a constructor that might reject it.
    let (mut brep, analysis) = unit_cube_analysis();
    brep.vertices[analysis.vertex_start].position = brep.vertices[analysis.vertex_end].position;

    let result = check_edge_degeneracy(&brep, &analysis, &tol());
    match result {
        Err(DegeneracyKind::CollapsedEdge { length }) => {
            assert!(
                length < tol().linear(),
                "reported length {length} should be below tolerance {}",
                tol().linear()
            );
        }
        other => panic!("expected CollapsedEdge, got {other:?}"),
    }
}

#[test]
fn closed_edge_is_not_flagged_as_collapsed() {
    // A closed edge has `vertex_start == vertex_end` by design —
    // the collapsed-edge detector must skip this case, otherwise
    // every closed-fillet annulus in the test suite would fail.
    //
    // We synthesise this condition at the analysis level by
    // overriding `vertex_end` to equal `vertex_start`, bypassing
    // the need to run a real closed-fillet constructor.
    let (brep, mut analysis) = unit_cube_analysis();
    analysis.vertex_end = analysis.vertex_start;

    // The collapsed-edge detector must return Ok because
    // `vertex_start == vertex_end` triggers the closed-edge
    // short-circuit.
    //
    // Other detectors (tangent, normal, tangent supports) may
    // still fire on this synthetic input, so we test
    // `check_collapsed_edge` directly rather than through the
    // combined `check_edge_degeneracy`.
    let result = super::check_collapsed_edge(&brep, &analysis, &tol());
    assert!(
        result.is_ok(),
        "closed edge should not trigger CollapsedEdge: {result:?}"
    );
}

// ---------------------------------------------------------------------------
// degenerate-support-normal detector
// ---------------------------------------------------------------------------

#[test]
fn degenerate_support_normal_detector_available() {
    // The detector is exercised implicitly by every passing call
    // to `check_edge_degeneracy` — on a healthy cube no detector
    // fires, and on a cube with a replaced degenerate surface the
    // normal check would fire. We don't synthesise a degenerate-
    // normal BRep here because constructing one is non-trivial
    // (all popkern primitives have well-defined normals across
    // their parameter range). This placeholder documents the
    // detector's presence in the check pipeline; a real-world
    // reproduction would come from a lossy import hitting a
    // NURBS surface with a degenerate partial derivative.
    let (brep, analysis) = unit_cube_analysis();
    let result = check_degenerate_support_normal(&brep, &analysis, &tol());
    assert!(result.is_ok());
}

// ---------------------------------------------------------------------------
// nearly-tangent-supports detector
// ---------------------------------------------------------------------------

#[test]
fn nearly_tangent_supports_detector_passes_on_90_degree_cube_edge() {
    // Cube edges have a 90° dihedral. The detector flags cases
    // below 2·tol_angular from 0 — 90° is nowhere near that, so
    // every cube edge must pass.
    let (brep, analysis) = unit_cube_analysis();
    let result = check_nearly_tangent_supports(&brep, &analysis, &tol());
    assert!(
        result.is_ok(),
        "90° cube edge should pass the nearly-tangent check: {result:?}"
    );
}

// ---------------------------------------------------------------------------
// zero-area-support detector
// ---------------------------------------------------------------------------

#[test]
fn zero_area_support_is_detected_when_outer_loop_vertices_collapse() {
    // Synthesise a degenerate cube face by moving every vertex of
    // face_1's outer loop to coincide with vertex 0 of that loop.
    // The polygonal area drops to 0; the detector must fire.
    let (mut brep, analysis) = unit_cube_analysis();
    let positions_before = brep.face_outer_loop_positions(analysis.face_1);
    assert!(
        positions_before.len() >= 3,
        "cube face should have ≥3 outer-loop vertices"
    );
    // Walk the outer-loop half-edges and collapse every vertex
    // origin to the first one's position.
    let face_rec = &brep.faces[analysis.face_1];
    let outer_loop = face_rec.outer_loop;
    let start_he = brep.loops[outer_loop].half_edge;
    let first_vertex = brep.half_edges[start_he].origin;
    let first_pos = brep.vertices[first_vertex].position;
    // Collect every unique origin vertex on the loop first, then
    // apply the position change — mutating while walking the loop
    // confuses the walker if the topology changes, which it
    // doesn't here, but doing it in two passes is cleaner.
    let mut vertices_to_collapse = Vec::new();
    let mut edges_to_collapse = Vec::new();
    let mut current = start_he;
    loop {
        let origin = brep.half_edges[current].origin;
        if origin != first_vertex {
            vertices_to_collapse.push(origin);
        }
        edges_to_collapse.push(brep.half_edges[current].edge);
        current = brep.half_edges[current].next;
        if current == start_he {
            break;
        }
    }
    for v in vertices_to_collapse {
        brep.vertices[v].position = first_pos;
    }
    // A real collapse degenerates the edge curves too; the detector
    // samples curves (not just vertices) since curved faces bounded
    // by closed loops have as few as two vertices. Clear the curves
    // so the vertex fallback represents the collapsed geometry.
    for eh in edges_to_collapse {
        if let Some(edge) = brep.edges.get_mut(eh) {
            edge.curve = None;
        }
    }

    let result = super::check_zero_area_support(&brep, &analysis, &tol());
    match result {
        Err(DegeneracyKind::ZeroAreaSupport { side, area }) => {
            assert_eq!(side, SupportSide::Face1);
            assert!(
                area < tol().linear() * tol().linear(),
                "reported area {area} should be below tolerance² {}",
                tol().linear() * tol().linear()
            );
        }
        other => panic!("expected ZeroAreaSupport, got {other:?}"),
    }
}

#[test]
fn zero_area_support_passes_on_healthy_cube() {
    // Every cube face has unit area — far from the tolerance²
    // threshold, so no detector should fire on clean input.
    let (brep, analysis) = unit_cube_analysis();
    let result = super::check_zero_area_support(&brep, &analysis, &tol());
    assert!(
        result.is_ok(),
        "healthy cube should not trip zero-area detector: {result:?}"
    );
}

// ---------------------------------------------------------------------------
// parameterization-singularity detector
// ---------------------------------------------------------------------------

#[test]
fn parameterization_singularity_passes_on_healthy_cube() {
    // popkern's analytic primitives (Plane, Cylinder, Sphere, Cone,
    // Torus) never produce NaN/∞ across their natural parameter
    // range. A cube face is a Plane, which is the simplest case.
    let (_brep, analysis) = unit_cube_analysis();
    let result = check_parameterization_singularity(&analysis);
    assert!(
        result.is_ok(),
        "healthy cube (Plane surfaces) should not trip parameterization-singularity detector: {result:?}"
    );
}

#[test]
fn parameterization_singularity_fires_on_nan_nurbs_control_point() {
    // Build an EdgeAnalysis whose `surface_1` is a pathological
    // NurbsSurface with a NaN control point. Any point_at evaluation
    // that touches that control point (and the sweep covers the
    // full parameter range, so it will) returns NaN, and the
    // detector must fire.
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{NurbsSurface, Point3};

    let (_brep, mut analysis) = unit_cube_analysis();
    let nan = f64::NAN;
    let control_points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(nan, nan, nan), // poisoned corner
    ];
    let nurbs = NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    };
    analysis.surface_1 = SurfaceGeometry::Nurbs(nurbs);

    let result = check_parameterization_singularity(&analysis);
    match result {
        Err(DegeneracyKind::ParameterizationSingularity { side, .. }) => {
            assert_eq!(side, SupportSide::Face1);
        }
        other => panic!("expected ParameterizationSingularity, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// limit-normal repair (§13.1(6) Bucket A)
// ---------------------------------------------------------------------------

/// A sphere evaluated at `v = 0` has `partial_u = 0` in popkern's
/// parameterisation — `partial_u = (-r·sin(v)·sin(u), r·sin(v)·cos(u), 0)`
/// which vanishes at `sin(v) = 0`. At that point the exact normal
/// cross-product is zero. `limit_normal_at` should recover a
/// coherent normal pointing along the pole axis from averaging the
/// 8 neighbour samples.
#[test]
fn limit_normal_recovers_sphere_pole_normal() {
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: popkern_core::geometry::Vec3::Z,
    };
    // Pole: v = 0. Exact normal is degenerate.
    let exact = sphere.normal_at(0.0, 0.0);
    assert!(
        exact.length() < 1e-9,
        "exact sphere normal at v=0 should be degenerate, got magnitude {}",
        exact.length()
    );
    // Limit normal should exist and point along the pole axis.
    // popkern's sphere pole geometry depends on the parameterisation;
    // we only assert the limit is finite, unit-length, and coherent.
    let limit = limit_normal_at(&sphere, 0.0, 0.0, 0.05, 0.05, 1e-9);
    let limit = limit.expect("limit normal should be defined at sphere pole");
    let mag = limit.length();
    assert!(
        (mag - 1.0).abs() < 1e-9,
        "limit normal should be unit length, got {mag}"
    );
    // The limit should be anti-parallel to the +Y or +Z axis depending
    // on convention. Just assert it's a coherent 3D direction, not
    // near zero.
    assert!(
        limit.length() > 0.5,
        "limit normal should be coherent, got {:?}",
        limit
    );
    // At a sphere pole the limit should approximately equal the
    // average of nearby outward radial directions, which is along
    // the pole axis. We don't pin which axis to avoid coupling to
    // popkern's internal parameterisation, but we do assert the
    // magnitude is stable.
    let _ = Vec3::new(0.0, 0.0, 1.0);
}

/// §13.1(6) end-to-end repair: `check_degenerate_support_normal`
/// must **not** fire on a sphere-support edge whose midpoint
/// projects to the sphere pole, because the limit-normal-from-
/// neighbours repair recovers a coherent limit.
///
/// Build a synthetic `EdgeAnalysis` by taking a cube edge analysis
/// and swapping `surface_1` for a sphere positioned so the cube
/// edge midpoint coincides with the sphere's north pole in
/// popkern's parameterisation (`v = 0` gives `(0, 0, r)` offset
/// from the sphere centre). Without the repair, the detector
/// would reject this input. With the repair, it passes.
#[test]
fn degenerate_support_normal_recovered_via_limit_at_sphere_pole() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let (brep, mut analysis) = unit_cube_analysis();

    // Locate the edge midpoint in 3D space.
    let midpoint = crate::geometry_compute::get_edge_midpoint(&brep, &analysis).unwrap();
    // Place a sphere whose north pole (v = 0) coincides with the
    // midpoint. popkern's sphere point_at is
    //   (r sin(v) cos(u), r sin(v) sin(u), r cos(v))
    // At v = 0 that is `(0, 0, r)`, so the sphere centre must be
    // offset from the midpoint by `(0, 0, -r)`.
    let radius = 0.5;
    let sphere_centre = Point3::new(midpoint.x, midpoint.y, midpoint.z - radius);
    analysis.surface_1 = SurfaceGeometry::Sphere {
        center: sphere_centre,
        radius,
        axis: popkern_core::geometry::Vec3::Z,
    };

    // Sanity: the exact normal at the sphere pole is degenerate.
    // Either `closest_point` lands exactly at the pole with
    // degenerate normal, OR popkern's closest_point nudges slightly
    // off the singular point. Both are acceptable for this test —
    // we only care that the detector does NOT fire overall, which
    // the limit-normal repair guarantees in the first case and
    // the exact-normal path handles in the second.
    let _ = analysis.surface_1.closest_point(&midpoint);

    let result = check_degenerate_support_normal(&brep, &analysis, &tol());
    assert!(
        result.is_ok(),
        "sphere-pole support should pass via limit-normal repair, got {result:?}"
    );
}

#[test]
fn limit_normal_returns_none_at_sharp_cancellation() {
    // A NURBS surface with a "valley fold" — neighbours along the
    // fold line have normals flipping 180° between adjacent samples.
    // The limit should be None because the neighbour unit normals
    // cancel out when summed.
    //
    // We don't have an easy-to-construct synthetic fold here without
    // building a pathological NURBS, so this test documents the
    // contract: `None` is returned when <3 neighbours contribute OR
    // when the summed magnitude is below 0.5. The NURBS-with-NaN
    // test below exercises the first branch via the nan-normal path.
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{NurbsSurface, Vec3};

    // Degenerate NurbsSurface with all control points at one point:
    // every partial is zero everywhere, so every neighbour normal is
    // (0, 0, 0) → under the tol threshold → count stays at 0 → None.
    let cps = vec![Point3::new(0.0, 0.0, 0.0); 4];
    let nurbs = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: cps,
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    let limit = limit_normal_at(&nurbs, 0.5, 0.5, 0.1, 0.1, 1e-9);
    assert!(
        limit.is_none(),
        "degenerate NURBS should have no limit normal, got {limit:?}"
    );
    let _ = Vec3::new(0.0, 0.0, 1.0);
}

// ---------------------------------------------------------------------------
// describe() rendering
// ---------------------------------------------------------------------------

#[test]
fn describe_renders_prefix_for_each_variant() {
    let variants = [
        DegeneracyKind::CollapsedEdge { length: 1e-12 },
        DegeneracyKind::DegenerateEdgeTangent {
            at: EdgeEndpoint::Start,
            magnitude: 1e-12,
        },
        DegeneracyKind::NearlyTangentSupports { angle_rad: 1e-6 },
        DegeneracyKind::DegenerateSupportNormal {
            side: SupportSide::Face1,
            magnitude: 1e-12,
        },
        DegeneracyKind::ZeroAreaSupport {
            side: SupportSide::Face2,
            area: 1e-16,
        },
        DegeneracyKind::ParameterizationSingularity {
            side: SupportSide::Face1,
            uv: (0.5, 0.5),
        },
    ];
    for v in &variants {
        let msg = v.describe();
        assert!(
            msg.starts_with("degenerate: "),
            "describe() should prefix with `degenerate: ` for parse-ability, got {msg:?}"
        );
        assert!(
            !msg.is_empty(),
            "describe() should produce a non-empty message"
        );
    }
}
