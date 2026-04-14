//! Unit tests for the §13.1 degeneracy detectors.

use super::*;
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
    let result = super::check_degenerate_support_normal(&brep, &analysis, &tol());
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
    let result = super::check_nearly_tangent_supports(&brep, &analysis, &tol());
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
    let mut current = start_he;
    loop {
        let origin = brep.half_edges[current].origin;
        if origin != first_vertex {
            vertices_to_collapse.push(origin);
        }
        current = brep.half_edges[current].next;
        if current == start_he {
            break;
        }
    }
    for v in vertices_to_collapse {
        brep.vertices[v].position = first_pos;
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
    let result = super::check_parameterization_singularity(&analysis);
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

    let result = super::check_parameterization_singularity(&analysis);
    match result {
        Err(DegeneracyKind::ParameterizationSingularity { side, .. }) => {
            assert_eq!(side, SupportSide::Face1);
        }
        other => panic!("expected ParameterizationSingularity, got {other:?}"),
    }
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
