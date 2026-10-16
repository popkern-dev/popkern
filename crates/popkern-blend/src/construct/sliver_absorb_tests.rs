//! Tests for the general §13.1(2) sliver absorption pass.
//!
//! Structure:
//! - **Negative paths** (clean input → no absorptions)
//! - **Class A (zero-area)**: detector + absorber on a synthesised
//!   zero-area face
//! - **Class C (redundant-edge)**: `surfaces_equal` helper unit test
//! - **Regression guard**: the full 12-edge cube fillet still
//!   produces exactly the 26-face count pinned by
//!   `tests/probe_multi_chord.rs`

use super::*;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// NOTE: the collapsed-face fixtures below move VERTICES and leave
// the cube's stored Line curves stale. That is deliberate and pinned:
// mid-pipeline surgery does exactly this (PK-159 interference seam),
// so vertex positions are the classifiers' live truth for open edges
// — stored curves arbitrate only on topologically closed edges,
// where the vertex skeleton is structurally silent (PK-296).

// ---------------------------------------------------------------------------
// Negative path — clean cube has no slivers at all.
// ---------------------------------------------------------------------------

#[test]
fn absorb_slivers_is_noop_on_pristine_cube() {
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let original_faces = brep.faces.len();
    let warnings = absorb_slivers(&mut brep, &tol());
    assert_eq!(
        brep.faces.len(),
        original_faces,
        "pristine cube should not lose faces to sliver absorption"
    );
    assert!(
        warnings.is_empty(),
        "pristine cube should emit no retained-sliver warnings"
    );
}

// ---------------------------------------------------------------------------
// Class D — corner fragment regression (the existing pattern)
// ---------------------------------------------------------------------------
//
// The end-to-end test `twelve_edge_cube_fillet_face_type_count`
// in `tests/probe_multi_chord.rs` is the canonical regression
// guard for Class D absorption. It asserts `6 + 12 + 8 = 26`
// faces after a 12-edge cube fillet, which only holds if the
// corner-fragment absorber fires correctly. Re-asserting the
// invariant here would be redundant; the integration test stays
// as the load-bearing check.

// ---------------------------------------------------------------------------
// Class A — zero-area face
// ---------------------------------------------------------------------------

#[test]
fn find_zero_area_target_detects_collapsed_cube_face() {
    // Synthesise a zero-area face by collapsing every vertex of
    // one cube face's outer loop to a single position. The face
    // retains its 4-edge topology but its polygonal area drops
    // to (nearly) zero. The Class A classifier should return
    // that face.
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let target = brep.faces.keys().next().unwrap();
    let face_data = &brep.faces[target];
    let start_he = brep.loops[face_data.outer_loop].half_edge;
    let pivot_vertex = brep.half_edges[start_he].origin;
    let pivot_pos = brep.vertices[pivot_vertex].position;

    let mut to_collapse = Vec::new();
    let mut cursor = start_he;
    loop {
        let origin = brep.half_edges[cursor].origin;
        if origin != pivot_vertex {
            to_collapse.push(origin);
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == start_he {
            break;
        }
    }
    for v in to_collapse {
        brep.vertices[v].position = pivot_pos;
    }

    let retained = std::collections::HashSet::new();
    let result = find_zero_area_target(&brep, &tol(), &retained);
    let (face, _he, area) = result.expect("zero-area target should be detected");
    assert_eq!(
        face, target,
        "the detected zero-area face should match the one we collapsed"
    );
    let tol_sq = tol().linear() * tol().linear();
    assert!(
        area < tol_sq,
        "reported area {area} should be below tolerance² {tol_sq}"
    );
}

#[test]
fn zero_area_not_detected_on_pristine_cube() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let retained = std::collections::HashSet::new();
    let result = find_zero_area_target(&brep, &tol(), &retained);
    assert!(
        result.is_none(),
        "pristine cube face areas are unit-scale, not zero"
    );
}

// ---------------------------------------------------------------------------
// Class B — short edge
// ---------------------------------------------------------------------------

#[test]
fn find_short_edge_target_detects_collapsed_edge() {
    // Collapse one edge of a cube face by moving one vertex to
    // coincide with its successor. The Class B classifier should
    // flag the face as having a short edge.
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let target = brep.faces.keys().next().unwrap();
    let face_data = &brep.faces[target];
    let start_he = brep.loops[face_data.outer_loop].half_edge;
    let v_origin = brep.half_edges[start_he].origin;
    let v_next = {
        let next_he = brep.half_edges[start_he].next;
        brep.half_edges[next_he].origin
    };
    brep.vertices[v_origin].position = brep.vertices[v_next].position;

    let retained = std::collections::HashSet::new();
    let result = find_short_edge_target(&brep, &tol(), &retained);
    let (_face, _he, min_length) = result.expect("short-edge target should be detected");
    assert!(
        min_length < tol().linear(),
        "reported min_length {min_length} should be below tolerance"
    );
}

#[test]
fn short_edge_not_detected_on_pristine_cube() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let retained = std::collections::HashSet::new();
    let result = find_short_edge_target(&brep, &tol(), &retained);
    assert!(
        result.is_none(),
        "pristine cube has unit-length edges — no short-edge detection"
    );
}

/// PK-296: vertex-skeleton measurements lie about loops carrying
/// closed circle edges. A cylinder primitive's lateral loop has the
/// vertex skeleton [bot, bot, top, top] (vertex-polygon area exactly
/// 0) and its rim edges are closed circles (endpoint distance exactly
/// 0 despite 2πr of arc). Both classifiers previously flagged these
/// macroscopic faces as slivers — on the flagship boss plate the
/// closed-fillet strip and both supports shipped three spurious
/// SliverRetained warnings per operation.
#[test]
fn classifiers_do_not_flag_closed_circle_loops() {
    let brep = BRep::cylinder(0.35, 1.0, tol()).unwrap();
    let retained = std::collections::HashSet::new();
    assert!(
        find_zero_area_target(&brep, &tol(), &retained).is_none(),
        "cylinder faces are macroscopic — the curve-sampled boundary area must see them"
    );
    assert!(
        find_short_edge_target(&brep, &tol(), &retained).is_none(),
        "closed rim circles have 2πr of arc — endpoint distance must not classify them short"
    );
}

// ---------------------------------------------------------------------------
// Class C — redundant edge via surfaces_equal
// ---------------------------------------------------------------------------

#[test]
fn surfaces_equal_recognises_identical_planes() {
    let p1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let p2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.5, 0.3, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // Same plane (both are z = 0) even though origins differ.
    assert!(surfaces_equal(&p1, &p2, &tol()));
}

#[test]
fn surfaces_equal_rejects_parallel_but_offset_planes() {
    let p1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let p2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0), // offset along the normal
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // Parallel normals but different offsets → different planes.
    assert!(!surfaces_equal(&p1, &p2, &tol()));
}

#[test]
fn surfaces_equal_rejects_rotated_planes() {
    let p1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let p2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0), // rotated
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    assert!(!surfaces_equal(&p1, &p2, &tol()));
}

#[test]
fn surfaces_equal_recognises_identical_cylinders() {
    let c1 = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let c2 = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 5.0), // same axis line
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert!(surfaces_equal(&c1, &c2, &tol()));
}

#[test]
fn surfaces_equal_rejects_cylinders_with_different_radii() {
    let c1 = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let c2 = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    assert!(!surfaces_equal(&c1, &c2, &tol()));
}

#[test]
fn surfaces_equal_recognises_identical_spheres() {
    let s1 = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 2.0, 3.0),
        radius: 0.5,
        axis: popkern_core::geometry::Vec3::Z,
    };
    let s2 = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 2.0, 3.0),
        radius: 0.5,
        axis: popkern_core::geometry::Vec3::Z,
    };
    assert!(surfaces_equal(&s1, &s2, &tol()));
}

#[test]
fn surfaces_equal_rejects_mismatched_surface_types() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: popkern_core::geometry::Vec3::Z,
    };
    assert!(!surfaces_equal(&plane, &sphere, &tol()));
}
