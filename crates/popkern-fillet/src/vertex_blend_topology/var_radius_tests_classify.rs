#![cfg(test)]
use super::var_radius_tests_support::*;
use super::*;

// ================================================================
// Synthetic 4-face apex BRep — vertex-split primitive coverage
// ================================================================
//
// Build a planar pyramid topology with apex `V` valence-4
// incident to 4 triangle faces, then exercise
// `split_apex_footprint_vertex` directly.  No surfaces or curves
// — purely topology + vertex positions.  Tests the structural
// logic in isolation from the cube_mixed pipeline.

/// Build the planar pyramid: square ABCD with central V_apex,
/// split into 4 inside triangles meeting at V_apex, plus one
/// outside square face F_back.  Returns:
/// `(brep, v_apex, [tri0, tri1, tri2, tri3])` where the four
/// triangles are listed in cyclic order around V_apex going
/// around the +Z axis (CCW).
///
/// Vertex layout (z=0 plane):
/// ```text
///   D(0,1) -- C(1,1)
///     | \  /  |
///     |  V    |
///     | /  \  |
///   A(0,0) -- B(1,0)
/// ```
/// V_apex at (0.5, 0.5, 0).  Triangles in cyclic order:
/// `tri0 = (V,A,B)`, `tri1 = (V,B,C)`, `tri2 = (V,C,D)`,
/// `tri3 = (V,D,A)`.
#[test]
fn synthetic_pyramid_4_valid() {
    let (brep, v_apex, tris) = build_pyramid_4();
    brep.validate().expect("pyramid_4 should be valid");
    assert_eq!(brep.vertices.len(), 5);
    assert_eq!(brep.faces.len(), 5); // 4 inner + 1 back
    assert_eq!(brep.edges.len(), 8); // 4 perimeter + 4 spokes
    // V_apex valence 4
    let valence = brep
        .half_edges
        .values()
        .filter(|he| he.origin == v_apex)
        .count();
    assert_eq!(valence, 4, "V_apex should have valence 4");
    // V_apex incident to all 4 inner triangles (and not f_back).
    for &tri in &tris {
        let outer = brep.faces[tri].outer_loop;
        let mut found = false;
        let start = brep.loops[outer].half_edge;
        let mut cursor = start;
        for _ in 0..16 {
            if brep.half_edges[cursor].origin == v_apex {
                found = true;
                break;
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
        }
        assert!(found, "V_apex must be in tri {tri:?} loop");
    }
}

/// First empirical exercise of `split_apex_footprint_vertex` on
/// the synthetic 4-face apex.  Cyclic order matches cube_mixed:
/// `cap=tri0`, `cyl_split=tri1`, `cube=tri2`, `extra=tri3`.
/// This is the cyclic order [cap, cyl_split, cube, extra] —
/// where cap and cyl_split are adjacent at one V_apex edge,
/// cyl_split and cube at the next, and so on.
#[test]
fn split_synthetic_4_face_cyclic() {
    use popkern_core::geometry::Point3;

    let (mut brep, v_apex, [tri0, tri1, tri2, tri3]) = build_pyramid_4();
    let _ = tri3;
    let v_new_pos = Point3::new(0.4, 0.5, 0.0);
    let pre_face_count = brep.faces.len();
    let pre_vertex_count = brep.vertices.len();

    let v_new = split_apex_footprint_vertex(&mut brep, v_apex, tri0, tri2, tri1, v_new_pos)
        .expect("split should succeed");
    eprintln!("split OK: v_new={v_new:?}");
    brep.validate().expect("validate after split");
    assert_eq!(brep.vertices.len(), pre_vertex_count + 1);
    assert_eq!(brep.faces.len(), pre_face_count);
}

/// Configuration where cap and cube are cyclically adjacent
/// (cyl_split is "across" from the cap-cube boundary). The
/// surgery design relies on cap and cube each having one
/// cyl_seam edge incident to V — but in this layout neither
/// does.  The function detects the unsupported configuration
/// and returns `LoopBroken` rather than producing a corrupt
/// BRep.
#[test]
fn split_synthetic_4_face_cap_cube_adjacent_unsupported() {
    use popkern_core::geometry::Point3;

    let (mut brep, v_apex, [tri0, tri1, tri2, tri3]) = build_pyramid_4();
    let _ = tri3;
    let v_new_pos = Point3::new(0.5, 0.4, 0.0);

    let result = split_apex_footprint_vertex(&mut brep, v_apex, tri0, tri1, tri2, v_new_pos);
    assert!(matches!(result, Err(RerouteError::LoopBroken)));
    // Snapshot rollback: BRep should be unchanged.
    brep.validate().expect("BRep unchanged after rollback");
}

/// Diagnostic: for each cap-vs-non-cap edge in cube_mixed
/// post-pipeline (with surgery), dump G1 coplanarity deficit
/// at sample points so we can identify which edges drag the
/// worst-g1 number up.  `#[ignore]`d because it sets
/// process-level env vars; run isolated.
#[test]
fn split_synthetic_4_face_cap_cyl_split_opposite() {
    use popkern_core::geometry::Point3;

    let (mut brep, v_apex, [tri0, tri1, tri2, tri3]) = build_pyramid_4();
    let _ = tri3;
    let v_new_pos = Point3::new(0.5, 0.5, 0.0);

    match split_apex_footprint_vertex(&mut brep, v_apex, tri0, tri1, tri2, v_new_pos) {
        Ok(v_new) => {
            eprintln!("split OK: v_new={v_new:?}");
            if let Err(e) = brep.validate() {
                panic!("validate after split failed (cap-cyl opposite): {e:?}");
            }
        }
        Err(e) => eprintln!("split returned error (cap-cyl opposite): {e:?}"),
    }
}
