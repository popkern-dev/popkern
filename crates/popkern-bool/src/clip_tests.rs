// Structural lint marker: this file is test code, included via #[cfg(test)] in clip.rs.
#[cfg(test)]
use super::*;
use popkern_core::{
    geometry::{Point3, nurbs_curve::NurbsCurve},
    tolerance::Tolerance,
    topology::BRep,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn nurbs_line(p0: Point3, p1: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

fn cos_line(p0: Point3, p1: Point3) -> CurveOnSurface {
    let c = nurbs_line(p0, p1);
    CurveOnSurface {
        curve_3d: c.clone(),
        curve_uv1: c.clone(),
        curve_uv2: c,
        is_closed: false,
    }
}

// ── face_vertex_positions: invalid half_edge returns empty (line 22) ────

#[test]
fn face_vertex_positions_loop_not_in_brep_returns_empty() {
    use popkern_core::topology::{Face, LoopHandle, ShellHandle};
    let mut brep = BRep::new(tol());
    // Use a default LoopHandle that is NOT inserted into brep.loops.
    let face_h = brep.faces.insert(Face {
        outer_loop: LoopHandle::default(), // not in brep.loops → None branch
        inner_loops: Vec::new(),
        shell: ShellHandle::default(),
        surface: None,
        name: None,
    });
    let positions = face_vertex_positions(&brep, face_h);
    assert!(
        positions.is_empty(),
        "face with no loop in brep should return empty positions (line 22)"
    );
}

// ── clip_line_to_face: < 3 positions → None (line 54) ──────────────────

#[test]
fn clip_curve_to_single_face_degenerate_face_no_panic() {
    // Face with only 2 vertices (no polygon — spur edge).
    let mut brep = BRep::new(tol());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let _mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev");
    let curve = cos_line(Point3::new(-0.5, 0.0, 0.0), Point3::new(1.5, 0.0, 0.0));
    // Must not panic — result may be None or Some.
    let _result = clip_curve_to_single_face(&brep, mvfs.face, &curve);
}

// ── clip_curve_to_single_face: off-plane → None (line 189) ─────────────

#[test]
fn clip_curve_to_single_face_off_plane_returns_none() {
    // Square face in z=0 plane.
    let mut brep = BRep::new(tol());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).expect("mef");
    // Curve offset by z=1 — not on the face's z=0 plane.
    let curve = cos_line(Point3::new(0.2, 0.5, 1.0), Point3::new(0.8, 0.5, 1.0));
    let result = clip_curve_to_single_face(&brep, mef.face, &curve);
    assert!(result.is_none(), "off-plane curve should return None");
}

// ── clip_curve_to_faces: curve not in face_b → None (line 293 `?`) ─────

#[test]
fn clip_curve_to_faces_curve_misses_face_b_returns_none() {
    // face_a: square [0,1]×[0,1] in z=0
    let mut brep_a = BRep::new(tol());
    let mvfs_a = brep_a.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs_a");
    let mev1 = brep_a
        .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep_a
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep_a
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_a = brep_a
        .mef(mev3.half_edge_twin, mev1.half_edge)
        .expect("mef_a");
    let face_a = mef_a.face;

    // face_b: square [3,4]×[0,1] — far from face_a
    let mut brep_b = BRep::new(tol());
    let mvfs_b = brep_b.mvfs(Point3::new(3.0, 0.0, 0.0)).expect("mvfs_b");
    let mev4 = brep_b
        .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(4.0, 0.0, 0.0))
        .expect("mev4");
    let mev5 = brep_b
        .mev(mev4.half_edge_twin, Point3::new(4.0, 1.0, 0.0))
        .expect("mev5");
    let mev6 = brep_b
        .mev(mev5.half_edge_twin, Point3::new(3.0, 1.0, 0.0))
        .expect("mev6");
    let mef_b = brep_b
        .mef(mev6.half_edge_twin, mev4.half_edge)
        .expect("mef_b");
    let face_b = mef_b.face;

    // Curve goes through face_a [0,1] but not face_b [3,4].
    let curve = cos_line(Point3::new(0.5, 0.0, 0.0), Point3::new(0.5, 1.0, 0.0));
    let result = clip_curve_to_faces(&brep_a, face_a, &brep_b, face_b, &curve);
    assert!(
        result.is_none(),
        "curve in face_a but not face_b should return None"
    );
}

// ── clip_line_to_face: s_enter > s_exit mid-loop → None (line 117) ─────

/// A diagonal curve to the right of the square face triggers the early-exit
/// branch (line 117) inside the Cyrus-Beck loop when the exit boundary is
/// pushed to a negative parameter before s_enter is reached.
#[test]
fn clip_line_to_face_diagonal_outside_polygon_hits_early_exit() {
    // face: [0,1]×[0,1] square in z=0
    let mut brep_a = BRep::new(tol());
    let mvfs = brep_a.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep_a
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep_a
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep_a
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_a = brep_a
        .mef(mev3.half_edge_twin, mev1.half_edge)
        .expect("mef_a");
    let face_a = mef_a.face;

    // face_b: duplicate structure so clip_curve_to_faces has a valid second face.
    let mut brep_b = BRep::new(tol());
    let mvfs_b = brep_b.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs_b");
    let mev4 = brep_b
        .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev4");
    let mev5 = brep_b
        .mev(mev4.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev5");
    let mev6 = brep_b
        .mev(mev5.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev6");
    let mef_b = brep_b
        .mef(mev6.half_edge_twin, mev4.half_edge)
        .expect("mef_b");
    let face_b = mef_b.face;

    // Diagonal curve from (2,0,0) to (3,1,0): entirely outside [0,1]×[0,1],
    // not parallel to any edge. The right edge sets s_exit < 0 before
    // s_enter is updated, triggering line 117.
    let curve = cos_line(Point3::new(2.0, 0.0, 0.0), Point3::new(3.0, 1.0, 0.0));
    let result = clip_curve_to_faces(&brep_a, face_a, &brep_b, face_b, &curve);
    assert!(
        result.is_none(),
        "diagonal curve outside polygon should return None"
    );
}

// ── clip_line_to_face: s_enter >= s_exit after loop → None (line 122) ───

/// A curve tangent to the polygon corner triggers the post-loop check
/// (line 122) when s_enter equals s_exit at the corner parameter.
#[test]
fn clip_line_to_face_tangent_to_corner_hits_post_loop_check() {
    // face: [0,1]×[0,1] square in z=0
    let mut brep_a = BRep::new(tol());
    let mvfs = brep_a.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep_a
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep_a
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep_a
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_a = brep_a
        .mef(mev3.half_edge_twin, mev1.half_edge)
        .expect("mef_a");
    let face_a = mef_a.face;

    let mut brep_b = BRep::new(tol());
    let mvfs_b = brep_b.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs_b");
    let mev4 = brep_b
        .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev4");
    let mev5 = brep_b
        .mev(mev4.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev5");
    let mev6 = brep_b
        .mev(mev5.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev6");
    let mef_b = brep_b
        .mef(mev6.half_edge_twin, mev4.half_edge)
        .expect("mef_b");
    let face_b = mef_b.face;

    // Curve from (-1,-1,0) to (0,0,0): approaches the origin corner of the
    // [0,1]×[0,1] square. Cyrus-Beck clamps s_enter == s_exit == 1.0 (the
    // corner), so the post-loop check at line 122 returns None.
    let curve = cos_line(Point3::new(-1.0, -1.0, 0.0), Point3::new(0.0, 0.0, 0.0));
    let result = clip_curve_to_faces(&brep_a, face_a, &brep_b, face_b, &curve);
    assert!(result.is_none(), "corner-tangent curve should return None");
}

// ── clip_curve_to_faces: disjoint clipped ranges → None (line 301) ──────

/// Two faces that cover disjoint x-ranges of the same line cause
/// t_start >= t_end after intersecting the two individual clip ranges,
/// triggering the None return at line 301.
#[test]
fn clip_curve_to_faces_disjoint_ranges_returns_none() {
    // face_a: [0, 0.4]×[0,1] in z=0
    let mut brep_a = BRep::new(tol());
    let mvfs_a = brep_a.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs_a");
    let mev1 = brep_a
        .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(0.4, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep_a
        .mev(mev1.half_edge_twin, Point3::new(0.4, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep_a
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_a = brep_a
        .mef(mev3.half_edge_twin, mev1.half_edge)
        .expect("mef_a");
    let face_a = mef_a.face;

    // face_b: [0.6, 1.0]×[0,1] in z=0
    let mut brep_b = BRep::new(tol());
    let mvfs_b = brep_b.mvfs(Point3::new(0.6, 0.0, 0.0)).expect("mvfs_b");
    let mev4 = brep_b
        .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev4");
    let mev5 = brep_b
        .mev(mev4.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev5");
    let mev6 = brep_b
        .mev(mev5.half_edge_twin, Point3::new(0.6, 1.0, 0.0))
        .expect("mev6");
    let mef_b = brep_b
        .mef(mev6.half_edge_twin, mev4.half_edge)
        .expect("mef_b");
    let face_b = mef_b.face;

    // Curve from (0,0.5,0) to (1,0.5,0): passes through both faces but in
    // disjoint parameter ranges — [0,0.4] for face_a and [0.6,1.0] for face_b.
    // After intersecting: t_start=0.6 >= t_end=0.4 → line 301.
    let curve = cos_line(Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0));
    let result = clip_curve_to_faces(&brep_a, face_a, &brep_b, face_b, &curve);
    assert!(
        result.is_none(),
        "disjoint clipped ranges should return None at line 301"
    );
}

// ── face_uv_boundary tests ──────────────────────────────────────────────────

#[test]
fn face_uv_boundary_returns_pcurves_for_cube() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let face = cube.faces.keys().next().expect("has faces");
    let boundary = super::face_uv_boundary(&cube, face);
    assert!(
        boundary.is_some(),
        "cube faces should have pcurves on all half-edges"
    );
    let edges = boundary.expect("verified");
    assert_eq!(edges.len(), 4, "cube face has 4 boundary edges");
}

#[test]
fn face_uv_boundary_returns_none_without_pcurves() {
    // Build a face with no pcurves (manual Euler ops, no surface geometry).
    let mut brep = BRep::new(tol());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).expect("mef");
    let boundary = super::face_uv_boundary(&brep, mef_r.face);
    assert!(
        boundary.is_none(),
        "faces without pcurves should return None"
    );
}

// ── find_uv_crossings tests ────────────────────────────────────────────────

#[test]
fn find_uv_crossings_detects_diagonal_through_cube_face() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let face = cube.faces.keys().next().expect("has faces");

    // Get the face's UV boundary to understand the UV domain.
    let boundary = super::face_uv_boundary(&cube, face).expect("has pcurves");
    assert!(!boundary.is_empty());

    // Build a diagonal UV curve that should cross two boundary edges.
    // Cube faces have pcurves in the [0,1] UV domain typically.
    let diag_uv = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(-0.5, 0.5, 0.0), Point3::new(1.5, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };

    let crossings = super::find_uv_crossings(&cube, face, &diag_uv, &tol());
    assert!(crossings.is_some(), "should find crossings");
    let crossings = crossings.expect("verified");
    assert!(
        crossings.len() >= 2,
        "diagonal through face should have >= 2 crossings, got {}",
        crossings.len()
    );
    // Crossings should be sorted by t_curve.
    for w in crossings.windows(2) {
        assert!(
            w[0].t_curve <= w[1].t_curve,
            "crossings should be sorted by t_curve"
        );
    }
}

#[test]
fn find_uv_crossings_outside_face_returns_empty() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let face = cube.faces.keys().next().expect("has faces");

    // UV curve far outside the face's UV domain.
    let far = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(10.0, 10.0, 0.0), Point3::new(11.0, 11.0, 0.0)],
        weights: vec![1.0, 1.0],
    };

    let crossings = super::find_uv_crossings(&cube, face, &far, &tol());
    assert!(crossings.is_some());
    let crossings = crossings.expect("verified");
    assert_eq!(
        crossings.len(),
        0,
        "curve outside face should have 0 crossings"
    );
}

// ── clip_curve_to_crossings tests ───────────────────────────────────────────

#[test]
fn clip_curve_to_crossings_degree1_clips_correctly() {
    let curve = nurbs_line(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0));
    // Clip to [0.2, 0.8] in parameter space (original range [0, 1]).
    let clipped = super::clip_curve_to_crossings(&curve, 0.2, 0.8);
    assert!(clipped.is_some());
    let c = clipped.expect("verified");
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = c.parameter_range();
    let p0 = c.point_at(t0);
    let p1 = c.point_at(t1);
    assert!((p0.x - 2.0).abs() < 1e-6, "start should be at x=2");
    assert!((p1.x - 8.0).abs() < 1e-6, "end should be at x=8");
}

#[test]
fn clip_curve_to_crossings_zero_length_returns_none() {
    let curve = nurbs_line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    let result = super::clip_curve_to_crossings(&curve, 0.5, 0.5);
    assert!(result.is_none(), "zero-length clip should return None");
}

// UV-space clipping tests are in clip_uv_tests.rs (included from clip.rs)
