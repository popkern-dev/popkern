//! The UV clip's crossings on a face whose boundary crosses the seam
//! (split from `clip_uv_tests.rs`).

use super::*;

// ── find_uv_crossings with seam-crossing boundary ───────────────────────

/// Helper: build a quad face and assign pcurves to its half-edges.
pub(in super::super) fn build_face_with_pcurves(
    uv_corners: [(Point3, Point3); 4],
) -> (BRep, popkern_core::topology::FaceHandle) {
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);
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

    let face_data = &brep.faces[mef.face];
    let loop_data = brep.loops.get(face_data.outer_loop).expect("loop");
    let start = loop_data.half_edge;
    let mut cursor = start;
    for (i, (p0, p1)) in uv_corners.iter().enumerate() {
        let pcurve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![*p0, *p1],
            weights: vec![1.0, 1.0],
        };
        brep.half_edges[cursor].pcurve = Some(pcurve.into());
        cursor = brep.half_edges[cursor].next;
        if cursor == start && i < 3 {
            break;
        }
    }
    (brep, mef.face)
}

#[test]
fn find_uv_crossings_near_seam_finds_crossings() {
    // Boundary spans the seam: u=[5.5, 0.8]. The seam-shift dance
    // is gated on the face carrying a PERIODIC surface (the raw
    // u-value thresholds fired spuriously on planar faces whose UV
    // coordinates are plain lengths), so this fixture must attach
    // the cylinder its pcurves imply.
    //
    // The face is the strip through the seam, u from 5.5 to 0.8 + 2π: its
    // vertices and edges lie there in space (the bottom and top edges the
    // short arcs across the seam), while its horizontal pcurves are stored
    // the long way round (5.5 → 0.8), the class the window reads from the
    // edges (chart-service.md §6b step 2).
    let (mut brep, face) = build_face_with_pcurves([
        (Point3::new(5.5, 0.0, 0.0), Point3::new(0.8, 0.0, 0.0)),
        (Point3::new(0.8, 0.0, 0.0), Point3::new(0.8, 1.0, 0.0)),
        (Point3::new(0.8, 1.0, 0.0), Point3::new(5.5, 1.0, 0.0)),
        (Point3::new(5.5, 1.0, 0.0), Point3::new(5.5, 0.0, 0.0)),
    ]);
    let cylinder = popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    // Each half-edge's edge in space, from (u0, v0) to (u1, v1) across the
    // seam the short way.
    let tau = std::f64::consts::TAU;
    let walks = [
        ((5.5, 0.0), (0.8 + tau, 0.0)),
        ((0.8 + tau, 0.0), (0.8 + tau, 1.0)),
        ((0.8 + tau, 1.0), (5.5, 1.0)),
        ((5.5, 1.0), (5.5, 0.0)),
    ];
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut cursor = start;
    for (from, to) in walks {
        use popkern_core::geometry::surface::SurfaceEval as _;
        let points: Vec<Point3> = (0..=64)
            .map(|k| {
                let f = k as f64 / 64.0;
                cylinder.point_at(from.0 + (to.0 - from.0) * f, from.1 + (to.1 - from.1) * f)
            })
            .collect();
        let he = brep.half_edges[cursor].clone();
        let canonical = brep.edges[he.edge].half_edge == cursor;
        let mut ordered = points.clone();
        if !canonical {
            ordered.reverse();
        }
        let n = ordered.len();
        let mut knots = vec![0.0, 0.0];
        knots.extend((1..n - 1).map(|k| k as f64 / (n - 1) as f64));
        knots.extend([1.0, 1.0]);
        brep.edges[he.edge].curve = Some(popkern_core::geometry::curve::CurveOnEdge {
            geometry: popkern_core::geometry::curve::CurveGeometry::Nurbs(NurbsCurve {
                degree: 1,
                knots,
                control_points: ordered,
                weights: vec![1.0; n],
            }),
            t_start: 0.0,
            t_end: 1.0,
        });
        brep.vertices[he.origin].position = points[0];
        cursor = he.next;
    }
    brep.faces[face].surface = Some(popkern_core::geometry::SurfaceOnFace {
        geometry: cylinder,
        orientation: true,
    });

    // Curve near u=0.1 crossing the face in v direction.
    let curve_uv = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.1, -0.5, 0.0), Point3::new(0.1, 1.5, 0.0)],
        weights: vec![1.0, 1.0],
    };

    let crossings = super::super::find_uv_crossings(&brep, face, &curve_uv, &tol());
    assert!(crossings.is_some(), "should produce crossings");
    let crossings = crossings.expect("verified");
    assert!(
        crossings.len() >= 2,
        "curve crossing seam-spanning face should have >= 2 crossings, got {}",
        crossings.len()
    );
}
