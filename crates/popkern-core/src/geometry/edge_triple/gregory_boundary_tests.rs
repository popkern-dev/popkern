//! Tests for the Gregory-boundary dispatcher (stage 5).
//!
//! Strategy: build a synthetic `RectangularGregory` patch whose
//! `u=0` boundary lies on a known cylinder, request a triple via
//! `compute_edge_geometry` with
//! `EdgeHint::SurfaceBoundary { side: A, boundary: UMin }`, and
//! verify the G1–G6 invariants.
//!
//! The tests deliberately stay close to fillet-pipeline shapes
//! (cylinder fillet × Gregory cap) since that's the load-bearing
//! case for stage 6's user-visible fix.

use crate::geometry::edge_triple::compute_edge_geometry;
use crate::geometry::edge_triple::inputs::{BoundarySide, EdgeAnchors, EdgeHint, SurfaceSide};
use crate::geometry::edge_triple::types::{EdgeTripleProvenance, GregorySide, ParameterMap};
use crate::geometry::gregory::{GregoryPatch, RectangularGregory};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Build a synthetic 4×4 RectangularGregory whose `u=0` boundary
/// is the straight chord from `chord_a` to `chord_b` (a degenerate
/// "boundary" on a flat patch — sufficient to exercise the dispatch
/// logic without requiring a full fillet pipeline run).
fn flat_gregory_with_u0_chord(chord_a: Point3, chord_b: Point3) -> GregoryPatch {
    // u=0 boundary CPs at j=0..4: linearly interpolate a→b at
    // (0, 1/3, 2/3, 1). Place the rest of the 4×4 grid as a
    // bilinear sheet for a valid (but trivial) patch.
    let lerp = |s: f64| {
        Point3::new(
            chord_a.x + s * (chord_b.x - chord_a.x),
            chord_a.y + s * (chord_b.y - chord_a.y),
            chord_a.z + s * (chord_b.z - chord_a.z),
        )
    };
    let u0_b = [lerp(0.0), lerp(1.0 / 3.0), lerp(2.0 / 3.0), lerp(1.0)];
    // Place u=1 boundary one unit "below" (negative-z) the chord.
    let off = Vec3::new(1.0, 0.0, 0.0);
    let u1_b = [
        Point3::new(u0_b[0].x + off.x, u0_b[0].y + off.y, u0_b[0].z + off.z),
        Point3::new(u0_b[1].x + off.x, u0_b[1].y + off.y, u0_b[1].z + off.z),
        Point3::new(u0_b[2].x + off.x, u0_b[2].y + off.y, u0_b[2].z + off.z),
        Point3::new(u0_b[3].x + off.x, u0_b[3].y + off.y, u0_b[3].z + off.z),
    ];
    // Bilinear interior: row 1 = u0 + 1/3·off, row 2 = u0 + 2/3·off
    let row1 = [
        Point3::new(u0_b[0].x + off.x / 3.0, u0_b[0].y, u0_b[0].z),
        Point3::new(u0_b[1].x + off.x / 3.0, u0_b[1].y, u0_b[1].z),
        Point3::new(u0_b[2].x + off.x / 3.0, u0_b[2].y, u0_b[2].z),
        Point3::new(u0_b[3].x + off.x / 3.0, u0_b[3].y, u0_b[3].z),
    ];
    let row2 = [
        Point3::new(u0_b[0].x + off.x * 2.0 / 3.0, u0_b[0].y, u0_b[0].z),
        Point3::new(u0_b[1].x + off.x * 2.0 / 3.0, u0_b[1].y, u0_b[1].z),
        Point3::new(u0_b[2].x + off.x * 2.0 / 3.0, u0_b[2].y, u0_b[2].z),
        Point3::new(u0_b[3].x + off.x * 2.0 / 3.0, u0_b[3].y, u0_b[3].z),
    ];
    let boundary_cps = [u0_b, row1, row2, u1_b];
    // Inner CPs: bilinear (no twist). Match boundary_cps interior.
    let inner_cps = [
        (boundary_cps[1][1], boundary_cps[1][1]),
        (boundary_cps[1][2], boundary_cps[1][2]),
        (boundary_cps[2][2], boundary_cps[2][2]),
        (boundary_cps[2][1], boundary_cps[2][1]),
    ];
    GregoryPatch::Rectangular(Box::new(RectangularGregory {
        boundary_cps,
        inner_cps,
    }))
}

/// Build a synthetic Gregory whose `u=0` boundary is an *axial*
/// segment on a known cylinder. CPs are placed at constant theta
/// with varying height, so the cubic Bezier is exact on the
/// cylinder (a straight line at fixed angle, parameterized by
/// height).
fn gregory_with_u0_arc_on_cylinder(
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
    theta_a: f64,
    theta_b: f64,
    height: f64,
) -> GregoryPatch {
    use crate::geometry::analytic_surfaces::orthonormal_basis;
    let (perp1, perp2) = orthonormal_basis(cyl_axis);
    // Use the average angle for the boundary (so we get a vertical
    // line on the cylinder at fixed theta). theta_a, theta_b
    // describe the parametric span we want — pass them in only for
    // documentation; the actual boundary is axial.
    let theta = theta_a + 0.5 * (theta_b - theta_a);
    let radial = perp1 * (cyl_radius * theta.cos()) + perp2 * (cyl_radius * theta.sin());
    let axial_pt = |s: f64| {
        let h = height + s * 0.5; // segment of height [h, h+0.5]
        Point3::new(
            cyl_origin.x + radial.x + cyl_axis.x * h,
            cyl_origin.y + radial.y + cyl_axis.y * h,
            cyl_origin.z + radial.z + cyl_axis.z * h,
        )
    };
    let u0 = [
        axial_pt(0.0),
        axial_pt(1.0 / 3.0),
        axial_pt(2.0 / 3.0),
        axial_pt(1.0),
    ];
    let off = cyl_axis * 0.5;
    let u1 = [
        Point3::new(u0[0].x + off.x, u0[0].y + off.y, u0[0].z + off.z),
        Point3::new(u0[1].x + off.x, u0[1].y + off.y, u0[1].z + off.z),
        Point3::new(u0[2].x + off.x, u0[2].y + off.y, u0[2].z + off.z),
        Point3::new(u0[3].x + off.x, u0[3].y + off.y, u0[3].z + off.z),
    ];
    let row1 = [
        Point3::new(
            u0[0].x + off.x / 3.0,
            u0[0].y + off.y / 3.0,
            u0[0].z + off.z / 3.0,
        ),
        Point3::new(
            u0[1].x + off.x / 3.0,
            u0[1].y + off.y / 3.0,
            u0[1].z + off.z / 3.0,
        ),
        Point3::new(
            u0[2].x + off.x / 3.0,
            u0[2].y + off.y / 3.0,
            u0[2].z + off.z / 3.0,
        ),
        Point3::new(
            u0[3].x + off.x / 3.0,
            u0[3].y + off.y / 3.0,
            u0[3].z + off.z / 3.0,
        ),
    ];
    let row2 = [
        Point3::new(
            u0[0].x + off.x * 2.0 / 3.0,
            u0[0].y + off.y * 2.0 / 3.0,
            u0[0].z + off.z * 2.0 / 3.0,
        ),
        Point3::new(
            u0[1].x + off.x * 2.0 / 3.0,
            u0[1].y + off.y * 2.0 / 3.0,
            u0[1].z + off.z * 2.0 / 3.0,
        ),
        Point3::new(
            u0[2].x + off.x * 2.0 / 3.0,
            u0[2].y + off.y * 2.0 / 3.0,
            u0[2].z + off.z * 2.0 / 3.0,
        ),
        Point3::new(
            u0[3].x + off.x * 2.0 / 3.0,
            u0[3].y + off.y * 2.0 / 3.0,
            u0[3].z + off.z * 2.0 / 3.0,
        ),
    ];
    let boundary_cps = [u0, row1, row2, u1];
    let inner_cps = [
        (boundary_cps[1][1], boundary_cps[1][1]),
        (boundary_cps[1][2], boundary_cps[1][2]),
        (boundary_cps[2][2], boundary_cps[2][2]),
        (boundary_cps[2][1], boundary_cps[2][1]),
    ];
    GregoryPatch::Rectangular(Box::new(RectangularGregory {
        boundary_cps,
        inner_cps,
    }))
}

/// Gregory × Cylinder: the load-bearing case for the "messy top"
/// fix. Set up a Gregory whose u=0 boundary is a Bezier arc with
/// CPs on a known cylinder — note that a cubic Bezier through 4
/// points on a circle is a polynomial *approximation* of the
/// circular arc, so we don't expect bit-exact reconstruction; the
/// test verifies G2 holds at the dispatcher's stated tolerance.
#[test]
fn gregory_x_cylinder_umin_constructs_triple() {
    let cyl_origin = Point3::new(0.0, 0.0, 0.0);
    let cyl_axis = Vec3::new(0.0, 0.0, 1.0);
    let cyl_radius = 1.0;
    let theta_a = 0.0;
    let theta_b = std::f64::consts::FRAC_PI_4; // 45°
    let height = 0.5;
    let greg =
        gregory_with_u0_arc_on_cylinder(cyl_origin, cyl_axis, cyl_radius, theta_a, theta_b, height);
    let surf_a = SurfaceGeometry::Gregory(greg);
    let surf_b = SurfaceGeometry::Cylinder {
        origin: cyl_origin,
        axis: cyl_axis,
        radius: cyl_radius,
    };
    // Anchors = the cubic Bezier's endpoint CPs (Bezier through CPs
    // hits the first/last CPs at t=0/t=1 by construction).
    let (anchor_a, anchor_b) =
        if let SurfaceGeometry::Gregory(GregoryPatch::Rectangular(r)) = &surf_a {
            (r.boundary_cps[0][0], r.boundary_cps[0][3])
        } else {
            panic!("test setup: expected Rectangular Gregory")
        };
    let anchors = EdgeAnchors::Open {
        start: anchor_a,
        end: anchor_b,
    };
    let hint = EdgeHint::SurfaceBoundary {
        side: SurfaceSide::A,
        boundary: BoundarySide::UMin,
    };
    let opts = crate::geometry::edge_triple::EdgeTripleOptions::default();
    let triple = compute_edge_geometry(
        &surf_a,
        &surf_b,
        anchors,
        Some(hint),
        &Tolerance::default(),
        &opts,
    )
    .expect("dispatch should succeed");
    use crate::geometry::curve::CurveEval;
    use crate::geometry::pcurve::Pcurve;
    let n_samples = 16;
    for i in 0..n_samples {
        let s = (i as f64) / ((n_samples - 1) as f64);
        let t = triple.t_start + s * (triple.t_end - triple.t_start);
        let p3 = triple.curve_3d.point_at(t);
        let uv_a = triple.pcurve_a.point_at(t);
        let p_a = surf_a.point_at(uv_a.x, uv_a.y);
        assert!(
            p_a.distance_to(&p3) < 1e-9,
            "G2 (greg) violated at s={s}: surf_a(pcurve_a) = {p_a:?}, curve_3d = {p3:?}, gap = {}",
            p_a.distance_to(&p3)
        );
        let uv_b = triple.pcurve_b.point_at(t);
        let p_b = surf_b.point_at(uv_b.x, uv_b.y);
        // Cylinder side is fitted, so allow a slightly larger
        // tolerance — but well below user-visible (1e-4 vs 1e-1
        // baseline measured in fillet_zoo).
        assert!(
            p_b.distance_to(&p3) < 1e-4,
            "G2 (cyl) violated at s={s}: surf_b(pcurve_b) = {p_b:?}, curve_3d = {p3:?}, gap = {}",
            p_b.distance_to(&p3)
        );
        let _: &Pcurve = &triple.pcurve_b;
    }
}

/// Plane × Gregory: dispatch returns an EdgeTriple whose `curve_3d`
/// is a degree-3 NurbsCurve interpolating the Gregory's `u=0` chord
/// (exact — this is a flat boundary on a flat patch). The plane-side
/// pcurve is a degree-3 fit through 32 projected samples.
#[test]
fn gregory_x_plane_umin_constructs_triple() {
    let chord_a = Point3::new(0.0, 0.0, 0.0);
    let chord_b = Point3::new(1.0, 0.0, 0.0);
    let greg = flat_gregory_with_u0_chord(chord_a, chord_b);
    let surf_a = SurfaceGeometry::Gregory(greg);
    // Plane z=0 (chord lies on it).
    let surf_b = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let anchors = EdgeAnchors::Open {
        start: chord_a,
        end: chord_b,
    };
    let hint = EdgeHint::SurfaceBoundary {
        side: SurfaceSide::A,
        boundary: BoundarySide::UMin,
    };
    let opts = crate::geometry::edge_triple::EdgeTripleOptions::default();
    let triple = compute_edge_geometry(
        &surf_a,
        &surf_b,
        anchors,
        Some(hint),
        &Tolerance::default(),
        &opts,
    )
    .expect("dispatch should succeed");
    // Provenance: GregoryBoundary on side A.
    match triple.provenance {
        EdgeTripleProvenance::GregoryBoundary { side, .. } => {
            assert_eq!(side, GregorySide::A);
        }
        ref other => panic!("expected GregoryBoundary provenance, got {other:?}"),
    }
    // ParameterMap: Identity.
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    // Anchor CP equality: curve_3d's first/last CPs should match the
    // anchors at machine precision.
    use crate::geometry::curve::CurveEval;
    let p_t0 = triple.curve_3d.point_at(triple.t_start);
    let p_t1 = triple.curve_3d.point_at(triple.t_end);
    assert!(
        p_t0.distance_to(&chord_a) < 1e-12,
        "curve_3d(t_start) = {p_t0:?}, expected {chord_a:?}"
    );
    assert!(
        p_t1.distance_to(&chord_b) < 1e-12,
        "curve_3d(t_end) = {p_t1:?}, expected {chord_b:?}"
    );
    // G2: surf_a(pcurve_a(t)) ≡ curve_3d(t) AND surf_b(pcurve_b(t))
    // ≡ curve_3d(t) at sampled points.
    use crate::geometry::pcurve::Pcurve;
    let n_samples = 16;
    for i in 0..n_samples {
        let s = (i as f64) / ((n_samples - 1) as f64);
        let t = triple.t_start + s * (triple.t_end - triple.t_start);
        let p3 = triple.curve_3d.point_at(t);
        let uv_a = triple.pcurve_a.point_at(t);
        let p_a = surf_a.point_at(uv_a.x, uv_a.y);
        assert!(
            p_a.distance_to(&p3) < 1e-9,
            "G2 (a) violated at s={s}: surf_a(pcurve_a) = {p_a:?}, curve_3d = {p3:?}"
        );
        let uv_b = triple.pcurve_b.point_at(t);
        let p_b = surf_b.point_at(uv_b.x, uv_b.y);
        // Plane is fit via samples; tolerance loosened to 1e-6.
        assert!(
            p_b.distance_to(&p3) < 1e-6,
            "G2 (b) violated at s={s}: surf_b(pcurve_b) = {p_b:?}, curve_3d = {p3:?}"
        );
        // unused, but make sure pcurves typecheck
        let _: &Pcurve = &triple.pcurve_a;
        let _: &Pcurve = &triple.pcurve_b;
    }
}

#[test]
fn gregory_boundary_dispatch_trims_to_anchor_subspan() {
    let chord_a = Point3::new(0.0, 0.0, 0.0);
    let chord_b = Point3::new(1.0, 0.0, 0.0);
    let surf_a = SurfaceGeometry::Gregory(flat_gregory_with_u0_chord(chord_a, chord_b));
    let surf_b = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let start = surf_a.point_at(0.0, 0.25);
    let end = surf_a.point_at(0.0, 0.75);
    let triple = compute_edge_geometry(
        &surf_a,
        &surf_b,
        EdgeAnchors::Open { start, end },
        Some(EdgeHint::SurfaceBoundary {
            side: SurfaceSide::A,
            boundary: BoundarySide::UMin,
        }),
        &Tolerance::default(),
        &crate::geometry::edge_triple::EdgeTripleOptions::default(),
    )
    .expect("boundary-owned construction should trim to the anchor subspan");

    use crate::geometry::curve::CurveEval;
    let repaired_start = triple.curve_3d.point_at(triple.t_start);
    let repaired_end = triple.curve_3d.point_at(triple.t_end);
    assert!(repaired_start.distance_to(&start) < 1.0e-7);
    assert!(repaired_end.distance_to(&end) < 1.0e-7);
    let uv_start = triple.pcurve_a.point_at(triple.t_start);
    let uv_end = triple.pcurve_a.point_at(triple.t_end);
    assert!(uv_start.x.abs() < 1.0e-9);
    assert!(uv_end.x.abs() < 1.0e-9);
    assert!((uv_start.y - 0.25).abs() < 1.0e-7);
    assert!((uv_end.y - 0.75).abs() < 1.0e-7);
}
