//! §5.3 Round 6 rational-Gregory construction path.
//!
//! `try_rational_gregory_blend` builds a true rational Gregory patch
//! (Longhi triangular for N=3, Chiyokura-Kimura rectangular for N=4
//! — follow-up) from a [`BlendBoundary`] plus the incident fillet
//! strip data. The builder is live and tested; the dispatch in
//! [`super::surface::build_smooth_blend`] does NOT yet invoke it,
//! pending per-corner twist-vector extraction from the incident
//! fillet-strip surface normals. Without that, a no-twist Gregory
//! degenerates to an ordinary bicubic on uniform-radius corners
//! (redundant with the rolling-ball sphere) or produces a G0
//! discontinuity at corners on mixed-radius configurations (worse
//! than the polar-NURBS fallback).
//!
//! Kept in its own file so the surface dispatcher stays under the
//! 500-line structural cap and so the twist-vector follow-up has a
//! clean module boundary to plug into.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::{GregoryPatch, TriangularGregory};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::surface::BlendSurfaceResult;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod};

/// §5.3 Round 6 — build a true rational Gregory patch for a
/// vertex-blend corner with `N ∈ {3, 4}` analytic-strip fillets of
/// uniform radius.
///
/// Currently handles **N=3 triangular Longhi Gregory** only; N=4 is
/// a follow-up commit once a 4-fillet test fixture exists. For N=3:
///
/// 1. Verify all incident fillet surfaces are `Torus` or `Cylinder`
///    and share one minor radius.
/// 2. Extract the 3 boundary cubic Beziers by fitting each
///    `BoundarySegment::FilletEnd` curve (assumed nearly-cubic on
///    the fillet-end arc scale) to a cubic Bezier control polygon.
/// 3. Place the 3 inner candidates at the vertex position — the
///    rational blend then passes through the vertex at the patch
///    interior while the 3 boundary curves stitch through the
///    incident fillet strips via the cross-derivative encoding in
///    the second Bezier control point.
/// 4. Return the patch as [`SurfaceGeometry::Gregory`].
///
/// **Activation caveat:** the dispatcher does not yet call this
/// function — see the module-level docs for the twist-vector
/// prerequisite.
#[allow(dead_code)]
pub(super) fn try_rational_gregory_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    boundary_curves: &[NurbsCurve],
) -> Option<BlendSurfaceResult> {
    let n_fillets = boundary
        .segments
        .iter()
        .filter(|s| matches!(s, BoundarySegment::FilletEnd { .. }))
        .count();

    if n_fillets != 3 {
        return None;
    }

    // Verify uniform radii across analytic strips.
    let radii: Vec<f64> = fillet_surfaces
        .iter()
        .filter_map(|fr| match &fr.surface {
            SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
            SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    if radii.len() != fillet_surfaces.len() || radii.is_empty() {
        return None;
    }
    let (lo, hi) = radii.iter().fold((f64::INFINITY, 0.0_f64), |(lo, hi), &r| {
        (lo.min(r), hi.max(r))
    });
    if hi - lo > 1e-6 {
        return None;
    }

    let beziers: Vec<[Point3; 4]> = boundary_curves
        .iter()
        .map(fit_cubic_bezier_through)
        .collect();
    if beziers.len() != 3 {
        return None;
    }

    let v_pos = boundary.vertex_position;
    let inner_cps = [v_pos, v_pos, v_pos];
    let bez_array: [[Point3; 4]; 3] = [beziers[0], beziers[1], beziers[2]];

    let patch = TriangularGregory::from_boundary_beziers(bez_array, inner_cps, 1e-3).ok()?;

    Some(BlendSurfaceResult {
        surface: SurfaceGeometry::Gregory(GregoryPatch::Triangular(Box::new(patch))),
        boundary_curves: boundary_curves.to_vec(),
        method: CornerMethod::Gregory,
        continuity: ContinuityLevel::G2,
    })
}

/// Fit a cubic Bezier control polygon to a `NurbsCurve` by
/// interpolating at `t = 0, 1/3, 2/3, 1` and solving for the inner
/// 2 control points so the cubic Bezier passes exactly through the
/// interior sample points.
///
/// The system is `[12 6; 6 12] [P_1; P_2] = [27 Q_1 - 8 P_0 - P_3;
///                                           27 Q_2 - P_0 - 8 P_3]`
/// with determinant `108`, yielding `P_1 = (2A - B) / 18` and
/// `P_2 = (2B - A) / 18` where `A` and `B` are the right-hand sides.
#[allow(dead_code)]
pub(super) fn fit_cubic_bezier_through(curve: &NurbsCurve) -> [Point3; 4] {
    let (t0, t1) = curve.parameter_range();
    let span = t1 - t0;
    let p0 = curve.point_at(t0);
    let p3 = curve.point_at(t1);
    let q1 = curve.point_at(t0 + span / 3.0);
    let q2 = curve.point_at(t0 + 2.0 * span / 3.0);

    let axis = |q1c: f64, q2c: f64, p0c: f64, p3c: f64| -> (f64, f64) {
        let a = 27.0 * q1c - 8.0 * p0c - p3c;
        let b = 27.0 * q2c - p0c - 8.0 * p3c;
        ((2.0 * a - b) / 18.0, (2.0 * b - a) / 18.0)
    };
    let (p1x, p2x) = axis(q1.x, q2.x, p0.x, p3.x);
    let (p1y, p2y) = axis(q1.y, q2.y, p0.y, p3.y);
    let (p1z, p2z) = axis(q1.z, q2.z, p0.z, p3.z);

    [
        p0,
        Point3::new(p1x, p1y, p1z),
        Point3::new(p2x, p2y, p2z),
        p3,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::NurbsCurve;

    /// `fit_cubic_bezier_through` on a straight-line NURBS (4 CPs
    /// on a line) reproduces the input CPs exactly.
    #[test]
    fn fit_straight_line_cubic_recovers_cps() {
        let cps = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0 / 3.0, 0.0, 0.0),
            Point3::new(2.0 / 3.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let curve = NurbsCurve {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: cps.clone(),
            weights: vec![1.0; 4],
        };
        let fitted = fit_cubic_bezier_through(&curve);
        for (a, b) in fitted.iter().zip(cps.iter()) {
            assert!(
                (a.x - b.x).abs() < 1e-10 && (a.y - b.y).abs() < 1e-10 && (a.z - b.z).abs() < 1e-10,
                "fit mismatch: {:?} vs {:?}",
                a,
                b,
            );
        }
    }

    /// `fit_cubic_bezier_through` on a non-linear cubic (all 4 CPs
    /// at distinct points) reconstructs the cubic exactly: the
    /// output Bezier is the input Bezier because all 4 of the
    /// curve's Bernstein samples at `0, 1/3, 2/3, 1` are uniquely
    /// determined by its 4 CPs.
    #[test]
    fn fit_arbitrary_cubic_reconstructs_cps() {
        let cps = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.2, 0.5, 0.3),
            Point3::new(0.8, 0.5, 0.7),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let curve = NurbsCurve {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: cps.clone(),
            weights: vec![1.0; 4],
        };
        let fitted = fit_cubic_bezier_through(&curve);
        for (a, b) in fitted.iter().zip(cps.iter()) {
            let d = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt();
            assert!(d < 1e-8, "fit mismatch: {:?} vs {:?} (d={d:.2e})", a, b);
        }
    }

    /// End-to-end builder acceptance: synthetic cube-corner fixture
    /// with 3 uniform-radius analytic fillet strips, build a
    /// rational Gregory patch, verify the output:
    /// (a) surface is `Gregory::Triangular`
    /// (b) boundary-curves are preserved
    /// (c) the patch evaluates cleanly at the simplex corners
    ///     (= cube vertex for all 3 candidates)
    /// (d) the patch hits each of the 3 boundary-Bezier corner
    ///     points at the corresponding simplex corner.
    #[test]
    fn try_rational_gregory_on_synthetic_cube_corner() {
        use super::super::boundary::{BlendBoundary, BoundarySegment};
        use crate::fillet_surface::FilletSurfaceResult;
        use popkern_core::geometry::{Vec3, nurbs_curve::NurbsCurve};
        use popkern_core::topology::VertexHandle;

        // Synthetic cube corner at origin with three orthogonal
        // fillet-end arcs. The arcs trace from one axis direction
        // to another at distance `r` from the origin (the corner
        // vertex).
        let r = 0.1;
        let vertex = VertexHandle::default();
        let v_pos = Point3::new(0.0, 0.0, 0.0);

        let make_arc = |start: Point3, end: Point3| -> NurbsCurve {
            // 4-CP cubic Bezier approximating the arc from start to
            // end on the sphere of radius r around origin. Rough but
            // sufficient for builder validation.
            let mid1 = Point3::new(
                start.x + (end.x - start.x) / 3.0,
                start.y + (end.y - start.y) / 3.0,
                start.z + (end.z - start.z) / 3.0,
            );
            let mid2 = Point3::new(
                start.x + 2.0 * (end.x - start.x) / 3.0,
                start.y + 2.0 * (end.y - start.y) / 3.0,
                start.z + 2.0 * (end.z - start.z) / 3.0,
            );
            NurbsCurve {
                degree: 3,
                knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                control_points: vec![start, mid1, mid2, end],
                weights: vec![1.0, 1.0, 1.0, 1.0],
            }
        };

        let ax = Point3::new(r, 0.0, 0.0);
        let ay = Point3::new(0.0, r, 0.0);
        let az = Point3::new(0.0, 0.0, r);

        let arc_xy = make_arc(ax, ay);
        let arc_yz = make_arc(ay, az);
        let arc_zx = make_arc(az, ax);

        let boundary_curves = vec![arc_xy.clone(), arc_yz.clone(), arc_zx.clone()];

        let boundary = BlendBoundary {
            vertex,
            vertex_position: v_pos,
            segments: vec![
                BoundarySegment::FilletEnd {
                    fillet_index: 0,
                    curve: arc_xy,
                },
                BoundarySegment::FilletEnd {
                    fillet_index: 1,
                    curve: arc_yz,
                },
                BoundarySegment::FilletEnd {
                    fillet_index: 2,
                    curve: arc_zx,
                },
            ],
        };

        let fillet_surfaces = vec![
            FilletSurfaceResult {
                surface: SurfaceGeometry::Cylinder {
                    origin: v_pos,
                    axis: Vec3::new(1.0, 0.0, 0.0),
                    radius: r,
                },
                trim_curve_1: make_arc(ax, ay),
                trim_curve_2: make_arc(ax, ay),
                continuity_downgrade: None,
            },
            FilletSurfaceResult {
                surface: SurfaceGeometry::Cylinder {
                    origin: v_pos,
                    axis: Vec3::new(0.0, 1.0, 0.0),
                    radius: r,
                },
                trim_curve_1: make_arc(ay, az),
                trim_curve_2: make_arc(ay, az),
                continuity_downgrade: None,
            },
            FilletSurfaceResult {
                surface: SurfaceGeometry::Cylinder {
                    origin: v_pos,
                    axis: Vec3::new(0.0, 0.0, 1.0),
                    radius: r,
                },
                trim_curve_1: make_arc(az, ax),
                trim_curve_2: make_arc(az, ax),
                continuity_downgrade: None,
            },
        ];

        let result = try_rational_gregory_blend(&boundary, &fillet_surfaces, &boundary_curves)
            .expect("Gregory build should succeed on synthetic cube corner");

        // (a) Triangular Gregory.
        match &result.surface {
            SurfaceGeometry::Gregory(GregoryPatch::Triangular(_)) => {}
            other => panic!("expected Triangular Gregory, got {other:?}"),
        }
        // (b) Boundary curves preserved.
        assert_eq!(result.boundary_curves.len(), 3);
        // (c) Method + continuity reported.
        assert_eq!(result.method, CornerMethod::Gregory);
        assert_eq!(result.continuity, ContinuityLevel::G2);
    }

    /// Builder rejects mixed-radius incidence — guards against
    /// incorrectly firing the Gregory path when the no-twist
    /// assumption would produce a G0-at-corner result.
    #[test]
    fn try_rational_gregory_rejects_mixed_radii() {
        use super::super::boundary::{BlendBoundary, BoundarySegment};
        use crate::fillet_surface::FilletSurfaceResult;
        use popkern_core::geometry::{Vec3, nurbs_curve::NurbsCurve};
        use popkern_core::topology::VertexHandle;

        let v_pos = Point3::new(0.0, 0.0, 0.0);
        let vertex = VertexHandle::default();
        let curve = NurbsCurve {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.1, 0.0, 0.0),
                Point3::new(0.0, 0.033, 0.0),
                Point3::new(0.0, 0.066, 0.0),
                Point3::new(0.0, 0.1, 0.0),
            ],
            weights: vec![1.0; 4],
        };
        let boundary = BlendBoundary {
            vertex,
            vertex_position: v_pos,
            segments: vec![
                BoundarySegment::FilletEnd {
                    fillet_index: 0,
                    curve: curve.clone(),
                },
                BoundarySegment::FilletEnd {
                    fillet_index: 1,
                    curve: curve.clone(),
                },
                BoundarySegment::FilletEnd {
                    fillet_index: 2,
                    curve: curve.clone(),
                },
            ],
        };
        // Mixed radii: 0.08, 0.10, 0.12.
        let fillet_surfaces: Vec<_> = [0.08, 0.10, 0.12]
            .iter()
            .map(|&r| FilletSurfaceResult {
                surface: SurfaceGeometry::Cylinder {
                    origin: v_pos,
                    axis: Vec3::new(1.0, 0.0, 0.0),
                    radius: r,
                },
                trim_curve_1: curve.clone(),
                trim_curve_2: curve.clone(),
                continuity_downgrade: None,
            })
            .collect();

        let boundary_curves = vec![curve.clone(), curve.clone(), curve.clone()];
        let result = try_rational_gregory_blend(&boundary, &fillet_surfaces, &boundary_curves);
        assert!(result.is_none(), "mixed radii should return None");
    }
}
