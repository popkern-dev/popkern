//! Validity gates for the record-side witness (PK-494) — split from
//! `record_side_witness.rs` under the 500-line budget: the
//! periodic-honest chart window and the segment obstruction test.
//! Rationale and the measured failure each gate answers live on the
//! code, verbatim from the parent.

use popkern_core::geometry::intersection::surface_curve::intersect_surface_curve_dense;
use popkern_core::geometry::intersection::{SurfaceCurveIntersection, surface_subdivision::Aabb};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// The face's chart window from its boundary samples, PERIODIC-HONEST:
/// aperiodic axes take the plain min/max box; periodic axes take the
/// minimal covering arc of the sample angles (unwrapped around their
/// circular mean). `None` when a periodic axis has no arc smaller
/// than 90% of its period — such a window says nothing and every
/// consumer must fail CLOSED (the naive box on a seam-straddling cone
/// accepted off-trim projections and flipped a chamfered-nut serve).
pub(super) fn chart_window(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    positions: &[Point3],
) -> Option<((f64, f64), (f64, f64))> {
    let ((bu0, bu1), (bv0, bv1)) =
        popkern_core::geometry::uv_projection::exact_bounded_parameter_range(positions, surface);
    let (pu, pv) = surface.periodicity();
    let axis = |lo: f64, hi: f64, period: Option<f64>| -> Option<(f64, f64)> {
        let Some(p) = period.filter(|p| p.is_finite() && *p > 0.0) else {
            return Some((lo.min(hi), lo.max(hi)));
        };
        if (hi - lo).abs() <= 0.9 * p {
            return Some((lo.min(hi), lo.max(hi)));
        }
        None
    };
    Some((axis(bu0, bu1, pu)?, axis(bv0, bv1, pv)?))
}

/// True when the OPEN segment `s → q` crosses any face of `brep`
/// strictly before `q` — or when any hit is too ambiguous to rule
/// out (conservative: ambiguity counts as obstruction).
pub(super) fn segment_is_obstructed(
    s: &Point3,
    q: &Point3,
    brep: &BRep,
    aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> bool {
    let dir = Vec3::new(q.x - s.x, q.y - s.y, q.z - s.z);
    let len = dir.length();
    if len < tolerance.numerical_floor() {
        return true;
    }
    let ray = popkern_core::geometry::curve::CurveGeometry::Line {
        origin: *s,
        direction: dir,
    };
    // Segment box, tolerance-inflated, for the face cull.
    let m = tolerance.linear() * 10.0;
    let seg_box = Aabb {
        min: Point3::new(s.x.min(q.x) - m, s.y.min(q.y) - m, s.z.min(q.z) - m),
        max: Point3::new(s.x.max(q.x) + m, s.y.max(q.y) + m, s.z.max(q.z) + m),
    };
    // Parameter bands: ignore hits at the very start (the sample can
    // legitimately sit near its OWN operand geometry, not this brep,
    // so the start band is pure tolerance) and hits that ARE q (the
    // record's wall, by position, not by parameter).
    let t_eps = (tolerance.linear() / len).min(0.4);
    let q_gate = tolerance.linear() * 10.0;
    for (fh, fbox) in aabbs {
        if !fbox.overlaps(&seg_box, m) {
            continue;
        }
        let Some(sof) = brep.faces.get(*fh).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        // Face chart window from its own boundary samples (native).
        let positions = brep.face_outer_loop_sample_positions(*fh, 8);
        if positions.len() < 3 {
            return true; // unqueryable face in range — abstain
        }
        let Some((wu, wv)) = chart_window(&sof.geometry, &positions) else {
            return true; // seam-straddling window — fail closed (obstructed)
        };
        let pad_u = 0.05 * (wu.1 - wu.0).abs().max(1e-9);
        let pad_v = 0.05 * (wv.1 - wv.0).abs().max(1e-9);
        let window_u = (wu.0 - pad_u, wu.1 + pad_u);
        let window_v = (wv.0 - pad_v, wv.1 + pad_v);
        let Ok(result) = intersect_surface_curve_dense(
            &sof.geometry,
            window_u,
            window_v,
            &ray,
            (0.0, 1.0),
            tolerance,
        ) else {
            return true; // solver failure in range — abstain
        };
        for isect in &result.intersections {
            match isect {
                SurfaceCurveIntersection::Point { t, .. } => {
                    if *t <= t_eps || *t >= 1.0 {
                        continue;
                    }
                    let hp = Point3::new(s.x + dir.x * *t, s.y + dir.y * *t, s.z + dir.z * *t);
                    if hp.distance_to(q) <= q_gate {
                        continue; // the record's own wall
                    }
                    // In-box check keeps basis-extension phantoms out;
                    // anything surviving is an obstruction.
                    if hp.x >= fbox.min.x - m
                        && hp.x <= fbox.max.x + m
                        && hp.y >= fbox.min.y - m
                        && hp.y <= fbox.max.y + m
                        && hp.z >= fbox.min.z - m
                        && hp.z <= fbox.max.z + m
                    {
                        return true;
                    }
                }
                SurfaceCurveIntersection::Overlap { .. } => return true,
            }
        }
    }
    false
}
