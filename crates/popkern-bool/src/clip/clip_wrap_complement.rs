//! A closed record's crossing hull taken as its complement when the
//! in-face span wraps through the record's parameter seam (split out of
//! `clip_uv` with PK-1190's boundary-riding samples).
//!
//! The naive crossing hull of a closed record can be the WRONG SIDE:
//! when the inside span wraps through the parameter seam, the hull is
//! exactly its complement (measured: a bore-mouth ellipse breaching a
//! wall kept only the 8% sliver beyond the wall, t[0.46, 0.54], and
//! discarded the 92% real mouth arc — the wall face's opening then
//! assembled from the sliver and failed closure by 27.1). The complement
//! is a WRAPPED window (exit beyond the range end; the extraction
//! re-anchors through the seam).
//!
//! It is taken only on UNANIMOUS multi-sample evidence — every interior
//! sample of the hull outside AND every interior sample of the
//! complement inside. A single-midpoint test misfired on flush fixtures
//! whose closed curves hug face boundaries: on-boundary containment
//! verdicts are coin flips there, and one flipped midpoint inverted a
//! correct hull (measured: 8 tuned union/chamfer gates broke). A sample
//! riding the face's boundary edge counts as on the face: a record lying
//! along that edge (a slot's round end on a plate, its arc across its
//! cylinder's seam and the chart's circle start) kept the 180 degrees
//! off the face instead of the 180 along its edge (PK-1190). A correct
//! hull hugging the boundary keeps its samples inside, so it is never
//! inverted.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

/// The wrapped complement `(t_exit, t_enter + period)` of the hull
/// `(t_enter, t_exit)` of a closed record on `face`, when the hull lies
/// off the face and the complement on it; `None` keeps the hull.
pub(super) fn wrapped_complement(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    is_closed: bool,
    (t_enter, t_exit): (f64, f64),
    trace: bool,
) -> Option<(f64, f64)> {
    let (t_min, t_max) = curve_uv.parameter_range();
    if !(is_closed && t_exit - t_enter > 1e-12 && (t_exit - t_enter) < (t_max - t_min)) {
        return None;
    }
    let period = t_max - t_min;
    // Riding means coincident with the boundary edge (the slot end's
    // arc lies on its edge to rounding), a thousandth of the clip's
    // boundary band: the band itself read the complement of two balls'
    // near-tangent circles 1e-3 into each other as riding, and swapped
    // their correct hull.
    let band = super::clip_uv_band::UvBand::of(brep, face, curve_uv).scaled(1e-3);
    let inside = |t: f64| {
        let tt = if t > t_max { t - period } else { t };
        let q = curve_uv.point_at(tt);
        crate::split_inner_loop::point_inside_face_uv(brep, face, q, &brep.tolerance)
            || super::clip_boundary_riding::point_rides_boundary(brep, face, q, band)
    };
    let samples = |a: f64, b: f64| -> [f64; 3] {
        [a + (b - a) * 0.25, a + (b - a) * 0.5, a + (b - a) * 0.75]
    };
    let hull_all_out = samples(t_enter, t_exit).iter().all(|&t| !inside(t));
    let comp_all_in = samples(t_exit, t_enter + period).iter().all(|&t| inside(t));
    if !(hull_all_out && comp_all_in) {
        return None;
    }
    if trace {
        eprintln!(
            "    [uvclip-wrap] {face:?} hull ({t_enter:.3},{t_exit:.3}) outside — complement ({t_exit:.3},{:.3})",
            t_enter + period
        );
    }
    Some((t_exit, t_enter + period))
}
