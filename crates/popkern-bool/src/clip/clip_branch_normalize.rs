//! Period-branch normalization for the UV clip query (split from
//! `clip_uv_crossings.rs` under the 500-line file cap).

use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

/// The v-branch counterpart for a v-BOUNDED face of a v-periodic chart
/// (a torus band): translate a record's pcurve rigidly by whole
/// v-periods so its v-centre lands nearest the boundary's (PK-775).
/// Applied once, before the clip, so crossings, containment and the
/// extracted pieces all see one branch. The torus×plane
/// marcher authored a cube-wall record on v ∈ [5.24, 7.33], the 2π
/// branch of a band on v ∈ [0, π/2]; with u alone normalized no crossing
/// was found and the record dropped as Outside. A face wrapping v keeps
/// its record as authored.
pub(crate) fn normalize_v_branch(brep: &BRep, face: FaceHandle, normalized: &mut NurbsCurve) {
    let Some(period) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .and_then(|s| {
            popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(&s.geometry).1
        })
        .filter(|p| *p > 0.0)
    else {
        return;
    };
    let Some(boundary) = crate::clip::face_uv_boundary(brep, face) else {
        return;
    };
    let (lo, hi, sum, n) = boundary
        .iter()
        .flat_map(|(_, pc)| {
            let (t0, t1) = pc.parameter_range();
            (0..=8).map(move |i| pc.point_at(t0 + (t1 - t0) * f64::from(i) / 8.0).y)
        })
        .fold(
            (f64::INFINITY, f64::NEG_INFINITY, 0.0, 0usize),
            |(lo, hi, s, n), v| (lo.min(v), hi.max(v), s + v, n + 1),
        );
    if n == 0 || hi - lo >= period * (1.0 - 1e-6) {
        return;
    }
    let q = &normalized.control_points;
    let q_center = q.iter().map(|p| p.y).sum::<f64>() / q.len().max(1) as f64;
    let k = ((sum / n as f64 - q_center) / period).round();
    if k != 0.0 {
        for p in &mut normalized.control_points {
            p.y += k * period;
        }
        if popkern_core::trace::trace_enabled("clip") {
            eprintln!("    [uvclip-branch] {face:?} query v shifted by {k:+.0} period(s)");
        }
    }
}

/// A split containment vote on a record with no boundary crossings —
/// some samples inside the face, some outside — means the crossing
/// search missed real crossings: the query lies on another period of
/// a single-period chart, wholly OFF the boundary's u-extent (touching
/// it at most at the seam), where no boundary pcurve can meet it.
/// The query landed rigidly (no per-point fold) on the period it
/// overlaps most, for the caller to classify again; `None` when it is
/// not off-branch or the chart is not periodic. A rigid shift keeps
/// the curve parameter, so windows found on the landed query are the
/// original's. Measured (PK-732): the x0 = −20 equal-radius tee's
/// record from the singular point S+ — on the branch cylinder's seam —
/// ran u ∈ [−π, 0]; the 5-sample vote read 2 inside and the record
/// dropped whole although both its ends lie in the face (its mirror
/// twin at u ∈ [0, π] kept its two end windows). Landing every
/// off-branch query instead re-classified a unanimous-outside record
/// whose ends merely touch the boundary (the x0 = 0 tee's x < 0
/// half-ellipse) as inside.
pub(super) fn land_split_vote_query(
    brep: &BRep,
    face: FaceHandle,
    query: &NurbsCurve,
) -> Option<NurbsCurve> {
    let period = popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(
        &brep.faces.get(face)?.surface.as_ref()?.geometry,
    )
    .0?;
    let extent = |points: &mut dyn Iterator<Item = f64>| {
        points.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
            (lo.min(u), hi.max(u))
        })
    };
    let ring = crate::clip::face_uv_boundary(brep, face)?;
    let (b_lo, b_hi) = extent(&mut ring.iter().flat_map(|(_, pc)| {
        let (a, b) = pc.parameter_range();
        (0..=16).map(move |i| pc.point_at(a + (b - a) * i as f64 / 16.0).x)
    }));
    let (q_lo, q_hi) = extent(&mut query.control_points.iter().map(|p| p.x));
    if !(q_lo.is_finite() && q_hi.is_finite() && b_lo.is_finite() && b_hi.is_finite())
        || b_hi - b_lo > period * (1.0 + 1e-9)
    {
        return None;
    }
    let slack = 1e-9 * period;
    let overlap = |k: f64| (q_hi + k * period).min(b_hi) - (q_lo + k * period).max(b_lo);
    if overlap(0.0) > slack {
        return None;
    }
    let mut best = (0.0, overlap(0.0));
    let mut k = ((b_lo - q_hi) / period).floor();
    while k <= ((b_hi - q_lo) / period).ceil() {
        if overlap(k) > best.1 + slack {
            best = (k, overlap(k));
        }
        k += 1.0;
    }
    if best.0 == 0.0 || best.1 <= slack {
        return None;
    }
    if popkern_core::trace::trace_enabled("clip") {
        eprintln!(
            "    [uvclip-branch] {face:?} split vote: off-branch query landed by {:+.0} period(s)",
            best.0
        );
    }
    let mut landed = query.clone();
    for p in &mut landed.control_points {
        p.x += best.0 * period;
    }
    Some(landed)
}

/// Whether the uv point `q` lies on the box `lo..hi` (widened by
/// `band`) once shifted by whole periods of `face`'s periodic chart
/// directions. A record traced on the branch one period below the
/// face's chart (u in [-pi, 0] on a [0, 2pi] cylinder face) failed the
/// literal test, and the clip kept a hull lying outside the face
/// (PK-765: an equal-radius tee whose branch seam runs through the
/// singular point served a phantom half-ellipse at v = -40, then
/// declined on its curveless edge).
pub(super) fn on_periodic_box(
    brep: &BRep,
    face: FaceHandle,
    q: popkern_core::geometry::Point3,
    lo: (f64, f64),
    hi: (f64, f64),
    (band_u, band_v): (f64, f64),
) -> bool {
    let (u_period, v_period) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map_or((None, None), |s| {
            popkern_core::geometry::surface::SurfaceEval::parameter_periodicity(&s.geometry)
        });
    let in_span = |x: f64, lo: f64, hi: f64, period: Option<f64>, band: f64| match period {
        Some(p) if p > 0.0 => x + ((lo - band - x) / p).ceil() * p <= hi + band,
        _ => x >= lo - band && x <= hi + band,
    };
    in_span(q.x, lo.0, hi.0, u_period, band_u) && in_span(q.y, lo.1, hi.1, v_period, band_v)
}
