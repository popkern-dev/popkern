//! Profile analysis for sweep operations.
//!
//! Detects tangent discontinuities, classifies open/closed profiles,
//! and centers profiles at the origin for frame-based sweeping.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_math::PortableMath;

use crate::error::SweepError;
use crate::types::SweepProfile;

/// Angle threshold (radians) for detecting tangent discontinuities.
/// Approximately 15 degrees.
const DISCONTINUITY_ANGLE: f64 = 0.2618; // ~15°

/// Number of sample points for centroid and discontinuity detection.
const PROFILE_SAMPLES: usize = 64;

/// Analyzed profile ready for sweeping.
#[derive(Debug)]
pub struct AnalyzedProfile {
    /// Smooth spans of the profile (split at tangent discontinuities).
    /// Every span's knot domain is normalized to [0, 1]: the sweep
    /// surface and topology stages evaluate span/surface corners at
    /// exactly 0 and 1, and a span carrying its source domain (an
    /// adopted rim arc's angles) put u = 0 mid-wire — ring vertices
    /// landed 90° off every cap curve (PK-479). Sub-[0,1] split
    /// domains only ever worked by clamped evaluation.
    pub spans: Vec<NurbsCurve>,
    /// Span count per wire, in `spans` order. A single-wire profile is
    /// one entry; a multi-wire profile lists the OUTER wire first and
    /// each inner (hole) wire after it. The topology builder wires one
    /// side-face ring per wire and one cap per path end whose outer
    /// loop is wire 0 and whose inner loops are the rest (PK-477 —
    /// flattening all wires into one ring cross-stitched the walls).
    pub wires: Vec<usize>,
    /// Whether the profile is closed.
    pub is_closed: bool,
    /// Centroid offset (profile was translated by -centroid).
    pub centroid: Point3,
}

/// Analyze a sweep profile: detect discontinuities, center, classify.
pub fn analyze_profile(
    profile: &SweepProfile,
    tolerance: &Tolerance,
) -> Result<AnalyzedProfile, SweepError> {
    match profile {
        SweepProfile::Wire(curve) => analyze_wire(curve, tolerance),
        SweepProfile::MultiWire(curves) => analyze_multi_wire(curves, tolerance),
        SweepProfile::Face { brep, face } => analyze_face(brep, *face, tolerance),
    }
}

/// Analyze a face profile by extracting its boundary loops as
/// curve-faithful wire chains (PK-479): each loop half-edge adopts
/// its edge's stored curve over its stored range, oriented to the
/// traversal, so circular holes stay exact circles instead of
/// degenerating to their single vertex. Edges without stored curves
/// fall back to the vertex chord. The profile-plane winding logic
/// projects onto XY (the wire contract shared with `analyze_wire`).
fn analyze_face(
    brep: &BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
) -> Result<AnalyzedProfile, SweepError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| SweepError::InvalidProfile("face handle not found in BRep".into()))?;

    let mut chains = vec![loop_to_wire_chain(brep, face_data.outer_loop)?];
    for &inner_loop in &face_data.inner_loops {
        chains.push(loop_to_wire_chain(brep, inner_loop)?);
    }
    analyze_wire_chains(chains, tolerance)
}

/// Analyze a single wire profile.
fn analyze_wire(curve: &NurbsCurve, tolerance: &Tolerance) -> Result<AnalyzedProfile, SweepError> {
    if curve.control_points.is_empty() {
        return Err(SweepError::InvalidProfile("empty wire".into()));
    }

    let centroid = compute_centroid(curve);
    let centered = translate_curve(curve, centroid);
    let is_closed = check_closed(&centered, tolerance);
    // Closed profiles are normalized to counter-clockwise winding in
    // the profile plane (positive shoelace area). The sweep topology
    // authors side-face loops against the (u = profile, v = path)
    // frame and takes du×dv as each face's outward direction, so a
    // clockwise profile yields a coherently inside-out solid — every
    // flag-folded normal INTO the material, which per-face winding
    // honesty cannot see. Measured (PK-475): the external-thread
    // groove tool, swept from a clockwise-wound profile, entered the
    // subtraction folded-inward; every tool fragment walked its
    // shared helix edges reversed, the shank wall's crest-band cells
    // were handed wrong-twin seam ladders, and the served loops gaped
    // by exactly one pitch. Winding is a construction invariant here,
    // not a caller contract.
    let centered = if is_closed && signed_area_2d(&centered) < 0.0 {
        reversed_curve(&centered)
    } else {
        centered
    };
    let spans: Vec<NurbsCurve> = detect_and_split(&centered)?
        .into_iter()
        .map(normalized_domain)
        .collect();

    let wires = vec![spans.len()];
    Ok(AnalyzedProfile {
        spans,
        wires,
        is_closed,
        centroid,
    })
}

/// Compute centroid as the average of sampled points along a curve.
fn compute_centroid(curve: &NurbsCurve) -> Point3 {
    let (t_min, t_max) = curve.parameter_range();
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;
    let n = PROFILE_SAMPLES;

    // Half-open sampling: for a CLOSED wire, t_max evaluates to the
    // same point as t_min — including both double-counts the seam
    // and shifts the centroid by r/(n+1) toward it. Every swept
    // body was rigidly mistranslated by that bias (measured:
    // 0.003846 = 0.25/65 on an r = 0.25 circle profile, visible as
    // an off-centre recognized cylinder). Half-open sampling is
    // exact for symmetric profiles and unbiased for all closed
    // wires; for open wires the missing endpoint is one sample in
    // PROFILE_SAMPLES on an anchor whose only job is consistency.
    for i in 0..n {
        let t = t_min + (i as f64 / n as f64) * (t_max - t_min);
        let p = curve.point_at(t);
        sum_x += p.x;
        sum_y += p.y;
        sum_z += p.z;
    }

    let count = n as f64;
    Point3::new(sum_x / count, sum_y / count, sum_z / count)
}

/// Translate a curve's control points by -offset (centering at origin).
fn translate_curve(curve: &NurbsCurve, offset: Point3) -> NurbsCurve {
    let shift = Vec3::new(-offset.x, -offset.y, -offset.z);
    let control_points = curve.control_points.iter().map(|p| *p + shift).collect();
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points,
        weights: curve.weights.clone(),
    }
}

/// Check if a curve is closed (start point ≈ end point).
fn check_closed(curve: &NurbsCurve, tolerance: &Tolerance) -> bool {
    let (t_min, t_max) = curve.parameter_range();
    let start = curve.point_at(t_min);
    let end = curve.point_at(t_max);
    tolerance.points_coincident(&start, &end)
}

/// Detect tangent discontinuities and split the curve at those parameters.
pub(crate) fn detect_and_split(curve: &NurbsCurve) -> Result<Vec<NurbsCurve>, SweepError> {
    let (t_min, t_max) = curve.parameter_range();
    let n = PROFILE_SAMPLES;

    // Tangent discontinuities of a NURBS can only occur at interior
    // knots of multiplicity ≥ degree: test one-sided tangents AT
    // those knots and split exactly there. (The previous sampled
    // scan split at sample-interval midpoints, chipping polygon
    // corners and double-splitting kinks that landed near a sample —
    // a closed rectangle produced 7 spans instead of 4.)
    let mut split_params = Vec::new();
    let mut candidate_knots: Vec<f64> = Vec::new();
    let mut i = 0;
    while i < curve.knots.len() {
        let k = curve.knots[i];
        let mut mult = 1;
        while i + mult < curve.knots.len() && (curve.knots[i + mult] - k).abs() < 1e-12 {
            mult += 1;
        }
        if k > t_min + 1e-10 && k < t_max - 1e-10 && mult >= curve.degree {
            candidate_knots.push(k);
        }
        i += mult;
    }
    let eps = (t_max - t_min) * 1.0e-7;
    // A POLYLINE's interior vertices are authored corners — the
    // angular threshold exists for sampled smooth curves where
    // numeric noise fakes kinks. A 24-gon arc's corners sit exactly
    // AT the 15° threshold, merged into one "smooth" span, and the
    // sweep emitted a single C0-kinked NURBS wall (pcurve-less, and
    // the boolean's marching stage hung on it). Degree-1 profiles
    // split at every genuine corner; only collinear joints merge.
    let threshold = if curve.degree == 1 {
        1.0e-6
    } else {
        DISCONTINUITY_ANGLE
    };
    // Curvature jump across a knot marks an authored TYPE boundary
    // (line↔arc) even when tangents agree: a mixed sketch path joins
    // straight and circular spans tangentially, and merging them
    // yields an unrecognizable NURBS wall. A full circle's quadrant
    // joints have EQUAL curvature on both sides and stay merged.
    let kappa = |t: f64| -> f64 {
        let d1 = curve.derivative_at(t, 1);
        let d2 = curve.derivative_at(t, 2);
        let num = d1.cross(&d2).length();
        let den = d1.length().powi(3);
        if den < 1.0e-30 { 0.0 } else { num / den }
    };
    for &k in &candidate_knots {
        let before = curve.tangent_at(k - eps);
        let after = curve.tangent_at(k + eps);
        if angle_between(&before, &after) > threshold {
            split_params.push(k);
            continue;
        }
        if curve.degree >= 2 {
            let k_l = kappa(k - eps);
            let k_r = kappa(k + eps);
            if (k_l - k_r).abs() > 0.25 * (k_l + k_r) + 1.0e-9 {
                split_params.push(k);
            }
        }
    }
    let _ = n;

    if split_params.is_empty() {
        return Ok(vec![curve.clone()]);
    }

    // Split curve at each discontinuity parameter, working from right to left.
    split_params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut spans = Vec::new();
    let mut remaining = curve.clone();

    for &t in split_params.iter().rev() {
        let (t_lo, t_hi) = remaining.parameter_range();
        // Only split if parameter is strictly inside the range.
        if t > t_lo + 1e-10 && t < t_hi - 1e-10 {
            let (left, right) = remaining.split_at(t).map_err(|e| {
                SweepError::InvalidProfile(format!("split failed at t={t:.6}: {e}"))
            })?;
            spans.push(right);
            remaining = left;
        }
    }
    spans.push(remaining);
    spans.reverse();

    // Normalize every span to the [0, 1] parameter range: split
    // sub-curves keep their sub-range knot vectors, and downstream
    // consumers (surface generation corner evaluation, pcurve
    // authoring) address spans at u ∈ [0, 1] — evaluating a
    // [0.25, 0.5]-ranged span at u = 1 extrapolates garbage.
    for span in &mut spans {
        let (lo, hi) = span.parameter_range();
        let width = hi - lo;
        if width > 1e-12 && ((lo).abs() > 1e-12 || (hi - 1.0).abs() > 1e-12) {
            for k in &mut span.knots {
                *k = (*k - lo) / width;
            }
        }
    }

    Ok(spans)
}

/// Angle between two vectors in radians.
fn angle_between(a: &Vec3, b: &Vec3) -> f64 {
    let la = a.length();
    let lb = b.length();
    if la < 1e-30 || lb < 1e-30 {
        return 0.0;
    }
    let cos_theta = a.dot(b) / (la * lb);
    cos_theta.clamp(-1.0, 1.0).pacos()
}

/// Reverse a NURBS curve's direction: mirrored knot vector, reversed
/// control points and weights. The parameter range is preserved.
fn reversed_curve(curve: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = curve.parameter_range();
    let mut knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots,
        control_points,
        weights,
    }
}

/// Approximate signed area of a curve projected onto XY plane (shoelace formula).
fn signed_area_2d(curve: &NurbsCurve) -> f64 {
    let (t_min, t_max) = curve.parameter_range();
    let n = PROFILE_SAMPLES;
    let mut area = 0.0;

    for i in 0..n {
        let t0 = t_min + (i as f64 / n as f64) * (t_max - t_min);
        let t1 = t_min + ((i + 1) as f64 / n as f64) * (t_max - t_min);
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);
        area += (p0.x * p1.y) - (p1.x * p0.y);
    }

    area * 0.5
}

#[cfg(test)]
#[path = "profile_tests.rs"]
mod tests;

#[path = "profile_multi_wire.rs"]
mod profile_multi_wire;
use profile_multi_wire::{
    analyze_multi_wire, analyze_wire_chains, loop_to_wire_chain, normalized_domain,
};
