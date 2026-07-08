//! Profile analysis for sweep operations.
//!
//! Detects tangent discontinuities, classifies open/closed profiles,
//! and centers profiles at the origin for frame-based sweeping.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, LoopHandle};

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
    pub spans: Vec<NurbsCurve>,
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

/// Analyze a face profile by extracting its boundary loops as polyline wires.
fn analyze_face(
    brep: &BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
) -> Result<AnalyzedProfile, SweepError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| SweepError::InvalidProfile("face handle not found in BRep".into()))?;

    let outer_curve = loop_to_polyline(brep, face_data.outer_loop)?;

    if face_data.inner_loops.is_empty() {
        analyze_wire(&outer_curve, tolerance)
    } else {
        let mut curves = vec![outer_curve];
        for &inner_loop in &face_data.inner_loops {
            curves.push(loop_to_polyline(brep, inner_loop)?);
        }
        analyze_multi_wire(&curves, tolerance)
    }
}

/// Walk a loop's half-edges and build a degree-1 polyline NurbsCurve through the vertex positions.
///
/// The curve is closed: the last vertex duplicates the first.
fn loop_to_polyline(brep: &BRep, loop_handle: LoopHandle) -> Result<NurbsCurve, SweepError> {
    let loop_data = brep
        .loops
        .get(loop_handle)
        .ok_or_else(|| SweepError::InvalidProfile("loop handle not found in BRep".into()))?;

    let start_he = loop_data.half_edge;
    let mut positions: Vec<Point3> = Vec::new();
    let mut current_he = start_he;

    loop {
        let he = brep.half_edges.get(current_he).ok_or_else(|| {
            SweepError::InvalidProfile("half-edge handle not found in BRep".into())
        })?;
        let vertex = brep
            .vertices
            .get(he.origin)
            .ok_or_else(|| SweepError::InvalidProfile("vertex handle not found in BRep".into()))?;
        positions.push(vertex.position);
        current_he = he.next;
        if current_he == start_he {
            break;
        }
    }

    if positions.is_empty() {
        return Err(SweepError::InvalidProfile("loop has no half-edges".into()));
    }

    // Close the polyline by repeating the first vertex.
    positions.push(positions[0]);

    let n = positions.len(); // number of control points
    // Degree-1 polyline knots: clamped, uniform interior knots.
    // knot count = n + degree + 1 = n + 2
    let mut knots = Vec::with_capacity(n + 2);
    knots.push(0.0);
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);

    let weights = vec![1.0; n];

    Ok(NurbsCurve {
        degree: 1,
        knots,
        control_points: positions,
        weights,
    })
}

/// Analyze a single wire profile.
fn analyze_wire(curve: &NurbsCurve, tolerance: &Tolerance) -> Result<AnalyzedProfile, SweepError> {
    if curve.control_points.is_empty() {
        return Err(SweepError::InvalidProfile("empty wire".into()));
    }

    let centroid = compute_centroid(curve);
    let centered = translate_curve(curve, centroid);
    let is_closed = check_closed(&centered, tolerance);
    let spans = detect_and_split(&centered)?;

    Ok(AnalyzedProfile {
        spans,
        is_closed,
        centroid,
    })
}

/// Analyze multiple wires: use the largest (by area) as the outer profile.
fn analyze_multi_wire(
    curves: &[NurbsCurve],
    tolerance: &Tolerance,
) -> Result<AnalyzedProfile, SweepError> {
    if curves.is_empty() {
        return Err(SweepError::InvalidProfile("empty multi-wire".into()));
    }

    // Find the curve with the largest absolute signed area (the outer wire).
    let mut best_idx = 0;
    let mut best_area = 0.0_f64;
    for (i, curve) in curves.iter().enumerate() {
        let area = signed_area_2d(curve).abs();
        if area > best_area {
            best_area = area;
            best_idx = i;
        }
    }

    let centroid = compute_centroid(&curves[best_idx]);

    let mut all_spans = Vec::new();
    let mut any_open = false;

    for curve in curves {
        let centered = translate_curve(curve, centroid);
        if !check_closed(&centered, tolerance) {
            any_open = true;
        }
        let spans = detect_and_split(&centered)?;
        all_spans.extend(spans);
    }

    Ok(AnalyzedProfile {
        spans: all_spans,
        is_closed: !any_open,
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
fn detect_and_split(curve: &NurbsCurve) -> Result<Vec<NurbsCurve>, SweepError> {
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
    for &k in &candidate_knots {
        let before = curve.tangent_at(k - eps);
        let after = curve.tangent_at(k + eps);
        if angle_between(&before, &after) > threshold {
            split_params.push(k);
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
    cos_theta.clamp(-1.0, 1.0).acos()
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
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Create a unit circle as a rational NURBS (degree 2, 9 control points).
    fn unit_circle() -> NurbsCurve {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsCurve {
            degree: 2,
            knots: vec![
                0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
            ],
            control_points: vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(-1.0, 1.0, 0.0),
                Point3::new(-1.0, 0.0, 0.0),
                Point3::new(-1.0, -1.0, 0.0),
                Point3::new(0.0, -1.0, 0.0),
                Point3::new(1.0, -1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
        }
    }

    /// Create a simple line segment from (0,0,0) to (5,0,0).
    fn line_segment() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(5.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    /// Create a rectangle as a degree-1 polyline (4 corners + close).
    fn rectangle_polyline() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
            control_points: vec![
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(-1.0, 1.0, 0.0),
                Point3::new(-1.0, -1.0, 0.0),
                Point3::new(1.0, -1.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            weights: vec![1.0; 5],
        }
    }

    #[test]
    fn circle_is_closed_one_span() {
        let tol = Tolerance::default().with_linear(1e-3);
        let profile = SweepProfile::Wire(unit_circle());
        let analyzed = analyze_profile(&profile, &tol).ok();
        assert!(analyzed.is_some());
        let analyzed = analyzed.unwrap(); // ok in tests
        assert!(analyzed.is_closed);
        // Circle is smooth, so should remain 1 span.
        assert_eq!(analyzed.spans.len(), 1);
    }

    #[test]
    fn line_segment_is_open() {
        let tol = Tolerance::default();
        let profile = SweepProfile::Wire(line_segment());
        let analyzed = analyze_profile(&profile, &tol).ok();
        assert!(analyzed.is_some());
        let analyzed = analyzed.unwrap(); // ok in tests
        assert!(!analyzed.is_closed);
        assert_eq!(analyzed.spans.len(), 1);
    }

    #[test]
    fn rectangle_is_closed_and_splits() {
        let tol = Tolerance::default().with_linear(1e-3);
        let profile = SweepProfile::Wire(rectangle_polyline());
        let analyzed = analyze_profile(&profile, &tol).ok();
        assert!(analyzed.is_some());
        let analyzed = analyzed.unwrap(); // ok in tests
        assert!(analyzed.is_closed);
        // A rectangle polyline has 90° corners → 4 spans (one per side).
        assert!(
            analyzed.spans.len() >= 3,
            "expected >=3 spans for rectangle, got {}",
            analyzed.spans.len()
        );
    }

    #[test]
    fn centroid_near_origin_for_centered_circle() {
        let tol = Tolerance::default().with_linear(1e-3);
        let profile = SweepProfile::Wire(unit_circle());
        let analyzed = analyze_profile(&profile, &tol).ok();
        assert!(analyzed.is_some());
        let analyzed = analyzed.unwrap(); // ok in tests
        // Circle is centered at origin, so centroid should be near (0,0,0).
        let c = analyzed.centroid;
        assert!(
            c.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 0.1,
            "centroid should be near origin, got ({}, {}, {})",
            c.x,
            c.y,
            c.z
        );
    }

    #[test]
    fn empty_wire_returns_error() {
        let tol = Tolerance::default();
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![],
            weights: vec![],
        };
        let result = analyze_profile(&SweepProfile::Wire(curve), &tol);
        assert!(result.is_err());
    }

    #[test]
    fn multi_wire_uses_largest_centroid() {
        let tol = Tolerance::default().with_linear(1e-3);
        // Outer: large circle, Inner: small line segment
        let outer = unit_circle();
        let inner = line_segment();
        let profile = SweepProfile::MultiWire(vec![inner, outer]);
        let analyzed = analyze_profile(&profile, &tol).ok();
        assert!(analyzed.is_some());
        let analyzed = analyzed.unwrap(); // ok in tests
        // Centroid should be near the circle's center, not the line's center.
        let c = analyzed.centroid;
        assert!(
            c.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 0.2,
            "centroid should be near circle center"
        );
    }

    #[test]
    fn signed_area_circle_positive() {
        let circle = unit_circle();
        let area = signed_area_2d(&circle);
        // Should be approximately π for unit circle.
        assert!(
            (area.abs() - PI).abs() < 0.2,
            "area should be ~π, got {}",
            area.abs()
        );
    }

    #[test]
    fn face_profile_from_cube() {
        use popkern_core::topology::BRep;
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let profile = SweepProfile::Face {
            brep: Box::new(brep.clone()),
            face,
        };
        let result = analyze_profile(&profile, &Tolerance::default());
        assert!(
            result.is_ok(),
            "face profile should analyze: {:?}",
            result.err()
        );
        let analyzed = result.unwrap();
        assert!(analyzed.is_closed);
        assert!(!analyzed.spans.is_empty());
    }
}
