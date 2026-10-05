//! Sampling budgets for sweep-path discretization: how many frames a
//! path earns and how finely its arc-length table resolves it. Both
//! budgets must track the path's OWN complexity — absolute caps here
//! were the silent quality cliffs of the multi-wrap fixed-cap family
//! (PK-461 catalogue, PK-495 instances): correct up to some helix
//! turn count, then phase-luck degradation with no diagnostic.

use crate::arc_length;
use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_math::PortableMath;

/// Default number of arc-length samples for path parameterization.
pub(crate) const ARC_LENGTH_SAMPLES: usize = 100;

/// Minimum number of frames for sweep operations.
pub(crate) const MIN_FRAMES: usize = 10;

/// Arc-length table resolution for a path: the floor is
/// [`ARC_LENGTH_SAMPLES`], raised to 32 samples per knot span so the
/// piecewise-linear `s ↔ t` interpolation resolves the path's own
/// parameterization *within* each span. Two per span only located the
/// spans: on the degree-5 helix arcs (PK-675) the speed varies inside
/// each arc, the table's chord error misplaced frames along the coil,
/// and the skinned thread tool walls left their helicoids by 8.7e-5
/// to 2.3e-4 axially (M3..M10). At 8 per span the error is 2.3e-5 to
/// 6.3e-5, at 32 it is 2.0e-6 to 5.5e-6, and 128 buys about a factor
/// of three more (measured, `helicoid_fidelity_tests`). A fixed 100 was a silent aliasing defect on
/// long coils: on a 102-turn thread helix each table cell spanned
/// 1.02 turns, the cell boundaries' phase drifted 2% of a turn per
/// cell against the helix period, and the interpolated frame
/// placement wobbled with a 51-turn beat — the skinned flank dipped
/// ≥ 0.03 radially at the beat maxima (measured, PK-495: the radial
/// level-curve SSI arm's station scans found zero crossings at
/// exactly the beat phases, period 25.5 turns per flank, the two
/// flanks in antiphase). The cap bounds table cost on pathologically
/// knotted paths: it binds past 512 helix turns (4 spans per turn),
/// beyond the served thread classes, and past it the in-span
/// resolution thins smoothly; placement does not fail.
pub(crate) fn arc_length_samples(path: &NurbsCurve) -> usize {
    ARC_LENGTH_SAMPLES.max(32 * knot_spans(path)).min(65_536)
}

/// Compute sweep frame count from path arc length and complexity.
pub(crate) fn compute_sample_count(table: &arc_length::ArcLengthTable, path: &NurbsCurve) -> usize {
    // Roughly 10 samples per unit length — PLUS a knot-span floor:
    // frame density must track the path's own complexity, not its
    // absolute size. A 5-turn helix authored at millimetre-in-SI
    // scale has total length ~0.1, the length rule handed it
    // MIN_FRAMES (~2 frames per turn), and the skinned tool surface
    // genuinely folded (PK-460 — the sweep's self-intersection gate
    // caught it, but the construction, not the gate, was wrong). Six
    // samples per distinct knot span keeps every quarter-arc of a
    // stitched helix resolved (15° per frame) at any scale. The
    // length rule keeps its long-tuned [MIN_FRAMES, 200] band —
    // large-scale sweeps are bit-identical to before; the span floor
    // only RAISES counts the length rule under-serves.
    //
    // Guard band: the span floor is allowed past 1024 frames only as
    // far as the path's MEASURED total turning justifies — spans
    // alone are not evidence of geometric complexity (a fitted
    // straight line can carry thousands of knots), but an absolute
    // cap here was a silent quality cliff: a constant-pitch thread
    // helix stitches 4 quarter-arc spans per turn, so a fixed 1024
    // bit at ~43 turns and starved a 100-turn tool wall down to 10
    // frames/turn — the skinned wall's radial ripple grew from its
    // 2.4e-3 baseline to 7e-2 (measured,
    // `probe_wall_radial_ripple`), large enough that the radial
    // level-curve SSI arm's structure gates declined the flank×shank
    // pairs and the thread serve frontier ended in a strand-gate
    // decline (PK-495). Turning-derived frames keep every genuinely
    // coiled path fully resolved (15°/frame, the same density
    // six-per-span gives a quarter arc) while zero-turning knot
    // pathology still clamps.
    let by_length = ((table.total_length * 10.0).ceil() as usize).clamp(MIN_FRAMES, 200);
    let by_spans = 6 * knot_spans(path);
    let span_frames = if by_spans > 1024 {
        by_spans.min(1024.max(turning_frames(path)))
    } else {
        by_spans
    };
    by_length.max(span_frames)
}

/// Distinct (non-degenerate) knot spans of a path.
fn knot_spans(path: &NurbsCurve) -> usize {
    path.knots
        .windows(2)
        .filter(|w| w[1] > w[0] + 1e-12)
        .count()
}

/// Frames justified by the path's total turning: accumulated tangent
/// direction change over a span-resolving sampling, at 15° per frame
/// — the density `6 * spans` gives a stitched quarter arc, so the
/// two measures agree exactly on the paths the span floor was tuned
/// for. Degenerate tangents are skipped (they carry no direction);
/// the sampling itself is capped, which UNDER-estimates turning on
/// pathologically knotted paths and therefore only ever tightens the
/// guard, never loosens it.
fn turning_frames(path: &NurbsCurve) -> usize {
    let (t0, t1) = path.parameter_range();
    if t1 <= t0 || !t1.is_finite() || !t0.is_finite() {
        return 0;
    }
    let n = (4 * knot_spans(path)).clamp(64, 16_384);
    let mut turning = 0.0_f64;
    let mut prev: Option<popkern_core::geometry::Vec3> = None;
    for i in 0..=n {
        let t = t0 + (t1 - t0) * i as f64 / n as f64;
        let d = path.derivative_at(t, 1);
        let Some(dir) = d.normalized() else { continue };
        if let Some(p) = prev {
            turning += p.dot(&dir).clamp(-1.0, 1.0).pacos();
        }
        prev = Some(dir);
    }
    let per_frame = 15.0_f64.to_radians();
    (turning / per_frame).ceil() as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{HelixPitch, HelixSpec};
    use popkern_core::geometry::{Point3, Vec3};

    fn helix_path(turns: f64) -> NurbsCurve {
        let spec = HelixSpec {
            axis_origin: Point3::new(0.0, 0.0, 0.0),
            axis_dir: Vec3::new(0.0, 0.0, 1.0),
            start_radius: 3.0,
            end_radius: 3.0,
            pitch: HelixPitch::Constant(1.0),
            turns,
            start_angle: 0.0,
            left_handed: false,
        };
        crate::helix::build_helix_path(&spec).expect("helix path").0
    }

    /// PK-495 regression: a 100-turn constant-pitch helix earns its
    /// full turning-justified frame budget (~24/turn) — the old
    /// absolute 1024 cap starved it to 10/turn and the skinned tool
    /// wall rippled 7e-2 radially.
    #[test]
    fn hundred_turn_helix_frames_track_turning() {
        let path = helix_path(100.0);
        let table = arc_length::ArcLengthTable::build(&path, arc_length_samples(&path));
        let n = compute_sample_count(&table, &path);
        assert!(n >= 2000, "100-turn helix starved to {n} frames");
        assert!(n <= 3200, "100-turn helix frame budget unbounded: {n}");
    }

    /// Zero-turning knot pathology still clamps: a straight polyline
    /// with hundreds of knot spans carries no geometric complexity
    /// and must not earn frames past the guard band.
    #[test]
    fn zero_turning_knot_pathology_still_clamps() {
        let spans = 300usize;
        let mut knots = vec![0.0];
        knots.extend((0..=spans).map(|i| i as f64));
        knots.push(spans as f64);
        let control_points = (0..=spans)
            .map(|i| Point3::new(0.0, 0.0, i as f64 * 0.01))
            .collect::<Vec<_>>();
        let weights = vec![1.0; control_points.len()];
        let path = NurbsCurve {
            degree: 1,
            knots,
            control_points,
            weights,
        };
        let table = arc_length::ArcLengthTable::build(&path, arc_length_samples(&path));
        let n = compute_sample_count(&table, &path);
        assert!(n <= 1024, "zero-turning path escaped the guard band: {n}");
    }

    /// The arc-length table resolves every knot span (two samples per
    /// span) so its piecewise-linear `s ↔ t` map cannot beat against
    /// a long coil's period; short paths keep the tuned floor.
    #[test]
    fn arc_length_table_is_span_resolving() {
        let long = helix_path(100.0);
        let spans = long
            .knots
            .windows(2)
            .filter(|w| w[1] > w[0] + 1e-12)
            .count();
        assert!(arc_length_samples(&long) >= 32 * spans);
        // A short coil still resolves within its spans (PK-675): the
        // span rule, not the fixed floor, sets its table.
        let short = helix_path(2.0);
        let short_spans = knot_spans(&short);
        assert_eq!(
            arc_length_samples(&short),
            ARC_LENGTH_SAMPLES.max(32 * short_spans)
        );
        // A path with few spans keeps the floor.
        let line = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        assert_eq!(arc_length_samples(&line), ARC_LENGTH_SAMPLES);
    }
}
