//! Record-adjacent material-side classification (PK-494 layer 3).
//!
//! A fragment bounded by a Shared intersection record carries a LOCAL,
//! exact in/out witness against the other operand: at the record the
//! other operand's wall passes through the fragment's own surface, so
//! which side of that wall the fragment lies on IS the containment
//! verdict there — no rays, no global parity, no dependence on coil
//! scale. Point-in-solid ray voting misread exactly this class at the
//! thread coil's end region (the void-backed corner sliver above the
//! last groove record classified Outside from a sample 0.08 inside
//! the tool — every ray threading the 47-turn coil's cut-open end;
//! the kept sliver then triple-claimed its record downstream and the
//! serve declined or gapped by claim-order luck).
//!
//! Method, per Shared record on the fragment's outer loop whose two
//! source faces live on OPPOSITE operands:
//! 1. Take the record's curve-true midpoint `m` and the OTHER
//!    operand's face surface there. The Y-side graph pcurve gives the
//!    chart point exactly when present; otherwise a local projection
//!    seeded at the record end helps `closest_point` stay on the
//!    local sheet.
//! 2. Locally project the fragment's classification sample onto Y's
//!    surface (Gauss-Newton from the midpoint's chart point — LOCAL
//!    by construction: the verdict wanted is the side of THIS wall
//!    piece, not of whatever sheet is globally nearest).
//! 3. The sign of `dot(sample − q, N_Y(q))` (Y's outward normal) is
//!    the side: positive = outside Y, negative = inside Y. Near-zero
//!    (tangency, grazing contact) abstains.
//! 4. Votes from every qualifying record must agree, else abstain —
//!    the ray verdict then stands.
//!
//! VALIDITY MODEL (v2): a record's side verdict equals the solid
//! containment verdict exactly when two measurable conditions hold —
//! (a) the sample's local projection q lands ON the other face's
//! TRIMMED extent (its own chart window and spatial AABB, the same
//! gates the ray classifier applies to accepted crossings: an
//! off-trim projection is a basis-extension phantom, the naive
//! witness's overlapping-cubes failure), and (b) the open segment
//! sample→q crosses NO other-operand boundary strictly before q (an
//! obstruction means another wall lies between the sample and this
//! record's wall and the local side says nothing about containment —
//! the far-field failure). Any ambiguity in either test abstains and
//! the ray verdict stands: the witness only ever REPLACES a ray
//! verdict with a locally-proven one.
//!
//! Three further validity gates, each landed against a measured
//! failure (PK-494 layer 3):
//! - FOOT-POINT stationarity: (sample − q) must be orthogonal to both
//!   surface tangents at q — a projection clamped to the chart corner
//!   produced six identical bogus votes across a cone ring.
//! - PERIODIC-HONEST chart windows: minimal covering arc on periodic
//!   axes, failing CLOSED when the boundary samples straddle the seam.
//! - LOOP-HONEST trim, fail-closed: the outer loop must contain q
//!   with a REAL pcurve verdict and every inner loop must really
//!   exclude it — the plain box gate passed a projection into a
//!   cap's bore HOLE, and the bool wrappers fold missing-pcurve
//!   ambiguity OPEN (a rotated keyed-bore cap's fresh inner loop
//!   lacked pcurves and let a q-in-hole vote through).
//!
//! The stored per-face `orientation` flag LIES on some construction
//! paths ("extrude caps lie in the flag, chamfer cones in the loop
//! winding"), so votes carry the RAW chart-normal cosine and
//! `calibrate_face_signs` recovers each face's honest outward sign
//! from the ray-verdict consensus (≥75% majority of ≥2 votes); only
//! minority mismatches are overridden, via `calibrated_verdict`.
//!
//! STATUS: DEFAULT ON. `POPKERN_RECORD_SIDE_CLASSIFY=0` pins the
//! ray-only path for A/B (PK-239 kill-switch precedent; per-call
//! read). Enabled after: bool suite 531/531 both ways, fastener
//! thread 8/8 + generators 7/7 + holes 4/4 with the arm on, and the
//! multiwrap probe ladder watertight through 61 tool turns.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_core::topology::FaceHandle;

use crate::error::PointClassification;
use crate::merge::graph::{EdgeSource, IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

#[path = "record_side_witness_gates.rs"]
mod record_side_witness_gates;
use record_side_witness_gates::{chart_window, segment_is_obstructed};

/// Decisiveness gate on the normalized side dot: below it a record is
/// tangential/grazing at the sample's offset and abstains.
const SIDE_COS_GATE: f64 = 0.02;
/// Records examined per fragment (longest chords first).
const MAX_RECORDS: usize = 8;

/// A validity-gated raw side vote: the other operand's face this
/// record belongs to, and the UNCALIBRATED side cosine of the
/// fragment sample against that face's RAW chart normal (no
/// orientation flag — the stored flags lie on some construction
/// paths; `calibrate_face_signs` recovers the honest sign from the
/// ray-verdict consensus).
#[derive(Debug, Clone, Copy)]
pub(crate) struct RecordVote {
    /// The other operand's face at the voting record.
    pub(crate) other_face: FaceHandle,
    /// Whether that face belongs to operand A.
    pub(crate) other_is_a: bool,
    /// Side cosine against the RAW (unflagged) chart normal.
    pub(crate) cosine: f64,
}

/// Collect the validity-gated raw votes for `fragment`'s sample
/// against the other operand — empty when the arm is disabled or no
/// record qualifies.
pub(crate) fn record_votes(
    graph: &IntersectionGraph,
    fragment: &Fragment,
    sample: &Point3,
    brep_a: &BRep,
    brep_b: &BRep,
    other_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Vec<RecordVote> {
    let enabled = std::env::var_os("POPKERN_RECORD_SIDE_CLASSIFY")
        .map(|v| v != "0")
        .unwrap_or(true);
    if !enabled {
        return Vec::new();
    }
    let own = fragment.face.operand;
    let other_brep = match own {
        Operand::A => brep_b,
        Operand::B => brep_a,
    };

    // Qualifying records, longest chord first: long records carry the
    // most reliable local frame; slop-scale junction slivers go last.
    let mut candidates: Vec<(usize, f64)> = fragment
        .outer_loop
        .iter()
        .copied()
        .filter(|&ghe| matches!(graph.half_edges[ghe].source, EdgeSource::Shared { .. }))
        .map(|ghe| {
            let twin = graph.half_edges[ghe].twin;
            let a = graph.nodes[graph.half_edges[ghe].origin].position;
            let b = graph.nodes[graph.half_edges[twin].origin].position;
            (ghe, a.distance_to(&b))
        })
        .collect();
    candidates.sort_by(|x, y| y.1.partial_cmp(&x.1).unwrap_or(std::cmp::Ordering::Equal));
    candidates.truncate(MAX_RECORDS);

    let mut votes: Vec<RecordVote> = Vec::new();
    for (ghe, chord) in candidates {
        if chord <= tolerance.linear() * 10.0 {
            continue;
        }
        let EdgeSource::Shared { face_a, face_b } = graph.half_edges[ghe].source else {
            continue;
        };
        // The other operand's face at this record. A residue pair
        // registered twice on ONE operand has no opposite-side wall
        // here — skip it.
        let other_face = match own {
            Operand::A => face_b,
            Operand::B => face_a,
        };
        let Some(surface) = other_brep
            .faces
            .get(other_face)
            .and_then(|f| f.surface.as_ref())
        else {
            continue;
        };

        // Record midpoint (curve-true when a curve exists).
        let twin = graph.half_edges[ghe].twin;
        let mid = {
            use popkern_core::geometry::curve::CurveEval;
            match graph.half_edges[ghe]
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[twin].curve.as_ref())
            {
                Some(c) => c.geometry.point_at(0.5 * (c.t_start + c.t_end)),
                None => {
                    let a = graph.nodes[graph.half_edges[ghe].origin].position;
                    let b = graph.nodes[graph.half_edges[twin].origin].position;
                    Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
                }
            }
        };

        // Chart seed on Y at the midpoint: the Y-side pcurve is the
        // chart-true witness when either half of the pair carries it
        // for the other face; otherwise fall back to a global
        // projection of the midpoint (it lies ON the surface, so the
        // seeded multisheet machinery converges to the local sheet).
        let seed_uv = [ghe, twin]
            .into_iter()
            .filter_map(|h| {
                let he = &graph.half_edges[h];
                let pc = he.pcurve.as_ref()?;
                let fr = he.face.as_ref()?;
                (fr.operand != own && fr.face == other_face).then(|| {
                    let q = pc.point_at(0.5 * (pc.t_start + pc.t_end));
                    (q.x, q.y)
                })
            })
            .next()
            .or_else(|| surface.geometry.closest_point(&mid));
        let Some((su, sv)) = seed_uv else { continue };

        // Local projection of the sample onto Y's sheet at the record.
        let Some((qu, qv)) = popkern_core::geometry::uv_projection::project_point_to_uv_near(
            &surface.geometry,
            sample,
            su,
            sv,
        ) else {
            continue;
        };
        let q = surface.geometry.point_at(qu, qv);
        // Validity (a): q must land ON the other face's trimmed
        // extent — its chart window (from the face's own boundary
        // samples; native window, small pad) and its spatial AABB.
        {
            let positions = other_brep.face_outer_loop_sample_positions(other_face, 8);
            if positions.len() < 3 {
                continue;
            }
            let Some(((wu0, wu1), (wv0, wv1))) = chart_window(&surface.geometry, &positions) else {
                continue; // seam-straddling window — fail closed
            };
            let pad_u = 0.02 * (wu1 - wu0).abs().max(1e-9);
            let pad_v = 0.02 * (wv1 - wv0).abs().max(1e-9);
            // Periodic-honest containment: compare on the principal
            // branch nearest the window.
            let (qpu, qpv) = {
                let (pu, pv) = surface.geometry.periodicity();
                let fold = |x: f64, lo: f64, hi: f64, period: Option<f64>| -> f64 {
                    match period.filter(|p| p.is_finite() && *p > 0.0) {
                        Some(p) => {
                            let c = 0.5 * (lo + hi);
                            x - p * ((x - c) / p).round()
                        }
                        None => x,
                    }
                };
                (fold(qu, wu0, wu1, pu), fold(qv, wv0, wv1, pv))
            };
            if qpu < wu0 - pad_u || qpu > wu1 + pad_u || qpv < wv0 - pad_v || qpv > wv1 + pad_v {
                continue;
            }
            // LOOP-HONEST trim, FAIL-CLOSED: the box passes a
            // projection landing in a HOLE (the chamfered nut's cap
            // fragment projected into the bore's inner loop and voted
            // on a half-space that is not face material). The
            // pcurve-sampled even-odd tests are the ones ray
            // classification trusts — but their bool wrappers fold a
            // missing-pcurve "no verdict" OPEN, and a rotated
            // keyed-bore cap whose fresh inner bore loop lacked
            // pcurves let a q-in-hole vote through. Every loop must
            // give a REAL verdict here: outer Some(true), each inner
            // Some(false); any None abstains.
            {
                let Some(face_data) = other_brep.faces.get(other_face) else {
                    continue;
                };
                if crate::classify_uv_loops::uv_point_in_loop(
                    other_brep,
                    face_data.outer_loop,
                    qpu,
                    qpv,
                ) != Some(true)
                {
                    continue;
                }
                if !face_data.inner_loops.iter().all(|&lh| {
                    crate::classify_uv_loops::uv_point_in_loop(other_brep, lh, qpu, qpv)
                        == Some(false)
                }) {
                    continue;
                }
            }
            let m = tolerance.linear() * 10.0;
            let (mut lo, mut hi) = (positions[0], positions[0]);
            for p in &positions {
                lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
                hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
            }
            // Boundary samples bound the box up to curve sag; a
            // conservative diagonal-relative margin covers it.
            let diag = lo.distance_to(&hi);
            let m2 = m + 0.05 * diag;
            if q.x < lo.x - m2
                || q.x > hi.x + m2
                || q.y < lo.y - m2
                || q.y > hi.y + m2
                || q.z < lo.z - m2
                || q.z > hi.z + m2
            {
                continue;
            }
        }
        // Validity (b): the open segment sample→q must cross NO
        // other-operand boundary strictly before q. Any hit, and any
        // ambiguity, abstains.
        if segment_is_obstructed(sample, &q, other_brep, other_aabbs, tolerance) {
            continue;
        }
        let n = surface.geometry.normal_at(qu, qv);
        let n_len = n.length();
        let v = Vec3::new(sample.x - q.x, sample.y - q.y, sample.z - q.z);
        let v_len = v.length();
        if n_len < tolerance.numerical_floor() || v_len <= tolerance.linear() * 10.0 {
            continue; // degenerate normal, or sample ON the wall — abstain
        }
        // FOOT-POINT stationarity: a projection clamped to the chart
        // window (or stalled at a corner) is not a foot point, and a
        // side dot against it is meaningless — the chamfered-nut cap
        // fragment collected SIX identical votes from one clamped
        // corner of the cone ring at |s−q| = 0.74. A genuine foot has
        // (sample − q) orthogonal to both tangents.
        {
            let su = surface.geometry.partial_u(qu, qv);
            let sv_t = surface.geometry.partial_v(qu, qv);
            let ortho = |t: Vec3| -> bool {
                let tl = t.length();
                tl < tolerance.numerical_floor() || (v.dot(&t) / (v_len * tl)).abs() < 0.05
            };
            if !ortho(su) || !ortho(sv_t) {
                continue;
            }
        }
        let cosine = n.dot(&v) / (n_len * v_len);
        if cosine.abs() <= SIDE_COS_GATE {
            continue; // grazing — abstain
        }
        if popkern_core::trace::trace_enabled("uv-scramble") {
            eprintln!(
                "  [vote] frag {:?} other_face {:?} (A:{}) cos {cosine:+.3} q ({:.3},{:.3},{:.3}) |s-q| {:.4} chord {chord:.4}",
                fragment.face,
                other_face,
                own == Operand::B,
                q.x,
                q.y,
                q.z,
                sample.distance_to(&q)
            );
        }
        votes.push(RecordVote {
            other_face,
            other_is_a: own == Operand::B,
            cosine,
        });
    }
    votes
}

/// Per-face outward-sign calibration from the ray-verdict consensus:
/// for each other-operand face, choose the sign that makes the
/// majority of (vote cosine, ray verdict) pairs consistent
/// (`Outside` ⟺ oriented cosine > 0). Faces with fewer than
/// MIN_CALIBRATION_VOTES votes, or without a MIN_CALIBRATION_RATIO
/// majority, stay uncalibrated (no overrides there). Keyed by
/// (other_is_a, face).
pub(crate) fn calibrate_face_signs(
    all: &[(PointClassification, Vec<RecordVote>)],
) -> std::collections::HashMap<(bool, FaceHandle), f64> {
    const MIN_CALIBRATION_VOTES: usize = 2;
    const MIN_CALIBRATION_RATIO: f64 = 0.75;
    let mut tally: std::collections::HashMap<(bool, FaceHandle), (usize, usize)> =
        std::collections::HashMap::new();
    for (class, votes) in all {
        let ray_outside = match class {
            PointClassification::Outside => true,
            PointClassification::Inside => false,
            _ => continue,
        };
        for v in votes {
            let entry = tally.entry((v.other_is_a, v.other_face)).or_insert((0, 0));
            if (v.cosine > 0.0) == ray_outside {
                entry.0 += 1; // consistent with sign +1
            } else {
                entry.1 += 1; // consistent with sign −1
            }
        }
    }
    if popkern_core::trace::trace_enabled("uv-scramble") {
        for (key, (plus, minus)) in &tally {
            eprintln!("  [calib] face {key:?} +{plus} -{minus}");
        }
    }
    tally
        .into_iter()
        .filter_map(|(key, (plus, minus))| {
            let total = plus + minus;
            if total < MIN_CALIBRATION_VOTES {
                return None;
            }
            let (winner, share) = if plus >= minus {
                (1.0, plus as f64 / total as f64)
            } else {
                (-1.0, minus as f64 / total as f64)
            };
            (share >= MIN_CALIBRATION_RATIO).then_some((key, winner))
        })
        .collect()
}

/// The calibrated witness verdict: every vote's face must be
/// calibrated and every oriented vote must agree, else `None` (the
/// ray verdict stands).
pub(crate) fn calibrated_verdict(
    votes: &[RecordVote],
    signs: &std::collections::HashMap<(bool, FaceHandle), f64>,
) -> Option<PointClassification> {
    let mut verdict: Option<PointClassification> = None;
    if votes.is_empty() {
        return None;
    }
    for v in votes {
        let sign = signs.get(&(v.other_is_a, v.other_face))?;
        let this = if v.cosine * sign > 0.0 {
            PointClassification::Outside
        } else {
            PointClassification::Inside
        };
        match verdict {
            None => verdict = Some(this),
            Some(prev) if prev == this => {}
            Some(_) => return None,
        }
    }
    verdict
}
