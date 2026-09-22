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
//! STATUS (PK-494 layer 3, exploration arm): DEFAULT OFF —
//! `POPKERN_RECORD_SIDE_CLASSIFY=1` enables it for frontier work.
//! The near-field validity model is not yet commercial-grade: the
//! naive witness regressed the overlapping-cubes suite (side-of-wall
//! is not side-of-solid beyond the wall's trimmed extent), and the
//! first gating attempt (chord-relative distance + chart-window
//! containment) shifted which records vote and produced
//! unanimous-WRONG sets on 25 fixtures — chart-window tests are
//! themselves unreliable on periodic charts. The design fork this
//! arm needs: a validity model where a record votes only when the
//! sample verifiably lies in the record's own near-field (distance
//! bounded by local feature size, projection within the TRIMMED
//! face, path from record to sample free of other records), plus
//! defined multi-record voting semantics. Until then the ray path
//! stands and this arm exists for measured experiments only.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::PointClassification;
use crate::merge::graph::{EdgeSource, IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// Decisiveness gate on the normalized side dot: below it a record is
/// tangential/grazing at the sample's offset and abstains.
const SIDE_COS_GATE: f64 = 0.02;
/// Records examined per fragment (longest chords first).
const MAX_RECORDS: usize = 8;

/// The record-side containment verdict for `fragment`'s sample point
/// against the OTHER operand, or `None` when no bounding record gives
/// a decisive, unanimous answer.
pub(crate) fn record_side_witness(
    graph: &IntersectionGraph,
    fragment: &Fragment,
    sample: &Point3,
    brep_a: &BRep,
    brep_b: &BRep,
    tolerance: &Tolerance,
) -> Option<PointClassification> {
    let enabled = std::env::var_os("POPKERN_RECORD_SIDE_CLASSIFY")
        .map(|v| v == "1")
        .unwrap_or(false);
    if !enabled {
        return None;
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

    let mut verdict: Option<PointClassification> = None;
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
        // NEAR-FIELD validity gates: side-of-wall equals side-of-solid
        // only close to the record and on the wall's own trimmed
        // extent — a cube fragment's face-scale sample sits behind the
        // other face's PLANE while fully outside the CUBE (measured:
        // the overlapping-cubes suite regressed under the ungated
        // witness). Gate A: the sample's wall distance stays within a
        // quarter of the record chord. Gate B: the projection lands
        // inside the other face's own chart window.
        if sample.distance_to(&q) > 0.25 * chord {
            continue;
        }
        {
            let positions = other_brep.face_outer_loop_sample_positions(other_face, 8);
            if positions.is_empty() {
                continue;
            }
            let ((wu0, wu1), (wv0, wv1)) =
                popkern_core::geometry::uv_projection::exact_bounded_parameter_range(
                    &positions,
                    &surface.geometry,
                );
            let du = (wu1 - wu0).abs().max(1e-12);
            let dv = (wv1 - wv0).abs().max(1e-12);
            let pad_u = 0.05 * du;
            let pad_v = 0.05 * dv;
            if qu < wu0.min(wu1) - pad_u
                || qu > wu0.max(wu1) + pad_u
                || qv < wv0.min(wv1) - pad_v
                || qv > wv0.max(wv1) + pad_v
            {
                continue;
            }
        }
        let mut n = surface.geometry.normal_at(qu, qv);
        if !surface.orientation {
            n = n * -1.0;
        }
        let n_len = n.length();
        let v = Vec3::new(sample.x - q.x, sample.y - q.y, sample.z - q.z);
        let v_len = v.length();
        if n_len < tolerance.numerical_floor() || v_len <= tolerance.linear() * 10.0 {
            continue; // degenerate normal, or sample ON the wall — abstain
        }
        let cosine = n.dot(&v) / (n_len * v_len);
        let this = if cosine > SIDE_COS_GATE {
            PointClassification::Outside
        } else if cosine < -SIDE_COS_GATE {
            PointClassification::Inside
        } else {
            continue; // grazing — abstain
        };
        match verdict {
            None => verdict = Some(this),
            Some(prev) if prev == this => {}
            Some(_) => return None, // records disagree — rays stand
        }
    }
    verdict
}
