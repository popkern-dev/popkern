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
//! STATUS: exploration arm, DEFAULT OFF —
//! `POPKERN_RECORD_SIDE_CLASSIFY=1` enables it for frontier work
//! (PK-239 kill-switch precedent; per-call read).
//!
//! REMAINING DEFECT (measured, overlapping-cubes suite with the arm
//! on): the validity conditions hold but the VERDICT SIGN is wrong on
//! faces whose stored `orientation` flag lies about the outward
//! normal — a documented kernel condition ("extrude caps lie in the
//! flag, chamfer cones in the loop winding", see the cellular
//! coorientation self-calibration). Worked cube example: fragment
//! sample (0.75,1,1) inside B=[0.5,1.5]³, q=(0.5,1,1) on B's x-wall,
//! honest N=(-1,0,0) gives dot<0 → Inside (correct); the stored flag
//! yields Outside. NEXT MECHANISM: per-face outward-sign calibration
//! — either the majority vote of (ray verdict × witness dot) across
//! all fragments touching each face's records (self-consistent, no
//! new probes; overrides only minority mismatches), or a shared
//! honest-normal service (probe/mesh-calibrated once per face per
//! boolean) that the cellular coorientation logic could also adopt.
//! Enable by default only after the full corpus and the multiwrap
//! probe ladder are green under it.

use popkern_core::geometry::intersection::surface_curve::intersect_surface_curve_dense;
use popkern_core::geometry::intersection::{SurfaceCurveIntersection, surface_subdivision::Aabb};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_core::topology::FaceHandle;

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
    other_aabbs: &[(FaceHandle, Aabb)],
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
        // Validity (a): q must land ON the other face's trimmed
        // extent — its chart window (from the face's own boundary
        // samples; native window, small pad) and its spatial AABB.
        {
            let positions = other_brep.face_outer_loop_sample_positions(other_face, 8);
            if positions.len() < 3 {
                continue;
            }
            let ((wu0, wu1), (wv0, wv1)) =
                popkern_core::geometry::uv_projection::exact_bounded_parameter_range(
                    &positions,
                    &surface.geometry,
                );
            let pad_u = 0.02 * (wu1 - wu0).abs().max(1e-9);
            let pad_v = 0.02 * (wv1 - wv0).abs().max(1e-9);
            if qu < wu0.min(wu1) - pad_u
                || qu > wu0.max(wu1) + pad_u
                || qv < wv0.min(wv1) - pad_v
                || qv > wv0.max(wv1) + pad_v
            {
                continue;
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

/// True when the OPEN segment `s → q` crosses any face of `brep`
/// strictly before `q` — or when any hit is too ambiguous to rule
/// out (conservative: ambiguity counts as obstruction).
fn segment_is_obstructed(
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
        let (wu, wv) = {
            let r = popkern_core::geometry::uv_projection::exact_bounded_parameter_range(
                &positions,
                &sof.geometry,
            );
            (r.0, r.1)
        };
        let pad_u = 0.05 * (wu.1 - wu.0).abs().max(1e-9);
        let pad_v = 0.05 * (wv.1 - wv.0).abs().max(1e-9);
        let window_u = (wu.0.min(wu.1) - pad_u, wu.0.max(wu.1) + pad_u);
        let window_v = (wv.0.min(wv.1) - pad_v, wv.0.max(wv.1) + pad_v);
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
