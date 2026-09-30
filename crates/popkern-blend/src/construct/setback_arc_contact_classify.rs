//! Contact classification for the setback arc/boundary crossing
//! detector: turns clustered UV hits into the surgery's two trim
//! events, or an honest decline.
//!
//! On iso-section setback corners the pre-carve apex cut lies AT the
//! setback station, so the arc is UV-COINCIDENT with the strip's
//! apex boundary BY CONSTRUCTION and projection noise reads the
//! coincidence as dense pseudo-crossing runs (measured on the squat
//! box 2.4x1.6x0.4 r=0.1: 21-hit runs per chord; the unit cube read
//! a clean 2 on identical local corner geometry — a numeric
//! knife-edge, not an aspect effect). Heuristics over run shape
//! (hit count, endpoint adjacency) fragment on sampling holes; the
//! decisive discriminator is TRANSVERSALITY: a genuine trim crossing
//! actually crosses — the arc's signed UV side of the edge polyline
//! flips with finite extent on both sides — while a coincident run
//! wiggles at projection-noise amplitude (~1e-9 of the chart) and
//! never leaves the edge.
//!
//! Nested inside `setback_vertex_topology` (shares the crossing
//! detector's chart conventions).

use super::*;

/// Arc-t adjacency band chaining raw UV hits into one contact
/// cluster (shared with the detector's chain clustering).
pub(super) const CLUSTER_ADJACENCY: f64 = 0.05;

/// Side extent, as a fraction of the arc's UV diameter, that a mid-arc
/// contact must reach on both sides to count as a crossing. Over the
/// cluster window (±[`CLUSTER_ADJACENCY`] in arc t) a transversal
/// crossing leaves the edge by about `0.05 · diam · sin(angle)`, above
/// 1e-2 for any crossing steeper than ~12° (measured genuine crossings:
/// 4e-2); a graze stays at projection-scale excursions (measured
/// 2.3e-4, PK-591).
const MID_ARC_CROSSING_EXTENT: f64 = 1.0e-2;

/// One contact event: the best (min-margin) representative hit plus
/// the arc-t span and population of its cluster.
pub(super) struct ContactCluster {
    /// `(arc_t, edge, edge_point, margin)` of the closest contact.
    pub rep: (f64, EdgeHandle, Point3, f64),
    /// First arc-t of the cluster.
    pub t_first: f64,
    /// Last arc-t of the cluster.
    pub t_last: f64,
    /// Raw hits chained into this cluster.
    pub count: usize,
}

impl ContactCluster {
    fn endpoint(&self) -> bool {
        self.rep.0 < 1.0e-3 || self.rep.0 > 1.0 - 1.0e-3
    }
}

/// Classify clustered UV hits into exactly two trim events.
#[allow(clippy::too_many_arguments)]
pub(super) fn classify_contacts(
    brep: &BRep,
    face: FaceHandle,
    arc_eval: &dyn CurveEval,
    arc_uv: &[(f64, f64, f64)],
    edge_uvs: &std::collections::BTreeMap<EdgeHandle, Vec<(f64, f64, f64)>>,
    mut clusters: Vec<ContactCluster>,
    trace: bool,
) -> Result<[(f64, EdgeHandle, Point3); 2], String> {
    let chord = arc_eval.point_at(0.0).distance_to(&arc_eval.point_at(1.0));
    // Same-edge graze merge across cluster gaps: an arc that
    // TERMINATES on a rail approaches it near-tangentially, and the
    // sampled polyline reads the tangent wiggle as an extra
    // transversal crossing beyond the adjacency band (measured on
    // the N=3 cube corner: exact endpoint contacts at 1e-17 next to
    // same-edge grazes at t≈0.03/0.97 with 2.1e-6 margins).
    {
        let graze_margin = (1.0e-3 * chord).max(1.0e-9);
        loop {
            let mut drop_index = None;
            'search: for i in 0..clusters.len() {
                for j in (i + 1)..clusters.len() {
                    if clusters[i].rep.1 == clusters[j].rep.1
                        && (clusters[i].rep.0 - clusters[j].rep.0).abs() <= 0.1
                        && clusters[i].rep.3.min(clusters[j].rep.3) <= graze_margin
                    {
                        drop_index = Some(if clusters[i].rep.3 > clusters[j].rep.3 {
                            i
                        } else {
                            j
                        });
                        break 'search;
                    }
                }
            }
            match drop_index {
                Some(index) => {
                    clusters.remove(index);
                }
                None => break,
            }
        }
    }
    // Dedup near-identical arc parameters (a crossing at a shared
    // vertex between consecutive edges registers on both).
    clusters.dedup_by(|a, b| (a.rep.0 - b.rep.0).abs() < 1.0e-3);
    // Exact-corner endpoint contact — synthesized WHENEVER missing,
    // not only on an under-count: with uniform radii the trim
    // intersections are exact (measured 4e-8 corner gaps on
    // pentagon5_uniform), so the arc ENDS ON the strip boundary
    // instead of crossing it and there is no sign flip for the
    // polyline test to find.
    {
        let band = (1.0e-3 * chord).max(1.0e-9);
        for t_missing in [0.0, 1.0] {
            if clusters
                .iter()
                .any(|c| (c.rep.0 - t_missing).abs() < 1.0e-3)
            {
                continue;
            }
            let p_end = arc_eval.point_at(t_missing);
            let mut best: Option<(f64, EdgeHandle, Point3)> = None;
            for he in face_loop_half_edges(brep, face)? {
                let edge = brep.half_edges[he].edge;
                let Some(curve) = brep.edges[edge].curve.as_ref() else {
                    continue;
                };
                let (e0, e1) = (curve.t_start, curve.t_end);
                if !(e0.is_finite() && e1.is_finite()) {
                    continue;
                }
                let edge_eval: &dyn CurveEval = &curve.geometry;
                let mut prev = edge_eval.point_at(e0);
                for i in 1..=64 {
                    let q = edge_eval.point_at(e0 + (e1 - e0) * i as f64 / 64.0);
                    let (d, foot) = point_segment_distance(p_end, prev, q);
                    if best.as_ref().map(|(bd, _, _)| d < *bd).unwrap_or(true) {
                        best = Some((d, edge, foot));
                    }
                    prev = q;
                }
            }
            if let Some((d, edge, point)) = best
                && d <= band
            {
                if trace {
                    eprintln!(
                        "    [crossings] endpoint contact at arc t {t_missing:.1} \
                         d {d:.3e} — synthesized (exact-corner class)"
                    );
                }
                clusters.push(ContactCluster {
                    rep: (t_missing, edge, point, d),
                    t_first: t_missing,
                    t_last: t_missing,
                    count: 1,
                });
            }
        }
        clusters.sort_by(|a, b| {
            a.rep
                .0
                .partial_cmp(&b.rep.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    // Transversality: every non-endpoint candidate must actually
    // LEAVE its edge in UV. Signed side extents of the arc polyline
    // around the candidate, measured against the edge polyline: a
    // boundary-COINCIDENT pseudo-crossing never leaves projection
    // noise on EITHER side (measured ±1e-16 of the chart on the
    // squat-box runs), while every genuine contact reaches finite
    // extent on at least one side — including the blunt-apex class
    // where the arc TERMINATES just past the boundary (measured
    // t=0.016/0.984 crossings with extents +2.6e-14/−4.5e-2: one
    // side only, and demanding both sides declined the pentagon5/
    // hexagon6 blunt fans that have served since PK-187), so the
    // one-sided allowance is confined to the arc's terminal band.
    // Endpoint contacts are the arc's designed corner landings and
    // are exempt.
    {
        let diam = {
            let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            );
            for &(u, v, _) in arc_uv {
                lo_u = lo_u.min(u);
                hi_u = hi_u.max(u);
                lo_v = lo_v.min(v);
                hi_v = hi_v.max(v);
            }
            ((hi_u - lo_u).powi(2) + (hi_v - lo_v).powi(2)).sqrt()
        };
        let extent_gate = (1.0e-4 * diam).max(1.0e-12);
        clusters.retain(|c| {
            if c.endpoint() {
                return true;
            }
            let Some(edge_poly) = edge_uvs.get(&c.rep.1) else {
                return true;
            };
            let (u0, v0, _) = edge_poly[0];
            let (u1, v1, _) = edge_poly[edge_poly.len() - 1];
            let len = ((u1 - u0).powi(2) + (v1 - v0).powi(2)).sqrt();
            if len < 1.0e-12 {
                return true;
            }
            let normal = ((v0 - v1) / len, (u1 - u0) / len);
            let t_lo = (c.t_first - CLUSTER_ADJACENCY).max(0.0);
            let t_hi = (c.t_last + CLUSTER_ADJACENCY).min(1.0);
            let (mut pos, mut neg) = (0.0_f64, 0.0_f64);
            for &(au, av, s) in arc_uv {
                if s < t_lo || s > t_hi {
                    continue;
                }
                let mut nearest: Option<(f64, (f64, f64))> = None;
                for w in edge_poly.windows(2) {
                    let (d, foot) = uv_segment_foot((au, av), (w[0].0, w[0].1), (w[1].0, w[1].1));
                    if nearest.map(|(bd, _)| d < bd).unwrap_or(true) {
                        nearest = Some((d, foot));
                    }
                }
                let Some((_, foot)) = nearest else { continue };
                let side = (au - foot.0) * normal.0 + (av - foot.1) * normal.1;
                pos = pos.max(side);
                neg = neg.max(-side);
            }
            // One side suffices only where the arc TERMINATES (the
            // blunt-apex landings at t 0.016/0.984). Mid-arc, a contact
            // must genuinely cross: finite extent on BOTH sides, at the
            // scale a transversal crossing reaches over the cluster
            // window. The S-shaped mixed-radius arc grazes a strip
            // boundary (t 0.218, extents ±2-4e-4 of the chart) and read as
            // a third crossing, declining the carve; genuine mid-arc
            // crossings reach ±7e-2 (PK-591).
            let terminal = c.t_first < CLUSTER_ADJACENCY || c.t_last > 1.0 - CLUSTER_ADJACENCY;
            let transversal = if terminal {
                pos >= extent_gate || neg >= extent_gate
            } else {
                let crossing_gate = MID_ARC_CROSSING_EXTENT * diam;
                pos >= crossing_gate && neg >= crossing_gate
            };
            if !transversal && trace {
                eprintln!(
                    "    [crossings] dropped non-transversal contact on {:?}: t {:.3}..{:.3} \
                     ({} hits, min margin {:.3e}, side extents +{pos:.3e}/-{neg:.3e} < {extent_gate:.3e})",
                    c.rep.1, c.t_first, c.t_last, c.count, c.rep.3
                );
            }
            transversal
        });
    }
    if clusters.len() != 2 {
        return Err(format!(
            "expected 2 {}, found {} ({:?})",
            crate::construct::setback_decline::ARC_BOUNDARY_CROSSINGS,
            clusters.len(),
            clusters
                .iter()
                .map(|c| (c.rep.0, c.rep.3, c.count))
                .collect::<Vec<_>>()
        ));
    }
    Ok([
        (clusters[0].rep.0, clusters[0].rep.1, clusters[0].rep.2),
        (clusters[1].rep.0, clusters[1].rep.1, clusters[1].rep.2),
    ])
}

/// Min distance from UV point `p` to segment `a..b`, with the foot.
fn uv_segment_foot(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> (f64, (f64, f64)) {
    let ab = (b.0 - a.0, b.1 - a.1);
    let len2 = ab.0 * ab.0 + ab.1 * ab.1;
    let s = if len2 > 0.0 {
        (((p.0 - a.0) * ab.0 + (p.1 - a.1) * ab.1) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let foot = (a.0 + s * ab.0, a.1 + s * ab.1);
    let d = ((p.0 - foot.0).powi(2) + (p.1 - foot.1).powi(2)).sqrt();
    (d, foot)
}

/// Min distance from `p` to 3-D segment `a..b`, with the foot point.
fn point_segment_distance(p: Point3, a: Point3, b: Point3) -> (f64, Point3) {
    let ab = b - a;
    let len2 = ab.dot(&ab);
    let s = if len2 > 0.0 {
        ((p - a).dot(&ab) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let foot = Point3::new(a.x + s * ab.x, a.y + s * ab.y, a.z + s * ab.z);
    (p.distance_to(&foot), foot)
}
