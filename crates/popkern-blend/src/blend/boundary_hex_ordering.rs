//! Ring-ordering candidates for `reorient_arcs_for_hex_closure` —
//! the successor walk, the angular sort, and the bridge-validity
//! score that arbitrates between them at valence >= 4. Split from
//! `boundary_hex.rs` (500-line structural cap).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::FaceHandle;

use crate::blend::surface_bridge::point_on_surface_distance;

/// Bridge validity score for a candidate ring ordering: the number of
/// consecutive end→start bridges with NO common support face within
/// `tol_slack` (primary), then the summed bridge length (tiebreak).
/// Lower is better under tuple ordering.
pub(super) fn ring_bridge_score(
    candidate: &(Vec<usize>, Vec<bool>, Vec<NurbsCurve>),
    fillet_endpoints: &[[Point3; 2]],
    face_supports: &[(FaceHandle, SurfaceGeometry)],
    tol_slack: f64,
) -> (usize, f64) {
    let (order, flip, _) = candidate;
    let n = order.len();
    let mut unsupported = 0usize;
    let mut total_len = 0.0f64;
    for k in 0..n {
        let idx = order[k];
        let end = if flip[k] {
            fillet_endpoints[idx][0]
        } else {
            fillet_endpoints[idx][1]
        };
        let next = (k + 1) % n;
        let nidx = order[next];
        let start = if flip[next] {
            fillet_endpoints[nidx][1]
        } else {
            fillet_endpoints[nidx][0]
        };
        total_len += end.distance_to(&start);
        let supported = face_supports.iter().any(|(_, surface)| {
            point_on_surface_distance(&end, surface) <= tol_slack
                && point_on_surface_distance(&start, surface) <= tol_slack
        });
        if !supported {
            unsupported += 1;
        }
    }
    (unsupported, total_len)
}

/// The original successor walk: chain arcs end-to-start by shared
/// support faces (strict then slack tolerance), with a nearest-
/// endpoint fallback for saddle trims that sit off every support.
pub(super) fn successor_walk_ordering(
    fillet_arcs: &[NurbsCurve],
    fillet_endpoints: &[[Point3; 2]],
    face_supports: &[(FaceHandle, SurfaceGeometry)],
    tol_linear: f64,
) -> Option<(Vec<usize>, Vec<bool>, Vec<NurbsCurve>)> {
    let n = fillet_arcs.len();
    let tol_slack = tol_linear.max(1e-5) * 10.0;
    let mut order = vec![0usize; n];
    let mut flip = vec![false; n];
    let mut used = vec![false; n];
    order[0] = 0;
    used[0] = true;

    for k in 1..n {
        let prev_idx = order[k - 1];
        let prev_endpoints = if flip[k - 1] {
            [fillet_endpoints[prev_idx][1], fillet_endpoints[prev_idx][0]]
        } else {
            fillet_endpoints[prev_idx]
        };
        let prev_end = prev_endpoints[1];

        // Pass 1: strict tolerance — well-formed BReps where both
        // endpoints lie exactly on a shared support face.
        // Pass 2: slack tolerance — saddle / mixed-convexity cases
        // where the boolean result leaves trim endpoints a few
        // multiples of `tol_linear` off their logical face.
        // Pass 3: nearest-endpoint fallback — §5.3 item #5 saddle
        // path. When no face-support shares both endpoints within
        // tolerance (a notched cube's trim endpoints may sit on
        // the fillet surface, not the support face), pick the
        // candidate whose start (or end, post-flip) minimises the
        // straight-line distance to `prev_end`. Hex-adjacency is
        // a topological property; geometrically-close endpoints
        // are correct up to a ridge discontinuity the bridge
        // builder handles via its straight-line fallback.
        let mut best: Option<(usize, bool)> = None;
        for &tol_try in &[tol_linear, tol_slack] {
            if best.is_some() {
                break;
            }
            'search: for (cand_idx, ep) in fillet_endpoints.iter().enumerate() {
                if used[cand_idx] {
                    continue;
                }
                for cand_flip in [false, true] {
                    let cand_start = if cand_flip { ep[1] } else { ep[0] };
                    for (_, surface) in face_supports {
                        let d1 = point_on_surface_distance(&prev_end, surface);
                        let d2 = point_on_surface_distance(&cand_start, surface);
                        if d1 <= tol_try && d2 <= tol_try {
                            best = Some((cand_idx, cand_flip));
                            break 'search;
                        }
                    }
                }
            }
        }
        if best.is_none() {
            // §5.3 item #5 saddle fallback: nearest-endpoint
            // matching. On mixed-convexity BReps the trim endpoints
            // can sit on the fillet surface with ~10⁻² offsets from
            // any support face, so the face-adjacency test above
            // rejects every pair. Hex-adjacency is ultimately a
            // topological property — picking the candidate whose
            // start is geometrically closest to `prev_end`
            // preserves the correct cycle order and lets the
            // bridge builder fill in a straight-line segment on
            // the best-effort face.
            let mut nearest: Option<(f64, usize, bool)> = None;
            for (cand_idx, ep) in fillet_endpoints.iter().enumerate() {
                if used[cand_idx] {
                    continue;
                }
                for cand_flip in [false, true] {
                    let cand_start = if cand_flip { ep[1] } else { ep[0] };
                    let d = prev_end.distance_to(&cand_start);
                    if nearest.is_none_or(|(best_d, _, _)| d < best_d) {
                        nearest = Some((d, cand_idx, cand_flip));
                    }
                }
            }
            if let Some((_, idx, fl)) = nearest {
                best = Some((idx, fl));
            }
        }
        let (next_idx, next_flip) = best?;
        order[k] = next_idx;
        flip[k] = next_flip;
        used[next_idx] = true;
    }

    let oriented_arcs: Vec<NurbsCurve> = (0..n)
        .map(|k| {
            let arc = &fillet_arcs[order[k]];
            if flip[k] {
                reverse_nurbs_curve(arc)
            } else {
                arc.clone()
            }
        })
        .collect();

    Some((order, flip, oriented_arcs))
}

/// Angular ordering for valence >= 4: sort arcs by the angle of
/// their midpoint direction about the endpoint centroid's best-fit
/// axis, and orient each arc counterclockwise in that frame.
pub(super) fn angular_arc_ordering(
    fillet_arcs: &[NurbsCurve],
    fillet_endpoints: &[[Point3; 2]],
) -> Option<(Vec<usize>, Vec<bool>, Vec<NurbsCurve>)> {
    use popkern_core::geometry::curve::CurveEval;
    let n = fillet_arcs.len();
    // Centroid of all endpoints.
    let mut c = [0.0_f64; 3];
    for ep in fillet_endpoints {
        for p in ep {
            c[0] += p.x;
            c[1] += p.y;
            c[2] += p.z;
        }
    }
    let count = (2 * n) as f64;
    let centroid = Point3::new(c[0] / count, c[1] / count, c[2] / count);
    // Best-fit axis from summed cross products of endpoint spokes.
    // Each arc's [a, b] endpoint order is ARBITRARY, so the raw
    // per-arc cross product carries an arbitrary sign — a raw sum
    // partially cancels and can tip the axis sideways, folding the
    // angular sort into a non-adjacent order (measured on
    // pyramid4_mixed: order [1,0,2,3] paired diagonal fillets, the
    // FaceEdge bridges grabbed a diagonal face, and every sphere-fit
    // read amplitude 0.56 where its adjacent faces give 1.0).
    // Sign-align every term against the largest-magnitude one before
    // summing; the overall axis SIGN stays irrelevant (a reversed
    // axis reverses the cycle, which is still an adjacency).
    let spokes: Vec<popkern_core::geometry::Vec3> = fillet_endpoints
        .iter()
        .map(|ep| {
            let a = ep[0] - centroid;
            let b = ep[1] - centroid;
            a.cross(&b)
        })
        .collect();
    let reference = spokes.iter().copied().max_by(|x, y| {
        x.length()
            .partial_cmp(&y.length())
            .unwrap_or(std::cmp::Ordering::Equal)
    })?;
    let mut axis = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
    for x in &spokes {
        if x.dot(&reference) < 0.0 {
            axis = axis - *x;
        } else {
            axis = axis + *x;
        }
    }
    let axis = axis.normalized()?;
    // In-plane frame.
    let seed = fillet_endpoints[0][0] - centroid;
    let u0 = (seed - axis * seed.dot(&axis)).normalized()?;
    let v0 = axis.cross(&u0);
    let angle_of = |p: Point3| -> f64 {
        let d = p - centroid;
        d.dot(&v0).atan2(d.dot(&u0))
    };
    let mut keyed: Vec<(f64, usize)> = (0..n)
        .map(|k| {
            let (t0, t1) = fillet_arcs[k].parameter_range();
            let mid = fillet_arcs[k].point_at(0.5 * (t0 + t1));
            (angle_of(mid), k)
        })
        .collect();
    keyed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let order: Vec<usize> = keyed.iter().map(|&(_, k)| k).collect();
    if popkern_core::trace::trace_enabled("boundary-hex") {
        eprintln!(
            "[boundary-hex] centroid=({:.4},{:.4},{:.4}) axis=({:.4},{:.4},{:.4}) order={order:?}",
            centroid.x, centroid.y, centroid.z, axis.x, axis.y, axis.z
        );
        for &(angle, k) in &keyed {
            let (t0, t1) = fillet_arcs[k].parameter_range();
            let mid = fillet_arcs[k].point_at(0.5 * (t0 + t1));
            eprintln!(
                "[boundary-hex]   arc {k}: mid=({:.4},{:.4},{:.4}) angle={angle:.4} \
                 ep0=({:.4},{:.4},{:.4}) ep1=({:.4},{:.4},{:.4})",
                mid.x,
                mid.y,
                mid.z,
                fillet_endpoints[k][0].x,
                fillet_endpoints[k][0].y,
                fillet_endpoints[k][0].z,
                fillet_endpoints[k][1].x,
                fillet_endpoints[k][1].y,
                fillet_endpoints[k][1].z
            );
        }
    }
    // Orient each arc counterclockwise: the signed angular step from
    // start to end (wrapped to (-pi, pi]) must be positive.
    let mut flip = vec![false; n];
    for (slot, &k) in order.iter().enumerate() {
        let a0 = angle_of(fillet_endpoints[k][0]);
        let a1 = angle_of(fillet_endpoints[k][1]);
        let mut step = a1 - a0;
        while step > core::f64::consts::PI {
            step -= core::f64::consts::TAU;
        }
        while step <= -core::f64::consts::PI {
            step += core::f64::consts::TAU;
        }
        flip[slot] = step < 0.0;
    }
    let oriented_arcs: Vec<NurbsCurve> = order
        .iter()
        .zip(&flip)
        .map(|(&k, &fl)| {
            if fl {
                reverse_nurbs_curve(&fillet_arcs[k])
            } else {
                fillet_arcs[k].clone()
            }
        })
        .collect();
    Some((order, flip, oriented_arcs))
}

/// Reverse a NURBS curve's parameter direction by reversing the
/// control-point and weight order and reflecting the knot vector.
pub(super) fn reverse_nurbs_curve(curve: &NurbsCurve) -> NurbsCurve {
    let (t0, t1) = curve.parameter_range();
    let mut new_knots: Vec<f64> = curve.knots.iter().rev().map(|&k| t0 + (t1 - k)).collect();
    for k in &mut new_knots {
        *k = k.clamp(t0, t1);
    }
    let new_cps: Vec<Point3> = curve.control_points.iter().rev().copied().collect();
    let new_weights: Vec<f64> = curve.weights.iter().rev().copied().collect();
    new_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    NurbsCurve {
        degree: curve.degree,
        knots: new_knots,
        control_points: new_cps,
        weights: new_weights,
    }
}
