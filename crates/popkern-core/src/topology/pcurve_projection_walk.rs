//! Anchored UV projection walks for dense pcurve recomputation.
//!
//! Projecting a 3D edge curve into a surface's UV space by continuation
//! must pick a starting preimage, and no single anchor is safe:
//!
//! * Anchoring at an ENDPOINT is branch-ambiguous at degenerate patch
//!   corners — where two boundary curves of a patch meet at one 3D
//!   point (e.g. the rim corner of an apex-cap sector, whose exact
//!   preimages are both (1,0) and (0,1)), the global projection may
//!   return the preimage belonging to the OTHER boundary, and
//!   continuation then walks the entire pcurve along the wrong side of
//!   the patch (measured on pyramid4_mixed: 1.8e-2 cap/support gap
//!   from a spoke pcurve hugging the rim).
//! * Anchoring at the MIDPOINT resolves that case (interior curve
//!   points generically have a unique nearest preimage) but trusts a
//!   single global projection that can land on a folded or extended
//!   region of shaped cap patches (measured on cube_mixed: 5.0e-1).
//!
//! So both walks run and the result with the smaller measured 3D
//! deviation from the sampled curve wins. The endpoint-anchored walk
//! is always a candidate, so surfaces the historical walk handled
//! correctly cannot regress.

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};

/// Sample `curve_on_edge` at the curve `fractions` of `[t_lo, t_hi]`
/// (ascending from 0 to 1; the end samples pinned to the half-edge's
/// vertex positions) and project each sample into
/// `surface` UV space with anchored continuation and seam unwrapping.
/// Runs an endpoint-anchored and a midpoint-anchored walk and returns
/// the UV points (z = 0, curve order) of whichever walk's surface
/// image deviates least from the 3D samples. Samples whose projection
/// fails are dropped (PK-095 reseeding included).
pub(super) fn anchored_uv_samples(
    surface: &SurfaceGeometry,
    curve_on_edge: &CurveOnEdge,
    (t_lo, t_hi): (f64, f64),
    origin_pos: Point3,
    dest_pos: Point3,
    endpoint_uv_seeds: ((f64, f64), (f64, f64)),
    fractions: &[f64],
) -> Vec<Point3> {
    let n_samples = fractions.len();
    if n_samples < 3 {
        return Vec::new();
    }
    let mut points_3d = Vec::with_capacity(n_samples);
    for (i, &frac) in fractions.iter().enumerate() {
        let t = t_lo + (t_hi - t_lo) * frac;
        let p = if i == 0 {
            origin_pos
        } else if i + 1 == n_samples {
            dest_pos
        } else {
            curve_on_edge.geometry.point_at(t)
        };
        points_3d.push(p);
    }

    let start_walk = anchored_walk(surface, &points_3d, fractions, endpoint_uv_seeds, 0);
    let mid_walk = anchored_walk(
        surface,
        &points_3d,
        fractions,
        endpoint_uv_seeds,
        n_samples / 2,
    );
    let start_deviation = walk_deviation(surface, &points_3d, &start_walk);
    let mid_deviation = walk_deviation(surface, &points_3d, &mid_walk);
    let winner = if mid_deviation < start_deviation {
        mid_walk
    } else {
        start_walk
    };
    winner
        .into_iter()
        .flatten()
        .map(|(u, v)| Point3::new(u, v, 0.0))
        .collect()
}

/// One continuation walk anchored at `anchor_index`: the anchor is
/// projected globally (endpoint-lerp seed as the PK-095 divergence
/// fallback), then the walk continues outward to both ends, each step
/// seeded from the last successful projection on its side and
/// seam-unwrapped against it (OCCT style: half-period threshold).
fn anchored_walk(
    surface: &SurfaceGeometry,
    points_3d: &[Point3],
    fractions: &[f64],
    endpoint_uv_seeds: ((f64, f64), (f64, f64)),
    anchor_index: usize,
) -> Vec<Option<(f64, f64)>> {
    let n_samples = points_3d.len();
    let ((u_s, v_s), (u_e, v_e)) = endpoint_uv_seeds;
    let seed_at = |i: usize| {
        let frac = fractions[i];
        (u_s + (u_e - u_s) * frac, v_s + (v_e - v_s) * frac)
    };

    // Analytic inversion for cylinder charts: the generic Newton
    // projection wanders across period branches and, on rulings,
    // converges to a wrong nearby location (measured: a strip-end
    // ruling's 23-sample cubic image 1.84 off its own line on
    // offset_bracket while the line lies ON the cylinder to machine
    // zero — the walk authored it and the survival rule then kept it
    // forever, since every re-walk reproduced the same wandering).
    // The same cure the boundary-first mesher uses; the BRANCH comes
    // from the endpoint seeds (which carry the face's chart window),
    // never from the principal value — a principal-branch candidate
    // was measured to break rotated charts and is annotated as
    // falsified in pcurve_recompute_analytic's history.
    let cylinder_invert = |p: &Point3| -> Option<(f64, f64)> {
        let SurfaceGeometry::Cylinder { origin, axis, .. } = surface else {
            return None;
        };
        let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(*axis);
        let d = *p - *origin;
        Some((d.dot(&perp2).atan2(d.dot(&perp1)), d.dot(axis)))
    };
    let is_cylinder = matches!(surface, SurfaceGeometry::Cylinder { .. });

    let mut uv_by_index: Vec<Option<(f64, f64)>> = vec![None; n_samples];
    let anchor = if is_cylinder {
        cylinder_invert(&points_3d[anchor_index]).map(|(mut u, v)| {
            // Unwrap the anchor into the branch the endpoint seeds
            // carry (the face's chart window).
            let (su, _sv) = seed_at(anchor_index);
            if let (Some(pu), _) = surface.periodicity() {
                while u - su > pu * 0.5 {
                    u -= pu;
                }
                while u - su < -pu * 0.5 {
                    u += pu;
                }
            }
            (u, v)
        })
    } else {
        project_point_to_uv(surface, &points_3d[anchor_index]).or_else(|| {
            let (su, sv) = seed_at(anchor_index);
            project_point_to_uv_near(surface, &points_3d[anchor_index], su, sv)
        })
    };
    let Some(anchor_uv) = anchor else {
        return uv_by_index;
    };
    uv_by_index[anchor_index] = Some(anchor_uv);

    let (period_u, period_v) = surface.periodicity();
    let mut walk = |indices: &mut dyn Iterator<Item = usize>| {
        let mut prev_uv = anchor_uv;
        for i in indices {
            let p = &points_3d[i];
            let (su, sv) = seed_at(i);
            // Continuation prevents branch jumping on multi-valued
            // inverse maps; the reseed fallbacks keep sample loss from
            // collapsing the pcurve to a straight 2-point line
            // (PK-095: 5e-2 UV error on fillet strip trims).
            let projected = if is_cylinder {
                cylinder_invert(p)
            } else {
                project_point_to_uv_near(surface, p, prev_uv.0, prev_uv.1)
                    .or_else(|| project_point_to_uv(surface, p))
                    .or_else(|| project_point_to_uv_near(surface, p, su, sv))
            };
            let Some((mut u, mut v)) = projected else {
                continue;
            };
            if let Some(pu) = period_u {
                while u - prev_uv.0 > pu * 0.5 {
                    u -= pu;
                }
                while u - prev_uv.0 < -pu * 0.5 {
                    u += pu;
                }
            }
            if let Some(pv) = period_v {
                while v - prev_uv.1 > pv * 0.5 {
                    v -= pv;
                }
                while v - prev_uv.1 < -pv * 0.5 {
                    v += pv;
                }
            }
            prev_uv = (u, v);
            uv_by_index[i] = Some((u, v));
        }
    };
    if anchor_index > 0 {
        walk(&mut (0..anchor_index).rev());
    }
    if anchor_index + 1 < n_samples {
        walk(&mut ((anchor_index + 1)..n_samples));
    }
    uv_by_index
}

/// Max 3D distance between each projected sample's surface image and
/// the curve sample it came from. Dropped samples count as the worst
/// deviation seen for a fully-projected walk would be optimistic —
/// they are penalized at infinity so a walk that keeps its samples on
/// the surface always beats one that loses them.
fn walk_deviation(
    surface: &SurfaceGeometry,
    points_3d: &[Point3],
    walk: &[Option<(f64, f64)>],
) -> f64 {
    let mut worst: f64 = 0.0;
    for (p, uv) in points_3d.iter().zip(walk.iter()) {
        match uv {
            Some((u, v)) => {
                worst = worst.max(p.distance_to(&surface.point_at(*u, *v)));
            }
            None => return f64::INFINITY,
        }
    }
    worst
}
