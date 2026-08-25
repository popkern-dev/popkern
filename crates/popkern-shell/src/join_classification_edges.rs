//! Per-edge measured classification for join planning (PK-325): the
//! station-sampled wedge test and the measured into-face direction.

use super::EdgeJoinClass;
use crate::curved_recompute::project_to_surface;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::STATION_FRACTIONS;

/// Classify one edge from measured displacements. Returns the class and
/// the measured ball radius (mean displacement magnitude).
#[allow(clippy::items_after_statements)]
pub(super) fn classify_edge(
    skin: &BRep,
    eh: EdgeHandle,
    off1: &SurfaceGeometry,
    off2: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Option<(EdgeJoinClass, f64)> {
    let trace = std::env::var_os("PK_JOIN_TRACE").is_some();
    let edge = skin.edges.get(eh)?;
    let Some(curve) = edge.curve.as_ref() else {
        if trace {
            eprintln!("[join] {eh:?}: no curve");
        }
        return None;
    };
    let he_handle = edge.half_edge;
    let he = skin.half_edges.get(he_handle)?;
    let twin_handle = he.twin;

    let smooth_gate = tolerance.linear() * 0.5;
    let mut votes: Vec<EdgeJoinClass> = Vec::new();
    let mut radii: Vec<f64> = Vec::new();

    for f in STATION_FRACTIONS {
        let t = curve.t_start + f * (curve.t_end - curve.t_start);
        let p = curve.geometry.point_at(t);
        let (Some(q1), Some(q2)) = (project_to_surface(off1, &p), project_to_surface(off2, &p))
        else {
            if trace {
                eprintln!("[join] {eh:?}: projection failed at f={f}");
            }
            return None; // projection ambiguity: leave to the ladder
        };
        let d1v = q1 - p;
        let d2v = q2 - p;
        let d1 = d1v.length();
        let d2 = d2v.length();
        if (d1v - d2v).length() <= smooth_gate {
            votes.push(EdgeJoinClass::Smooth);
            radii.push(0.5 * (d1 + d2));
            continue;
        }
        if d1 <= smooth_gate || d2 <= smooth_gate {
            // One side barely moves: no meaningful wedge model.
            votes.push(EdgeJoinClass::Resolvable);
            radii.push(0.5 * (d1 + d2));
            continue;
        }
        // Resolution floor: a rounding radius within ~100x tolerance is
        // indistinguishable from the sharp join at model resolution —
        // constructing micro-strips there only manufactures degenerate
        // faces whose orientation cannot even be measured. The sharp
        // intersection IS the Minkowski result to tolerance.
        if d1.min(d2) <= tolerance.linear() * 100.0 {
            votes.push(EdgeJoinClass::Resolvable);
            radii.push(0.5 * (d1 + d2));
            continue;
        }
        let n1 = d1v * (1.0 / d1);
        let n2 = d2v * (1.0 / d2);
        let g = n1.dot(&n2);
        let denom = 1.0 - g * g;
        if denom < 1.0e-12 {
            // Parallel displacement directions with different
            // magnitudes: a tangential step, not a wedge.
            votes.push(EdgeJoinClass::Irregular);
            radii.push(0.5 * (d1 + d2));
            continue;
        }
        let a = (d1 - g * d2) / denom;
        let b = (d2 - g * d1) / denom;
        let x_star = p + n1 * a + n2 * b;

        let (Some(w1), Some(w2)) = (
            into_face_direction(skin, he_handle, t),
            into_face_direction(skin, twin_handle, t),
        ) else {
            if trace {
                eprintln!("[join] {eh:?}: into-face failed at f={f}");
            }
            return None;
        };
        let eps = (0.02 * d1.min(d2)).max(tolerance.linear());
        let s1 = (x_star - q1).dot(&w1);
        let s2 = (x_star - q2).dot(&w2);
        if trace {
            eprintln!(
                "[join] {eh:?} f={f}: s1={s1:.3e} s2={s2:.3e} n1={n1:?} n2={n2:?} w1={w1:?} w2={w2:?} xstar={x_star:?} q1={q1:?} q2={q2:?}"
            );
        }
        if s1 < -eps && s2 < -eps {
            // A single rolling ball must touch BOTH offset surfaces:
            // unequal per-side displacement (variable per-face
            // distance) has no circular fill (design §10 E4).
            let side_gate = (0.5 * (d1 + d2) * 1.0e-3).max(tolerance.linear() * 100.0);
            if (d1 - d2).abs() > side_gate {
                votes.push(EdgeJoinClass::Irregular);
            } else {
                votes.push(EdgeJoinClass::Gap);
            }
        } else {
            votes.push(EdgeJoinClass::Resolvable);
        }
        radii.push(0.5 * (d1 + d2));
    }

    let first = votes.first().copied()?;
    let uniform = votes.iter().all(|v| *v == first);
    let mean_r = radii.iter().sum::<f64>() / radii.len() as f64;
    if !uniform || first == EdgeJoinClass::Irregular {
        return Some((EdgeJoinClass::Irregular, mean_r));
    }
    if first == EdgeJoinClass::Gap {
        // A single rolling ball must fit the whole edge: displacement
        // magnitudes must agree across sides and stations.
        let spread = radii
            .iter()
            .fold(0.0f64, |acc, r| acc.max((r - mean_r).abs()));
        let gate = (mean_r * 1.0e-3).max(tolerance.linear() * 100.0);
        if spread > gate {
            return Some((EdgeJoinClass::Irregular, mean_r));
        }
    }
    Some((first, mean_r))
}

/// Unit direction pointing from the edge into the face that owns
/// `he_handle`'s loop, measured at edge-curve parameter `t`.
///
/// The interior side of a loop traversal is MEASURED, never assumed:
/// the face interior is LEFT of the traversal iff the loop winds CCW
/// in UV (outer loops) or CW (inner loops — the interior is outside
/// the hole). Pole-touching charts can be legally wound against their
/// orientation flag (the cone lateral is flag-true with a measured CW
/// loop; the winding audit defers such charts to the tangency tier),
/// so the flag is not a winding witness here.
fn into_face_direction(skin: &BRep, he_handle: HalfEdgeHandle, t: f64) -> Option<Vec3> {
    let he = skin.half_edges.get(he_handle)?;
    let pcurve = he.pcurve.as_ref()?;
    let loop_handle = he.loop_;
    let face = skin.loops.get(loop_handle)?.face;
    let sof = skin.faces.get(face)?.surface.as_ref()?;
    let surface = &sof.geometry;
    let interior_side = loop_interior_side(skin, loop_handle, face, sof)?;

    // The pcurve parameter runs with the half-edge; the edge curve
    // parameter may run against it. Map `t` into the pcurve range by a
    // POSITIONAL witness: evaluate both direction candidates through
    // the surface and keep the one that lands on the station point.
    // (Endpoint matching is vacuous on closed edges — both curve ends
    // are the same vertex — and rim circles are exactly the gap edges
    // this classifier must serve.)
    let edge = skin.edges.get(he.edge)?;
    let curve = edge.curve.as_ref()?;
    let p_station = curve.geometry.point_at(t);
    let frac = (t - curve.t_start) / (curve.t_end - curve.t_start);
    let candidate = |pc_frac: f64| -> (f64, f64) {
        let pt = pcurve.t_start + pc_frac * (pcurve.t_end - pcurve.t_start);
        let uv = pcurve.point_at(pt);
        let q = surface.point_at(uv.x, uv.y);
        (pt, q.distance_to(&p_station))
    };
    let (pt_fwd, d_fwd) = candidate(frac);
    let (pt_rev, d_rev) = candidate(1.0 - frac);
    let pt = if d_fwd <= d_rev { pt_fwd } else { pt_rev };

    let h = (pcurve.t_end - pcurve.t_start).abs().max(1.0e-9) * 1.0e-6;
    let uv0 = pcurve.point_at(pt - h);
    let uv1 = pcurve.point_at(pt + h);
    let du = uv1.x - uv0.x;
    let dv = uv1.y - uv0.y;
    let len = (du * du + dv * dv).sqrt();
    if len < 1.0e-30 {
        return None;
    }
    // Left of the UV direction (CCW loops), right for CW loops.
    let (iu, iv) = (-dv / len * interior_side, du / len * interior_side);
    let uv = pcurve.point_at(pt);
    let su = surface.partial_u(uv.x, uv.y);
    let sv = surface.partial_v(uv.x, uv.y);
    let w = su * iu + sv * iv;
    let wl = w.length();
    if wl < 1.0e-30 {
        return None;
    }
    Some(w * (1.0 / wl))
}

/// Measured interior side for a loop traversal: `+1.0` when the face
/// interior is LEFT of the traversal, `-1.0` when RIGHT.
///
/// Ladder: core's period-unwrapped, pole-guarded loop area first; a
/// local pcurve shoelace for the charts the guards decline
/// (pole-touching rectangles like the full cone lateral, where the
/// seam edges cancel and the naive area is still sign-correct); the
/// winding-honesty flag convention as the last resort.
fn loop_interior_side(
    skin: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
    face: FaceHandle,
    sof: &popkern_core::geometry::surface::SurfaceOnFace,
) -> Option<f64> {
    let is_outer = skin.faces.get(face)?.outer_loop == loop_handle;
    let outer_sign = if is_outer { 1.0 } else { -1.0 };
    if let Ok(area) = skin.loop_signed_uv_area(loop_handle, &sof.geometry)
        && area.signed.abs() > area.perimeter.powi(2) * 1.0e-12
    {
        return Some(area.signed.signum() * outer_sign);
    }
    if let Some(area) = naive_loop_shoelace(skin, loop_handle)
        && area.0.abs() > area.1.powi(2) * 1.0e-12
    {
        return Some(area.0.signum() * outer_sign);
    }
    Some(if sof.orientation { 1.0 } else { -1.0 } * outer_sign)
}

/// Local (non-unwrapped) pcurve shoelace: `(signed area, perimeter)`.
fn naive_loop_shoelace(
    skin: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<(f64, f64)> {
    let start = skin.loops.get(loop_handle)?.half_edge;
    let mut poly: Vec<(f64, f64)> = Vec::new();
    let mut cur = start;
    for _ in 0..skin.half_edges.len() + 1 {
        let he = skin.half_edges.get(cur)?;
        let pcurve = he.pcurve.as_ref()?;
        for i in 0..16 {
            let t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * i as f64 / 16.0;
            let uv = pcurve.point_at(t);
            poly.push((uv.x, uv.y));
        }
        cur = he.next;
        if cur == start {
            break;
        }
    }
    if poly.len() < 3 {
        return None;
    }
    let mut area = 0.0;
    let mut perimeter = 0.0;
    for i in 0..poly.len() {
        let (x0, y0) = poly[i];
        let (x1, y1) = poly[(i + 1) % poly.len()];
        area += x0 * y1 - x1 * y0;
        perimeter += ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
    }
    Some((0.5 * area, perimeter))
}
