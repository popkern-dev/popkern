//! Setback-corner placement between angularly adjacent fillets.
//!
//! The corner `c_k` is where the two adjacent fillet strips meet on
//! their SHARED support face: the intersection of fillet k's contact
//! curve on that face with fillet k+1's contact curve on the same
//! face. Two structural facts drive the implementation:
//!
//! * **Pair selection is coplanarity, not distance.** Each fillet
//!   carries two trim curves (one per support face) and only the
//!   shared-face pair defines the corner. On mixed-radius corners a
//!   skew pair from the two NON-shared faces can pass closer in 3-D
//!   near the apex than the shared-face pair does inside its
//!   truncated window (measured 2.8e-3 on hexagon6_variable),
//!   dragging the min-distance corner tangentially off the strips by
//!   ~0.1 rad of sector angle. Coplanar pairs (both curves in one
//!   support plane) are therefore preferred whenever any exist.
//! * **The corner stays in the trims' parameter domain.** A
//!   tangent-ray extension past range-clamped minima was measured
//!   and REJECTED: contact lines can mathematically intersect far
//!   outside the blend region (pentagon5_uniform_g2, radius 0.16 —
//!   the extended crossing sits ABOVE the pyramid apex and breaks
//!   the paper assembly from machine-zero g2 to 2.1e-3 with 12x
//!   regularity loss). In-domain window-refined corners give every
//!   measured fixture machine-zero g2 closure, including the
//!   historically fold-rejected hexagon6_variable candidate
//!   (folds 6 -> 0, g2 2.3e-5 -> 8.2e-13).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use crate::fillet_surface::FilletSurfaceResult;

/// Corner between `order[k]` and `order[k + 1]` for every k, as
/// `(gap, corner_point)` — `gap` is the residual closest-approach
/// distance (a measured quality margin the caller reports).
///
/// `vertex_position` (the original apex) drives a MATERIAL-SIDE
/// gate: adjacent trim curves of symmetric (uniform-radius) fillets
/// intersect twice near the apex — once inside the blend region and
/// once mirrored ABOVE the apex (measured on pentagon5_uniform:
/// spurious corner at z = 0.811 over an apex at z = 0.800, killing
/// the surgery's crossing detection). Only candidates on the same
/// side of the apex as the trim-curve cloud qualify.
pub(super) fn setback_corners(
    fillet_surfaces: &[FilletSurfaceResult],
    order: &[usize],
    vertex_position: Point3,
) -> Vec<(f64, Point3)> {
    let n = order.len();
    let material_dir = {
        let mut acc = [0.0_f64; 3];
        let mut count = 0.0;
        for idx in order {
            for curve in [
                &fillet_surfaces[*idx].trim_curve_1,
                &fillet_surfaces[*idx].trim_curve_2,
            ] {
                let (t0, t1) = curve.parameter_range();
                for i in 0..=4 {
                    let p = curve.point_at(t0 + (t1 - t0) * i as f64 / 4.0);
                    acc[0] += p.x;
                    acc[1] += p.y;
                    acc[2] += p.z;
                    count += 1.0;
                }
            }
        }
        Vec3::new(
            acc[0] / count - vertex_position.x,
            acc[1] / count - vertex_position.y,
            acc[2] / count - vertex_position.z,
        )
    };
    let side_band = 1.0e-3 * material_dir.length();
    let material_side = |p: Point3| -> bool {
        match material_dir.normalized() {
            Some(dir) => (p - vertex_position).dot(&dir) >= -side_band,
            None => true,
        }
    };
    let mut corners = Vec::with_capacity(n);
    for k in 0..n {
        let fa = &fillet_surfaces[order[k]];
        let fb = &fillet_surfaces[order[(k + 1) % n]];
        // Among coplanar candidates, TANGENT COHERENCE breaks ties:
        // the corner must be where the two strips meet on their
        // SHARED support (surface normals agree there — the vertex-
        // enclosure data the cap assembly consumes is coherent). On
        // a trihedral cube corner every rail pair is coplanar
        // (straight lines) and two pairs intersect EXACTLY: the
        // shared-support contact corner (normals agree) and an
        // on-edge rail crossing whose strips meet at the full 90°
        // dihedral. Gap and coplanarity score both 1e-16 — only the
        // normal test separates them (measured: the creased pick fed
        // the assembly a 3.6e-2 corner band mismatch, min regularity
        // 8.7e-5, and three folded near-pole pieces).
        let mut best_coplanar: Option<(f64, Point3, bool)> = None;
        let mut best_any: Option<(f64, Point3)> = None;
        for ca in [&fa.trim_curve_1, &fa.trim_curve_2] {
            for cb in [&fb.trim_curve_1, &fb.trim_curve_2] {
                let candidate = closest_between_refined(ca, cb, &material_side);
                if pair_is_coplanar(ca, cb) {
                    let coherent = strips_normal_disagreement(fa, fb, candidate.1)
                        .is_some_and(|sin_angle| sin_angle < 1.0e-3);
                    let better = match best_coplanar {
                        None => true,
                        Some((d, _, best_coherent)) => match (coherent, best_coherent) {
                            (true, false) => true,
                            (false, true) => false,
                            _ => candidate.0 < d,
                        },
                    };
                    if better {
                        best_coplanar = Some((candidate.0, candidate.1, coherent));
                    }
                }
                if best_any.is_none_or(|(d, _)| candidate.0 < d) {
                    best_any = Some(candidate);
                }
            }
        }
        let chosen = best_coplanar
            .map(|(d, p, _)| (d, p))
            .or(best_any)
            .unwrap_or((f64::INFINITY, Point3::new(0.0, 0.0, 0.0)));
        if popkern_core::trace::trace_enabled("setback") {
            eprintln!(
                "  [corner-select] k={k} coplanar={:?} coherent={:?} any={:?} chosen=({:.4},{:.4},{:.4}) gap={:.3e}",
                best_coplanar.map(|(d, _, _)| d),
                best_coplanar.map(|(_, _, c)| c),
                best_any.map(|(d, _)| d),
                chosen.1.x,
                chosen.1.y,
                chosen.1.z,
                chosen.0
            );
        }
        corners.push(chosen);
    }
    corners
}

/// Unit bisector of the two strips' tangent planes at a corner.
///
/// A sweep's parametric normal carries its own u×v orientation, not
/// the face's: adjacent strips can face opposite ways (hexagon6_blunt
/// G2: one strip's v runs the other way, `n_b = −n_a` at two corners,
/// and the raw sum was a roundoff vector that rotated the corner jet
/// by up to 68° — measured 4.7e-3 corner-block mismatch and 1.1e-1
/// rim tangency, PK-399). The bisector of two PLANES is sign-
/// invariant: `n_b` is aligned with `n_a` first.
pub(super) fn tangent_plane_bisector(n_a: Vec3, n_b: Vec3) -> Vec3 {
    let n_b = if n_a.dot(&n_b) < 0.0 { n_b * -1.0 } else { n_b };
    (n_a + n_b).normalized().unwrap_or(n_a)
}

/// Rotate a corner jet's first difference into the plane normal to
/// `bisector`, preserving its length (curvature untouched).
pub(super) fn jet_into_plane(jet: super::CornerJet, bisector: Vec3) -> super::CornerJet {
    let len = jet.d1.length();
    let flat = jet.d1 - bisector * jet.d1.dot(&bisector);
    let d1 = flat.normalized().map(|u| u * len).unwrap_or(jet.d1);
    super::CornerJet { d1, d2: jet.d2 }
}

/// The cap boundary's signed turn at a corner (radians): the angle
/// from this rim's arriving direction (`rim_out`, its outgoing tangent
/// at `c`) to the neighbour rim's leaving direction (the negation of
/// `corner_march`, the neighbour's leg away from `c`), measured in the
/// strips' common tangent plane (normal `bisector`), positive towards
/// the vertex — where the χ̃ `station_march` (oriented away from the
/// vertex) says the ring is NOT. `None` when either direction has no
/// in-plane component.
///
/// With cross-sectional rim ends (PK-401) the turn equals the shared
/// face's apex angle at every corner. Straight-uv rim ends on
/// mixed-radius setbacks put the trim-intersection corner on a flat
/// or REFLEX turn (hexagon6_variable, radii 0.02/0.04/0.06: −25.9° at
/// the corner between the two 0.02 strips), where a corner-compatible
/// cross field sweeps through the rim direction mid-piece — measured
/// assembly regularity 2.6e-6 and live cap↔strip g1 2.4e-1 (PK-399).
pub(super) fn cap_corner_turn(
    bisector: Vec3,
    rim_out: Vec3,
    station_march: Vec3,
    corner_march: Vec3,
) -> Option<f64> {
    let t_hat = (rim_out - bisector * rim_out.dot(&bisector)).normalized()?;
    let side = (station_march - t_hat * station_march.dot(&t_hat)).normalized()?;
    let leaving = corner_march * -1.0;
    let leaving = leaving - bisector * leaving.dot(&bisector);
    if leaving.length() < 1.0e-15 {
        return None;
    }
    Some((-leaving.dot(&side)).atan2(leaving.dot(&t_hat)))
}

/// Whether the cap boundary turns CONVEXLY at a corner, as a
/// tensor-product corner requires: `sin` of the turn
/// (`cap_corner_turn`) at least `sin_floor`. Degenerate directions
/// count as convex.
pub(super) fn cap_corner_is_convex(
    bisector: Vec3,
    rim_out: Vec3,
    station_march: Vec3,
    corner_march: Vec3,
    sin_floor: f64,
) -> bool {
    cap_corner_turn(bisector, rim_out, station_march, corner_march)
        .is_none_or(|turn| turn.sin() > sin_floor)
}

/// Per corner `c_s`, the eq (14) corner-block disagreement between
/// the two pieces meeting there, as the assembly's overlap average
/// will see it (`kp20_paper_border_corner_mismatch` with the paper's
/// default `d¹₃`). Traced as `[corner-block]` under `setback`.
pub(super) fn corner_block_mismatches(
    grouped: &[[popkern_core::geometry::gregory::Kp20PaperGuideBorderPiece; 2]],
    trace: bool,
) -> Vec<f64> {
    let n = grouped.len();
    let d1_3 = 2.0 * (core::f64::consts::TAU / n as f64).cos() / 3.0;
    grouped
        .iter()
        .enumerate()
        .map(|(s, sector)| {
            let value = popkern_core::geometry::gregory::kp20_paper_border_corner_mismatch(
                n,
                &sector[0].rows,
                &sector[1].rows,
                d1_3,
            );
            if trace {
                eprintln!("  [corner-block] c_{s}: mismatch {value:.3e}");
            }
            value
        })
        .collect()
}

/// `sin` of the angle between the two strips' surface normals at
/// `p` — the tangent-coherence measure for corner-pair selection.
/// `None` when either projection fails (the caller treats that as
/// incoherent).
fn strips_normal_disagreement(
    fa: &FilletSurfaceResult,
    fb: &FilletSurfaceResult,
    p: Point3,
) -> Option<f64> {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let mut normals = Vec::with_capacity(2);
    for fillet in [fa, fb] {
        let (u, v) = project_point_to_uv(&fillet.surface, &p)?;
        let on_surface = fillet.surface.point_at(u, v);
        if on_surface.distance_to(&p) > 1.0e-3 {
            return None;
        }
        normals.push(fillet.surface.normal_at(u, v).normalized()?);
    }
    Some(normals[0].cross(&normals[1]).length())
}

/// Both curves lie in one plane (the shared planar support face).
/// The plane is fit over BOTH curves' samples via principal spans —
/// fitting from one nearly-straight NURBS trim alone is
/// ill-conditioned (its residual curvature wobble dominates the
/// normal). Collinear sample sets are trivially coplanar. Curved
/// shared faces fail the residual test and fall back to
/// min-distance selection.
pub(super) fn pair_is_coplanar(ca: &NurbsCurve, cb: &NurbsCurve) -> bool {
    const SAMPLES: usize = 8;
    let mut points: Vec<Point3> = Vec::with_capacity(2 * (SAMPLES + 1));
    for c in [ca, cb] {
        let (t0, t1) = c.parameter_range();
        for i in 0..=SAMPLES {
            points.push(c.point_at(t0 + (t1 - t0) * i as f64 / SAMPLES as f64));
        }
    }
    let inv = 1.0 / points.len() as f64;
    let centroid = points.iter().fold(Point3::new(0.0, 0.0, 0.0), |acc, p| {
        Point3::new(acc.x + p.x * inv, acc.y + p.y * inv, acc.z + p.z * inv)
    });
    // First principal span: farthest sample from the centroid.
    let d1 = points
        .iter()
        .map(|p| *p - centroid)
        .max_by(|a, b| {
            a.length()
                .partial_cmp(&b.length())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
    let scale = d1.length();
    if scale < 1e-12 {
        return true;
    }
    // Second span: sample maximizing the cross with the first.
    let (normal, cross_len) = points
        .iter()
        .map(|p| d1.cross(&(*p - centroid)))
        .map(|c| (c, c.length()))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or((Vec3::new(0.0, 0.0, 0.0), 0.0));
    if cross_len < 1e-9 * scale * scale {
        // Collinear set — trivially coplanar.
        return true;
    }
    let Some(normal) = normal.normalized() else {
        return true;
    };
    let worst = points
        .iter()
        .map(|p| (*p - centroid).dot(&normal).abs())
        .fold(0.0_f64, f64::max);
    worst < coplanar_band_scale(scale)
}

fn coplanar_band_scale(scale: f64) -> f64 {
    (1e-6 * scale).max(1e-9)
}

/// Closest approach between two trim curves, always within the
/// curves' parameter domains and restricted to candidates passing
/// `accept` (the material-side gate — see `setback_corners`).
///
/// A coarse grid scan seeds a damped Gauss–Newton solve of the pair
/// parameters (residual `A(s) − B(t)`, Jacobian `[A'(s), −B'(t)]`),
/// which lands on an exact crossing in one step for straight
/// contact lines and on the true local pair minimum for skew
/// curves. The former alternating 1-D window scan was measured to
/// STALL on shallow-angle pairs (PK-399): the blunt pentagon's
/// shared-face contact lines meet at ~21°, where coordinate descent
/// contracts by only ~cos² per round while the windows shrank by
/// 0.35, leaving a 3.4e-3 "gap" on lines that intersect exactly and
/// a corner 1.7e-3 off both strips (the same class read 4.1e-3 on
/// the rolling-ball control and 1e-8 on the sharp pyramid, whose
/// lines meet at ~55°). The window scan remains only as the fallback
/// when the Newton polish cannot improve on the grid seed.
fn closest_between_refined(
    ca: &NurbsCurve,
    cb: &NurbsCurve,
    accept: &dyn Fn(Point3) -> bool,
) -> (f64, Point3) {
    const S: usize = 96;
    let (a0, a1) = ca.parameter_range();
    let (b0, b1) = cb.parameter_range();
    let mut best = (f64::INFINITY, Point3::new(0.0, 0.0, 0.0));
    let mut best_params = (a0, b0);
    for i in 0..=S {
        let ta = a0 + (a1 - a0) * i as f64 / S as f64;
        let pa = ca.point_at(ta);
        for j in 0..=S {
            let tb = b0 + (b1 - b0) * j as f64 / S as f64;
            let pb = cb.point_at(tb);
            let d = pa.distance_to(&pb);
            if d < best.0 && accept(midpoint(pa, pb)) {
                best = (d, midpoint(pa, pb));
                best_params = (ta, tb);
            }
        }
    }
    if let Some(polished) = newton_closest_pair(ca, cb, best_params, accept)
        && polished.0 < best.0
    {
        return polished;
    }
    let step_a = (a1 - a0) / S as f64;
    let step_b = (b1 - b0) / S as f64;
    // Fallback: alternate 1-D window scans around the grid minimum,
    // shrinking the windows each round — monotone descent to the
    // local pair minimum.
    {
        let mut wa = step_a.abs();
        let mut wb = step_b.abs();
        for _ in 0..12 {
            let (ta_lo, ta_hi) = (
                (best_params.0 - wa).max(a0.min(a1)),
                (best_params.0 + wa).min(a0.max(a1)),
            );
            for i in 0..=16 {
                let ta = ta_lo + (ta_hi - ta_lo) * i as f64 / 16.0;
                let d = ca.point_at(ta).distance_to(&cb.point_at(best_params.1));
                if d < best.0 && accept(midpoint(ca.point_at(ta), cb.point_at(best_params.1))) {
                    best = (d, midpoint(ca.point_at(ta), cb.point_at(best_params.1)));
                    best_params.0 = ta;
                }
            }
            let (tb_lo, tb_hi) = (
                (best_params.1 - wb).max(b0.min(b1)),
                (best_params.1 + wb).min(b0.max(b1)),
            );
            for j in 0..=16 {
                let tb = tb_lo + (tb_hi - tb_lo) * j as f64 / 16.0;
                let d = ca.point_at(best_params.0).distance_to(&cb.point_at(tb));
                if d < best.0 && accept(midpoint(ca.point_at(best_params.0), cb.point_at(tb))) {
                    best = (d, midpoint(ca.point_at(best_params.0), cb.point_at(tb)));
                    best_params.1 = tb;
                }
            }
            wa *= 0.35;
            wb *= 0.35;
        }
    }
    best
}

/// Damped Gauss–Newton polish of a closest-pair seed `(s, t)`:
/// minimises `|A(s) − B(t)|²` with Levenberg damping on the 2×2
/// normal equations, domain-clamped steps, and monotone acceptance
/// (a step must reduce the distance AND keep its midpoint on the
/// material side). Returns the polished `(distance, midpoint)`, or
/// `None` when no accepted step exists (a singular seed).
fn newton_closest_pair(
    ca: &NurbsCurve,
    cb: &NurbsCurve,
    seed: (f64, f64),
    accept: &dyn Fn(Point3) -> bool,
) -> Option<(f64, Point3)> {
    const MAX_ITERATIONS: usize = 40;
    const MAX_DAMPING_TRIES: usize = 12;
    let (a0, a1) = ca.parameter_range();
    let (b0, b1) = cb.parameter_range();
    let (a_lo, a_hi) = (a0.min(a1), a0.max(a1));
    let (b_lo, b_hi) = (b0.min(b1), b0.max(b1));
    let span = (a_hi - a_lo).max(b_hi - b_lo).max(f64::MIN_POSITIVE);
    let (mut s, mut t) = (seed.0.clamp(a_lo, a_hi), seed.1.clamp(b_lo, b_hi));
    let mut pa = ca.point_at(s);
    let mut pb = cb.point_at(t);
    let mut d = pa.distance_to(&pb);
    if !d.is_finite() {
        return None;
    }
    let mut lambda = 0.0_f64;
    let mut improved = false;
    for _ in 0..MAX_ITERATIONS {
        if d < 1.0e-15 {
            break;
        }
        let r = pa - pb;
        let da = ca.derivative_at(s, 1);
        let db = cb.derivative_at(t, 1);
        let (aa, ab, bb) = (da.dot(&da), -da.dot(&db), db.dot(&db));
        let (g0, g1) = (da.dot(&r), -db.dot(&r));
        let mut stepped = false;
        for _ in 0..MAX_DAMPING_TRIES {
            let m00 = aa * (1.0 + lambda);
            let m11 = bb * (1.0 + lambda);
            let det = m00 * m11 - ab * ab;
            if det.abs() > f64::MIN_POSITIVE {
                let ds = -(m11 * g0 - ab * g1) / det;
                let dt = -(m00 * g1 - ab * g0) / det;
                let s_new = (s + ds).clamp(a_lo, a_hi);
                let t_new = (t + dt).clamp(b_lo, b_hi);
                let pa_new = ca.point_at(s_new);
                let pb_new = cb.point_at(t_new);
                let d_new = pa_new.distance_to(&pb_new);
                if d_new < d && accept(midpoint(pa_new, pb_new)) {
                    let moved = (s_new - s).abs().max((t_new - t).abs());
                    s = s_new;
                    t = t_new;
                    pa = pa_new;
                    pb = pb_new;
                    d = d_new;
                    lambda = if lambda < 1.0e-12 { 0.0 } else { lambda * 0.1 };
                    stepped = true;
                    improved = true;
                    if moved < 1.0e-15 * span {
                        return Some((d, midpoint(pa, pb)));
                    }
                    break;
                }
            }
            lambda = if lambda == 0.0 { 1.0e-3 } else { lambda * 10.0 };
        }
        if !stepped {
            break;
        }
    }
    improved.then(|| (d, midpoint(pa, pb)))
}

fn midpoint(a: Point3, b: Point3) -> Point3 {
    Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
}

#[cfg(test)]
#[path = "kp20_paper_border_corners_tests.rs"]
mod tests;
