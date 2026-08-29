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
fn pair_is_coplanar(ca: &NurbsCurve, cb: &NurbsCurve) -> bool {
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

/// Closest approach between two trim curves, grid-scanned then
/// window-refined, always within the curves' parameter domains and
/// restricted to candidates passing `accept` (the material-side
/// gate — see `setback_corners`).
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
    let step_a = (a1 - a0) / S as f64;
    let step_b = (b1 - b0) / S as f64;
    // Local refinement: the coarse grid step can exceed the entire
    // corner region (full-length rails ~4 units vs a ~2e-2 corner on
    // hexagon6_variable), leaving the grid minimum an artifact.
    // Alternate 1-D window scans around the grid minimum, shrinking
    // the windows each round — monotone descent to the local pair
    // minimum.
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

fn midpoint(a: Point3, b: Point3) -> Point3 {
    Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
}
