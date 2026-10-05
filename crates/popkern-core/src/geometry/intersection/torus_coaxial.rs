//! Analytic torus intersection solver for coaxial pairs (PK-279).
//!
//! A torus and a second surface of revolution sharing the torus axis
//! — cylinder, cone, sphere (centre on the axis), or another torus —
//! intersect in exact CIRCLES centred on that axis. The classes solved
//! here are the everyday fillet configurations: a rim fillet against
//! the bore/boss wall it blends into, two stacked rim fillets, a rim
//! fillet against a countersink cone or a corner sphere. Before this
//! arm every such pair took NURBS promotion → capped AABB seeding →
//! marching.
//!
//! The mechanism is profile intersection in the (ρ, h) half-plane of
//! the shared axis (ρ = distance from axis, h = axial offset from the
//! torus centre): the torus profile is the circle centred (R, 0) with
//! radius r; the partner contributes a vertical line (cylinder), an
//! apex ray pair (cone), or a circle (sphere, torus). Each transversal
//! profile intersection point revolves into one exact 3D circle.
//!
//! Contracts and boundaries:
//! - **Tangency-silent** (PK-252): profile contacts within the linear
//!   tolerance band of tangency emit nothing, matching the
//!   plane×cylinder precedent and the marching path's behaviour on
//!   grazing fillet contact.
//! - **Coincident pairs decline**: an identical coaxial torus pair
//!   returns `None` so the coincident detector downstream owns the
//!   `CoincidentRegion` verdict.
//! - **Fallback boundary**: skew or parallel-offset axes return
//!   `None` and keep the marching path (a non-coaxial cylinder×torus
//!   intersection is a degree-8 curve with no exact form).
//! - Emitted circles honour the caller's tube/axial parameter windows
//!   (periodic-aware), matching the cylinder×cone precedent.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::intersection::circle_of_revolution::{
    circle_curve_on, station_in_window, torus_tube_station,
};
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Try to analytically intersect a torus with a coaxial surface of
/// revolution (cylinder, cone, axis-centred sphere, or torus).
///
/// Returns `Some(intersections)` when the pair is claimed (including
/// the empty and tangency-silent outcomes), or `None` when the pair is
/// not a coaxial torus combination — the caller falls back to the
/// general pipeline. Plane×torus pairs are owned by the dedicated
/// [`super::plane_torus`] solver, not this arm.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_torus_coaxial(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    match (surface1, surface2) {
        (SurfaceGeometry::Torus { .. }, _) if partner_kind(surface2).is_some() => intersect_impl(
            surface1,
            (u1_range, v1_range),
            surface2,
            (u2_range, v2_range),
            tolerance,
            false,
        ),
        (_, SurfaceGeometry::Torus { .. }) if partner_kind(surface1).is_some() => intersect_impl(
            surface2,
            (u2_range, v2_range),
            surface1,
            (u1_range, v1_range),
            tolerance,
            true,
        ),
        _ => None,
    }
}

/// Surfaces this arm can pair with a torus. Planes are excluded — the
/// plane×torus solver owns that pair in all orientations.
fn partner_kind(surface: &SurfaceGeometry) -> Option<()> {
    match surface {
        SurfaceGeometry::Cylinder { .. }
        | SurfaceGeometry::Cone { .. }
        | SurfaceGeometry::Sphere { .. }
        | SurfaceGeometry::Torus { .. } => Some(()),
        _ => None,
    }
}

/// A transversal profile intersection point in the (ρ, h) half-plane:
/// `rho` is SIGNED distance from the axis (negative values arise only
/// on spindle-torus inner lobes and revolve into the antipodal side),
/// `h` the axial offset from the torus centre, and `partner_station`
/// the partner surface's v-parameter at the circle when it has a
/// meaningful axial station to honour (cylinder, cone, torus).
struct ProfilePoint {
    rho: f64,
    h: f64,
    partner_station: Option<(f64, Option<f64>)>,
}

/// `torus` is the anchor operand; `swapped` = the torus was operand 2.
fn intersect_impl(
    torus: &SurfaceGeometry,
    torus_windows: ((f64, f64), (f64, f64)),
    partner: &SurfaceGeometry,
    partner_windows: ((f64, f64), (f64, f64)),
    tolerance: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let torus_v_range = torus_windows.1;
    let partner_v_range = partner_windows.1;
    let SurfaceGeometry::Torus {
        center,
        axis,
        major_radius,
        minor_radius,
    } = torus
    else {
        return None;
    };
    let axis_n = axis.normalized()?;
    let lin = tolerance.linear();

    let points = match partner {
        SurfaceGeometry::Cylinder {
            origin: po,
            axis: pa,
            radius: pr,
        } => {
            let pa_n = pa.normalized()?;
            let h0 = (*po - *center).dot(&axis_n);
            let extent = major_radius.abs() + minor_radius + pr.abs() + h0.abs();
            coaxial_gates(*center, axis_n, *po, Some(pa_n), extent, tolerance)?;
            let sense = pa_n.dot(&axis_n).signum();
            profile_line(*major_radius, *minor_radius, *pr, h0, sense, lin)
        }
        SurfaceGeometry::Cone {
            apex,
            axis: pa,
            half_angle,
        } => {
            let pa_n = pa.normalized()?;
            let h_apex_gate = (*apex - *center).dot(&axis_n);
            let extent = major_radius.abs() + minor_radius + h_apex_gate.abs();
            coaxial_gates(*center, axis_n, *apex, Some(pa_n), extent, tolerance)?;
            let tan_h = half_angle.ptan();
            if tan_h.abs() < tolerance.angular() {
                // Degenerate cone (≈ cylinder about its own axis with
                // zero radius growth): no transversal circle class to
                // claim here.
                return None;
            }
            let sense = pa_n.dot(&axis_n).signum();
            let h_apex = (*apex - *center).dot(&axis_n);
            profile_cone(*major_radius, *minor_radius, tan_h, h_apex, sense, lin)
        }
        SurfaceGeometry::Sphere {
            center: sc,
            radius: sr,
        } => {
            let extent = major_radius.abs() + minor_radius + sr.abs();
            coaxial_gates(*center, axis_n, *sc, None, extent, tolerance)?;
            let h_s = (*sc - *center).dot(&axis_n);
            profile_circle(*major_radius, *minor_radius, (0.0, h_s), *sr, lin)?
                .into_iter()
                .map(|(rho, h)| ProfilePoint {
                    rho,
                    h,
                    partner_station: None,
                })
                .collect()
        }
        SurfaceGeometry::Torus {
            center: c2,
            axis: a2,
            major_radius: big_r2,
            minor_radius: r2,
        } => {
            let a2_n = a2.normalized()?;
            let dz_gate = (*c2 - *center).dot(&axis_n);
            let extent = major_radius.abs() + minor_radius + big_r2.abs() + r2 + dz_gate.abs();
            coaxial_gates(*center, axis_n, *c2, Some(a2_n), extent, tolerance)?;
            let dz = dz_gate;
            // The partner's tube chart measures axial offset along ITS
            // OWN axis: antiparallel coaxial pairs (join strips author
            // material-side axes) flip the station sign. Without the
            // sense the window clip rejected true crossings — the
            // cylinder/cone branches carry it; this one forgot.
            let sense2 = a2_n.dot(&axis_n).signum();
            if dz.abs() <= lin
                && (big_r2 - major_radius).abs() <= lin
                && (r2 - minor_radius).abs() <= lin
            {
                // Coincident coaxial tori: the coincident detector
                // owns the CoincidentRegion verdict.
                return None;
            }
            profile_circle(*major_radius, *minor_radius, (*big_r2, dz), *r2, lin)?
                .into_iter()
                .map(|(rho, h)| ProfilePoint {
                    rho,
                    h,
                    partner_station: Some((
                        torus_tube_station((rho).abs(), sense2 * (h - dz), *big_r2),
                        Some(std::f64::consts::TAU),
                    )),
                })
                .collect()
        }
        _ => return None,
    };

    let (a1, a2v) = orthonormal_basis(axis_n);
    let mut curves = Vec::new();
    for point in points {
        if point.rho.abs() <= lin {
            continue; // zero-radius circle (axis contact)
        }
        // Honour the caller's tube window on the torus side.
        let station = torus_tube_station(point.rho.abs(), point.h, *major_radius);
        if !station_in_window(station, torus_v_range, Some(std::f64::consts::TAU)) {
            continue;
        }
        if let Some((partner_station, period)) = point.partner_station
            && !station_in_window(partner_station, partner_v_range, period)
        {
            continue;
        }
        let circle_center = *center + axis_n * point.h;
        // A negative signed ρ (spindle inner lobe) revolves into the
        // antipodal phase: negate the whole frame (phase shift of π,
        // same rotation sense) so the circle radius stays positive.
        let (f1, f2) = if point.rho >= 0.0 {
            (a1, a2v)
        } else {
            (-a1, -a2v)
        };
        let (s1, w1, s2, w2) = if swapped {
            (partner, partner_windows, torus, torus_windows)
        } else {
            (torus, torus_windows, partner, partner_windows)
        };
        let Some(curve) = circle_curve_on(
            s1,
            w1,
            s2,
            w2,
            circle_center,
            f1,
            f2,
            point.rho.abs(),
            tolerance,
        ) else {
            return None; // chart failure → decline the arm entirely
        };
        curves.push(SurfaceSurfaceIntersection::Curve(curve));
    }
    Some(curves)
}

/// Coaxiality gates: partner axis parallel to the torus axis (when it
/// has one) and the partner's anchor point on the torus axis line.
/// Returns `None` (→ marching fallback) when either gate fails.
///
/// The parallelism gate is GEOMETRIC: `extent` is the configuration's
/// characteristic reach, so a pair snaps to the coaxial solution only
/// when the snap's worst-case error `sin(axis angle)·extent` stays
/// within the linear tolerance. Authored coaxial features agree to
/// floating-point precision; near-miss pairs keep the (slower but
/// honest) marching path.
fn coaxial_gates(
    center: Point3,
    axis_n: Vec3,
    partner_anchor: Point3,
    partner_axis: Option<Vec3>,
    extent: f64,
    tolerance: &Tolerance,
) -> Option<()> {
    if let Some(pa) = partner_axis
        && pa.cross(&axis_n).length() * extent.max(1.0e-9) > tolerance.linear()
    {
        return None;
    }
    let d = partner_anchor - center;
    let off_axis = d - axis_n * d.dot(&axis_n);
    if off_axis.length() > tolerance.linear() {
        return None;
    }
    Some(())
}

/// Torus profile circle (centre (R, 0), radius r) against the vertical
/// line ρ = rc: the coaxial-cylinder profile. `h0`/`sense` translate
/// axial offsets into the cylinder's own v parameter.
fn profile_line(big_r: f64, r: f64, rc: f64, h0: f64, sense: f64, lin: f64) -> Vec<ProfilePoint> {
    let d_rho = rc - big_r;
    if (d_rho.abs() - r).abs() <= lin {
        return Vec::new(); // tangency-silent band
    }
    if d_rho.abs() > r {
        return Vec::new(); // no contact
    }
    let h = (r * r - d_rho * d_rho).sqrt();
    [h, -h]
        .into_iter()
        .map(|h| ProfilePoint {
            rho: rc,
            h,
            partner_station: Some((sense * (h - h0), None)),
        })
        .collect()
}

/// Torus profile circle against the two nappe rays of a coaxial cone
/// (apex on the axis at height `h_apex`, radius growth `tan_h` per
/// unit of |cone v|, cone axis `sense`-aligned with the torus axis).
fn profile_cone(
    big_r: f64,
    r: f64,
    tan_h: f64,
    h_apex: f64,
    sense: f64,
    lin: f64,
) -> Vec<ProfilePoint> {
    let tan_a = tan_h.abs();
    let mut out = Vec::new();
    for nappe in [1.0f64, -1.0] {
        // Ray from (0, h_apex) with direction (tan_a, sense·nappe) in
        // (ρ, h), parameterized by τ = |cone v| ≥ 0.
        let dir_h = sense * nappe;
        let len = (tan_a * tan_a + 1.0).sqrt();
        // Distance from the profile-circle centre (R, 0) to the ray's
        // carrier line, for the tangency-silent band.
        let d_line = (tan_a * (0.0 - h_apex) - dir_h * big_r).abs() / len;
        if (d_line - r).abs() <= lin {
            continue; // tangency-silent band
        }
        if d_line > r {
            continue;
        }
        // (tan_a·τ − R)² + (h_apex + dir_h·τ)² = r², quadratic in τ.
        let a = tan_a * tan_a + 1.0;
        let b = 2.0 * (-big_r * tan_a + h_apex * dir_h);
        let c0 = big_r * big_r + h_apex * h_apex - r * r;
        let disc = b * b - 4.0 * a * c0;
        if disc <= 0.0 {
            continue;
        }
        let sq = disc.sqrt();
        for tau in [(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
            if tau <= lin {
                continue; // wrong nappe side or apex graze
            }
            out.push(ProfilePoint {
                rho: tan_a * tau,
                h: h_apex + dir_h * tau,
                partner_station: Some((nappe * tau, None)),
            });
        }
    }
    // The two carrier lines meet at the apex; a root pair straddling
    // the apex can duplicate a station. Dedupe by profile proximity.
    out.dedup_by(|a, b| (a.rho - b.rho).abs() <= lin && (a.h - b.h).abs() <= lin);
    out
}

/// Torus profile circle (centre (R, 0), radius r) against a second
/// profile circle at `c2 = (ρ, h)` with radius `r2` — the coaxial
/// sphere and torus×torus profiles. Returns the transversal
/// intersection points, an empty list for separated or
/// tangency-silent contact, and `None` for numerically degenerate
/// concentric same-radius profiles (caller declines the arm).
fn profile_circle(
    big_r: f64,
    r: f64,
    c2: (f64, f64),
    r2: f64,
    lin: f64,
) -> Option<Vec<(f64, f64)>> {
    let dx = c2.0 - big_r;
    let dy = c2.1;
    let dist = (dx * dx + dy * dy).sqrt();
    if dist <= lin {
        // Concentric profiles: either coincident (handled upstream) or
        // disjoint nested circles — never transversal.
        return Some(Vec::new());
    }
    if (dist - (r + r2)).abs() <= lin || (dist - (r - r2).abs()).abs() <= lin {
        return Some(Vec::new()); // tangency-silent band
    }
    if dist > r + r2 || dist < (r - r2).abs() {
        return Some(Vec::new());
    }
    let along = (dist * dist + r * r - r2 * r2) / (2.0 * dist);
    let half_chord_sq = r * r - along * along;
    if half_chord_sq <= 0.0 {
        return Some(Vec::new());
    }
    let half_chord = half_chord_sq.sqrt();
    let (ux, uy) = (dx / dist, dy / dist);
    let base = (big_r + ux * along, uy * along);
    Some(vec![
        (base.0 - uy * half_chord, base.1 + ux * half_chord),
        (base.0 + uy * half_chord, base.1 - ux * half_chord),
    ])
}

#[cfg(test)]
#[path = "torus_coaxial_tests.rs"]
mod tests;
