//! The display mesh's row shading contract (PK-914): how far a needle
//! row's interpolated normal field departs from the rate at which the
//! surface turns at the row's base.
//!
//! A renderer interpolates vertex normals linearly across each facet,
//! so across a needle (a facet much longer than it is high) its
//! shading normal turns at the secant rate from its base edge to its
//! apex. The surface turns at its own rate at the base. Reflection
//! lines (zebra stripes) bend with the shading rate, so where the two
//! differ, the stripes kink although the surface does not. The chord
//! height bounds where a facet sits and the angular facet bound bounds
//! how far its normal turns; neither bounds this.
//!
//! A G2 blend strip beside its support shows the gap. Its cross
//! curvature rises from zero at a rail, and its rail needles are as
//! long as the rail's samples are apart (0.306 on pentagon5_variable)
//! but only as high as the end arc's first sample (0.0064). That first
//! row turns at 3.5 rad/mm while the surface turns at zero at the
//! rail, and the zebra reads a crease of 1.06e-2 per pixel at a seam
//! that is G1 to 0 rad.
//!
//! [`row_shading_error`] measures the gap as an angle: the apex
//! normal's departure from its first-order prediction at the base (the
//! foot normal turned by the surface's shape operator over the apex
//! offset, in arc length). For curvature rising as `k' d` from the
//! base, a row of height `h` departs by `k' h² / 2`, and the crease
//! it shows is that angle over `h`. Holding the angle to a bound
//! narrows rows where the curvature varies fast and leaves them alone
//! where it does not: a constant-curvature round departs by nothing at
//! any row height (to third order), and so does a plane. The bound is
//! a fixed fraction of the angular facet bound
//! ([`ROW_SHADING_DIVISOR`]), so the contract is scale-free like that
//! bound. The display pass
//! (`angular_refine::refine_face_shading`) holds a face's seam band to
//! it: the zebra reads a seam over a few pixels on either side, so the
//! rows behind the first count as much as the first.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{PortableMath, Vec3};

use crate::mesh::MeshVertex;
use crate::vertex_normal::{chart_singular, singular_without_unique_normal};

/// The row shading bound is the angular facet bound over this. The
/// display bound of 12° gives 0.375° (6.5e-3 rad). The
/// pentagon5_variable strips' first rows departed by 2.2e-2 rad (k'
/// about 1070 /mm² over a 0.0064 row). Held to the bound, the worst
/// strip-support seam reads 3.5e-3 per pixel, the surface view's own
/// 4.6e-3 floor or below. A bound of a 64th read the same at 1.4x the
/// added triangles, and an 8th or a 16th left the worst seam over
/// 4.7e-3 (PK-914).
pub const ROW_SHADING_DIVISOR: f64 = 32.0;

/// A facet is a needle, held to the contract, when its base edge is
/// at least this many times its height over it. A compact facet's row
/// is as long as it is high, and the boundary sampling and the angular
/// bound already size both. A needle's row height is sized by nothing
/// along its base, so only splits across the row narrow it.
pub const NEEDLE_ASPECT: f64 = 4.0;

/// A needle under this fraction of the chord height high is below the
/// display's resolution, and the contract does not hold it. Near a
/// vanishing strip end the cross-section shrinks to nothing (a
/// pentagon corner's strip, 1e-4 across), so the normal turns over
/// no distance. No row height there keeps the bound, and splitting
/// such a row spent the whole growth budget (PK-914).
pub const MIN_ROW_FRACTION: f64 = 1.0 / 16.0;

/// The row shading error (radians) of a needle on the base edge
/// `a`–`b` with apex `c`: the apex normal's departure from the base's
/// first-order prediction. Returns `None` for a facet the contract
/// does not hold: one whose base is under [`NEEDLE_ASPECT`] times its
/// height, whose apex does not stand over the base, or whose height
/// is under [`MIN_ROW_FRACTION`] of `chord_height`; or one it cannot
/// measure (a chart-singular vertex or foot, where the normal is
/// undefined, or a degenerate facet).
pub fn row_shading_error(
    surface: &SurfaceGeometry,
    a: &MeshVertex,
    b: &MeshVertex,
    c: &MeshVertex,
    chord_height: f64,
) -> Option<f64> {
    if [a, b, c]
        .iter()
        .any(|x| singular_without_unique_normal(surface, x.uv.0, x.uv.1))
    {
        return None;
    }
    let ab = b.position - a.position;
    let l2 = ab.dot(&ab);
    if l2 <= 1e-300 {
        return None;
    }
    let ac = c.position - a.position;
    let height = ab.cross(&ac).length() / l2.sqrt();
    if l2.sqrt() < NEEDLE_ASPECT * height || height < MIN_ROW_FRACTION * chord_height {
        return None;
    }
    // The apex's foot on the base, in 3-D and in the chart. An apex
    // beyond the base's ends does not stand over it: it lies along
    // it, and narrowing its row narrows nothing.
    let t = ac.dot(&ab) / l2;
    if !(0.0..=1.0).contains(&t) {
        return None;
    }
    let lerp = |p: f64, q: f64| p + (q - p) * t;
    let foot_uv = (lerp(a.uv.0, b.uv.0), lerp(a.uv.1, b.uv.1));
    if chart_singular(surface, foot_uv.0, foot_uv.1) {
        return None;
    }
    let (na, nb, nc) = (unit(a.normal)?, unit(b.normal)?, unit(c.normal)?);
    let n_foot = unit(na * (1.0 - t) + nb * t)?;
    let (u0, v0) = foot_uv;
    // The apex offset in the tangent plane at the foot, as a chart
    // step, scaled to the chord's length: the prediction is in arc
    // length, blind to how the chart parametrizes the surface.
    let w = c.position - surface.point_at(u0, v0);
    let w_t = w - n_foot * w.dot(&n_foot);
    let (su, sv) = (surface.partial_u(u0, v0), surface.partial_v(u0, v0));
    let (g11, g12, g22) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let det = g11 * g22 - g12 * g12;
    if det <= 1e-24 * (g11 * g22).max(1e-300) || w_t.length() <= 1e-300 {
        return None;
    }
    let (r1, r2) = (su.dot(&w_t), sv.dot(&w_t));
    let scale = w.length() / w_t.length();
    let du = scale * (g22 * r1 - g12 * r2) / det;
    let dv = scale * (g11 * r2 - g12 * r1) / det;
    // The surface's unit-normal Jacobian at the foot, by central
    // differences at a millionth of the step.
    let hu = (1e-6 * du.abs()).max(1e-12);
    let hv = (1e-6 * dv.abs()).max(1e-12);
    let at = |u: f64, v: f64| unit(surface.normal_at(u, v));
    let nu = (at(u0 + hu, v0)? - at(u0 - hu, v0)?) * (0.5 / hu);
    let nv = (at(u0, v0 + hv)? - at(u0, v0 - hv)?) * (0.5 / hv);
    // The mesh normals may face the other way (a reversed face).
    let sign = if at(u0, v0)?.dot(&n_foot) < 0.0 {
        -1.0
    } else {
        1.0
    };
    // The foot normal turned by the predicted turn (an angle, along
    // the tangent the turn points in).
    let turn = (nu * du + nv * dv) * sign;
    let angle = turn.length();
    let predicted = match unit(turn) {
        Some(dir) => n_foot * angle.pcos() + dir * angle.psin(),
        None => n_foot,
    };
    Some(nc.cross(&predicted).length().patan2(nc.dot(&predicted)))
}

/// The row shading bound under the angular facet bound `max_angle`.
pub fn row_shading_bound(max_angle: f64) -> f64 {
    max_angle / ROW_SHADING_DIVISOR
}

fn unit(v: Vec3) -> Option<Vec3> {
    let l = v.length();
    (l > 1e-300).then(|| Vec3::new(v.x / l, v.y / l, v.z / l))
}

#[cfg(test)]
#[path = "shading_contract_tests.rs"]
mod tests;
