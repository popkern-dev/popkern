//! Variable-radius NURBS sweep builder for [`crate::modify_inplace`].
//! Split out so `modify_inplace_geometry.rs` stays under the
//! 500-line structural cap.

use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

use crate::fillet_surface::{build_clamped_uniform_knots, rational_arc_row};

/// Inputs for [`build_variable_plane_plane_sweep`].
pub(crate) struct VariableSweepInputs {
    pub p_start: Point3,
    pub p_end: Point3,
    pub plane_a_normal: Vec3,
    pub plane_b_normal: Vec3,
    pub signed_perp_a: Vec3,
    pub signed_perp_b: Vec3,
    pub r_start: f64,
    pub r_end: f64,
    pub n_sections: usize,
    pub n_arc_samples: usize,
}

/// Build a true variable-radius NURBS sweep surface for a
/// plane-plane fillet. Samples `n_sections` cross-section arcs
/// along the sharp edge from `(p_start, r_start)` to
/// `(p_end, r_end)`, each at the locally interpolated radius
/// (linear interpolation), and skins them into the same
/// section-exact rational loft class as the fresh-blend pipeline
/// (degree 1 × 2 rational rows — each cross-section an EXACT
/// circular arc; PK-179: the earlier degree-1×1 sampled grid's
/// faceted sections put a 2-6e-4 floor under terminal settlement
/// and the modify/recreate equivalence).
///
/// Used by to produce a geometrically exact variable-
/// radius fillet surface in place. For arbitrary radius profiles
/// (e.g. `FilletRadius::ByPoints`) use
/// [`build_parametric_plane_plane_sweep`] directly with a
/// piecewise-linear closure.
pub(crate) fn build_variable_plane_plane_sweep(inputs: &VariableSweepInputs) -> SurfaceGeometry {
    let r_start = inputs.r_start;
    let r_end = inputs.r_end;
    build_parametric_plane_plane_sweep(inputs, |t| r_start + t * (r_end - r_start))
}

/// Generalised variant of [`build_variable_plane_plane_sweep`] that
/// takes a `radius_at: Fn(f64) -> f64` closure for the local-radius
/// profile. The linear-radius variant calls this with a linear
/// closure; the `ByPoints` in-place path passes a
/// piecewise-linear closure over the full points list, so
/// intermediate radius points shape the swept surface.
pub(crate) fn build_parametric_plane_plane_sweep<F: Fn(f64) -> f64>(
    inputs: &VariableSweepInputs,
    radius_at: F,
) -> SurfaceGeometry {
    let VariableSweepInputs {
        p_start,
        p_end,
        plane_a_normal,
        plane_b_normal,
        signed_perp_a,
        signed_perp_b,
        r_start: _,
        r_end: _,
        n_sections,
        n_arc_samples,
    } = *inputs;
    let n_sec = n_sections.max(2);
    let n_arc = n_arc_samples.max(3);

    // Half-angle of the dihedral and bisector direction. The
    // rolling ball at radius R sits on the bisector at distance
    // R / sin(half_angle) from the sharp edge.
    let cos_angle = plane_a_normal.dot(&plane_b_normal).clamp(-1.0, 1.0);
    let half_angle = cos_angle.pacos() * 0.5;
    let sin_half = half_angle.psin().max(1e-12);

    let bisector_raw = Vec3::new(
        signed_perp_a.x + signed_perp_b.x,
        signed_perp_a.y + signed_perp_b.y,
        signed_perp_a.z + signed_perp_b.z,
    );
    let bisector_len = bisector_raw.length().max(1e-12);
    let bisector = Vec3::new(
        bisector_raw.x / bisector_len,
        bisector_raw.y / bisector_len,
        bisector_raw.z / bisector_len,
    );

    let mut control_points: Vec<Point3> = Vec::with_capacity(n_sec * n_arc);
    let mut rational_rows: Option<Vec<([Point3; 3], f64)>> = Some(Vec::with_capacity(n_sec));

    for i in 0..n_sec {
        let t = i as f64 / (n_sec - 1) as f64;
        let p_i = Point3::new(
            p_start.x + t * (p_end.x - p_start.x),
            p_start.y + t * (p_end.y - p_start.y),
            p_start.z + t * (p_end.z - p_start.z),
        );
        let r_i = radius_at(t);

        let footprint_a = Point3::new(
            p_i.x + signed_perp_a.x * r_i,
            p_i.y + signed_perp_a.y * r_i,
            p_i.z + signed_perp_a.z * r_i,
        );
        let footprint_b = Point3::new(
            p_i.x + signed_perp_b.x * r_i,
            p_i.y + signed_perp_b.y * r_i,
            p_i.z + signed_perp_b.z * r_i,
        );
        let center_i = Point3::new(
            p_i.x + bisector.x * (r_i / sin_half),
            p_i.y + bisector.y * (r_i / sin_half),
            p_i.z + bisector.z * (r_i / sin_half),
        );

        // Sample the cross-section arc from footprint_a through to
        // footprint_b, centred at center_i. Build an orthonormal
        // (u, w) basis in the arc plane, sweep angle from 0 to the
        // arc's full sweep, emit `n_arc` points.
        let v1 = Vec3::new(
            footprint_a.x - center_i.x,
            footprint_a.y - center_i.y,
            footprint_a.z - center_i.z,
        );
        let v2 = Vec3::new(
            footprint_b.x - center_i.x,
            footprint_b.y - center_i.y,
            footprint_b.z - center_i.z,
        );
        let r1 = v1.length().max(1e-15);
        let u_basis = Vec3::new(v1.x / r1, v1.y / r1, v1.z / r1);
        let dot = v2.dot(&u_basis);
        let w_raw = Vec3::new(
            v2.x - dot * u_basis.x,
            v2.y - dot * u_basis.y,
            v2.z - dot * u_basis.z,
        );
        let w_len = w_raw.length();
        if w_len < 1e-15 {
            rational_rows = None;
            // Degenerate cross-section — emit a linear ramp.
            for j in 0..n_arc {
                let s = j as f64 / (n_arc - 1) as f64;
                control_points.push(Point3::new(
                    footprint_a.x + s * (footprint_b.x - footprint_a.x),
                    footprint_a.y + s * (footprint_b.y - footprint_a.y),
                    footprint_a.z + s * (footprint_b.z - footprint_a.z),
                ));
            }
            continue;
        }
        if let Some(rows) = rational_rows.as_mut() {
            match rational_arc_row(&center_i, &footprint_a, &footprint_b) {
                Some(row) => rows.push(row),
                None => rational_rows = None,
            }
        }
        let w_basis = Vec3::new(w_raw.x / w_len, w_raw.y / w_len, w_raw.z / w_len);
        let cos_alpha = dot / r1;
        let sin_alpha = w_len / r1;
        let alpha = sin_alpha.patan2(cos_alpha);
        for j in 0..n_arc {
            let s = j as f64 / (n_arc - 1) as f64;
            let theta = s * alpha;
            let cos_t = theta.pcos();
            let sin_t = theta.psin();
            control_points.push(Point3::new(
                center_i.x + r_i * (cos_t * u_basis.x + sin_t * w_basis.x),
                center_i.y + r_i * (cos_t * u_basis.y + sin_t * w_basis.y),
                center_i.z + r_i * (cos_t * u_basis.z + sin_t * w_basis.z),
            ));
        }
    }

    // Section-exact rational loft (same class as the fresh-blend
    // pipeline's sweep): u = along edge (n_sec stations), v = one
    // rational-quadratic row per station — every cross-section is
    // an EXACT circular arc. Falls back to the sampled degree-1
    // grid only when a station cannot form a single-segment
    // rational row (degenerate or >~115 degree section).
    let degree_u = 1usize.min(n_sec - 1);
    if let Some(rows) = rational_rows
        && rows.len() == n_sec
    {
        let mut cps = Vec::with_capacity(n_sec * 3);
        let mut weights = Vec::with_capacity(n_sec * 3);
        for (row, w_mid) in &rows {
            cps.extend_from_slice(row);
            weights.extend_from_slice(&[1.0, *w_mid, 1.0]);
        }
        return SurfaceGeometry::Nurbs(NurbsSurface {
            degree_u,
            degree_v: 2,
            knots_u: build_clamped_uniform_knots(n_sec, degree_u),
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: cps,
            num_u: n_sec,
            num_v: 3,
            weights,
        });
    }
    let degree_v = 1usize.min(n_arc - 1);
    let knots_u = clamped_uniform_knots(n_sec, degree_u);
    let knots_v = clamped_uniform_knots(n_arc, degree_v);
    let weights = vec![1.0; n_sec * n_arc];

    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u: n_sec,
        num_v: n_arc,
        weights,
    })
}

fn clamped_uniform_knots(n_cps: usize, degree: usize) -> Vec<f64> {
    // n_cps + degree + 1 knots, clamped at both ends.
    let n_knots = n_cps + degree + 1;
    let interior = n_knots.saturating_sub(2 * (degree + 1));
    let mut knots: Vec<f64> = vec![0.0; degree + 1];
    knots.extend((1..=interior).map(|i| i as f64 / (interior + 1) as f64));
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    knots
}
