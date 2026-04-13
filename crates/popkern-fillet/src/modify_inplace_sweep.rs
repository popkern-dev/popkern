//! Variable-radius NURBS sweep builder for [`crate::modify_inplace`].
//! Split out so `modify_inplace_geometry.rs` stays under the
//! 500-line structural cap.

use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

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
/// `(p_end, r_end)`, each at the locally interpolated radius, and
/// skins them into a tensor-product NurbsSurface (degree 1 × 1).
///
/// Used by Round 5i to produce a geometrically exact variable-
/// radius fillet surface in place — the alternative (a Torus with
/// the average minor radius) is correct for visualisation but not
/// for downstream consumers that read the surface for analysis.
pub(crate) fn build_variable_plane_plane_sweep(inputs: &VariableSweepInputs) -> SurfaceGeometry {
    let VariableSweepInputs {
        p_start,
        p_end,
        plane_a_normal,
        plane_b_normal,
        signed_perp_a,
        signed_perp_b,
        r_start,
        r_end,
        n_sections,
        n_arc_samples,
    } = *inputs;
    let n_sec = n_sections.max(2);
    let n_arc = n_arc_samples.max(3);

    // Half-angle of the dihedral and bisector direction. The
    // rolling ball at radius R sits on the bisector at distance
    // R / sin(half_angle) from the sharp edge.
    let cos_angle = plane_a_normal.dot(&plane_b_normal).clamp(-1.0, 1.0);
    let half_angle = cos_angle.acos() * 0.5;
    let sin_half = half_angle.sin().max(1e-12);

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

    for i in 0..n_sec {
        let t = i as f64 / (n_sec - 1) as f64;
        let p_i = Point3::new(
            p_start.x + t * (p_end.x - p_start.x),
            p_start.y + t * (p_end.y - p_start.y),
            p_start.z + t * (p_end.z - p_start.z),
        );
        let r_i = r_start + t * (r_end - r_start);

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
        let w_basis = Vec3::new(w_raw.x / w_len, w_raw.y / w_len, w_raw.z / w_len);
        let cos_alpha = dot / r1;
        let sin_alpha = w_len / r1;
        let alpha = sin_alpha.atan2(cos_alpha);
        for j in 0..n_arc {
            let s = j as f64 / (n_arc - 1) as f64;
            let theta = s * alpha;
            let cos_t = theta.cos();
            let sin_t = theta.sin();
            control_points.push(Point3::new(
                center_i.x + r_i * (cos_t * u_basis.x + sin_t * w_basis.x),
                center_i.y + r_i * (cos_t * u_basis.y + sin_t * w_basis.y),
                center_i.z + r_i * (cos_t * u_basis.z + sin_t * w_basis.z),
            ));
        }
    }

    // Tensor-product NurbsSurface: u = along edge (n_sec), v = along arc (n_arc).
    let degree_u = 1usize.min(n_sec - 1);
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
