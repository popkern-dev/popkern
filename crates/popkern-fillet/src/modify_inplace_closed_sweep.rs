//! Closed-rim plane-cylinder varying-radius NURBS sweep builder.
//!
//! Round 5s: true around-the-rim in-place path for
//! `modify_closed_in_place_bypoints`. Given the support plane,
//! cylinder, and a `radius_at: Fn(f64) -> f64` profile (param 0..1
//! = normalized angle around the cylinder axis), emits:
//!
//! - a tensor-product `NurbsSurface` for the fillet face (u = angle
//!   around rim, v = cross-section arc)
//! - a closed-rim `NurbsCurve` for the plane-side contact loop
//!   (each control point at distance `cyl_radius ± r(θ)` from the
//!   cylinder axis in the plane)
//! - a closed-rim `NurbsCurve` for the cylinder-side contact loop
//!   (each control point on the cylinder at height `s_plane · r(θ)`
//!   above the plane)
//! - the seam-edge endpoint positions at `θ = 0` on each side
//!
//! The existing [`crate::modify_inplace::modify_closed_in_place`]
//! constant-radius path produces an analytic Torus surface; this
//! module is only used when the caller passes a `ByPoints` profile
//! with non-uniform minor radii.
//!
//! # Assumptions
//!
//! The plane support must be perpendicular to the cylinder axis
//! (the typical "cube face with cylindrical hole" case). An oblique
//! plane would require a more general rolling-ball centre curve;
//! the current callers always satisfy this.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{CurveOnEdge, NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::BlendError;
use crate::modify_inplace::{ClosedLayout, collect_closed_layout};
use crate::modify_inplace_geometry::{determine_fillet_signs, project_to_plane_along_axis};
use crate::unfillet_closed_helpers::{
    CylinderGeometry, PlaneGeometry, extract_cylinder_geometry, extract_plane_geometry,
    find_seam_edge_at_vertex_excluding, rebuild_line_curve,
};

/// Round 5s: around-the-rim varying-radius closed-fillet in-place
/// update. Unlike `modify_closed_in_place`, which assigns an
/// analytic Torus, this path assembles a true NURBS sweep whose
/// minor radius follows the `points` piecewise-linear profile
/// (first component = normalized angle around the cylinder axis,
/// second = local minor radius).
///
/// The two contact loops become closed NURBS curves threaded
/// through per-section points, so the plane-side curve's radial
/// distance from the cylinder axis and the cylinder-side curve's
/// height above the plane both vary around the rim.
///
/// The seam vertex is re-anchored at `θ = 0`, so `pts[0].1` is
/// the radius at the seam.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] for the same reasons as
///   `modify_closed_in_place`, plus when `points` is empty.
pub fn modify_closed_in_place_bypoints(
    brep: &BRep,
    fillet_face: FaceHandle,
    points: &[(f64, f64)],
) -> Result<BRep, BlendError> {
    if points.is_empty() {
        return Err(BlendError::InvalidParameter {
            parameter: "points".to_string(),
            value: "empty".to_string(),
            reason: "ByPoints radius profile must have at least one point".to_string(),
        });
    }

    let layout: ClosedLayout = collect_closed_layout(brep, fillet_face)?;
    let plane_geom = extract_plane_geometry(brep, layout.plane_support)?;
    let cyl_geom = extract_cylinder_geometry(brep, layout.cyl_support)?;

    let v_plane_pos_old = brep.vertices[layout.v_plane_side].position;
    let v_cyl_pos_old = brep.vertices[layout.v_cyl_side].position;
    let signs = determine_fillet_signs(&plane_geom, &cyl_geom, &v_plane_pos_old, &v_cyl_pos_old)?;
    let plane_circle_centre = project_to_plane_along_axis(&plane_geom, &cyl_geom);

    let pts: Vec<(f64, f64)> = points.to_vec();
    let radius_at = |t: f64| crate::center_curve::bypoints_local_radius(&pts, t);
    let sweep = build_closed_plane_cylinder_sweep(&ClosedSweepInputs {
        plane: plane_geom,
        cyl: cyl_geom,
        plane_circle_centre,
        s_plane: signs.s_plane,
        s_radial: signs.s_radial,
        n_sections: 36,
        n_arc_samples: 9,
        radius_at: &radius_at,
    });

    let mut out = brep.clone();
    out.vertices[layout.v_plane_side].position = sweep.seam_plane_pos;
    out.vertices[layout.v_cyl_side].position = sweep.seam_cyl_pos;

    out.edges[layout.plane_edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(sweep.plane_contact_curve),
        t_start: 0.0,
        t_end: 1.0,
    });
    out.edges[layout.cyl_edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(sweep.cyl_contact_curve),
        t_start: 0.0,
        t_end: 1.0,
    });

    let cyl_loop_start = out.loops[brep.half_edges[layout.cyl_twin_h].loop_].half_edge;
    if let Some(seam_edge) = find_seam_edge_at_vertex_excluding(
        &out,
        cyl_loop_start,
        layout.cyl_twin_h,
        layout.v_cyl_side,
    ) {
        rebuild_line_curve(&mut out, seam_edge);
    }

    out.faces[fillet_face].surface = Some(SurfaceOnFace {
        geometry: sweep.surface,
        orientation: true,
    });
    crate::construct::strip_orientation::reconcile_strip_orientation_boundary(
        &mut out,
        fillet_face,
    );

    Ok(out)
}

/// Inputs for [`build_closed_plane_cylinder_sweep`].
pub(crate) struct ClosedSweepInputs<'a> {
    pub plane: PlaneGeometry,
    pub cyl: CylinderGeometry,
    /// Cylinder axis projected onto the plane. This is the common
    /// centre of both contact loops.
    pub plane_circle_centre: Point3,
    /// +1 if the cylinder-side contact lies on the positive
    /// cylinder-axis side of the plane; -1 otherwise.
    pub s_plane: f64,
    /// +1 if the plane-side contact lies outside the cylinder
    /// (radial distance > cylinder radius); -1 if inside.
    pub s_radial: f64,
    /// Angular sections around the rim (u direction of the sweep).
    /// Typically 36 — finer than necessary but cheap at degree 1.
    pub n_sections: usize,
    /// Cross-section arc samples (v direction of the sweep).
    pub n_arc_samples: usize,
    /// Local minor radius at normalized angular parameter
    /// `t ∈ [0, 1]`. `t = 0` is the seam (θ = 0); `t = 1` wraps
    /// back to the same angular position.
    pub radius_at: &'a dyn Fn(f64) -> f64,
}

/// Output of [`build_closed_plane_cylinder_sweep`].
pub(crate) struct ClosedSweepResult {
    pub surface: SurfaceGeometry,
    pub plane_contact_curve: NurbsCurve,
    pub cyl_contact_curve: NurbsCurve,
    pub seam_plane_pos: Point3,
    pub seam_cyl_pos: Point3,
}

/// Build an around-the-rim NURBS sweep + matching contact curves
/// for a plane-cylinder varying-radius fillet.
pub(crate) fn build_closed_plane_cylinder_sweep(inputs: &ClosedSweepInputs) -> ClosedSweepResult {
    // Orthonormal basis (u, v) spanning the plane perpendicular to
    // the cylinder axis. Angle 0 points along u.
    let n = normalise(&inputs.cyl.axis);
    let u_seed = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = project_out_and_normalise(&u_seed, &n);
    let v = n.cross(&u);

    let n_sec = inputs.n_sections.max(8);
    let n_arc = inputs.n_arc_samples.max(3);

    let plane_n = normalise(&inputs.plane.normal);
    let cyl_radius = inputs.cyl.radius;
    let s_plane = inputs.s_plane;
    let s_radial = inputs.s_radial;
    let centre = inputs.plane_circle_centre;

    // n_sec+1 sections so the rim wraps cleanly: section 0 is
    // θ=0, section n_sec is θ=2π, duplicating the seam.
    let total_sections = n_sec + 1;
    let mut surface_cps: Vec<Point3> = Vec::with_capacity(total_sections * n_arc);
    let mut plane_curve_cps: Vec<Point3> = Vec::with_capacity(total_sections);
    let mut cyl_curve_cps: Vec<Point3> = Vec::with_capacity(total_sections);

    let mut seam_plane_pos = Point3::new(0.0, 0.0, 0.0);
    let mut seam_cyl_pos = Point3::new(0.0, 0.0, 0.0);

    for i in 0..total_sections {
        let t = i as f64 / n_sec as f64;
        let theta = t * std::f64::consts::TAU;
        let (cos_t, sin_t) = (theta.cos(), theta.sin());
        let dir = Vec3::new(
            u.x * cos_t + v.x * sin_t,
            u.y * cos_t + v.y * sin_t,
            u.z * cos_t + v.z * sin_t,
        );

        let r_i = (inputs.radius_at)(t);

        // Plane-side contact: on the plane at radial distance
        // `cyl_radius + s_radial * r_i` from the cyl axis.
        let plane_r = cyl_radius + s_radial * r_i;
        let plane_point = Point3::new(
            centre.x + dir.x * plane_r,
            centre.y + dir.y * plane_r,
            centre.z + dir.z * plane_r,
        );

        // Cylinder-side contact: on the cylinder at height
        // `s_plane * r_i` above the plane.
        let cyl_centre_at_t = Point3::new(
            centre.x + plane_n.x * s_plane * r_i,
            centre.y + plane_n.y * s_plane * r_i,
            centre.z + plane_n.z * s_plane * r_i,
        );
        let cyl_point = Point3::new(
            cyl_centre_at_t.x + dir.x * cyl_radius,
            cyl_centre_at_t.y + dir.y * cyl_radius,
            cyl_centre_at_t.z + dir.z * cyl_radius,
        );

        plane_curve_cps.push(plane_point);
        cyl_curve_cps.push(cyl_point);

        if i == 0 {
            seam_plane_pos = plane_point;
            seam_cyl_pos = cyl_point;
        }

        // Rolling-ball centre: in the "corner" between the plane
        // and the cylinder, at (radial, axial) = (cyl_radius +
        // s_radial*r_i, s_plane*r_i).
        let ball_centre = Point3::new(
            centre.x + dir.x * plane_r + plane_n.x * s_plane * r_i,
            centre.y + dir.y * plane_r + plane_n.y * s_plane * r_i,
            centre.z + dir.z * plane_r + plane_n.z * s_plane * r_i,
        );

        // Cross-section arc centre r_i from ball_centre, passing
        // through plane_point and cyl_point.
        emit_cross_section_arc(
            &ball_centre,
            &plane_point,
            &cyl_point,
            r_i,
            n_arc,
            &mut surface_cps,
        );
    }

    let degree_u = 1usize.min(total_sections - 1);
    let degree_v = 1usize.min(n_arc - 1);
    let knots_u = clamped_uniform_knots(total_sections, degree_u);
    let knots_v = clamped_uniform_knots(n_arc, degree_v);
    let weights = vec![1.0; total_sections * n_arc];
    let surface = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points: surface_cps,
        num_u: total_sections,
        num_v: n_arc,
        weights,
    });

    let plane_contact_curve = build_closed_nurbs_curve(plane_curve_cps);
    let cyl_contact_curve = build_closed_nurbs_curve(cyl_curve_cps);

    ClosedSweepResult {
        surface,
        plane_contact_curve,
        cyl_contact_curve,
        seam_plane_pos,
        seam_cyl_pos,
    }
}

/// Append `n_arc` control points tracing the cross-section arc
/// from `plane_point` to `cyl_point`, centred at `ball_centre`
/// with radius `r_i`. The arc sweeps through the angle between
/// (plane_point - ball_centre) and (cyl_point - ball_centre) as
/// measured in the plane spanned by those two vectors.
fn emit_cross_section_arc(
    ball_centre: &Point3,
    plane_point: &Point3,
    cyl_point: &Point3,
    r_i: f64,
    n_arc: usize,
    out: &mut Vec<Point3>,
) {
    let v_bp = Vec3::new(
        plane_point.x - ball_centre.x,
        plane_point.y - ball_centre.y,
        plane_point.z - ball_centre.z,
    );
    let v_bc = Vec3::new(
        cyl_point.x - ball_centre.x,
        cyl_point.y - ball_centre.y,
        cyl_point.z - ball_centre.z,
    );
    let r_bp = v_bp.length().max(1e-15);
    let u_basis = Vec3::new(v_bp.x / r_bp, v_bp.y / r_bp, v_bp.z / r_bp);
    let d = v_bc.dot(&u_basis);
    let w_raw = Vec3::new(
        v_bc.x - d * u_basis.x,
        v_bc.y - d * u_basis.y,
        v_bc.z - d * u_basis.z,
    );
    let w_len = w_raw.length();
    if w_len < 1e-15 {
        // Degenerate — emit a linear ramp so the sweep stays a
        // valid (if flat) NURBS.
        for j in 0..n_arc {
            let s = j as f64 / (n_arc - 1) as f64;
            out.push(Point3::new(
                plane_point.x + s * (cyl_point.x - plane_point.x),
                plane_point.y + s * (cyl_point.y - plane_point.y),
                plane_point.z + s * (cyl_point.z - plane_point.z),
            ));
        }
        return;
    }
    let w_basis = Vec3::new(w_raw.x / w_len, w_raw.y / w_len, w_raw.z / w_len);
    let cos_alpha = d / r_bp;
    let sin_alpha = w_len / r_bp;
    let alpha = sin_alpha.atan2(cos_alpha);
    for j in 0..n_arc {
        let s = j as f64 / (n_arc - 1) as f64;
        let phi = s * alpha;
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        out.push(Point3::new(
            ball_centre.x + r_i * (cos_phi * u_basis.x + sin_phi * w_basis.x),
            ball_centre.y + r_i * (cos_phi * u_basis.y + sin_phi * w_basis.y),
            ball_centre.z + r_i * (cos_phi * u_basis.z + sin_phi * w_basis.z),
        ));
    }
}

fn build_closed_nurbs_curve(cps: Vec<Point3>) -> NurbsCurve {
    let n = cps.len();
    let degree = 1;
    let weights = vec![1.0; n];
    let knots = clamped_uniform_knots(n, degree);
    NurbsCurve {
        degree,
        knots,
        control_points: cps,
        weights,
    }
}

fn clamped_uniform_knots(n_cps: usize, degree: usize) -> Vec<f64> {
    let n_knots = n_cps + degree + 1;
    let interior = n_knots.saturating_sub(2 * (degree + 1));
    let mut knots: Vec<f64> = vec![0.0; degree + 1];
    knots.extend((1..=interior).map(|i| i as f64 / (interior + 1) as f64));
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    knots
}

fn normalise(v: &Vec3) -> Vec3 {
    let len = v.length().max(1e-15);
    Vec3::new(v.x / len, v.y / len, v.z / len)
}

fn project_out_and_normalise(v: &Vec3, axis: &Vec3) -> Vec3 {
    let dot = v.dot(axis);
    let raw = Vec3::new(v.x - dot * axis.x, v.y - dot * axis.y, v.z - dot * axis.z);
    normalise(&raw)
}
