//! §5.3 C4b.2 parallelogram + Chiyokura–Kimura twist extraction for
//! rectangular Gregory sub-patches.
//!
//! Extracted from `gregory_three_quad.rs` to keep that file under the
//! 500-line structural cap. Holds both the cubic (bidegree-3) and
//! quintic (bidegree-5) variants of the parallelogram inner-CP
//! computation plus the degree-parameterised C-K twist extractor.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};

use crate::fillet_surface::FilletSurfaceResult;

/// Compute Chiyokura–Kimura inner-CP pairs for a cubic (bidegree-3)
/// Gregory sub-patch from its 4 boundary Beziers + the incident
/// fillet strips. See `gregory_three_quad.rs` theory docs.
pub(crate) fn compute_parallelogram_inner_cps(
    boundaries: [&[Point3; 4]; 4],
    corner_strips: [Option<&FilletSurfaceResult>; 4],
) -> [(Point3, Point3); 4] {
    let b0 = boundaries[0];
    let b1 = boundaries[1];
    let b2 = boundaries[2];
    let b3 = boundaries[3];

    let corners = [b0[0], b0[3], b1[3], b2[3]];
    let raw_parallelogram = [
        parallelogram_add(b3[2], b0[1], b0[0]),
        parallelogram_add(b0[2], b1[1], b0[3]),
        parallelogram_add(b2[1], b1[2], b1[3]),
        parallelogram_add(b3[1], b2[2], b2[3]),
    ];
    let tangents: [(Vec3, Vec3); 4] = [
        (scale3(sub3(b3[2], b0[0])), scale3(sub3(b0[1], b0[0]))),
        (scale3(sub3(b0[2], b0[3])), scale3(sub3(b1[1], b1[0]))),
        (scale3(sub3(b1[2], b1[3])), scale3(sub3(b2[1], b2[0]))),
        (scale3(sub3(b2[2], b2[3])), scale3(sub3(b3[1], b3[0]))),
    ];

    let mut inner_cps: [(Point3, Point3); 4] =
        [(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)); 4];
    for (i, slot) in inner_cps.iter_mut().enumerate() {
        let corner = corners[i];
        let p_para = raw_parallelogram[i];
        let (e_u, e_v) = tangents[i];
        let p_refined = if let Some(strip) = corner_strips[i] {
            chiyokura_twisted_inner_cp_with_degree_sq(&strip.surface, corner, p_para, e_u, e_v, 9.0)
                .unwrap_or(p_para)
        } else {
            p_para
        };
        *slot = (p_refined, p_refined);
    }
    inner_cps
}

/// Deg5 variant of [`compute_parallelogram_inner_cps`]. Uses
/// degree-5 boundary 1-ring indices and the `degree² = 25` C-K
/// twist factor.
pub(crate) fn compute_parallelogram_inner_cps_deg5(
    boundaries: [&[Point3; 6]; 4],
    corner_strips: [Option<&FilletSurfaceResult>; 4],
) -> [(Point3, Point3); 4] {
    let b0 = boundaries[0];
    let b1 = boundaries[1];
    let b2 = boundaries[2];
    let b3 = boundaries[3];

    let corners = [b0[0], b0[5], b1[5], b2[5]];
    let raw_parallelogram = [
        parallelogram_add(b3[4], b0[1], b0[0]),
        parallelogram_add(b0[4], b1[1], b0[5]),
        parallelogram_add(b2[1], b1[4], b1[5]),
        parallelogram_add(b3[1], b2[4], b2[5]),
    ];
    let scale5 = |v: Vec3| Vec3::new(5.0 * v.x, 5.0 * v.y, 5.0 * v.z);
    let tangents: [(Vec3, Vec3); 4] = [
        (scale5(sub3(b3[4], b0[0])), scale5(sub3(b0[1], b0[0]))),
        (scale5(sub3(b0[4], b0[5])), scale5(sub3(b1[1], b1[0]))),
        (scale5(sub3(b1[4], b1[5])), scale5(sub3(b2[1], b2[0]))),
        (scale5(sub3(b2[4], b2[5])), scale5(sub3(b3[1], b3[0]))),
    ];

    let mut inner_cps: [(Point3, Point3); 4] =
        [(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)); 4];
    for (i, slot) in inner_cps.iter_mut().enumerate() {
        let corner = corners[i];
        let p_para = raw_parallelogram[i];
        let (e_u, e_v) = tangents[i];
        let p_refined = if let Some(strip) = corner_strips[i] {
            chiyokura_twisted_inner_cp_with_degree_sq(
                &strip.surface,
                corner,
                p_para,
                e_u,
                e_v,
                25.0,
            )
            .unwrap_or(p_para)
        } else {
            p_para
        };
        *slot = (p_refined, p_refined);
    }
    inner_cps
}

/// Chiyokura–Kimura twist-corrected inner CP, parameterised on
/// the sub-patch Bezier degree squared: `S_uv(0, 0) = degree²·(P_11
/// − P_parallelogram)`, so `P_11 = P_parallelogram + τ/degree²`.
/// Pass `degree_sq = 9` for cubic, `25` for quintic.
fn chiyokura_twisted_inner_cp_with_degree_sq(
    surface: &SurfaceGeometry,
    corner: Point3,
    p_parallelogram: Point3,
    e_u: Vec3,
    e_v: Vec3,
    degree_sq: f64,
) -> Option<Point3> {
    let (u, v) = project_point_to_uv(surface, &corner)?;
    let t_u = surface.partial_u(u, v);
    let t_v = surface.partial_v(u, v);
    let t_uu = surface.second_partial_uu(u, v);
    let t_uv = surface.second_partial_uv(u, v);
    let t_vv = surface.second_partial_vv(u, v);

    let g00 = t_u.dot(&t_u);
    let g01 = t_u.dot(&t_v);
    let g11 = t_v.dot(&t_v);
    let det = g00 * g11 - g01 * g01;
    if det.abs() < 1e-18 {
        return None;
    }
    let solve = |e: &Vec3| -> (f64, f64) {
        let r0 = t_u.dot(e);
        let r1 = t_v.dot(e);
        ((g11 * r0 - g01 * r1) / det, (g00 * r1 - g01 * r0) / det)
    };
    let (alpha_u, beta_u) = solve(&e_u);
    let (alpha_v, beta_v) = solve(&e_v);

    let cross = alpha_u * beta_v + beta_u * alpha_v;
    let tau = Vec3::new(
        alpha_u * alpha_v * t_uu.x + cross * t_uv.x + beta_u * beta_v * t_vv.x,
        alpha_u * alpha_v * t_uu.y + cross * t_uv.y + beta_u * beta_v * t_vv.y,
        alpha_u * alpha_v * t_uu.z + cross * t_uv.z + beta_u * beta_v * t_vv.z,
    );

    Some(Point3::new(
        p_parallelogram.x + tau.x / degree_sq,
        p_parallelogram.y + tau.y / degree_sq,
        p_parallelogram.z + tau.z / degree_sq,
    ))
}

fn parallelogram_add(p_a: Point3, p_b: Point3, base: Point3) -> Point3 {
    Point3::new(
        p_a.x + p_b.x - base.x,
        p_a.y + p_b.y - base.y,
        p_a.z + p_b.z - base.z,
    )
}

pub(crate) fn sub3(a: Point3, b: Point3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

fn scale3(v: Vec3) -> Vec3 {
    Vec3::new(3.0 * v.x, 3.0 * v.y, 3.0 * v.z)
}
