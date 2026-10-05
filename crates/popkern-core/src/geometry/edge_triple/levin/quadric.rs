//! Quadric form: `SurfaceGeometry` → `4×4` symmetric matrix.
//!
//! For each analytic surface that is algebraically a quadric (degree
//! ≤ 2), produce the homogeneous matrix `A` such that
//! `[x, y, z, w] · A · [x, y, z, w]ᵀ = 0` defines the surface (with
//! `w = 1` selecting the affine slice).
//!
//! Surfaces handled:
//! - **Plane** (`(p - P) · n = 0`) — rank-1 quadric `A = v·vᵀ` with
//!   `v = (n_x, n_y, n_z, -n·P)`. Squared linear form.
//! - **Sphere** (`‖p − C‖² = ρ²`) — rank-4 generic.
//! - **Cylinder** (right circular, `(p − O)ᵀ·M·(p − O) = r²` with
//!   projector `M = I − a·aᵀ`) — rank-3 (one degenerate direction:
//!   the axis).
//! - **Cone** (right circular, `cos²θ·‖p − V‖² − ((p − V)·a)² = 0`)
//!   — rank-3 (apex is a singular point).
//!
//! Torus and NURBS surfaces are *not* quadrics (degree 4 / arbitrary)
//! and are out of scope here. The dispatch layer should call
//! `surface_to_quadric` only for `Plane | Sphere | Cylinder | Cone`.

use crate::geometry::Vec3;
use crate::geometry::surface::SurfaceGeometry;
use popkern_math::PortableMath;

/// 4×4 symmetric matrix in row-major layout.
///
/// Stored densely (16 entries) rather than as a packed upper triangle
/// — the 64-byte footprint stays cache-friendly and the symmetry is
/// maintained as an invariant by every constructor in this module.
#[derive(Copy, Clone, Debug, PartialEq)]
pub(crate) struct Mat4Sym(pub(crate) [[f64; 4]; 4]);

impl Mat4Sym {
    /// Evaluate the quadratic form `pᵀ · A · p` at affine point `(x, y, z)`
    /// (i.e., homogeneous `(x, y, z, 1)`).
    pub(crate) fn eval_affine(&self, x: f64, y: f64, z: f64) -> f64 {
        let p = [x, y, z, 1.0];
        let mut s = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                s += p[i] * self.0[i][j] * p[j];
            }
        }
        s
    }

    /// Maximum absolute deviation from symmetry (`A_ij - A_ji`).
    /// Test-only sanity check.
    #[cfg(test)]
    pub(crate) fn symmetry_residual(&self) -> f64 {
        let mut max = 0.0_f64;
        for i in 0..4 {
            for j in (i + 1)..4 {
                let d = (self.0[i][j] - self.0[j][i]).abs();
                if d > max {
                    max = d;
                }
            }
        }
        max
    }
}

/// Convert an analytic surface to its homogeneous quadric form.
///
/// Returns `None` for surfaces that are not quadrics (Torus, NURBS,
/// Gregory).
pub(crate) fn surface_to_quadric(surf: &SurfaceGeometry) -> Option<Mat4Sym> {
    match surf {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis: _,
            v_axis: _,
        } => Some(plane_quadric(*origin, *normal)),
        SurfaceGeometry::Sphere { center, radius } => Some(sphere_quadric(*center, *radius)),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some(cylinder_quadric(*origin, *axis, *radius)),
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => Some(cone_quadric(*apex, *axis, *half_angle)),
        _ => None,
    }
}

fn plane_quadric(origin: crate::geometry::Point3, normal: Vec3) -> Mat4Sym {
    // Plane equation: n·(p − P) = 0  →  n·p − n·P = 0.
    // Squared form: (n·p − n·P)² = pᵀ·(v·vᵀ)·p with
    //   v = (n_x, n_y, n_z, −n·P)
    let n = normal;
    let np = n.x * origin.x + n.y * origin.y + n.z * origin.z;
    let v = [n.x, n.y, n.z, -np];
    let mut a = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            a[i][j] = v[i] * v[j];
        }
    }
    Mat4Sym(a)
}

fn sphere_quadric(center: crate::geometry::Point3, radius: f64) -> Mat4Sym {
    // ‖p − C·w‖² − ρ²·w² = 0
    // = pᵀ I p − 2w·Cᵀp + w²·(CᵀC − ρ²)
    let c = [center.x, center.y, center.z];
    let cdot = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    let mut a = [[0.0; 4]; 4];
    for i in 0..3 {
        a[i][i] = 1.0;
        a[i][3] = -c[i];
        a[3][i] = -c[i];
    }
    a[3][3] = cdot - radius * radius;
    Mat4Sym(a)
}

fn cylinder_quadric(origin: crate::geometry::Point3, axis: Vec3, radius: f64) -> Mat4Sym {
    // M = I − a·aᵀ  (projector onto plane perp to axis)
    // (p − O·w)ᵀ·M·(p − O·w) − r²·w² = 0
    // = pᵀ M p − 2w·OᵀMp + w²·(OᵀMO − r²)
    let m = projector_perp_to(axis);
    let o = [origin.x, origin.y, origin.z];
    let mo = mat3_vec3(&m, &o);
    let omo = o[0] * mo[0] + o[1] * mo[1] + o[2] * mo[2];
    let mut a = [[0.0; 4]; 4];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = m[i][j];
        }
        a[i][3] = -mo[i];
        a[3][i] = -mo[i];
    }
    a[3][3] = omo - radius * radius;
    Mat4Sym(a)
}

fn cone_quadric(apex: crate::geometry::Point3, axis: Vec3, half_angle: f64) -> Mat4Sym {
    // Q3 = cos²θ · I − a·aᵀ
    // (p − V·w)ᵀ·Q3·(p − V·w) = 0
    // = pᵀ Q3 p − 2w·VᵀQ3p + w²·VᵀQ3V
    let cos2 = half_angle.pcos().powi(2);
    let mut q3 = [[0.0; 3]; 3];
    let a_arr = [axis.x, axis.y, axis.z];
    for i in 0..3 {
        for j in 0..3 {
            q3[i][j] = -a_arr[i] * a_arr[j];
            if i == j {
                q3[i][j] += cos2;
            }
        }
    }
    let v = [apex.x, apex.y, apex.z];
    let qv = mat3_vec3(&q3, &v);
    let vqv = v[0] * qv[0] + v[1] * qv[1] + v[2] * qv[2];
    let mut a = [[0.0; 4]; 4];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = q3[i][j];
        }
        a[i][3] = -qv[i];
        a[3][i] = -qv[i];
    }
    a[3][3] = vqv;
    Mat4Sym(a)
}

fn projector_perp_to(axis: Vec3) -> [[f64; 3]; 3] {
    let a = [axis.x, axis.y, axis.z];
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = -a[i] * a[j];
            if i == j {
                m[i][j] += 1.0;
            }
        }
    }
    m
}

fn mat3_vec3(m: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i] += m[i][j] * v[j];
        }
    }
    out
}
