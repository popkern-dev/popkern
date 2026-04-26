//! Rank decomposition of a degenerate quadric fiber from the pencil.
//!
//! After `pencil::find_degenerate_fiber` returns a rank-deficient
//! `D = λA + μB`, this module classifies `D` by rank and extracts
//! the geometric primitives:
//!
//! - **Rank 1** (`D = α · v · vᵀ`) — the squared form of a plane
//!   `v · p = 0`. Always real.
//! - **Rank 2 with mixed-sign eigenvalues** — `D = a · aᵀ − b · bᵀ
//!   = ½ ((a+b)(a−b)ᵀ + (a−b)(a+b)ᵀ)`, factoring into two distinct
//!   real planes `(a + b) · p = 0` and `(a − b) · p = 0`.
//! - **Rank 2 with same-sign eigenvalues** — `D` is semidefinite,
//!   zero set is at most a line. No real plane factorisation.
//! - **Rank 3** — projective cone with apex at the kernel
//!   eigenvector. Parameterisation deferred to stage 2b.34.
//! - **Rank 4** — not actually degenerate (should never arise here).
//! - **Zero** — the whole matrix is below tolerance.
//!
//! Each plane is returned as a 4-vector `(n_x, n_y, n_z, d)` in
//! standard form: `(n_x, n_y, n_z)` unit-length, with `d` such that
//! `n · p + d = 0` defines the plane (affine `p`, homogeneous `w = 1`).

use super::eigen::jacobi_eigh_4x4;
use super::quadric::Mat4Sym;

/// Geometric shape of a degenerate quadric, with the linear forms
/// (or apex + eigenstructure) already extracted.
#[derive(Clone, Debug)]
pub(crate) enum DegenerateShape {
    /// Squared linear form `(v · p)² = 0` → single plane.
    SinglePlane([f64; 4]),
    /// Two distinct real planes whose product is `D`.
    TwoPlanes([f64; 4], [f64; 4]),
    /// Rank-2 semidefinite — no real plane factorisation. The zero
    /// set is at most a line/point in real space.
    Rank2SemiDefinite,
    /// Projective cone with apex at the rank-1 kernel direction.
    /// `eigenvalues` and `eigenvectors` are the *non-zero* part of
    /// the spectrum (3 entries each), aligned by index. Stage 2b.34
    /// uses these to parameterise the cone's ruling lines.
    Rank3Cone {
        apex: [f64; 4],
        eigenvalues: [f64; 3],
        eigenvectors: [[f64; 4]; 3],
    },
    /// Not actually degenerate (shouldn't occur from pencil search).
    NonDegenerate,
    /// All entries below tolerance — pencil search shouldn't return
    /// this either.
    Zero,
}

/// Classify a 4×4 symmetric matrix by its rank structure.
///
/// `scale` is a length scale used to set the eigenvalue magnitude
/// threshold relative to the matrix's natural size. Pass the
/// max-absolute entry of the original `A` and `B` (or 1.0 if the
/// quadrics are already normalised).
pub(crate) fn classify_degenerate(d: &Mat4Sym, scale: f64) -> DegenerateShape {
    let (eig, q) = jacobi_eigh_4x4(d);
    let entry_max =
        d.0.iter()
            .flat_map(|row| row.iter())
            .map(|v| v.abs())
            .fold(0.0_f64, f64::max);
    let bound = entry_max.max(scale).max(1.0);
    let tol = 1e-10 * bound;

    let mut idx = [0_usize, 1, 2, 3];
    idx.sort_by(|&i, &j| {
        eig[j]
            .abs()
            .partial_cmp(&eig[i].abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let nonzero_count = idx.iter().filter(|&&i| eig[i].abs() > tol).count();

    match nonzero_count {
        0 => DegenerateShape::Zero,
        1 => {
            let v = column(&q, idx[0]);
            DegenerateShape::SinglePlane(normalise_plane(v))
        }
        2 => {
            let l1 = eig[idx[0]];
            let l2 = eig[idx[1]];
            let v1 = column(&q, idx[0]);
            let v2 = column(&q, idx[1]);
            if l1 * l2 < 0.0 {
                // Mixed signs: D = a·aᵀ − b·bᵀ with a along positive
                // eigenvalue, b along negative.
                let (a, b) = if l1 > 0.0 {
                    (scale_vec(&v1, l1.sqrt()), scale_vec(&v2, (-l2).sqrt()))
                } else {
                    (scale_vec(&v2, l2.sqrt()), scale_vec(&v1, (-l1).sqrt()))
                };
                let p1 = add_vec(&a, &b);
                let p2 = sub_vec(&a, &b);
                DegenerateShape::TwoPlanes(normalise_plane(p1), normalise_plane(p2))
            } else {
                DegenerateShape::Rank2SemiDefinite
            }
        }
        3 => {
            let apex = column(&q, idx[3]);
            let eigenvalues = [eig[idx[0]], eig[idx[1]], eig[idx[2]]];
            let eigenvectors = [column(&q, idx[0]), column(&q, idx[1]), column(&q, idx[2])];
            DegenerateShape::Rank3Cone {
                apex,
                eigenvalues,
                eigenvectors,
            }
        }
        4 => DegenerateShape::NonDegenerate,
        _ => unreachable!(),
    }
}

fn column(q: &[[f64; 4]; 4], j: usize) -> [f64; 4] {
    [q[0][j], q[1][j], q[2][j], q[3][j]]
}

fn scale_vec(v: &[f64; 4], s: f64) -> [f64; 4] {
    [v[0] * s, v[1] * s, v[2] * s, v[3] * s]
}

fn add_vec(a: &[f64; 4], b: &[f64; 4]) -> [f64; 4] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]
}

fn sub_vec(a: &[f64; 4], b: &[f64; 4]) -> [f64; 4] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]]
}

/// Normalise a plane 4-vector so the spatial part `(n_x, n_y, n_z)`
/// has unit length. If the spatial part is zero (the "plane at
/// infinity"), leave the vector untouched.
fn normalise_plane(v: [f64; 4]) -> [f64; 4] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if n < 1e-15 {
        return v;
    }
    [v[0] / n, v[1] / n, v[2] / n, v[3] / n]
}

/// Convert a normalised plane 4-vector to a `SurfaceGeometry::Plane`,
/// picking convenient `u_axis` / `v_axis` directions (the existing
/// helper `analytic_surfaces::orthonormal_basis` produces these).
pub(crate) fn plane_to_surface(p: [f64; 4]) -> crate::geometry::surface::SurfaceGeometry {
    let n = crate::geometry::Vec3::new(p[0], p[1], p[2]);
    // n · origin + d = 0 → origin = -d · n (closest point to world origin
    // on the plane, since n is unit).
    let origin = crate::geometry::Point3::new(-p[3] * p[0], -p[3] * p[1], -p[3] * p[2]);
    let (u_axis, v_axis) = crate::geometry::analytic_surfaces::orthonormal_basis(n);
    crate::geometry::surface::SurfaceGeometry::Plane {
        origin,
        normal: n,
        u_axis,
        v_axis,
    }
}
