//! Rigid-mode detection (cross-cutting Thread 2): the six rigid motions
//! are linear fields, so they lie exactly in the aggregated space. The
//! elastic operator annihilates them; R_ij = r_iᵀ A r_j therefore
//! measures only what the supports restrain. Its near-null eigenvectors
//! are the free motions, measured against a support-independent
//! stiffness scale Σ_k r_kᵀ diag(A) r_k (the trace of R itself is pure
//! round-off when a body has no support at all).

use serde::Serialize;

use crate::linalg::sparse::{Csr, SparseError};

/// A motion a solid is free to make.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub enum RigidMotion {
    /// Translation along a unit direction.
    Translation {
        /// Direction.
        direction: [f64; 3],
    },
    /// Rotation about the axis through `point` along unit `direction`.
    Rotation {
        /// A point on the axis.
        point: [f64; 3],
        /// Axis direction.
        direction: [f64; 3],
    },
}

pub use crate::linalg::dense_eigen::symmetric_eigen;

/// The six rigid modes about `centre` at the given node `positions`,
/// as vectors over dofs 3·node + component.
pub fn rigid_modes(positions: &[[f64; 3]], centre: [f64; 3]) -> [Vec<f64>; 6] {
    std::array::from_fn(|m| {
        let mut r = vec![0.0; 3 * positions.len()];
        for (node, x) in positions.iter().enumerate() {
            let d = [x[0] - centre[0], x[1] - centre[1], x[2] - centre[2]];
            let u = match m {
                0..=2 => std::array::from_fn(|i| if i == m { 1.0 } else { 0.0 }),
                _ => {
                    let mut w = [0.0; 3];
                    w[m - 3] = 1.0;
                    [
                        w[1] * d[2] - w[2] * d[1],
                        w[2] * d[0] - w[0] * d[2],
                        w[0] * d[1] - w[1] * d[0],
                    ]
                }
            };
            r[3 * node..3 * node + 3].copy_from_slice(&u);
        }
        r
    })
}

/// The rigid-mode spectrum of one solid: eigenvalues of R over the
/// diagonal stiffness scale, ascending, with the motion of each
/// eigenvector.
#[derive(Debug, Clone, PartialEq)]
pub struct RigidSpectrum {
    /// λ_k / scale, ascending (0 when the scale is 0).
    pub ratios: [f64; 6],
    /// The motion of each eigenvector.
    pub motions: [RigidMotion; 6],
}

fn motion(t: [f64; 3], w: [f64; 3], centre: [f64; 3]) -> RigidMotion {
    let wn = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    let tn = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
    if wn <= 1e-6 * tn {
        RigidMotion::Translation {
            direction: t.map(|x| x / tn),
        }
    } else {
        // Chasles: the axis passes through centre + (ω × t)/|ω|².
        let wxt = [
            w[1] * t[2] - w[2] * t[1],
            w[2] * t[0] - w[0] * t[2],
            w[0] * t[1] - w[1] * t[0],
        ];
        RigidMotion::Rotation {
            point: std::array::from_fn(|i| centre[i] + wxt[i] / (wn * wn)),
            direction: w.map(|x| x / wn),
        }
    }
}

/// The rigid-mode spectrum of a solid whose reduced operator is `a` (on
/// dofs whose node positions are `positions`), about `centre`.
pub fn rigid_spectrum(
    a: &Csr,
    positions: &[[f64; 3]],
    centre: [f64; 3],
) -> Result<RigidSpectrum, SparseError> {
    let modes = rigid_modes(positions, centre);
    let mut applied = Vec::with_capacity(6);
    for r in &modes {
        let mut y = vec![0.0; r.len()];
        a.mul_vec(r, &mut y)?;
        applied.push(y);
    }
    let r: [[f64; 6]; 6] = std::array::from_fn(|i| {
        std::array::from_fn(|j| modes[i].iter().zip(&applied[j]).map(|(x, y)| x * y).sum())
    });
    let diagonal = a.diagonal();
    let scale: f64 = modes
        .iter()
        .map(|m| m.iter().zip(&diagonal).map(|(x, d)| d * x * x).sum::<f64>())
        .sum();
    let (values, vectors) = symmetric_eigen(r);
    Ok(RigidSpectrum {
        ratios: values.map(|v| if scale > 0.0 { v / scale } else { 0.0 }),
        motions: std::array::from_fn(|k| {
            motion(
                [vectors[0][k], vectors[1][k], vectors[2][k]],
                [vectors[3][k], vectors[4][k], vectors[5][k]],
                centre,
            )
        }),
    })
}

/// The motions of `spectrum` whose ratio is at most `tol`.
pub fn motions_below(spectrum: &RigidSpectrum, tol: f64) -> Vec<RigidMotion> {
    spectrum
        .ratios
        .iter()
        .zip(&spectrum.motions)
        .filter(|(r, _)| **r <= tol)
        .map(|(_, m)| *m)
        .collect()
}
