//! Point-like supports (PH-033). A support's rotational restraint is
//! ∫|ω × d|² dA over its area (d from the support centroid): the
//! smallest eigenvalue of J = ∫(|d|² I − d dᵀ) dA gives the weakest
//! axis, and span = √(12 · λ_min / A) (the side of an equivalent square
//! patch). Measured on a cantilever clamped on a central patch of its
//! end face (support_tests.rs): a patch spanning 2.2 cells gives the
//! tip within 2% of the converged value, 2.5 cells within 5%, 1.3 cells
//! 32% off. Below `SUPPORT_SPAN_CELLS` the answer depends on the grid,
//! so the support is reported as point-like, naming the rotations it
//! barely holds.

use crate::clip::ClippedSolid;

use super::assemble::FaceMechanics;
use super::rigid::{RigidMotion, symmetric_eigen};

/// Supports spanning fewer cells than this are point-like.
pub const SUPPORT_SPAN_CELLS: f64 = 2.0;

/// The span of the supported area of one solid along its weakest axis,
/// and the rotations (about the support centroid) whose span is under
/// `min_span`. `None` when the solid has no supported area.
pub fn weak_rotations(
    clipped: &ClippedSolid,
    faces: &[FaceMechanics],
    min_span: f64,
) -> Option<(f64, Vec<RigidMotion>)> {
    let supported = |f: usize| {
        matches!(
            faces.get(f),
            Some(
                FaceMechanics::Support(_)
                    | FaceMechanics::Prescribed
                    | FaceMechanics::PrescribedTraction
            )
        )
    };
    let pieces: Vec<_> = clipped
        .pieces
        .iter()
        .filter(|p| supported(p.face.0 as usize) && p.area > 0.0)
        .collect();
    let area: f64 = pieces.iter().map(|p| p.area).sum();
    if area <= 0.0 {
        return None;
    }
    // Area centroid of the support, from the pieces' fan triangles.
    let mut c = [0.0; 3];
    for p in &pieces {
        let v = p.polygon.vertices();
        for k in 1..v.len().saturating_sub(1) {
            let t = [v[0], v[k], v[k + 1]];
            let e1 = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
            let e2 = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
            let n = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0],
            ];
            let a = 0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            for i in 0..3 {
                c[i] += a * (t[0][i] + t[1][i] + t[2][i]) / (3.0 * area);
            }
        }
    }
    // J about the centroid, exactly: each piece is fanned into
    // triangles, and ∫_T d dᵀ dA = A/12 (Σ a aᵀ + s sᵀ) with a the
    // corners relative to c and s their sum.
    let mut m = [[0.0; 3]; 3];
    for p in &pieces {
        let v = p.polygon.vertices();
        for k in 1..v.len().saturating_sub(1) {
            let t = [v[0], v[k], v[k + 1]].map(|x| [x[0] - c[0], x[1] - c[1], x[2] - c[2]]);
            let e1 = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
            let e2 = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
            let n = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0],
            ];
            let a = 0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            let sum = [0, 1, 2].map(|i| t[0][i] + t[1][i] + t[2][i]);
            for r in 0..3 {
                for q in 0..3 {
                    let corners: f64 = t.iter().map(|x| x[r] * x[q]).sum();
                    m[r][q] += a / 12.0 * (corners + sum[r] * sum[q]);
                }
            }
        }
    }
    let trace = m[0][0] + m[1][1] + m[2][2];
    let mut j = [[0.0; 6]; 6];
    for r in 0..3 {
        for q in 0..3 {
            let delta = if r == q { trace } else { 0.0 };
            j[r][q] = delta - m[r][q];
        }
    }
    // Pad to 6×6 with a block well above the 3×3 one, so the three
    // smallest eigenpairs are J's.
    let big = 10.0 * (j[0][0] + j[1][1] + j[2][2]) + 1.0;
    for (k, row) in j.iter_mut().enumerate().skip(3) {
        row[k] = big;
    }
    let (values, vectors) = symmetric_eigen(j);
    let span = |l: f64| (12.0 * l.max(0.0) / area).sqrt();
    let weakest = span(values[0]);
    let motions = (0..3)
        .filter(|&k| span(values[k]) < min_span)
        .map(|k| {
            let w = [vectors[0][k], vectors[1][k], vectors[2][k]];
            let n = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
            RigidMotion::Rotation {
                point: c,
                direction: w.map(|x| x / n),
            }
        })
        .collect();
    Some((weakest, motions))
}
