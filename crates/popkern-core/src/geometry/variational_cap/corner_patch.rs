//! The corner-patched cap (PK-619): a corner patch over each thin
//! side's corner region, and the main cap over the rest.
//!
//! A star cap gives every sector half of each of two adjacent sides, so
//! a side much shorter than its neighbours makes two lopsided sectors
//! (rim halves 3.5–6x apart on the steep trihedral apexes), and the cap
//! bends there: 6.3–50 / r_min against 0.4–2.7 in its balanced sectors.
//! Here each thin side T, with the tails of its neighbours A and B, gets
//! a four-sided patch [A tail, T, B head, K] cut off by the seam K, and
//! the main cap sees K in place of T. Every piece is an ordinary
//! [`fit_variational_cap`]; both sides of K fit G1 to the same ribbon of
//! tangent planes ([`corner_patch_seam`]), so the pieces meet G1.
//!
//! See `docs/architecture/variational-corner-fit.md`, "Corner-patched
//! cap", for the attribution, the measured constants and the sweep.

use super::corner_patch_seam::{corner_seam, length, restrict, reversed};
use super::layout;
use super::{
    VariationalCap, VariationalCapError, VariationalCapOptions, VariationalCapSide,
    VariationalCapTargets, fit_variational_cap,
};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;

/// Shape of the corner patches.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CornerPatchOptions {
    /// A side is thin when shorter than this times the mean of its two
    /// neighbours' lengths.
    pub thin_side_ratio: f64,
    /// The fraction of each neighbour the patch takes beside the thin
    /// side (< 0.5, so a side between two thin sides keeps its middle).
    pub tail_fraction: f64,
    /// The seam's handles, as a fraction of the chord between its ends.
    pub seam_bulge: f64,
}

impl Default for CornerPatchOptions {
    /// Measured 2026-10-02 on five steep apexes: tail 0.15 with bulge 0.45
    /// read 1.0–2.0 / r_min on every row's best border; bulge 0.25 is
    /// worse on every row, and tail 0.35 lets two patches pinch the
    /// bridge between them.
    fn default() -> Self {
        CornerPatchOptions {
            thin_side_ratio: 0.6,
            tail_fraction: 0.15,
            seam_bulge: 0.45,
        }
    }
}

/// One fitted piece of a corner-patched cap.
#[derive(Debug, Clone)]
pub struct CapPiece {
    /// The piece's sides, in loop order.
    pub sides: Vec<VariationalCapSide>,
    /// The piece's fit.
    pub cap: VariationalCap,
    /// Largest sampled interior principal curvature (model units).
    pub interior_curvature: f64,
}

/// The patch over one thin side's corner region.
#[derive(Debug, Clone)]
pub struct CornerPatch {
    /// The thin side's index in the opening.
    pub side: usize,
    /// The patch: sides [previous side's tail, the thin side, next
    /// side's head, the seam].
    pub piece: CapPiece,
    /// The seam K from p_B (on the next side) to p_A (on the previous
    /// side): the patch's fourth side, and the main cap's side `side`
    /// reversed.
    pub seam: NurbsCurve,
}

/// A main cap and its corner patches.
#[derive(Debug, Clone)]
pub struct CornerPatchedCap {
    /// The main cap: the opening's sides with each thin side replaced by
    /// its seam (reversed) and its neighbours shortened.
    pub main: CapPiece,
    /// One patch per thin side, in side order.
    pub corners: Vec<CornerPatch>,
}

impl CornerPatchedCap {
    /// The largest interior reading over every piece (model units).
    pub fn interior_curvature(&self) -> f64 {
        self.corners
            .iter()
            .map(|c| c.piece.interior_curvature)
            .fold(self.main.interior_curvature, f64::max)
    }
}

/// The thin sides of an opening, thinnest first, no two adjacent (a
/// thin side's neighbours lend it their tails). Each is shorter than
/// `ratio` times its neighbours' mean length.
///
/// # Errors
///
/// [`VariationalCapError::NoThinSide`] when no side qualifies.
pub fn thin_sides(
    sides: &[VariationalCapSide],
    ratio: f64,
) -> Result<Vec<usize>, VariationalCapError> {
    let n = sides.len();
    let lengths: Vec<f64> = sides.iter().map(|s| length(&s.curve)).collect();
    let mut ratios: Vec<(usize, f64)> = (0..n)
        .map(|s| {
            let mean = 0.5 * (lengths[(s + n - 1) % n] + lengths[(s + 1) % n]);
            (s, lengths[s] / mean)
        })
        .collect();
    ratios.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let mut chosen: Vec<usize> = Vec::new();
    for &(s, r) in &ratios {
        if r >= ratio {
            break;
        }
        let touches = chosen.iter().any(|&c| c == (s + 1) % n || (c + 1) % n == s);
        if !touches {
            chosen.push(s);
        }
    }
    if chosen.is_empty() {
        return Err(VariationalCapError::NoThinSide {
            shortest_ratio: ratios.first().map_or(f64::INFINITY, |r| r.1),
            threshold: ratio,
        });
    }
    chosen.sort_unstable();
    Ok(chosen)
}

/// Fit a corner-patched cap over `sides`: one corner patch per thin side
/// ([`thin_sides`]) and the main cap over the rest, each piece an
/// ordinary [`fit_variational_cap`] under `targets` and `options`.
///
/// # Errors
///
/// [`VariationalCapError::NoThinSide`]; [`VariationalCapError::CornerPatchSeam`]
/// when a seam cannot be built; [`VariationalCapError::PatchPiece`]
/// wrapping a piece's own decline.
pub fn fit_corner_patched_cap(
    sides: &[VariationalCapSide],
    targets: &VariationalCapTargets,
    options: &VariationalCapOptions,
    patch: &CornerPatchOptions,
    tolerance: f64,
) -> Result<CornerPatchedCap, VariationalCapError> {
    let n = sides.len();
    let thin = thin_sides(sides, patch.thin_side_ratio)?;
    let tail = patch.tail_fraction;
    let mut seams: Vec<(usize, NurbsCurve, SurfaceGeometry)> = Vec::with_capacity(thin.len());
    for &t in &thin {
        let (seam, ribbon) = corner_seam(sides, t, tail, patch.seam_bulge)
            .map_err(|reason| VariationalCapError::CornerPatchSeam { side: t, reason })?;
        seams.push((t, seam, ribbon));
    }
    let seam_error =
        |side: usize| move |reason: String| VariationalCapError::CornerPatchSeam { side, reason };
    let fit = |piece_sides: Vec<VariationalCapSide>, side: Option<usize>| {
        let cap = fit_variational_cap(&piece_sides, targets, options, tolerance).map_err(|e| {
            VariationalCapError::PatchPiece {
                side,
                source: Box::new(e),
            }
        })?;
        let interior_curvature = if cap.report.interior_curvature > 0.0 {
            cap.report.interior_curvature
        } else {
            layout::interior_curvature(&cap)
        };
        Ok::<_, VariationalCapError>(CapPiece {
            sides: piece_sides,
            cap,
            interior_curvature,
        })
    };
    // Every piece's sides: the corner patches', then the main cap's.
    let mut pieces: Vec<Result<Vec<VariationalCapSide>, VariationalCapError>> = seams
        .iter()
        .map(|(t, seam, ribbon)| {
            let (a, b) = (&sides[(t + n - 1) % n], &sides[(t + 1) % n]);
            Ok(vec![
                VariationalCapSide {
                    curve: restrict(&a.curve, 1.0 - tail, 1.0).map_err(seam_error(*t))?,
                    ..a.clone()
                },
                sides[*t].clone(),
                VariationalCapSide {
                    curve: restrict(&b.curve, 0.0, tail).map_err(seam_error(*t))?,
                    ..b.clone()
                },
                VariationalCapSide {
                    curve: seam.clone(),
                    neighbour: ribbon.clone(),
                    neighbour_reversed: false,
                },
            ])
        })
        .collect();
    pieces.push(main_sides(sides, &seams, &thin, tail));
    // The pieces share nothing but their seams, so they fit
    // concurrently (PK-620: the octagon's five pieces took 5-7 s one
    // after another). Each fit is a pure function of its sides, and the
    // results are read in the serial order, so the first error is the
    // serial fit's.
    let fitted = super::linalg::par_map(pieces.len(), |k| {
        let side = seams.get(k).map(|(t, _, _)| *t);
        match &pieces[k] {
            Ok(piece_sides) => Some(fit(piece_sides.clone(), side)),
            Err(_) => None,
        }
    });
    let mut corners = Vec::with_capacity(seams.len());
    let mut main = None;
    for (k, (piece_sides, fitted)) in pieces.into_iter().zip(fitted).enumerate() {
        piece_sides?;
        let piece = fitted.ok_or(VariationalCapError::Singular {
            stage: "corner patch piece",
        })??;
        match seams.get(k) {
            Some((t, seam, _)) => corners.push(CornerPatch {
                side: *t,
                piece,
                seam: seam.clone(),
            }),
            None => main = Some(piece),
        }
    }
    Ok(CornerPatchedCap {
        main: main.ok_or(VariationalCapError::Singular {
            stage: "corner patch main cap",
        })?,
        corners,
    })
}

/// The main cap's sides: each thin side replaced by its seam, reversed,
/// and its neighbours shortened by the patch's tail.
fn main_sides(
    sides: &[VariationalCapSide],
    seams: &[(usize, NurbsCurve, SurfaceGeometry)],
    thin: &[usize],
    tail: f64,
) -> Result<Vec<VariationalCapSide>, VariationalCapError> {
    let n = sides.len();
    let mut main = Vec::with_capacity(n);
    for (s, side) in sides.iter().enumerate() {
        if let Some((_, seam, ribbon)) = seams.iter().find(|(t, _, _)| *t == s) {
            main.push(VariationalCapSide {
                curve: reversed(seam),
                neighbour: ribbon.clone(),
                neighbour_reversed: false,
            });
            continue;
        }
        let lo = if thin.contains(&((s + n - 1) % n)) {
            tail
        } else {
            0.0
        };
        let hi = if thin.contains(&((s + 1) % n)) {
            1.0 - tail
        } else {
            1.0
        };
        main.push(VariationalCapSide {
            curve: restrict(&side.curve, lo, hi)
                .map_err(|reason| VariationalCapError::CornerPatchSeam { side: s, reason })?,
            ..side.clone()
        });
    }
    Ok(main)
}

#[cfg(test)]
#[path = "corner_patch_tests.rs"]
mod tests;
