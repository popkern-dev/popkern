//! N-sided hole filling — a G¹ patch set over an arbitrary closed
//! boundary loop.
//!
//! The PnS bi-4 cap machinery ([`super::gregory`]) can cover
//! any `n`-sided opening, but until now it was reachable only from the
//! fillet crate's blend-corner path, which builds its input from a
//! `BlendBoundary` around a filleted vertex. A hole in an imported or
//! damaged body has no such structure — it is just a ring of edges,
//! each with one surviving neighbour face.
//!
//! This module is the general form: give it the boundary curves and
//! the tangent-plane data of whatever surrounds them, and it returns
//! the cap patches. It knows nothing about fillets, vertices, or
//! blends.
//!
//! # Input contract
//!
//! One [`FillBoundaryEdge`] per side of the opening, in order around
//! the loop, each carrying:
//!
//! - three points sampled along the edge (`start`, `mid`, `end`),
//!   consecutive edges sharing endpoints;
//! - the unit surface normal of the *surviving* neighbour face at
//!   those points, and the edge tangent there.
//!
//! From that pair the transversal direction — in the neighbour's
//! tangent plane, perpendicular to the boundary, pointing into the
//! hole — is derived. That is what makes the cap meet its
//! surroundings tangentially rather than merely touching them.
//!
//! # What it does not do
//!
//! It does not modify topology. Turning the returned patches into
//! faces is the caller's job, because the two callers do it
//! differently: the blend path pokes an existing n-gon face, while a
//! healing path has to build the fan against a boundary loop. Keeping
//! the geometry separate is what lets both share it.

use super::gregory::{
    BetaReparam, Bi2BorderStrip, Bi2TensorBorderRing, CapAssembleInput, assemble_pns_bi4_cap,
    pns_bi4_to_nurbs, pole_from_boundary, sigma_center_jet_default,
};
use super::surface::SurfaceGeometry;
use super::{Point3, Vec3};

/// Lowest and highest valence the bi-4 cap covers.
const MIN_SIDES: usize = 3;
const MAX_SIDES: usize = 12;

/// F_κ order used for the fill. Matches the blend-corner path.
const FILL_KAPPA: usize = 3;

/// Fraction of the local boundary length used as the transversal
/// derivative magnitude.
///
/// The derivative rows set how far the cap reaches before it turns —
/// too short and the cap is flat with a kink at the rim, too long and
/// it billows. A third of the side length is the value the
/// blend-corner adapter has carried in production.
const DERIV_MAG_FRAC: f64 = 1.0 / 3.0;

/// One side of the opening, with the tangent-plane data of the face
/// that survived on the other side of it.
#[derive(Debug, Clone)]
pub struct FillBoundaryEdge {
    /// Boundary curve sampled at its start, midpoint, and end.
    /// `end` must coincide with the next edge's `start`.
    pub start: Point3,
    /// Midpoint of the boundary curve.
    pub mid: Point3,
    /// End of the boundary curve.
    pub end: Point3,
    /// Unit normal of the neighbouring face at `start`.
    pub normal_start: Vec3,
    /// Unit normal of the neighbouring face at `mid`.
    pub normal_mid: Vec3,
    /// Unit normal of the neighbouring face at `end`.
    pub normal_end: Vec3,
    /// Boundary tangent at `start`, pointing along the loop.
    pub tangent_start: Vec3,
    /// Boundary tangent at `mid`.
    pub tangent_mid: Vec3,
    /// Boundary tangent at `end`.
    pub tangent_end: Vec3,
    /// Extra `(t, point)` samples along the side, `t` in `(0, 1)`.
    ///
    /// Empty means "use `mid` alone", which reproduces the
    /// interpolating bi-2 exactly. Supplying more lets the position
    /// row be fitted instead, which matters when the side is a
    /// circular arc concatenated with a line and no quadratic can
    /// pass through all of it.
    pub extra_samples: Vec<(f64, Point3)>,
}

impl FillBoundaryEdge {
    /// The samples the position row is fitted against.
    fn interior_samples(&self) -> Vec<(f64, Point3)> {
        if self.extra_samples.is_empty() {
            return vec![(0.5, self.mid)];
        }
        let mut all = vec![(0.5, self.mid)];
        all.extend(
            self.extra_samples
                .iter()
                .filter(|(t, _)| *t > 0.0 && *t < 1.0)
                .copied(),
        );
        all
    }
}

/// Why a fill could not be produced.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum NSidedFillError {
    /// The opening has too few or too many sides for the bi-4 cap.
    #[error(
        "n-sided fill: boundary has {got} side(s); supported range is [{MIN_SIDES}, {MAX_SIDES}]"
    )]
    UnsupportedSideCount {
        /// Sides supplied.
        got: usize,
    },
    /// Consecutive edges do not share an endpoint, so the boundary is
    /// not a closed loop.
    #[error(
        "n-sided fill: boundary is not closed — side {side}'s end and side {next}'s start are {gap:.3e} apart (tolerance {tolerance:.3e})"
    )]
    BoundaryNotClosed {
        /// Index of the side whose end failed to meet its successor.
        side: usize,
        /// Index of the successor.
        next: usize,
        /// Measured gap.
        gap: f64,
        /// The tolerance it had to clear.
        tolerance: f64,
    },
    /// A normal or tangent was degenerate, so no transversal
    /// direction exists there.
    #[error(
        "n-sided fill: side {side} sample {sample} has a degenerate frame (|normal x tangent| = {cross:.3e}); the neighbour surface is singular or the boundary doubles back"
    )]
    DegenerateFrame {
        /// Index of the offending side.
        side: usize,
        /// Which of the three samples (0 = start, 1 = mid, 2 = end).
        sample: usize,
        /// The measured cross-product length.
        cross: f64,
    },
    /// A side had no length, so it cannot orient anything.
    #[error("n-sided fill: side {side} is degenerate (length {length:.3e})")]
    DegenerateSide {
        /// Index of the offending side.
        side: usize,
        /// Measured chord length.
        length: f64,
    },
    /// The cap assembly declined.
    #[error("n-sided fill: cap assembly failed: {0}")]
    Assembly(String),
}

/// The patch set covering an opening.
#[derive(Debug, Clone)]
pub struct NSidedFill {
    /// One surface per side of the opening, in boundary order.
    pub surfaces: Vec<SurfaceGeometry>,
    /// The central pole all patches meet at.
    pub pole: Point3,
    /// Worst deviation between a patch corner and the boundary point
    /// it was built from. A fill that does not interpolate its own
    /// boundary is not usable, so callers gate on this.
    pub worst_boundary_gap: f64,
}

/// Build a G¹ patch set covering the opening described by `boundary`.
///
/// The patches meet each other G¹ by construction (that is what the
/// bi-4 cap's eq (2) enforces) and meet the surrounding surfaces
/// tangentially to the extent the supplied normals are accurate.
///
/// `tolerance` is the linear tolerance used to check that consecutive
/// sides actually meet — it is an explicit parameter, never global
/// state.
///
/// # Errors
///
/// See [`NSidedFillError`]. Every rejection names the offending side
/// and, where relevant, the measured quantity: a fill that cannot be
/// trusted declines rather than emitting a folded patch.
pub fn fill_n_sided_boundary(
    boundary: &[FillBoundaryEdge],
    tolerance: f64,
) -> Result<NSidedFill, NSidedFillError> {
    let n = boundary.len();
    if !(MIN_SIDES..=MAX_SIDES).contains(&n) {
        return Err(NSidedFillError::UnsupportedSideCount { got: n });
    }
    check_closed(boundary, tolerance)?;

    let mut sectors = Vec::with_capacity(n);
    for (side, edge) in boundary.iter().enumerate() {
        sectors.push(strip_for_edge(side, edge)?);
    }
    let ring = Bi2TensorBorderRing { sectors };

    let pole = pole_from_boundary(&ring);
    let sigma =
        sigma_center_jet_default(n).map_err(|e| NSidedFillError::Assembly(e.to_string()))?;
    let cap = assemble_pns_bi4_cap(&CapAssembleInput {
        boundary: ring,
        beta: BetaReparam::bi4(n),
        sigma,
        pole,
        kappa: FILL_KAPPA,
    })
    .map_err(|e| NSidedFillError::Assembly(e.to_string()))?;

    // Measure how well the cap actually interpolates the boundary it
    // was handed. The assembly is closed-form, so this is cheap and
    // catches a mis-ordered or mis-oriented boundary immediately.
    let mut worst_boundary_gap = 0.0_f64;
    for (k, edge) in boundary.iter().enumerate() {
        let patch = &cap.patches[k];
        // Sector k's rim is the column i = 4: `patch[0][4]` sits at
        // that side's start and `patch[4][4]` at its end, while
        // `patch[0][0]` is the pole. (Measured, not assumed — the
        // first cut of this check read the wrong row and reported a
        // 1.1-unit "gap" on a cap that interpolates exactly.)
        worst_boundary_gap = worst_boundary_gap
            .max(patch[0][4].distance_to(&edge.start))
            .max(patch[4][4].distance_to(&edge.end));
    }

    let surfaces = cap
        .patches
        .iter()
        .map(|patch| SurfaceGeometry::Nurbs(pns_bi4_to_nurbs(patch)))
        .collect();

    Ok(NSidedFill {
        surfaces,
        pole,
        worst_boundary_gap,
    })
}

/// Least-squares middle BB coefficient for the side's bi-2 position
/// row.
///
/// For a degree-2 Bézier with fixed ends, `P(t) = (1-t)^2 * start +
/// 2t(1-t) * b1 + t^2 * end`, so each sample gives a linear equation
/// in `b1` with weight `w(t) = 2t(1-t)`. Minimising the squared
/// residual over the samples gives
///
/// ```text
///   b1 = sum_k w_k * r_k / sum_k w_k^2,
///   r_k = sample_k - (1-t_k)^2 * start - t_k^2 * end
/// ```
///
/// With the single midpoint sample this returns the interpolating
/// coefficient exactly, so behaviour is unchanged unless a caller
/// supplies more.
fn fit_middle_coefficient(edge: &FillBoundaryEdge) -> Point3 {
    let samples = edge.interior_samples();
    let mut num = [0.0_f64; 3];
    let mut den = 0.0_f64;
    for (t, point) in samples {
        let w = 2.0 * t * (1.0 - t);
        if w <= 0.0 {
            continue;
        }
        let a = (1.0 - t) * (1.0 - t);
        let c = t * t;
        let r = [
            point.x - a * edge.start.x - c * edge.end.x,
            point.y - a * edge.start.y - c * edge.end.y,
            point.z - a * edge.start.z - c * edge.end.z,
        ];
        num[0] += w * r[0];
        num[1] += w * r[1];
        num[2] += w * r[2];
        den += w * w;
    }
    if den <= 0.0 {
        return edge.mid;
    }
    Point3::new(num[0] / den, num[1] / den, num[2] / den)
}

fn check_closed(boundary: &[FillBoundaryEdge], tolerance: f64) -> Result<(), NSidedFillError> {
    let n = boundary.len();
    for side in 0..n {
        let next = (side + 1) % n;
        let gap = boundary[side].end.distance_to(&boundary[next].start);
        if gap > tolerance {
            return Err(NSidedFillError::BoundaryNotClosed {
                side,
                next,
                gap,
                tolerance,
            });
        }
    }
    Ok(())
}

/// Build one bi-2 border strip from a boundary side.
///
/// The position row is the degree-2 Bézier through `start`, `mid`,
/// `end`; the derivative row is the transversal direction at each
/// sample, scaled by a fraction of the side length.
fn strip_for_edge(side: usize, edge: &FillBoundaryEdge) -> Result<Bi2BorderStrip, NSidedFillError> {
    let length = edge.start.distance_to(&edge.mid) + edge.mid.distance_to(&edge.end);
    if length < 1e-12 {
        return Err(NSidedFillError::DegenerateSide { side, length });
    }
    let mag = length * DERIV_MAG_FRAC;

    // Bi-2 Bézier through the side. Interpolating start / mid / end
    // exactly gives b1 = 2*mid - (start + end)/2, but that pins the
    // curve to one interior sample and lets it drift elsewhere — on a
    // corner rim (a circular arc concatenated with a line) the drift
    // measured 1.5e-2 against a 5e-3 acceptance bar.
    //
    // When the caller supplies extra samples, the middle coefficient
    // is instead the least-squares best fit over all of them, which
    // spreads the error rather than parking it away from the midpoint.
    // With no extra samples this reduces exactly to the interpolating
    // form, so the simple case is unchanged.
    let mid_bb = fit_middle_coefficient(edge);

    let samples = [
        (0_usize, edge.normal_start, edge.tangent_start),
        (1, edge.normal_mid, edge.tangent_mid),
        (2, edge.normal_end, edge.tangent_end),
    ];
    let mut deriv = [Point3::new(0.0, 0.0, 0.0); 3];
    for (sample, normal, tangent) in samples {
        // Transversal = normal x tangent: lies in the neighbour's
        // tangent plane and is perpendicular to the boundary, which is
        // exactly the direction a G1 join needs.
        let transversal = normal.cross(&tangent);
        let len = transversal.length();
        if len < 1e-9 {
            return Err(NSidedFillError::DegenerateFrame {
                side,
                sample,
                cross: len,
            });
        }
        let unit = transversal * (1.0 / len);
        deriv[sample] = Point3::new(unit.x * mag, unit.y * mag, unit.z * mag);
    }

    Ok(Bi2BorderStrip {
        pos: [edge.start, mid_bb, edge.end],
        deriv,
    })
}

#[cfg(test)]
#[path = "n_sided_fill_tests.rs"]
mod tests;
