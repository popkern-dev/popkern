//! Contact interfaces between solids (plan 4.12 item 3, PK-1073).
//!
//! Two solids touching over a patch of surface stay separate bodies,
//! each with its own dofs, and are coupled across the patch: a thermal
//! contact conducts heat in proportion to the temperature jump (a
//! contact conductance), a sliding contact holds the two surfaces
//! together along the patch normal and leaves them free to slide
//! (frictionless, no separation). Both are a penalty on the jump of a
//! field across the patch, ∫_Γ (Δu)ᵀ P (Δv) dA with Δu = u_a − u_b:
//! P = h_c for heat, P = β n nᵀ for sliding.
//!
//! The patch comes from the host, which has the bodies' cell complex:
//! the faces two bodies share, tessellated, wound with their normal out
//! of the first solid into the second. Both solids live on one grid, so
//! at a point of the patch they share the cell's trilinear functions on
//! their own dof copies: each patch piece's surface moments give one
//! 8×8 facet mass M per piece, placed on the solids' fragments there
//! (each found by stepping off the patch into its own material) as
//! [[P M, −P M], [−P M, P M]], then reduced through each solid's
//! extension. A piece where either solid has no active cell (the patch
//! lies off one of the bodies) declines: the coupling would hold
//! nothing.

use crate::aggregate::Extension;
use crate::clip::{FacetPiece, clip_triangle};
use crate::element::{SurfaceMoments, surface_moments};
use crate::error::AnalysisError;
use crate::grid::Grid;
use crate::linalg::sparse::Csr;
use crate::model::{FaceId, SolidId};
use crate::moments::SolidMoments;

/// A patch two solids share.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactPatch {
    /// The two solids; the patch normal points out of the first.
    pub solids: (SolidId, SolidId),
    /// The patch, triangulated, wound with its normal out of
    /// `solids.0` (mm).
    pub triangles: Vec<[[f64; 3]; 3]>,
}

impl ContactPatch {
    /// The patch's area, mm².
    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let (u, v) = (sub(t[1], t[0]), sub(t[2], t[0]));
                let c = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
            })
            .sum()
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// One solid's side of a coupling: its moments, extension (reduced
/// through `reduce`, scalar or vector) and first free dof.
pub(crate) struct Side<'a> {
    pub moments: &'a SolidMoments,
    pub extension: &'a Extension,
    pub reduce: &'a Csr,
    pub offset: usize,
}

/// Each corner's active-node row in `side`, for the patch piece in
/// `cell` at `centre`: the side's fragment in that cell, found by
/// stepping `inward` (a patch on a grid plane lies in a face-only entry
/// of the body beyond it, measured at cell sizes dividing the offset).
/// A corner the entry lacks has `None`.
fn side_rows(
    grid: &Grid,
    side: &Side<'_>,
    cell: usize,
    centre: [f64; 3],
    inward: [f64; 3],
    corners: [usize; 8],
) -> Option<[Option<usize>; 8]> {
    let entry = side.moments.entry_at(grid, cell, centre, Some(inward))?;
    let keys = side.moments.node_keys(grid, entry)?;
    let copies = crate::moments::fragments::NODE_COPIES;
    Some(corners.map(|node| {
        keys.iter()
            .find(|&&k| k / copies == node)
            .and_then(|&k| side.extension.position(k))
    }))
}

/// The reduced coupling of `patch` between `sides` as global triplets,
/// with `dofs` dofs per node and per-piece node block `block(s, n)`
/// (8·dofs square, dof `dofs·a + i` of corner `a`).
pub(crate) fn coupling_triplets(
    grid: &Grid,
    patch: &ContactPatch,
    sides: (&Side<'_>, &Side<'_>),
    dofs: usize,
    block: impl Fn(&SurfaceMoments, [f64; 3]) -> Vec<Vec<f64>>,
) -> Result<Vec<(usize, usize, f64)>, AnalysisError> {
    let declined = |what| AnalysisError::Unsupported { what };
    let mut pieces: Vec<FacetPiece> = Vec::new();
    for (t, v) in patch.triangles.iter().enumerate() {
        clip_triangle(grid, (t as u32, FaceId(0)), *v, &mut pieces)
            .map_err(|_| declined("a contact patch reaches outside the grid"))?;
    }
    let h = grid.cell_size();
    let n_nodes = |s: &Side<'_>| s.extension.nodes.len();
    // Node-level blocks: [aa, ab, ba, bb].
    let mut blocks: [Vec<(usize, usize, f64)>; 4] = Default::default();
    for piece in &pieces {
        let ijk = grid
            .cell_ijk(piece.cell)
            .ok_or(declined("a contact patch piece lies in no cell"))?;
        let lo = grid.cell_min(ijk);
        let v = piece.polygon.vertices();
        let nv = v.len().max(1) as f64;
        let centre = [0, 1, 2].map(|a| v.iter().map(|x| x[a]).sum::<f64>() / nv);
        let n = piece.normal;
        let into_a = n.map(|c| -c);
        let corners = grid
            .cell_nodes(ijk)
            .ok_or(declined("a contact cell has no nodes"))?;
        let rows =
            |side: &Side<'_>, inward: [f64; 3]| -> Result<[Option<usize>; 8], AnalysisError> {
                side_rows(grid, side, piece.cell, centre, inward, corners)
                    .ok_or(declined("a contact patch lies off one of its bodies"))
            };
        let (ra, rb) = (rows(sides.0, into_a)?, rows(sides.1, n)?);
        let s = surface_moments(std::iter::once(piece), lo, h);
        let k = block(&s, n);
        // A corner off the side's cell lies off the patch's cell face
        // (the patch is on the face both cells share): its functions
        // vanish there, so its terms must too.
        let scale = k.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
        for a in 0..8 {
            for b in 0..8 {
                for i in 0..dofs {
                    for j in 0..dofs {
                        let value = k[dofs * a + i][dofs * b + j];
                        if value == 0.0 {
                            continue;
                        }
                        let (Some(ra_a), Some(rb_a), Some(ra_b), Some(rb_b)) =
                            (ra[a], rb[a], ra[b], rb[b])
                        else {
                            if value.abs() > 1e-9 * scale {
                                return Err(declined(
                                    "a contact patch piece leaves the cells of its bodies",
                                ));
                            }
                            continue;
                        };
                        let (pa, pb) = (dofs * ra_a + i, dofs * rb_a + i);
                        let (qa, qb) = (dofs * ra_b + j, dofs * rb_b + j);
                        blocks[0].push((pa, qa, value));
                        blocks[1].push((pa, qb, -value));
                        blocks[2].push((pb, qa, -value));
                        blocks[3].push((pb, qb, value));
                    }
                }
            }
        }
    }
    let internal = |e: crate::linalg::sparse::SparseError| AnalysisError::Internal {
        stage: "contact coupling",
        detail: e.to_string(),
    };
    let mut out = Vec::new();
    let pairs = [
        (sides.0, sides.0),
        (sides.0, sides.1),
        (sides.1, sides.0),
        (sides.1, sides.1),
    ];
    for ((row, col), trips) in pairs.into_iter().zip(&blocks) {
        let c = Csr::from_triplets(dofs * n_nodes(row), dofs * n_nodes(col), trips)
            .map_err(internal)?;
        let reduced = row
            .reduce
            .transpose()
            .mul(&c.mul(col.reduce).map_err(internal)?)
            .map_err(internal)?;
        for r in 0..reduced.rows() {
            for (c, value) in reduced.row(r) {
                out.push((row.offset + r, col.offset + c, value));
            }
        }
    }
    Ok(out)
}
