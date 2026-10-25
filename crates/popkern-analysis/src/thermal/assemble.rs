//! Per-solid assembly of the steady conduction system, in mm units.
//!
//! Volume terms come from the fragment moments; boundary terms from the
//! surface moments of each clipped facet piece, by the condition of its
//! face. The system is assembled on the active nodes and reduced by the
//! AgFEM extension.

use crate::aggregate::Aggregation;
use crate::clip::ClippedSolid;
use crate::element::poly::Poly27;
use crate::element::{ScalarTables, SurfaceMoments, surface_moments};
use crate::grid::Grid;
use crate::linalg::sparse::{Csr, SparseError};
use crate::moments::{SolidMoments, moment_index};

/// A boundary condition on one face, in mm units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MmCondition {
    /// No condition: insulated.
    Insulated,
    /// T = t0 + gradient·x (°C, °C/mm), imposed by Nitsche.
    Fixed {
        /// Value at the origin, °C.
        t0: f64,
        /// Gradient, °C/mm.
        gradient: [f64; 3],
    },
    /// Heat out = h (T − ambient), h in W/(mm²·K).
    Convection {
        /// Coefficient, W/(mm²·K).
        h: f64,
        /// Ambient, °C.
        ambient: f64,
    },
    /// Heat in = q, W/mm².
    Flux {
        /// Flux density, W/mm².
        q: f64,
    },
    /// Held to the prescribed field (`SolidInput::prescribed`), affine
    /// per piece; insulated without one.
    Prescribed,
    /// Convection h (W/(mm²·K)) plus grey radiation of emittance ε to an
    /// ambient (°C): per piece, one convection combining h with
    /// radiation's tangent at the piece's temperature in
    /// `SolidInput::previous` (the ambient without one).
    Radiating {
        /// Convective coefficient, W/(mm²·K).
        h: f64,
        /// Emittance.
        emittance: f64,
        /// Ambient, °C.
        ambient: f64,
    },
}

/// The node-level and reduced system of one solid.
#[derive(Debug, Clone)]
pub struct SolidSystem {
    /// Reduced operator Eᵀ A E.
    pub a: Csr,
    /// Reduced load Eᵀ f.
    pub b: Vec<f64>,
}

/// Everything the assembly reads for one solid.
pub struct SolidInput<'a> {
    /// The grid.
    pub grid: &'a Grid,
    /// Clipped facet pieces.
    pub clipped: &'a ClippedSolid,
    /// Fragment moments.
    pub moments: &'a SolidMoments,
    /// Aggregation and extension.
    pub aggregation: &'a Aggregation,
    /// Condition of each face, indexed by `FaceId`.
    pub conditions: &'a [MmCondition],
    /// The field `Prescribed` faces are held to.
    pub prescribed: Option<&'a super::TemperatureField>,
    /// The previous iterate, for `Radiating` faces' coefficients.
    pub previous: Option<&'a super::TemperatureField>,
    /// Geometry-only boundary data computed once for repeated passes.
    pub cache: Option<&'a BoundaryCache>,
}

/// One solid's boundary data that only geometry and the fixed faces
/// decide, kept across the passes of a radiating solve (PH-107): each
/// entry's Nitsche β, each convective or radiating piece's facet mass
/// and load (whose coefficient changes per pass), and each fixed or
/// prescribed piece's surface moments (a submodel's caps).
#[derive(Debug, Clone, Default)]
pub struct BoundaryCache {
    /// β per entry ([`super::stability::entry_betas`]).
    pub betas: Vec<f64>,
    /// Per piece (empty when not cached).
    pub facets: Vec<Option<CachedPiece>>,
}

/// A piece's cached boundary geometry.
#[derive(Debug, Clone)]
pub enum CachedPiece {
    /// ∫N_aN_b and ∫N_a over a convective or radiating piece.
    Facet(Box<([[f64; 8]; 8], [f64; 8])>),
    /// The surface moments of a fixed or prescribed piece.
    Moments(Box<SurfaceMoments>),
}

/// The [`BoundaryCache`] of `input`; facets only when `facets`.
pub fn boundary_cache(
    t: &ScalarTables,
    input: &SolidInput<'_>,
    physics: &Physics,
    facets: bool,
) -> BoundaryCache {
    let betas = super::stability::entry_betas(t, input, physics);
    let facets = if facets {
        let (grid, h) = (input.grid, input.grid.cell_size());
        input
            .clipped
            .pieces
            .iter()
            .map(|piece| {
                let condition = input.conditions.get(piece.face.0 as usize).copied()?;
                let lo = grid.cell_min(grid.cell_ijk(piece.cell).unwrap_or_default());
                let s = || surface_moments(std::iter::once(piece), lo, h);
                match condition {
                    MmCondition::Convection { .. } | MmCondition::Radiating { .. } => {
                        let s = s();
                        Some(CachedPiece::Facet(Box::new((
                            t.facet_mass(&s),
                            t.facet_load(&s),
                        ))))
                    }
                    MmCondition::Fixed { .. } | MmCondition::Prescribed => {
                        Some(CachedPiece::Moments(Box::new(s())))
                    }
                    MmCondition::Insulated | MmCondition::Flux { .. } => None,
                }
            })
            .collect()
    } else {
        Vec::new()
    };
    BoundaryCache { betas, facets }
}

/// One piece's boundary contributions: from the cache when it holds the
/// piece's facet terms (a convective condition), else computed.
fn piece_boundary(
    t: &ScalarTables,
    input: &SolidInput<'_>,
    (piece_index, piece): (usize, &crate::clip::FacetPiece),
    condition: MmCondition,
    physics: &Physics,
    (lo, h): ([f64; 3], f64),
) -> Option<([[f64; 8]; 8], [f64; 8])> {
    let cached = input
        .cache
        .and_then(|c| c.facets.get(piece_index))
        .and_then(|f| f.as_ref());
    match (condition, cached) {
        (MmCondition::Convection { h: coeff, ambient }, Some(CachedPiece::Facet(facet))) => {
            let (mass, load) = &**facet;
            Some((
                mass.map(|row| row.map(|v| coeff * v)),
                load.map(|v| coeff * ambient * v),
            ))
        }
        (_, Some(CachedPiece::Moments(s))) => boundary(t, condition, s, physics, lo, h),
        _ => {
            let s = surface_moments(std::iter::once(piece), lo, h);
            boundary(t, condition, &s, physics, lo, h)
        }
    }
}

/// The condition of `piece` (on a face with `condition`): a
/// `Prescribed` face's held to the field's tangent plane at the piece's
/// centroid (insulated where there is no field).
pub(super) fn piece_condition(
    condition: MmCondition,
    input: &SolidInput<'_>,
    piece: &crate::clip::FacetPiece,
) -> MmCondition {
    let v = piece.polygon.vertices();
    let nv = v.len().max(1) as f64;
    let c = [0, 1, 2].map(|a| v.iter().map(|x| x[a]).sum::<f64>() / nv);
    if let MmCondition::Radiating {
        h,
        emittance,
        ambient,
    } = condition
    {
        let t = input
            .previous
            .and_then(|f| f.temperature_at(c).ok().flatten())
            .unwrap_or(ambient);
        // Convection h toward the ambient plus radiation's tangent toward
        // its effective ambient, as one convection.
        let (slope, effective) = super::radiation::tangent(t, ambient, emittance);
        let slope = 1e-6 * slope;
        return MmCondition::Convection {
            h: h + slope,
            ambient: (h * ambient + slope * effective) / (h + slope),
        };
    }
    if condition != MmCondition::Prescribed {
        return condition;
    }
    match input
        .prescribed
        .and_then(|f| f.value_and_gradient_at(c).ok().flatten())
    {
        Some((t, g)) => MmCondition::Fixed {
            t0: t - (0..3).map(|a| g[a] * c[a]).sum::<f64>(),
            gradient: g,
        },
        None => MmCondition::Insulated,
    }
}

/// Physics constants of the solve, in mm units.
#[derive(Debug, Clone, Copy)]
pub struct Physics {
    /// Conductivity, W/(mm·K).
    pub k: f64,
    /// Nitsche β, W/(mm²·K).
    pub beta: f64,
    /// Volumetric source density, W/mm³.
    pub source: f64,
}

/// The linear field t0 + G·x as a polynomial in the local coordinates
/// of the cell with minimum corner `lo`.
pub fn linear_local(t0: f64, gradient: [f64; 3], lo: [f64; 3], h: f64) -> Poly27 {
    let mut p = [0.0; 27];
    let centre = lo.map(|c| c + 0.5 * h);
    p[0] = t0 + (0..3).map(|a| gradient[a] * centre[a]).sum::<f64>();
    p[moment_index(1, 0, 0)] = 0.5 * h * gradient[0];
    p[moment_index(0, 1, 0)] = 0.5 * h * gradient[1];
    p[moment_index(0, 0, 1)] = 0.5 * h * gradient[2];
    p
}

/// One cell's boundary contributions for one piece's condition.
fn boundary(
    t: &ScalarTables,
    condition: MmCondition,
    s: &SurfaceMoments,
    physics: &Physics,
    lo: [f64; 3],
    h: f64,
) -> Option<([[f64; 8]; 8], [f64; 8])> {
    match condition {
        MmCondition::Insulated | MmCondition::Prescribed | MmCondition::Radiating { .. } => None,
        MmCondition::Fixed { t0, gradient } => {
            let g = linear_local(t0, gradient, lo, h);
            Some((
                t.nitsche(s, physics.k, physics.beta, h),
                t.nitsche_load(s, &g, physics.k, physics.beta, h),
            ))
        }
        MmCondition::Convection { h: coeff, ambient } => {
            let m = t.facet_mass(s).map(|row| row.map(|v| coeff * v));
            let f = t.facet_load(s).map(|v| coeff * ambient * v);
            Some((m, f))
        }
        MmCondition::Flux { q } => Some(([[0.0; 8]; 8], t.facet_load(s).map(|v| q * v))),
    }
}

/// Assembles and reduces one solid's system.
pub fn assemble_solid(
    t: &ScalarTables,
    input: &SolidInput<'_>,
    physics: &Physics,
) -> Result<SolidSystem, SparseError> {
    let grid = input.grid;
    let h = grid.cell_size();
    let e = &input.aggregation.extension;
    let pieces = &input.clipped.pieces;
    let mut triplets = Vec::new();
    let mut load = vec![0.0; e.nodes.len()];
    let moments = input.moments;
    let computed;
    let betas = match input.cache {
        Some(c) => &c.betas,
        None => {
            computed = super::stability::entry_betas(t, input, physics);
            &computed
        }
    };
    for (entry, cell) in moments.cells.iter().enumerate() {
        let physics = &Physics {
            beta: betas.get(entry).copied().unwrap_or(physics.beta),
            ..*physics
        };
        let ijk = grid.cell_ijk(cell.cell).unwrap_or_default();
        let lo = grid.cell_min(ijk);
        let nodes = moments
            .node_keys(grid, entry)
            .unwrap_or_default()
            .map(|n| e.position(n).unwrap_or(usize::MAX));
        if nodes.contains(&usize::MAX) {
            return Err(SparseError::Shape("cell node missing from the extension"));
        }
        let mut k = t.conduction(&cell.m, physics.k, h);
        let mut f = t.volume_load(&cell.m).map(|v| physics.source * v);
        let start = pieces.partition_point(|p| p.cell < cell.cell);
        let end = pieces.partition_point(|p| p.cell <= cell.cell);
        for (index, piece) in pieces[start..end].iter().enumerate() {
            if moments.entry_of_piece(start + index, cell.cell) != Some(entry) {
                continue;
            }
            let condition = input
                .conditions
                .get(piece.face.0 as usize)
                .copied()
                .unwrap_or(MmCondition::Insulated);
            let condition = piece_condition(condition, input, piece);
            let terms = piece_boundary(
                t,
                input,
                (start + index, piece),
                condition,
                physics,
                (lo, h),
            );
            if let Some((bk, bf)) = terms {
                for a in 0..8 {
                    f[a] += bf[a];
                    for b in 0..8 {
                        k[a][b] += bk[a][b];
                    }
                }
            }
        }
        for a in 0..8 {
            load[nodes[a]] += f[a];
            for b in 0..8 {
                triplets.push((nodes[a], nodes[b], k[a][b]));
            }
        }
    }
    let a_nodes = Csr::from_triplets(e.nodes.len(), e.nodes.len(), &triplets)?;
    let ecsr = e.to_csr()?;
    let a = ecsr.transpose().mul(&a_nodes.mul(&ecsr)?)?;
    Ok(SolidSystem {
        a,
        b: e.restrict(&load),
    })
}

/// The reduced capacity matrix Eᵀ C E of one solid, C = ∫ ρc N_a N_b dV
/// with `rho_c` in J/(mm³·K).
pub fn assemble_capacity(
    t: &ScalarTables,
    input: &SolidInput<'_>,
    rho_c: f64,
) -> Result<Csr, SparseError> {
    let grid = input.grid;
    let e = &input.aggregation.extension;
    let moments = input.moments;
    let mut triplets = Vec::new();
    for (entry, cell) in moments.cells.iter().enumerate() {
        let nodes = moments
            .node_keys(grid, entry)
            .unwrap_or_default()
            .map(|n| e.position(n).unwrap_or(usize::MAX));
        if nodes.contains(&usize::MAX) {
            return Err(SparseError::Shape("cell node missing from the extension"));
        }
        let m = t.mass(&cell.m);
        for a in 0..8 {
            for b in 0..8 {
                triplets.push((nodes[a], nodes[b], rho_c * m[a][b]));
            }
        }
    }
    let c_nodes = Csr::from_triplets(e.nodes.len(), e.nodes.len(), &triplets)?;
    let ecsr = e.to_csr()?;
    ecsr.transpose().mul(&c_nodes.mul(&ecsr)?)
}

/// Heat leaving through each face (W, indexed by `FaceId`; negative
/// enters) for nodal temperatures `temps` (parallel to the extension's
/// nodes), added into `out`. On fixed faces this is the Nitsche flux
/// −k∇T·n + β(T − g); on convective faces h(T − T∞).
pub fn boundary_heat(
    t: &ScalarTables,
    input: &SolidInput<'_>,
    physics: &Physics,
    temps: &[f64],
    out: &mut [f64],
) {
    let grid = input.grid;
    let h = grid.cell_size();
    let e = &input.aggregation.extension;
    let pieces = &input.clipped.pieces;
    let moments = input.moments;
    let computed;
    let betas = match input.cache {
        Some(c) => &c.betas,
        None => {
            computed = super::stability::entry_betas(t, input, physics);
            &computed
        }
    };
    for (entry, cell) in moments.cells.iter().enumerate() {
        let physics = &Physics {
            beta: betas.get(entry).copied().unwrap_or(physics.beta),
            ..*physics
        };
        let ijk = grid.cell_ijk(cell.cell).unwrap_or_default();
        let lo = grid.cell_min(ijk);
        let tv: [f64; 8] = moments
            .node_keys(grid, entry)
            .unwrap_or_default()
            .map(|n| e.position(n).map_or(0.0, |p| temps[p]));
        let start = pieces.partition_point(|p| p.cell < cell.cell);
        let end = pieces.partition_point(|p| p.cell <= cell.cell);
        for (index, piece) in pieces[start..end].iter().enumerate() {
            if moments.entry_of_piece(start + index, cell.cell) != Some(entry) {
                continue;
            }
            let condition = input
                .conditions
                .get(piece.face.0 as usize)
                .copied()
                .unwrap_or(MmCondition::Insulated);
            let condition = piece_condition(condition, input, piece);
            let terms = piece_boundary(
                t,
                input,
                (start + index, piece),
                condition,
                physics,
                (lo, h),
            );
            let Some((bk, bf)) = terms else {
                continue;
            };
            if matches!(condition, MmCondition::Flux { .. }) {
                continue;
            }
            // Σ_a (K T − f)_a: the heat this piece passes out.
            let leaving: f64 = (0..8)
                .map(|a| (0..8).map(|b| bk[a][b] * tv[b]).sum::<f64>() - bf[a])
                .sum();
            if let Some(slot) = out.get_mut(piece.face.0 as usize) {
                *slot += leaving;
            }
        }
    }
}
