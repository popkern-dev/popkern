//! Stress recovery and error estimation (map §1.6, cross-cutting
//! Thread 4).
//!
//! 1. Each cell's stress field σ_h(ξ) is sampled at its fragment
//!    centroid.
//! 2. SPR fits a linear field per node over its 2×2×2 patch. Nodes
//!    whose patch has at least `MIN_PATCH` well-spread samples use their
//!    own fit; other nodes the fit of the nearest such node within two
//!    node steps (map §1.6).
//! 3. A display vertex interpolates nodal σ* trilinearly; where a node
//!    of its cell has no σ*, the vertex takes the cell's own EAS field
//!    σ_h at the vertex (exact in pure bending on full cells). Constant
//!    averaging is never used: it reports the neutral axis on a wall in
//!    bending.
//! 4. ZZ: e_K² = ∫ (σ* − σ_h) : C⁻¹ : (σ* − σ_h) per cell, exact by
//!    moments. In-cell vertices and thin walls get no estimate.

pub mod cell;
pub mod core;
pub mod interior;
mod locate;
mod normal;
pub mod sample;
pub mod spr;
pub mod tensor;
mod true_error;
pub mod wall;

use serde::Serialize;

use crate::aggregate::Extension;
use crate::clip::ClippedSolid;
use crate::element::VectorTables;
use crate::element::poly::{Poly27, integrate, mul, shape};
use crate::error::AnalysisError;
use crate::grid::Grid;
use crate::model::SolidMesh;
use crate::moments::SolidMoments;
use crate::moments::fragments::NODE_COPIES;

use cell::{DerivativeTables, StressPoly, centroid, stress_at, stress_poly};
use locate::{locate, vertex_inward};
use sample::{Sampler, SurfaceSample, surface_samples};
use spr::{LinearFit, Sample, fit};
use tensor::{Voigt, principal};

/// Fewest patch samples for a node's own fit (map §1.6: "at least 6
/// well-spread samples").
pub const MIN_PATCH: usize = 6;

/// Only cells at least this full sample a patch: smaller fragments have
/// no dofs of their own under AgFEM, and a patch made of slivers passes
/// a count test while carrying little information.
pub const PATCH_FILL: f64 = 0.25;

/// Walls thinner than this many cells get `UnderResolved` and no error
/// estimate (E4, PH-023).
pub const THIN_WALL_CELLS: f64 = 1.5;

/// How a vertex's stress was obtained.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Recovery {
    /// Trilinear interpolation of nodal SPR values.
    Patch,
    /// The cell's own EAS field at the vertex (the fallback).
    InCell,
}

/// Which stress the fallback reports (production uses `InCell`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fallback {
    /// The cell's EAS field at the vertex.
    InCell,
    /// Verification only: the in-cell field at every vertex, to measure
    /// the fallback where patches would otherwise be used.
    ForceInCell,
    /// Verification only: the cell's centroid sample at every vertex
    /// (constant averaging), to show what it gets wrong.
    ForceConstant,
}

/// The recovered stress at one display vertex.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct VertexStress {
    /// Stress, MPa (Voigt, tensor shears).
    pub stress: Voigt,
    /// Von Mises, MPa.
    pub von_mises: f64,
    /// Largest principal stress, MPa.
    pub max_principal: f64,
    /// How it was obtained.
    pub recovery: Recovery,
    /// Local relative ZZ error of the vertex's cell; `None` where no
    /// estimate is honest (in-cell fallback or a thin wall).
    pub error: Option<f64>,
}

/// Recovered stresses of one solid.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidStresses {
    /// One per mesh vertex.
    pub vertices: Vec<VertexStress>,
    /// The surface at cell resolution (every facet piece's centroid and
    /// corners), where peaks are searched: display triangles can be far
    /// larger than a cell.
    pub surface: Vec<SurfaceSample>,
    /// Global relative ZZ error over cells with an estimate.
    pub global_error: Option<f64>,
    /// Vertices that used the in-cell fallback.
    pub fallbacks: usize,
    /// Thinnest wall found below the threshold, mm, if any.
    pub thin_wall: Option<f64>,
    /// Largest principal stress over the raw fragment-centroid samples
    /// (before recovery), for comparison.
    pub raw_peak_principal: f64,
    /// Relative true error ‖σ_exact − σ_h‖ / ‖σ_exact‖ when an exact
    /// field was supplied (verification only).
    pub true_error: Option<f64>,
    /// The infill core's worst stress, when the solid has a core cell.
    pub core_peak: Option<core::CorePeak>,
    /// Wall samples at the core interface inside an infilled part,
    /// which the wall checks read alongside `surface` (empty without
    /// infill).
    pub interior: Vec<interior::InteriorSample>,
}

/// An exact stress field, for verification of the error estimate.
pub type ExactStress<'a> = &'a dyn Fn([f64; 3]) -> Voigt;

/// Everything recovery reads for one solid.
pub struct RecoveryInput<'a> {
    /// The grid.
    pub grid: &'a Grid,
    /// Fragment moments.
    pub moments: &'a SolidMoments,
    /// Active nodes and dof map.
    pub extension: &'a Extension,
    /// Nodal displacements, 3 per active node.
    pub u: &'a [f64],
    /// Triangles of the display mesh on faces whose tessellation does
    /// not match the face (no wall readings start or end on them),
    /// parallel to its triangles.
    pub unreliable: &'a [bool],
    /// The display mesh.
    pub mesh: &'a SolidMesh,
    /// Clipped pieces (for wall thickness).
    pub clipped: &'a ClippedSolid,
    /// Young's modulus and Poisson's ratio.
    pub material: (f64, f64),
    /// Fallback stress (production: `InCell`).
    pub fallback: Fallback,
    /// Exact field for computing the true error (verification only;
    /// `None` in production).
    pub exact: Option<ExactStress<'a>>,
    /// A thermal strain: expansion α and nodal temperature rises per
    /// entry; the stress subtracts it (`element::thermal_strain`).
    pub thermal: Option<(f64, &'a [[f64; 8]])>,
    /// An infill: each entry's core and the core's Lamé parameters; a
    /// cell mostly core recovers with the core's.
    pub core: Option<(&'a [crate::element::operators::Core], (f64, f64))>,
    /// The size of the cube the core's stress is averaged over before
    /// its peak is taken (the gyroid's unit cell), mm.
    pub core_window: f64,
}

fn internal(detail: String) -> AnalysisError {
    AnalysisError::Internal {
        stage: "recover",
        detail,
    }
}

/// Cells (or nodes) below which recovery stays on one thread.
const PARALLEL_MIN_CELLS: usize = 256;

/// Recovers stresses at the mesh vertices of one solid.
pub fn recover_solid(input: &RecoveryInput<'_>) -> Result<SolidStresses, AnalysisError> {
    let grid = input.grid;
    let h = grid.cell_size();
    let (e_mod, nu) = input.material;
    let (lambda, mu) = crate::element::vector::lame(e_mod, nu);
    let e = input.extension;
    let vt = VectorTables::new();
    let operators = crate::element::operators::CellOperators::new(&vt, h, lambda, mu)
        .and_then(|o| match input.core {
            Some((_, lame)) => o.with_core(lame),
            None => Ok(o),
        })
        .map_err(|e| internal(e.to_string()))?;
    let no_core = crate::element::operators::Core::None;
    let dt = DerivativeTables::new();
    let cells = &input.moments.cells;
    // 1. Cell fields and centroid samples, on threads, in entry order.
    let (moments, u_all, core_input, thermal) = (input.moments, input.u, input.core, input.thermal);
    let cell_field = |entry: usize| -> Result<(StressPoly, Sample, [usize; 8]), AnalysisError> {
        let c = &cells[entry];
        let ijk = grid
            .cell_ijk(c.cell)
            .ok_or_else(|| internal("cell outside grid".into()))?;
        let nodes = moments
            .node_keys(grid, entry)
            .ok_or_else(|| internal("cell outside grid".into()))?;
        let mut pos = [0usize; 8];
        for (slot, n) in pos.iter_mut().zip(nodes) {
            *slot = e
                .position(n)
                .ok_or_else(|| internal(format!("node {n} has no value")))?;
        }
        let u: [f64; 24] = std::array::from_fn(|r| u_all[3 * pos[r / 3] + r % 3]);
        let core = core_input
            .and_then(|(cores, _)| cores.get(entry))
            .unwrap_or(&no_core);
        let operator = operators
            .of_core(c, moments.is_fragment(entry), core)
            .map_err(|e| internal(e.to_string()))?;
        let (op, _) = operator.get();
        // A cell mostly core reads the core's stress.
        let mostly_core = match core {
            crate::element::operators::Core::None => false,
            crate::element::operators::Core::Full => true,
            crate::element::operators::Core::Partial(m) => m[0] > 0.5 * c.m[0],
        };
        let (cl, cm) = match core_input {
            Some((_, lame)) if mostly_core => lame,
            _ => (lambda, mu),
        };
        let mut poly = stress_poly(&dt, op, &c.m, &u, h, cl, cm);
        if let Some((alpha, theta)) = thermal
            && let Some(theta) = theta.get(entry)
        {
            let parts = core.parts(&c.m, (lambda, mu), core_input.map(|(_, l)| l));
            let parts: Vec<(&[f64; 27], f64, f64)> =
                parts.iter().map(|(m, l, u)| (m, *l, *u)).collect();
            cell::add_thermal(&mut poly, op, &c.m, &parts, h, (cl, cm), alpha, theta);
        }
        let xc = centroid(&c.m);
        let lo = grid.cell_min(ijk);
        let sample = Sample {
            x: [0, 1, 2].map(|a| lo[a] + 0.5 * h * (xc[a] + 1.0)),
            w: c.m[0].max(0.0),
            s: stress_at(&poly, xc),
        };
        Ok((poly, sample, pos))
    };
    let mut polys: Vec<StressPoly> = Vec::with_capacity(cells.len());
    let mut samples: Vec<Sample> = Vec::with_capacity(cells.len());
    let mut cell_nodes: Vec<[usize; 8]> = Vec::with_capacity(cells.len());
    for field in crate::parallel::map_range(0..cells.len(), PARALLEL_MIN_CELLS, cell_field) {
        let (poly, sample, pos) = field?;
        polys.push(poly);
        samples.push(sample);
        cell_nodes.push(pos);
    }
    // Raw peak over cells with a real fragment (η ≥ 1%).
    let raw_peak_principal = samples
        .iter()
        .zip(cells)
        .filter(|(_, c)| c.volume_fraction(h) >= 0.01)
        .map(|(s, _)| principal(&s.s)[0])
        .fold(f64::NEG_INFINITY, f64::max);
    let true_error = input
        .exact
        .map(|exact| true_error::true_error(input, &polys, exact))
        .transpose()?
        .flatten();
    // 2. Nodal σ*.
    let nodal = nodal_values(input.grid, input.moments, input.extension, &samples);
    // 3 and 4: ZZ per cell where all its nodes have σ*.
    let shapes: [Poly27; 8] = std::array::from_fn(shape);
    let mut cell_error: Vec<Option<f64>> = vec![None; cells.len()];
    let (mut err_sum, mut norm_sum) = (0.0, 0.0);
    // Each cell's energies on threads; summed in entry order.
    let energies = crate::parallel::map_range(0..cells.len(), PARALLEL_MIN_CELLS, |k| {
        let star = star_poly(&shapes, &cell_nodes[k], &nodal)?;
        let delta: StressPoly =
            std::array::from_fn(|q| std::array::from_fn(|i| star[q][i] - polys[k][q][i]));
        let m = &cells[k].m;
        Some((energy(&delta, m, e_mod, nu), energy(&star, m, e_mod, nu)))
    });
    for (k, terms) in energies.into_iter().enumerate() {
        let Some((e2, n2)) = terms else {
            continue;
        };
        err_sum += e2;
        norm_sum += n2;
        if n2 > 0.0 {
            cell_error[k] = Some((e2.max(0.0) / n2).sqrt());
        }
    }
    let thin = wall::thin_walls(grid, input.clipped, input.mesh, input.unreliable);
    let flat = normal::FlatCells::new(grid, input.clipped);
    let inside = |q: [f64; 3]| -> Option<(usize, [f64; 3])> {
        let ijk = grid.cell_of_point(q)?;
        if !crate::moments::lattice::contains(grid, input.clipped, ijk, q) {
            return None;
        }
        locate(grid, input.moments, q, None)
    };
    let sampler = Sampler {
        shapes: &shapes,
        cell_nodes: &cell_nodes,
        nodal: &nodal,
        polys: &polys,
        samples: &samples,
        cell_error: &cell_error,
        fallback: input.fallback,
        normal: Some(sample::Normal {
            flat: &flat,
            inside: &inside,
            h,
        }),
    };
    // Vertices.
    let mut vertices = Vec::with_capacity(input.mesh.vertices.len());
    let mut fallbacks = 0;
    let inward = vertex_inward(input.mesh);
    for (vi, v) in input.mesh.vertices.iter().enumerate() {
        let n = inward.get(vi).copied();
        let (k, xi) = locate(grid, input.moments, *v, n)
            .ok_or_else(|| internal(format!("vertex {vi} in no active cell")))?;
        let cell = cells.get(k).map_or(usize::MAX, |c| c.cell);
        let stress =
            sampler.at_surface((k, xi), (*v, cell, n.unwrap_or_default()), thin.vertex[vi]);
        if stress.recovery == Recovery::InCell && input.fallback != Fallback::ForceConstant {
            fallbacks += 1;
        }
        vertices.push(stress);
    }
    let in_cell = |piece: usize, cell: usize, p: [f64; 3]| -> Option<(usize, [f64; 3])> {
        let k = input.moments.entry_of_piece(piece, cell)?;
        let lo = grid.cell_min(grid.cell_ijk(cell)?);
        Some((
            k,
            [0, 1, 2].map(|a| (2.0 * (p[a] - lo[a]) / h - 1.0).clamp(-1.0, 1.0)),
        ))
    };
    let surface = surface_samples(&sampler, input.clipped, &in_cell, &thin.triangle);
    let interior = input.core.map_or_else(Vec::new, |(cores, _)| {
        interior::interior_wall_samples(grid, input.moments, cores, &polys, &cell_error)
    });
    let core_peak = input.core.and_then(|(cores, _)| {
        core::core_peak(
            grid,
            input.moments,
            cores,
            &polys,
            &cell_error,
            input.core_window,
        )
    });
    Ok(SolidStresses {
        vertices,
        surface,
        global_error: (norm_sum > 0.0).then(|| (err_sum.max(0.0) / norm_sum).sqrt()),
        fallbacks,
        thin_wall: thin.thinnest,
        raw_peak_principal,
        true_error,
        core_peak,
        interior,
    })
}

/// Nodal SPR values for every active node of `extension` from one
/// sample per entry of `moments` (`None` where no interior fit is
/// within reach). Shared by stress and heat-flux recovery.
pub(crate) fn nodal_values(
    grid: &Grid,
    moments: &SolidMoments,
    e: &Extension,
    samples: &[Sample],
) -> Vec<Option<Voigt>> {
    let h = grid.cell_size();
    let cells = &moments.cells;
    // Each node key's patch: the entries it is a corner of, at least
    // PATCH_FILL full.
    let mut patches: Vec<Vec<usize>> = vec![Vec::new(); e.nodes.len()];
    for (k, c) in cells.iter().enumerate() {
        if c.volume_fraction(h) < PATCH_FILL {
            continue;
        }
        for key in moments.node_keys(grid, k).unwrap_or_default() {
            if let Some(p) = e.position(key) {
                patches[p].push(k);
            }
        }
    }
    let node_ijk: Vec<[usize; 3]> = e
        .nodes
        .iter()
        .map(|&n| grid.node_ijk(n / NODE_COPIES).unwrap_or_default())
        .collect();
    // Interior fits.
    let fits: Vec<Option<LinearFit>> =
        crate::parallel::map_range(0..node_ijk.len(), PARALLEL_MIN_CELLS, |p| {
            let patch = &patches[p];
            if patch.len() < MIN_PATCH {
                return None;
            }
            let patch: Vec<Sample> = patch.iter().map(|&k| samples[k]).collect();
            fit(&patch, grid.node_position(node_ijk[p]), h)
        });
    let mut out = Vec::with_capacity(e.nodes.len());
    for (p, nijk) in node_ijk.iter().enumerate() {
        if let Some(f) = &fits[p] {
            out.push(Some(f.at(grid.node_position(*nijk))));
            continue;
        }
        // Nearest interior node within two steps, in fixed order.
        let mut best: Option<(usize, usize)> = None;
        for dz in -2i64..=2 {
            for dy in -2i64..=2 {
                for dx in -2i64..=2 {
                    let q = [
                        nijk[0] as i64 + dx,
                        nijk[1] as i64 + dy,
                        nijk[2] as i64 + dz,
                    ];
                    if q.iter().any(|&c| c < 0) {
                        continue;
                    }
                    let Some(n) = grid.node_index(q.map(|c| c as usize)) else {
                        continue;
                    };
                    for copy in 0..NODE_COPIES {
                        let Some(pos) = e.position(n * NODE_COPIES + copy) else {
                            continue;
                        };
                        if fits[pos].is_none() {
                            continue;
                        }
                        let d = (dx * dx + dy * dy + dz * dz) as usize;
                        if best.is_none_or(|(bd, bp)| (d, pos) < (bd, bp)) {
                            best = Some((d, pos));
                        }
                    }
                }
            }
        }
        out.push(
            best.and_then(|(_, pos)| fits[pos].as_ref())
                .map(|f| f.at(grid.node_position(*nijk))),
        );
    }
    out
}

/// σ* inside a cell as polynomials from its 8 nodal values.
fn star_poly(
    shapes: &[Poly27; 8],
    nodes: &[usize; 8],
    nodal: &[Option<Voigt>],
) -> Option<StressPoly> {
    let mut out: StressPoly = [[0.0; 27]; 6];
    for (a, &n) in nodes.iter().enumerate() {
        let v = nodal[n]?;
        for q in 0..6 {
            for (x, y) in out[q].iter_mut().zip(&shapes[a]) {
                *x += v[q] * y;
            }
        }
    }
    Some(out)
}

/// ∫ σ : C⁻¹ : σ over the fragment with moments `m`.
fn energy(s: &StressPoly, m: &[f64; 27], youngs: f64, poisson: f64) -> f64 {
    let sq = |p: &Poly27| integrate(&mul(p, p), m);
    let normal = sq(&s[0]) + sq(&s[1]) + sq(&s[2]);
    let shear = sq(&s[3]) + sq(&s[4]) + sq(&s[5]);
    let tr: Poly27 = std::array::from_fn(|i| s[0][i] + s[1][i] + s[2][i]);
    ((1.0 + poisson) * (normal + 2.0 * shear) - poisson * sq(&tr)) / youngs
}
