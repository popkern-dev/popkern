//! The steady system of a conduction problem, split so the steady and
//! transient solvers share it: [`prepare`] (grid, moments, aggregation,
//! assembly), [`hierarchy`] (multigrid) and [`finish`] (heat balance,
//! vertex values and the field from a solution).

use crate::aggregate::{Aggregation, aggregate};
use crate::clip::{ClippedSolid, clip_solid};
use crate::element::ScalarTables;
use crate::error::AnalysisError;
use crate::grid::Grid;
use crate::linalg::sparse::Csr;
use crate::model::AnalysisModel;
use crate::moments::{SolidMoments, rasterise_solid, volume_parity};
use crate::options::AnalysisOptions;
use crate::report::{SolveReport, VolumeParity, Warning};
use crate::solve::{Hierarchy, HierarchyOptions};

use super::assemble::{
    BoundaryCache, MmCondition, Physics, SolidInput, assemble_solid, boundary_cache, boundary_heat,
};
use super::{
    ConductionProblem, FaceCondition, HeatBalance, TemperatureField, ThermalResult,
    VOLUME_PARITY_TOL, bounds, grid_failure, internal, sample, to_mm,
};

/// The assembled steady system and what produced it.
pub(super) struct Prepared {
    /// The grid.
    pub grid: Grid,
    /// Per solid: clipped pieces, moments, aggregation.
    pub solids: Vec<(ClippedSolid, SolidMoments, Aggregation)>,
    /// Per solid physics.
    pub physics: Vec<Physics>,
    /// Face conditions, mm units.
    pub conditions: Vec<MmCondition>,
    /// Scalar tables.
    pub tables: ScalarTables,
    /// The report so far.
    pub report: SolveReport,
    /// First free dof of each solid.
    pub offsets: Vec<usize>,
    /// Lattice position of each free dof.
    pub positions: Vec<[usize; 3]>,
    /// Reduced operator.
    pub a: Csr,
    /// Reduced load.
    pub b: Vec<f64>,
    /// Per solid: boundary data kept across passes.
    pub caches: Vec<BoundaryCache>,
}

/// Grid, moments, aggregation and the reduced steady system;
/// `Radiating` faces take their coefficients from `previous`.
pub(super) fn prepare(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    options: &AnalysisOptions,
    prescribed: Option<&TemperatureField>,
    previous: Option<&TemperatureField>,
) -> Result<Prepared, AnalysisError> {
    model.validate()?;
    options.validate().map_err(crate::setup::SetupError::from)?;
    let (lo, hi) = bounds(model);
    let grid = Grid::for_model(lo, hi, options).map_err(|e| grid_failure(e, options))?;
    let h = grid.cell_size();
    let k = 1e-3 * problem.conductivity;
    let total_volume: f64 = model.solids.iter().map(|s| s.enclosed_volume()).sum();
    // Source density per solid, W/mm³: the part-wide share plus the
    // solid's own sources.
    let physics: Vec<Physics> = model
        .solids
        .iter()
        .map(|solid| {
            let own: f64 = problem
                .body_sources
                .iter()
                .filter(|(id, _)| *id == solid.id)
                .map(|(_, w)| w)
                .sum();
            Physics {
                k,
                beta: options.nitsche_factor * k / h,
                source: problem.source_w / total_volume + own / solid.enclosed_volume(),
            }
        })
        .collect();
    let conditions: Vec<MmCondition> = problem.faces.iter().map(|c| to_mm(*c)).collect();
    let tables = ScalarTables::new();
    let mut report = SolveReport::new(*options);
    report.cell_size = Some(h);
    let mut solids = Vec::new();
    for solid in &model.solids {
        let clipped = clip_solid(&grid, solid).map_err(|e| internal("clip", e))?;
        let moments = rasterise_solid(&grid, solid.id, &clipped, options.subcell_depth)?;
        volume_parity(&moments, solid.enclosed_volume(), VOLUME_PARITY_TOL)?;
        report.volume_parity.push(VolumeParity {
            solid: solid.id.0,
            rasterised: moments.volume,
            exact: solid.exact_volume,
        });
        report
            .warnings
            .extend(crate::model::tessellation::volume_warning(solid));
        let area = clipped.pieces.iter().map(|p| p.area).sum();
        let aggregation = aggregate(&grid, &moments, options.aggregation_threshold, area).map_err(
            |e| match e {
                crate::aggregate::AggregationError::Empty(solid) => {
                    AnalysisError::DegenerateGeometry {
                        solid: Some(solid),
                        detail: crate::error::Degeneracy::ZeroVolume,
                    }
                }
                other => internal("aggregate", other),
            },
        )?;
        report.cells_active += moments.cells.len();
        report.cells_cut += moments.cut;
        report.moment_fallback_cells += moments.fallback_cells;
        report.split_cells += moments.fragments.split.len();
        report.split_cells_lattice_moments += moments.fragments.lattice_moment_cells;
        report.cells_aggregated += aggregation.cells_aggregated;
        report.warnings.extend(aggregation.warnings.iter().cloned());
        solids.push((clipped, moments, aggregation));
    }
    let mut prepared = Prepared {
        grid,
        solids,
        physics,
        conditions,
        tables,
        report,
        offsets: Vec::new(),
        positions: Vec::new(),
        a: Csr::from_triplets(0, 0, &[]).map_err(|e| internal("assemble", e))?,
        b: Vec::new(),
        caches: Vec::new(),
    };
    // Facet terms are kept only for a radiating solve's repeated passes.
    let radiating = problem
        .faces
        .iter()
        .any(|f| matches!(f, FaceCondition::Radiating { .. }));
    prepared.caches = prepared
        .solids
        .iter()
        .enumerate()
        .map(|(s, (clipped, moments, aggregation))| {
            let input = SolidInput {
                grid: &prepared.grid,
                clipped,
                moments,
                aggregation,
                conditions: &prepared.conditions,
                prescribed,
                previous: None,
                cache: None,
            };
            boundary_cache(&prepared.tables, &input, &prepared.physics[s], radiating)
        })
        .collect();
    assemble(&mut prepared, prescribed, previous)?;
    Ok(prepared)
}

/// (Re)assembles `prepared`'s reduced system: `Radiating` faces take
/// their coefficients from `previous`.
pub(super) fn assemble(
    prepared: &mut Prepared,
    prescribed: Option<&TemperatureField>,
    previous: Option<&TemperatureField>,
) -> Result<(), AnalysisError> {
    // Block-diagonal global system.
    let mut triplets = Vec::new();
    let mut b = Vec::new();
    let mut positions = Vec::new();
    let mut offsets = Vec::new();
    for (s, (clipped, moments, aggregation)) in prepared.solids.iter().enumerate() {
        let input = SolidInput {
            grid: &prepared.grid,
            clipped,
            moments,
            aggregation,
            conditions: &prepared.conditions,
            prescribed,
            previous,
            cache: prepared.caches.get(s),
        };
        let system = assemble_solid(&prepared.tables, &input, &prepared.physics[s])
            .map_err(|e| internal("assemble", e))?;
        let offset = b.len();
        offsets.push(offset);
        for r in 0..system.a.rows() {
            for (c, v) in system.a.row(r) {
                triplets.push((offset + r, offset + c, v));
            }
        }
        b.extend(system.b);
        positions.extend(aggregation.extension.free_positions(&prepared.grid));
    }
    let n = b.len();
    prepared.report.dof = n;
    prepared.a = Csr::from_triplets(n, n, &triplets).map_err(|e| internal("assemble", e))?;
    prepared.b = b;
    prepared.positions = positions;
    prepared.offsets = offsets;
    Ok(())
}

/// The temperature field of free-dof solution `x` of `prepared`.
pub(super) fn field(prepared: &Prepared, x: &[f64]) -> TemperatureField {
    TemperatureField {
        grid: prepared.grid,
        solids: prepared
            .solids
            .iter()
            .enumerate()
            .map(|(s, (_, moments, aggregation))| {
                let e = &aggregation.extension;
                let start = prepared.offsets.get(s).copied().unwrap_or(0);
                let temps = e.expand(&x[start..start + e.free]);
                (moments.clone(), e.clone(), temps)
            })
            .collect(),
    }
}

/// The multigrid hierarchy of operator `a` (free-dof `positions`),
/// recording its levels in `report`.
pub(super) fn hierarchy(
    a: Csr,
    positions: Vec<[usize; 3]>,
    options: &AnalysisOptions,
    report: &mut SolveReport,
) -> Result<Hierarchy, AnalysisError> {
    let n = a.rows();
    let hierarchy_options = HierarchyOptions {
        coarsest_dofs: options.coarsest_dofs,
        split_components: true,
    };
    let hierarchy =
        Hierarchy::build(a, positions, vec![0; n], &hierarchy_options).map_err(|e| match e {
            crate::solve::MultigridError::CoarsestTooLarge { dofs, cap } => {
                AnalysisError::CoarseLevelTooLarge { dofs, cap }
            }
            crate::solve::MultigridError::Coarsest(c) => report.singular(c.to_string()),
            other => internal("multigrid", other),
        })?;
    report.levels = hierarchy.levels.len();
    report.level_dofs = hierarchy.levels.iter().map(|l| l.a.rows()).collect();
    Ok(hierarchy)
}

/// Of free-dof solution `x` of `prepared` (no prescribed faces;
/// radiating ones linearised about `previous`): (net heat out through the faces, heat in by flux and
/// sources, W; highest vertex temperature, °C).
pub(super) fn state(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    prepared: &Prepared,
    x: &[f64],
    previous: Option<&TemperatureField>,
) -> Result<(f64, f64, f64), AnalysisError> {
    let mut face_out = vec![0.0; model.faces.len()];
    let mut heat_in = problem.source_w + problem.body_sources.iter().map(|(_, w)| w).sum::<f64>();
    let mut max_c = f64::NEG_INFINITY;
    for (s, (clipped, moments, aggregation)) in prepared.solids.iter().enumerate() {
        let e = &aggregation.extension;
        let temps = e.expand(&x[prepared.offsets[s]..prepared.offsets[s] + e.free]);
        let input = SolidInput {
            grid: &prepared.grid,
            clipped,
            moments,
            aggregation,
            conditions: &prepared.conditions,
            prescribed: None,
            previous,
            cache: prepared.caches.get(s),
        };
        boundary_heat(
            &prepared.tables,
            &input,
            &prepared.physics[s],
            &temps,
            &mut face_out,
        );
        heat_in += clipped
            .pieces
            .iter()
            .map(|p| match prepared.conditions.get(p.face.0 as usize) {
                Some(MmCondition::Flux { q }) => q * p.area,
                _ => 0.0,
            })
            .sum::<f64>();
        let values = sample::vertex_values(&prepared.grid, &model.solids[s], moments, e, &temps)?;
        max_c = values.iter().map(|v| v.0).fold(max_c, f64::max);
    }
    Ok((face_out.iter().sum(), heat_in, max_c))
}

/// The pieces of a prepared system the heat balance reads.
type BalanceParts<'a> = (
    &'a Grid,
    &'a [(ClippedSolid, SolidMoments, Aggregation)],
    &'a [Physics],
    &'a [MmCondition],
    &'a ScalarTables,
    &'a [BoundaryCache],
);

/// The heat flows of free-dof solution `x` (radiating faces linearised
/// about `previous`), as the result reports them.
pub(super) fn heat_balance(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    (grid, solids, physics, conditions, tables, caches): BalanceParts<'_>,
    offsets: &[usize],
    x: &[f64],
    (prescribed, previous): (Option<&TemperatureField>, Option<&TemperatureField>),
) -> HeatBalance {
    let mut heat = HeatBalance {
        flux_in: 0.0,
        source: problem.source_w + problem.body_sources.iter().map(|(_, w)| w).sum::<f64>(),
        face_out: vec![0.0; model.faces.len()],
    };
    for (s, (clipped, moments, aggregation)) in solids.iter().enumerate() {
        let e = &aggregation.extension;
        let temps = e.expand(&x[offsets[s]..offsets[s] + e.free]);
        let input = SolidInput {
            grid,
            clipped,
            moments,
            aggregation,
            conditions,
            prescribed,
            previous,
            cache: caches.get(s),
        };
        boundary_heat(tables, &input, &physics[s], &temps, &mut heat.face_out);
        heat.flux_in += clipped
            .pieces
            .iter()
            .map(|p| match conditions.get(p.face.0 as usize) {
                Some(MmCondition::Flux { q }) => q * p.area,
                _ => 0.0,
            })
            .sum::<f64>();
    }
    heat
}

/// [`heat_balance`] of `prepared`.
pub(super) fn prepared_balance(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    prepared: &Prepared,
    x: &[f64],
    fields: (Option<&TemperatureField>, Option<&TemperatureField>),
) -> HeatBalance {
    heat_balance(
        model,
        problem,
        (
            &prepared.grid,
            &prepared.solids,
            &prepared.physics,
            &prepared.conditions,
            &prepared.tables,
            &prepared.caches,
        ),
        &prepared.offsets,
        x,
        fields,
    )
}

/// The result of free-dof solution `x` of `prepared`.
pub(super) fn finish(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    prepared: Prepared,
    x: &[f64],
    (prescribed, previous): (Option<&TemperatureField>, Option<&TemperatureField>),
) -> Result<ThermalResult, AnalysisError> {
    let Prepared {
        grid,
        solids,
        physics,
        conditions,
        tables,
        mut report,
        offsets,
        caches,
        ..
    } = prepared;
    let heat = heat_balance(
        model,
        problem,
        (&grid, &solids, &physics, &conditions, &tables, &caches),
        &offsets,
        x,
        (prescribed, previous),
    );
    let mut vertex_temperature = Vec::new();
    let mut vertex_heat_flux = Vec::new();
    let mut field_solids = Vec::new();
    for (s, (_, moments, aggregation)) in solids.iter().enumerate() {
        let e = &aggregation.extension;
        let free = &x[offsets[s]..offsets[s] + e.free];
        let temps = e.expand(free);
        let mesh = &model.solids[s];
        let values = sample::vertex_values(&grid, mesh, moments, e, &temps)?;
        let gradients = sample::recovered_gradients(&grid, mesh, moments, e, &temps)?;
        // q = −k∇T: k in W/(m·K), ∇T in K/mm, so ×1000 for W/m².
        vertex_heat_flux.push(
            gradients
                .iter()
                .map(|g| g.map(|c| -1e3 * problem.conductivity * c))
                .collect::<Vec<_>>(),
        );
        vertex_temperature.push(values.into_iter().map(|(t, _)| t).collect::<Vec<_>>());
        field_solids.push((moments.clone(), e.clone(), temps));
    }
    let max_temperature_c = vertex_temperature
        .iter()
        .flatten()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let max_heat_flux_w_per_m2 = vertex_heat_flux
        .iter()
        .flatten()
        .map(|q: &[f64; 3]| (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt())
        .fold(0.0, f64::max);
    if let Some(hdt) = problem.hdt_c
        && max_temperature_c > hdt - 10.0
    {
        report.warnings.push(Warning::ServiceTemperature {
            max_c: max_temperature_c,
            hdt_c: hdt,
        });
    }
    Ok(ThermalResult {
        vertex_temperature,
        max_temperature_c,
        vertex_heat_flux,
        max_heat_flux_w_per_m2,
        heat,
        field: TemperatureField {
            grid,
            solids: field_solids,
        },
        report,
        transient: Vec::new(),
    })
}
