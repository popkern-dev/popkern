//! The static solve: grid, moments, aggregation, assembly, multigrid
//! PCG, then face forces, vertex displacements and recovered stresses.

use crate::aggregate::aggregate;
use crate::clip::clip_solid;
use crate::element::VectorTables;
use crate::element::vector::lame;
use crate::element::vector_boundary::BoundaryTables;
use crate::error::{AnalysisError, Degeneracy};
use crate::grid::Grid;
use crate::linalg::sparse::Csr;
use crate::model::AnalysisModel;
use crate::moments::volume_parity;
use crate::options::AnalysisOptions;
use crate::report::{ForceBalance, SolveReport, VolumeParity};
use crate::thermal::sample::point_value;

use super::assemble::{
    Elastic, FaceMechanics, SolidInput, assemble_solid, face_forces, vector_extension,
};
use super::rigid::{motions_below, rigid_spectrum};
pub(super) use super::system::{solve_system, solve_system_reusing};
use super::{
    DisplacementField, ElasticProblem, POISSON_LIMIT, RIGID_MODE_TOL, StructuralResult, assemble,
    centroid, global_error, internal, support,
};

/// A solve's optional inputs.
#[derive(Clone, Copy, Default)]
pub(super) struct Extras<'a> {
    /// The field `Prescribed` faces are held to (submodels).
    pub prescribed: Option<assemble::Prescribed<'a>>,
    /// A thermal strain.
    pub thermal: Option<&'a super::thermoelastic::ThermalLoad<'a>>,
    /// An infill core.
    pub infill: Option<&'a super::infill::Infill>,
}

/// The assembled structural system and what produced it.
pub(super) struct Prepared {
    pub grid: Grid,
    pub elastic: Elastic,
    pub vt: VectorTables,
    pub bt: BoundaryTables,
    pub report: SolveReport,
    pub solids: Vec<(
        crate::model::SolidId,
        crate::clip::ClippedSolid,
        crate::moments::SolidMoments,
        crate::aggregate::Aggregation,
    )>,
    pub strains: Vec<Option<super::thermoelastic::ThermalStrains>>,
    /// Each solid's infill cores, per entry.
    pub cores: Vec<Option<Vec<crate::element::operators::Core>>>,
    /// The cube the core's stress is averaged over, mm (the gyroid's
    /// cell; 0 without infill).
    pub core_window: f64,
    /// Reduced stiffness per solid.
    pub blocks: Vec<Csr>,
    pub b: Vec<f64>,
    pub positions: Vec<[usize; 3]>,
    pub labels: Vec<u8>,
    pub offsets: Vec<usize>,
}

/// Grid, moments, aggregation and the reduced system of `problem`. A
/// `free` system (modal analysis) may have no support and rigid modes.
pub(super) fn prepare(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    options: &AnalysisOptions,
    (split, free): (crate::moments::Splitting, bool),
    extras: &Extras<'_>,
) -> Result<Prepared, AnalysisError> {
    let Extras {
        prescribed,
        thermal,
        infill,
    } = *extras;

    model.validate()?;
    options.validate().map_err(crate::setup::SetupError::from)?;
    if problem.poisson > POISSON_LIMIT {
        return Err(AnalysisError::NearlyIncompressible {
            poisson: problem.poisson,
            limit: POISSON_LIMIT,
        });
    }
    if !free
        && !problem.faces.iter().any(|f| {
            matches!(
                f,
                FaceMechanics::Support(_)
                    | FaceMechanics::Prescribed
                    | FaceMechanics::PrescribedTraction
            )
        })
    {
        return Err(AnalysisError::NoSupport);
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in model.solids.iter().flat_map(|s| &s.vertices) {
        for a in 0..3 {
            lo[a] = lo[a].min(v[a]);
            hi[a] = hi[a].max(v[a]);
        }
    }
    let grid = Grid::for_model(lo, hi, options).map_err(|e| internal("grid", e))?;
    let h = grid.cell_size();
    let (lambda, mu) = lame(problem.youngs, problem.poisson);
    // β must dominate the consistency term σ(u)·n, which scales with
    // the P-wave modulus λ + 2μ, not E: at ν = 0.49 it is ~17E, and
    // β = γ_D·E/h leaves the operator indefinite (measured, PH-023).
    let core_lame = infill.map(|f| {
        let (e, nu) = super::infill::gyroid(f.density, problem.youngs, problem.poisson);
        lame(e, nu)
    });
    let elastic = Elastic {
        core: core_lame,
        lambda,
        mu,
        beta: options.nitsche_factor * (lambda + 2.0 * mu) / h,
        enhanced_flux: true,
    };
    let (vt, bt) = (VectorTables::new(), BoundaryTables::new());
    let mut report = SolveReport::new(*options);
    report.cell_size = Some(h);
    let mut solids = Vec::new();
    for solid in &model.solids {
        let clipped = clip_solid(&grid, solid).map_err(|e| internal("clip", e))?;
        let moments = crate::moments::rasterise_with(
            &grid,
            solid.id,
            &clipped,
            options.subcell_depth,
            &|_| false,
            split,
        )?;
        volume_parity(
            &moments,
            solid.enclosed_volume(),
            crate::thermal::VOLUME_PARITY_TOL,
        )?;
        report.volume_parity.push(VolumeParity {
            solid: solid.id.0,
            rasterised: moments.volume,
            exact: solid.exact_volume,
        });
        report
            .warnings
            .extend(crate::model::tessellation::volume_warning(solid));
        let area = clipped.pieces.iter().map(|p| p.area).sum();
        let aggregation =
            aggregate(&grid, &moments, options.aggregation_threshold, area).map_err(|_| {
                AnalysisError::DegenerateGeometry {
                    solid: Some(solid.id),
                    detail: Degeneracy::ZeroVolume,
                }
            })?;
        report.cells_active += moments.cells.len();
        report.cells_cut += moments.cut;
        report.moment_fallback_cells += moments.fallback_cells;
        report.split_cells += moments.fragments.split.len();
        report.split_cells_lattice_moments += moments.fragments.lattice_moment_cells;
        report.cells_aggregated += aggregation.cells_aggregated;
        report.warnings.extend(aggregation.warnings.iter().cloned());
        solids.push((solid.id, clipped, moments, aggregation));
    }
    let cores: Vec<Option<Vec<crate::element::operators::Core>>> = solids
        .iter()
        .map(|(_, clipped, moments, _)| {
            infill.map(|f| super::infill::cores(&grid, clipped, moments, f.walls_mm))
        })
        .collect();
    report.infill_density = infill.map(|f| f.density);
    report.infill_core_volume = solids
        .iter()
        .zip(&cores)
        .map(|((_, _, moments, _), c)| super::infill::core_volume(moments, c.as_deref()))
        .sum();
    let strains: Vec<Option<super::thermoelastic::ThermalStrains>> = solids
        .iter()
        .map(|(_, _, moments, _)| {
            thermal
                .map(|t| {
                    super::thermoelastic::ThermalStrains::sample(
                        &grid,
                        moments,
                        t.alpha,
                        t.reference_c,
                        t.temperature,
                    )
                })
                .transpose()
        })
        .collect::<Result<_, _>>()?;
    let mut blocks = Vec::new();
    let mut b = Vec::new();
    let mut positions = Vec::new();
    let mut labels = Vec::new();
    let mut offsets = Vec::new();
    for (s, (id, clipped, moments, aggregation)) in solids.iter().enumerate() {
        let input = SolidInput {
            grid: &grid,
            clipped,
            moments,
            aggregation,
            faces: &problem.faces,
            prescribed,
            thermal: strains[s].as_ref(),
            cores: cores[s].as_deref(),
            scales: None,
        };
        let system =
            assemble_solid(&vt, &bt, &input, &elastic).map_err(|e| internal("assemble", e))?;
        report.eas_disabled_cells += system.eas_disabled;
        let free_xyz: Vec<[f64; 3]> = aggregation
            .extension
            .free_positions(&grid)
            .iter()
            .map(|&ijk| grid.node_position(ijk))
            .collect();
        let spectrum = rigid_spectrum(&system.a, &free_xyz, centroid(&grid, moments))
            .map_err(|e| internal("rigid modes", e))?;
        report.rigid_mode_ratios.push(spectrum.ratios);
        let motions = motions_below(&spectrum, RIGID_MODE_TOL);
        if !free && !motions.is_empty() {
            return Err(AnalysisError::FreeMotion {
                solid: *id,
                motions,
            });
        }
        let min_span = support::SUPPORT_SPAN_CELLS * h;
        if let Some((span, axes)) = support::weak_rotations(clipped, &problem.faces, min_span)
            && !axes.is_empty()
        {
            report
                .warnings
                .push(crate::report::Warning::PointLikeSupport {
                    solid: id.0,
                    span,
                    cell: h,
                    rotations: axes,
                });
        }
        offsets.push(b.len());
        blocks.push(system.a);
        b.extend(system.b);
        for ijk in aggregation.extension.free_positions(&grid) {
            positions.extend([ijk; 3]);
            labels.extend([0u8, 1, 2]);
        }
    }
    Ok(Prepared {
        grid,
        elastic,
        vt,
        bt,
        report,
        solids,
        strains,
        cores,
        core_window: infill.map_or(0.0, |f| f.cell_mm()),
        blocks,
        b,
        positions,
        labels,
        offsets,
    })
}

/// The solve, with `prescribed` the field `Prescribed` faces are held to.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve_core(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    options: &AnalysisOptions,
    fallback: crate::recover::Fallback,
    exact: Option<crate::recover::ExactStress<'_>>,
    split: crate::moments::Splitting,
    extras: &Extras<'_>,
) -> Result<StructuralResult, AnalysisError> {
    let mut prepared = prepare(model, problem, options, (split, false), extras)?;
    let blocks = std::mem::take(&mut prepared.blocks);
    let positions = std::mem::take(&mut prepared.positions);
    let labels = std::mem::take(&mut prepared.labels);
    let b = std::mem::take(&mut prepared.b);
    let report = &mut prepared.report;
    let solution = solve_system(blocks, positions, labels, &b, options, report)?;
    finish(
        model,
        problem,
        prepared,
        &solution.x,
        (fallback, exact),
        extras.prescribed,
    )
}

/// The result of free-dof displacements `x` of `prepared`.
pub(super) fn finish(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    prepared: Prepared,
    x: &[f64],
    (fallback, exact): (
        crate::recover::Fallback,
        Option<crate::recover::ExactStress<'_>>,
    ),
    prescribed: Option<assemble::Prescribed<'_>>,
) -> Result<StructuralResult, AnalysisError> {
    let Prepared {
        grid,
        elastic,
        vt,
        bt,
        mut report,
        solids,
        strains,
        cores,
        core_window,
        offsets,
        ..
    } = prepared;
    let mut face_reaction = vec![[0.0; 3]; model.faces.len()];
    let mut face_applied = vec![[0.0; 3]; model.faces.len()];
    let mut vertex_displacement = Vec::new();
    let mut field_solids = Vec::new();
    let mut stresses = Vec::new();
    for (s, (_, clipped, moments, aggregation)) in solids.iter().enumerate() {
        let e = &aggregation.extension;
        let ev = vector_extension(e).map_err(|e| internal("extension", e))?;
        let mut u = vec![0.0; ev.rows()];
        ev.mul_vec(&x[offsets[s]..offsets[s] + ev.cols()], &mut u)
            .map_err(|e| internal("extension", e))?;
        let input = SolidInput {
            grid: &grid,
            clipped,
            moments,
            aggregation,
            faces: &problem.faces,
            prescribed,
            thermal: strains[s].as_ref(),
            cores: cores[s].as_deref(),
            scales: None,
        };
        face_forces(
            &vt,
            &bt,
            &input,
            &elastic,
            &u,
            &mut face_reaction,
            &mut face_applied,
        )
        .map_err(|e| internal("forces", e))?;
        let components: [Vec<f64>; 3] =
            std::array::from_fn(|c| u.iter().skip(c).step_by(3).copied().collect());
        let mut disp = Vec::with_capacity(model.solids[s].vertices.len());
        for v in &model.solids[s].vertices {
            let mut d = [0.0; 3];
            for (slot, values) in d.iter_mut().zip(&components) {
                *slot = point_value(&grid, moments, e, values, *v)?.ok_or_else(|| {
                    AnalysisError::Internal {
                        stage: "sample",
                        detail: "vertex in no active cell".to_string(),
                    }
                })?;
            }
            disp.push(d);
        }
        vertex_displacement.push(disp);
        let (unreliable, mismatch) =
            crate::model::tessellation::mismatch_warning(&model.solids[s], &model.faces);
        report.warnings.extend(mismatch);
        let recovered = crate::recover::recover_solid(&crate::recover::RecoveryInput {
            grid: &grid,
            moments,
            extension: e,
            u: &u,
            unreliable: &unreliable,
            mesh: &model.solids[s],
            clipped,
            material: (problem.youngs, problem.poisson),
            fallback,
            exact,
            thermal: strains[s].as_ref().map(|t| (t.alpha, t.theta.as_slice())),
            core: cores[s].as_deref().zip(elastic.core),
            core_window,
        })?;
        report.spr_fallbacks += recovered.fallbacks;
        if let Some(wall) = recovered.thin_wall {
            report.warnings.push(crate::report::Warning::UnderResolved {
                solid: model.solids[s].id.0,
                wall,
                cell: grid.cell_size(),
            });
        }
        stresses.push(recovered);
        field_solids.push((moments.clone(), e.clone(), components));
    }
    let sum =
        |v: &[[f64; 3]]| -> [f64; 3] { std::array::from_fn(|i| v.iter().map(|f| f[i]).sum()) };
    let transmitted = sum(&face_reaction);
    report.force_balance = Some(ForceBalance {
        applied: sum(&face_applied),
        reaction: transmitted.map(|x| -x),
    });
    let field = DisplacementField {
        grid,
        solids: field_solids,
    };
    // Over the mesh vertices and the cell-resolution surface samples:
    // display triangles can be far larger than a cell, and the largest
    // displacement need not sit on a vertex (PH-050).
    let norm = |d: [f64; 3]| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let mut max_displacement = vertex_displacement
        .iter()
        .flatten()
        .map(|d| norm(*d))
        .fold(0.0, f64::max);
    for sample in stresses.iter().flat_map(|s| &s.surface) {
        if let Some(d) = field.displacement_at(sample.position)? {
            max_displacement = max_displacement.max(norm(d));
        }
    }
    Ok(StructuralResult {
        vertex_displacement,
        max_displacement,
        face_reaction,
        field,
        global_error: global_error(&stresses),
        stresses,
        report,
    })
}
