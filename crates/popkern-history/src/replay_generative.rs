//! Replay arm for [`Feature::Generative`](crate::feature::Feature::Generative)
//! (PK-956): the body is the design space; the stiffest design at the
//! recipe's volume fraction replaces it, through `popkern-generative`'s
//! cache (an unchanged design is not solved again).

use popkern_core::history::TrackingContext;
use popkern_core::topology::BRep;
use popkern_generative::{
    BuildDirection, GenerativeRequest, LoadGroup, SupportGroup, SupportKind, generate_cached,
};

use crate::error::ReplayError;
use crate::feature_generative::{GenerativeDesc, GenerativeSupportKindDesc};
use crate::parameter::ParameterTable;
use crate::replay::replay_fns::resolve_param;
use crate::replay_report::{FeatureEvidence, GenerativeRecord, KernelWarning, RecordSource};
use crate::replay_side_effects::ExecuteSideEffects;
use crate::resolve::resolve_faces;

fn declined(kind: &str, detail: String) -> ReplayError {
    ReplayError::FeatureFailed {
        label: "generative".into(),
        error: Box::new(crate::failure::FeatureFailure::GenerativeDeclined {
            generative_kind: kind.to_string(),
            detail,
        }),
    }
}

/// The solver's request for `spec` on `brep`.
fn request(
    spec: &GenerativeDesc,
    brep: &BRep,
    table: &ParameterTable,
) -> Result<GenerativeRequest, ReplayError> {
    let p = |id| resolve_param(table, id);
    let supports = spec
        .supports
        .iter()
        .map(|s| {
            Ok(SupportGroup {
                faces: resolve_faces(&s.faces, brep)?,
                kind: match s.kind {
                    GenerativeSupportKindDesc::Fixed => SupportKind::Fixed,
                    GenerativeSupportKindDesc::SurfaceRoller => SupportKind::SurfaceRoller,
                },
            })
        })
        .collect::<Result<Vec<_>, ReplayError>>()?;
    let loads = spec
        .loads
        .iter()
        .map(|l| {
            Ok(LoadGroup {
                faces: resolve_faces(&l.faces, brep)?,
                force_n: p(l.force)?,
                direction: l.direction,
            })
        })
        .collect::<Result<Vec<_>, ReplayError>>()?;
    Ok(GenerativeRequest {
        material: spec.material.clone(),
        supports,
        loads,
        keep: resolve_faces(&spec.keep.concat(), brep)?,
        volume_fraction: p(spec.volume_fraction)?,
        filter: p(spec.filter)?,
        keep_depth: spec.keep_depth.map(p).transpose()?,
        cell: spec.cell.map(p).transpose()?,
        build_direction: spec.build_direction.map(|b| BuildDirection {
            axis: b.axis as usize,
            negative: b.negative,
        }),
        max_iterations: spec.max_iterations as usize,
        robust: spec.robust,
        max_cells: spec.max_cells.map(|c| c as usize),
        mm_per_unit: spec.mm_per_unit,
    })
}

/// Replaces `brep` (the design space) with its design, every entity
/// named under `tracking`.
pub(crate) fn execute_generative(
    brep: &mut BRep,
    spec: &GenerativeDesc,
    tracking: &mut TrackingContext,
    table: &ParameterTable,
) -> Result<ExecuteSideEffects, ReplayError> {
    let request = request(spec, brep, table)?;
    let (result, cached) = generate_cached(brep, &request).map_err(|e| match e.kind() {
        "Cancelled" => ReplayError::FeatureFailed {
            label: "generative".into(),
            error: Box::new(crate::failure::FeatureFailure::Cancelled),
        },
        kind => declined(kind, e.to_string()),
    })?;
    let mut body = result.brep;
    popkern_generative::name_design(brep, &mut body, &result.sources, tracking);
    let r = &result.report;
    let mut evidence = FeatureEvidence::empty();
    evidence.warnings = r
        .warnings
        .iter()
        .map(|w| KernelWarning {
            source: RecordSource::Generative,
            message: w.clone(),
        })
        .collect();
    evidence.generative = Some(GenerativeRecord {
        compliance_nmm: r.compliance_nmm,
        volume_fraction: r.volume_fraction,
        body_fraction: r.body_fraction,
        iterations: r.iterations,
        converged: r.converged,
        grey_fraction: r.grey_fraction,
        member_thickness_mm: r.member_thickness_mm,
        unsupported_fraction: r.unsupported_fraction,
        cell_mm: r.cell_mm,
        faces: body.faces.len(),
        cached,
        design: result.design.clone(),
    });
    *brep = body;
    Ok(ExecuteSideEffects::reported(evidence))
}
