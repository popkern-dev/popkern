use super::atlas::{PolyhedralPatchAtlas, build_patch_atlas};
use super::types::{PatchEdge, PolyhedralPatch, PolyhedralSpline, PolyhedralSplineError};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PatchDomain {
    pub patch_index: usize,
    pub u_min: f64,
    pub u_max: f64,
    pub v_min: f64,
    pub v_max: f64,
}

impl PatchDomain {
    pub fn unit_square(patch_index: usize) -> Self {
        Self {
            patch_index,
            u_min: 0.0,
            u_max: 1.0,
            v_min: 0.0,
            v_max: 1.0,
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        let values = [self.u_min, self.u_max, self.v_min, self.v_max];
        if values.iter().any(|value| !value.is_finite()) {
            return Err("domain bounds must be finite".to_string());
        }
        if !(self.u_min < self.u_max && self.v_min < self.v_max) {
            return Err(format!(
                "expected increasing bounds, got u=[{},{}], v=[{},{}]",
                self.u_min, self.u_max, self.v_min, self.v_max
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PatchEdgeRole {
    Boundary,
    Seam { seam_index: usize },
    CollapsedSingularity { incident_seams: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PatchEdgeUse {
    pub patch_index: usize,
    pub edge: PatchEdge,
    pub role: PatchEdgeRole,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolyhedralSplineBrepReadinessReport {
    pub patch_count: usize,
    pub patch_domains: Vec<PatchDomain>,
    pub edge_uses: Vec<PatchEdgeUse>,
    pub atlas: PolyhedralPatchAtlas,
}

pub fn validate_brep_readiness(
    spline: &PolyhedralSpline,
    samples: usize,
    tolerance: f64,
) -> Result<PolyhedralSplineBrepReadinessReport, PolyhedralSplineError> {
    let atlas = build_patch_atlas(spline, samples, tolerance);
    let patch_domains = patch_domains(spline)?;
    let edge_uses = classify_patch_edges(spline, &atlas, tolerance)?;
    Ok(PolyhedralSplineBrepReadinessReport {
        patch_count: spline.patches.len(),
        patch_domains,
        edge_uses,
        atlas,
    })
}

fn patch_domains(spline: &PolyhedralSpline) -> Result<Vec<PatchDomain>, PolyhedralSplineError> {
    let mut domains = Vec::with_capacity(spline.patches.len());
    for patch_index in 0..spline.patches.len() {
        let domain = PatchDomain::unit_square(patch_index);
        domain
            .validate()
            .map_err(|reason| PolyhedralSplineError::InvalidPatchDomain {
                family: spline.family,
                patch_index,
                reason,
            })?;
        domains.push(domain);
    }
    Ok(domains)
}

fn classify_patch_edges(
    spline: &PolyhedralSpline,
    atlas: &PolyhedralPatchAtlas,
    tolerance: f64,
) -> Result<Vec<PatchEdgeUse>, PolyhedralSplineError> {
    if atlas.patch_count != spline.patches.len() {
        return Err(atlas_error(
            spline,
            format!(
                "atlas patch count {} does not match spline patch count {}",
                atlas.patch_count,
                spline.patches.len()
            ),
        ));
    }

    let mut edge_roles = vec![[None; 4]; spline.patches.len()];
    for (seam_index, seam) in atlas.seams.iter().enumerate() {
        assign_edge_role(
            spline,
            &mut edge_roles,
            seam.patch_a,
            seam.edge_a,
            PatchEdgeRole::Seam { seam_index },
            tolerance,
        )?;
        assign_edge_role(
            spline,
            &mut edge_roles,
            seam.patch_b,
            seam.edge_b,
            PatchEdgeRole::Seam { seam_index },
            tolerance,
        )?;
    }
    for &(patch_index, edge) in &atlas.boundary_edges {
        assign_edge_role(
            spline,
            &mut edge_roles,
            patch_index,
            edge,
            PatchEdgeRole::Boundary,
            tolerance,
        )?;
    }

    let mut edge_uses = Vec::with_capacity(spline.patches.len() * 4);
    for (patch_index, roles) in edge_roles.iter().enumerate() {
        for edge in all_edges() {
            let role = roles[edge_index(edge)].ok_or_else(|| {
                atlas_error(
                    spline,
                    format!("patch {patch_index} {edge:?} has no seam or boundary role"),
                )
            })?;
            edge_uses.push(PatchEdgeUse {
                patch_index,
                edge,
                role,
            });
        }
    }
    Ok(edge_uses)
}

fn assign_edge_role(
    spline: &PolyhedralSpline,
    edge_roles: &mut [[Option<PatchEdgeRole>; 4]],
    patch_index: usize,
    edge: PatchEdge,
    role: PatchEdgeRole,
    tolerance: f64,
) -> Result<(), PolyhedralSplineError> {
    let collapsed = spline
        .patches
        .get(patch_index)
        .is_some_and(|patch| edge_span(patch, edge) <= tolerance);
    let Some(roles) = edge_roles.get_mut(patch_index) else {
        return Err(atlas_error(
            spline,
            format!("atlas references missing patch {patch_index}"),
        ));
    };
    let slot = &mut roles[edge_index(edge)];
    if collapsed && matches!(role, PatchEdgeRole::Seam { .. }) {
        *slot = Some(match *slot {
            None => PatchEdgeRole::CollapsedSingularity { incident_seams: 1 },
            Some(PatchEdgeRole::CollapsedSingularity { incident_seams }) => {
                PatchEdgeRole::CollapsedSingularity {
                    incident_seams: incident_seams + 1,
                }
            }
            Some(existing) => {
                return Err(atlas_error(
                    spline,
                    format!(
                        "patch {patch_index} {edge:?} collapsed edge has incompatible roles: {existing:?} and {role:?}"
                    ),
                ));
            }
        });
        return Ok(());
    }
    if let Some(existing) = *slot {
        return Err(atlas_error(
            spline,
            format!("patch {patch_index} {edge:?} has multiple roles: {existing:?} and {role:?}"),
        ));
    }
    *slot = Some(role);
    Ok(())
}

fn edge_span(patch: &PolyhedralPatch, edge: PatchEdge) -> f64 {
    eval_edge(patch, edge, 0.0).distance_to(&eval_edge(patch, edge, 1.0))
}

fn eval_edge(patch: &PolyhedralPatch, edge: PatchEdge, t: f64) -> crate::geometry::Point3 {
    match edge {
        PatchEdge::UMin => patch.evaluate(0.0, t),
        PatchEdge::UMax => patch.evaluate(1.0, t),
        PatchEdge::VMin => patch.evaluate(t, 0.0),
        PatchEdge::VMax => patch.evaluate(t, 1.0),
    }
}

fn atlas_error(spline: &PolyhedralSpline, reason: String) -> PolyhedralSplineError {
    PolyhedralSplineError::AtlasTopology {
        family: spline.family,
        reason,
    }
}

fn all_edges() -> [PatchEdge; 4] {
    [
        PatchEdge::UMin,
        PatchEdge::UMax,
        PatchEdge::VMin,
        PatchEdge::VMax,
    ]
}

fn edge_index(edge: PatchEdge) -> usize {
    match edge {
        PatchEdge::UMin => 0,
        PatchEdge::UMax => 1,
        PatchEdge::VMin => 2,
        PatchEdge::VMax => 3,
    }
}
