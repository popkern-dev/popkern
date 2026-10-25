//! A generative run: the stiffest part at a volume fraction
//! (`[generative]`). The VTU is the design (one hexahedron per entry
//! with its density); the PNG field is the densest cell at each mesh
//! vertex.

use popkern_analysis::model::AnalysisModel;
use popkern_analysis::setup::StructuralSetup;
use popkern_analysis::structural::from_setup::solve_generative_setup;
use popkern_analysis::structural::generative::Design;
use popkern_analysis::structural::generative::surface::{SurfaceMesh, design_surface};
use popkern_analysis_brep::Resolution;
use serde_json::{Value, json};

use popkern_analysis::assess::{AssessOptions, assess};
use popkern_analysis::report::ResolvedFace;
use popkern_analysis::structural::from_setup::solve_structural;
use popkern_analysis::structural::generative::body::design_model;

use super::{Failure, analysis_error_kind, to_value, with_resolution};
use crate::vtu_design::write_design_vtu;

/// Optimises, fills `report`, and returns the design VTU and the
/// density at each mesh vertex, per solid.
pub(super) fn run_generative(
    model: &AnalysisModel,
    setup: &StructuralSetup,
    resolution: &Resolution,
    mesh_path: Option<&std::path::Path>,
    report: &mut Value,
) -> Result<(String, Vec<Vec<f64>>), Failure> {
    let (design, material) = solve_generative_setup(model, setup, &resolution.faces)?;
    let solve = with_resolution(design.report.clone(), &resolution.warnings);
    report["solve"] = to_value(&solve)?;
    let total: f64 = design.entries.iter().map(|e| e.volume).sum();
    let kept: f64 = design
        .entries
        .iter()
        .zip(&design.floors)
        .map(|(e, f)| e.volume * f)
        .sum();
    report["result"] = json!({
        "material": material.name,
        "compliance_nmm": design.compliance,
        "volume_fraction": design.volume_fraction,
        "kept_fraction": if total > 0.0 { kept / total } else { 0.0 },
        "iterations": design.history.len(),
        "converged": design.converged,
        "grey_fraction": design.grey_fraction,
        "unsupported_fraction": design.unsupported_fraction,
        "member_thickness_mm": design.member_thickness,
        "eroded_compliance_nmm": design.eroded_compliance,
        "cell_size": design.cell_size,
        "entries": design.entries.len(),
        "history": to_value(&design.history)?,
    });
    let surface = design_surface(&design, 0.5);
    report["result"]["design_mesh"] = json!({
        "triangles": surface.triangles.len(),
        "volume_mm3": surface.volume(),
        "closed": surface.is_closed(),
    });
    if let Some(path) = mesh_path {
        std::fs::write(path, binary_stl(&surface)).map_err(|e| Failure::new("DesignMesh", e))?;
        report["result"]["design_mesh"]["path"] = json!(path.display().to_string());
    }
    report["result"]["design_check"] = design_check(model, setup, resolution, &design, &surface);
    let field = vertex_densities(model, &design);
    Ok((write_design_vtu(&design), field))
}

/// The density at each mesh vertex: the densest design cell touching
/// it (fragments of a split cell included).
fn vertex_densities(model: &AnalysisModel, design: &Design) -> Vec<Vec<f64>> {
    let h = design.cell_size;
    let Some(origin) = design.entries.first().map(|e| e.corner) else {
        return model
            .solids
            .iter()
            .map(|s| vec![0.0; s.vertices.len()])
            .collect();
    };
    let key = |p: [f64; 3], round: bool| {
        [0, 1, 2].map(|a| {
            let t = (p[a] - origin[a]) / h;
            (if round { t.round() } else { t.floor() }) as i64
        })
    };
    let mut cells: Vec<([i64; 3], f64)> = design
        .entries
        .iter()
        .zip(&design.densities)
        .map(|(e, &d)| (key(e.corner, true), d))
        .collect();
    cells.sort_by(|a, b| a.0.cmp(&b.0).then(b.1.total_cmp(&a.1)));
    cells.dedup_by_key(|c| c.0);
    model
        .solids
        .iter()
        .map(|s| {
            s.vertices
                .iter()
                .map(|v| {
                    // A vertex on a cell boundary (every face of an
                    // axis-aligned part) borders up to eight cells:
                    // the densest of those in the design.
                    let k = key(*v, false);
                    (0..8)
                        .filter_map(|c: usize| {
                            let q = [0, 1, 2].map(|a| k[a] - ((c >> a) & 1) as i64);
                            cells.binary_search_by(|x| x.0.cmp(&q)).ok()
                        })
                        .map(|i| cells[i].1)
                        .fold(0.0, f64::max)
                })
                .collect()
        })
        .collect()
}

/// Binary STL of `mesh` (header, count, then per triangle its normal,
/// three vertices and a zero attribute; little-endian f32).
fn binary_stl(mesh: &SurfaceMesh) -> Vec<u8> {
    let mut out = Vec::with_capacity(84 + 50 * mesh.triangles.len());
    let mut header = [0u8; 80];
    let tag = b"popkern-analysis generative design surface";
    header[..tag.len()].copy_from_slice(tag);
    out.extend_from_slice(&header);
    out.extend_from_slice(&(mesh.triangles.len() as u32).to_le_bytes());
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.vertices[i as usize]);
        let u = [0, 1, 2].map(|d| b[d] - a[d]);
        let w = [0, 1, 2].map(|d| c[d] - a[d]);
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let n = if len > 0.0 {
            n.map(|x| x / len)
        } else {
            [0.0; 3]
        };
        for v in [n, a, b, c] {
            for x in v {
                out.extend_from_slice(&(x as f32).to_le_bytes());
            }
        }
        out.extend_from_slice(&0u16.to_le_bytes());
    }
    out
}

/// The strength preview on the design's surface as a part: its
/// triangles on the original faces keep them, so the setup's supports
/// and loads apply as named. A decline is reported, not fatal: the
/// design stands without it.
fn design_check(
    model: &AnalysisModel,
    setup: &StructuralSetup,
    resolution: &Resolution,
    design: &Design,
    surface: &SurfaceMesh,
) -> Value {
    if surface.triangles.is_empty() {
        return json!({
            "status": "declined",
            "kind": "EmptyDesignSurface",
            "message": "no part of the design reaches density 0.5; refine the cells or raise the volume fraction",
        });
    }
    let body = design_model(model, surface, design.cell_size);
    let resolved: Vec<ResolvedFace> = resolution
        .faces
        .iter()
        .map(|r| ResolvedFace {
            name: r.name.clone(),
            faces: r
                .faces
                .iter()
                .copied()
                .filter(|f| body.faces.get(f.0 as usize).is_some_and(|i| i.area > 0.0))
                .collect(),
        })
        .collect();
    let mut plain = setup.clone();
    plain.generative = None;
    let solved = solve_structural(&body, &plain, &resolved).and_then(|(result, material)| {
        let options = AssessOptions {
            build_direction: setup.build_direction,
            ..AssessOptions::default()
        };
        assess(&body, &result, material, &options)
            .map(|a| (result, a))
            .map_err(|e| popkern_analysis::AnalysisError::Internal {
                stage: "assess",
                detail: e.to_string(),
            })
    });
    match solved {
        Ok((result, assessment)) => json!({
            "status": "solved",
            "max_displacement_mm": result.max_displacement,
            "safety": assessment.safety,
        }),
        Err(e) => json!({
            "status": "declined",
            "kind": analysis_error_kind(&e),
            "message": e.to_string(),
        }),
    }
}
