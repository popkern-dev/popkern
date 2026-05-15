use super::atlas::build_patch_atlas;
use super::readiness::{PatchEdgeRole, validate_brep_readiness};
use super::types::{PatchEdge, PolyhedralSpline, PolyhedralSplineError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PolyhedralSplineWorkEstimate {
    pub patch_count: usize,
    pub patch_edge_count: usize,
    pub edge_pair_tests: usize,
    pub sampled_edge_evaluations: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyhedralSplineDebugArtifact {
    pub summary: String,
    pub obj: String,
}

pub fn estimate_polyhedral_spline_work(
    spline: &PolyhedralSpline,
    samples: usize,
) -> PolyhedralSplineWorkEstimate {
    let patch_count = spline.patches.len();
    let patch_edge_count = patch_count * 4;
    let patch_pairs = patch_count.saturating_sub(1) * patch_count / 2;
    let edge_pair_tests = patch_pairs * 16;
    PolyhedralSplineWorkEstimate {
        patch_count,
        patch_edge_count,
        edge_pair_tests,
        sampled_edge_evaluations: edge_pair_tests * samples.max(2),
    }
}

pub fn build_polyhedral_spline_debug_artifact(
    spline: &PolyhedralSpline,
    samples: usize,
    tolerance: f64,
) -> Result<PolyhedralSplineDebugArtifact, PolyhedralSplineError> {
    let readiness = validate_brep_readiness(spline, samples, tolerance)?;
    let work = estimate_polyhedral_spline_work(spline, samples);
    let atlas = build_patch_atlas(spline, samples, tolerance);
    let collapsed_count = readiness
        .edge_uses
        .iter()
        .filter(|edge_use| matches!(edge_use.role, PatchEdgeRole::CollapsedSingularity { .. }))
        .count();
    let summary = format!(
        "family={:?}\npatches={}\nseams={}\nboundary_edges={}\ncollapsed_singular_edges={}\nedge_pair_tests={}\nsampled_edge_evaluations={}\n",
        spline.family,
        spline.patches.len(),
        atlas.seams.len(),
        atlas.boundary_edges.len(),
        collapsed_count,
        work.edge_pair_tests,
        work.sampled_edge_evaluations
    );
    Ok(PolyhedralSplineDebugArtifact {
        summary,
        obj: build_debug_obj(spline, &atlas.boundary_edges),
    })
}

fn build_debug_obj(spline: &PolyhedralSpline, boundary_edges: &[(usize, PatchEdge)]) -> String {
    let mut out = String::from("# popkern polyhedral_spline debug artifact\n");
    let mut vertex_index = 1usize;
    for (patch_index, patch) in spline.patches.iter().enumerate() {
        out.push_str(&format!("o patch_{patch_index}\n"));
        let corners = [
            patch.evaluate(0.0, 0.0),
            patch.evaluate(1.0, 0.0),
            patch.evaluate(1.0, 1.0),
            patch.evaluate(0.0, 1.0),
        ];
        for point in corners {
            out.push_str(&format!("v {} {} {}\n", point.x, point.y, point.z));
        }
        out.push_str(&format!(
            "l {} {} {} {} {}\n",
            vertex_index,
            vertex_index + 1,
            vertex_index + 2,
            vertex_index + 3,
            vertex_index
        ));
        vertex_index += 4;
    }
    out.push_str("o boundary_edges\n");
    for &(patch_index, edge) in boundary_edges {
        out.push_str(&format!("# patch {patch_index} {edge:?}\n"));
    }
    out
}
