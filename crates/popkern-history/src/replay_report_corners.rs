//! Per-corner evidence for a served fillet (PK-603): which
//! construction closed each corner, so a document (and
//! popkern live) can show what a corner-cap preference actually got.
//! Derived from the kernel's corner patches and its mitre-fallback
//! warnings; the kernel records stay the single source of truth.

use std::collections::BTreeMap;

use popkern_blend::{
    BlendResult, BlendWarning, ContinuityLevel, CornerMethod, CornerPatchInfo, CornerSubPatchFamily,
};
use popkern_core::topology::VertexHandle;

use super::CornerEvidence;

/// One record per corner the fillet closed, in vertex-handle order
/// (deterministic across runs).
pub(super) fn corner_evidence(result: &BlendResult) -> Vec<CornerEvidence> {
    let mut by_vertex: BTreeMap<VertexHandle, Vec<&CornerPatchInfo>> = BTreeMap::new();
    for patch in &result.corner_patches {
        by_vertex.entry(patch.vertex).or_default().push(patch);
    }
    let mut corners: Vec<(VertexHandle, CornerEvidence)> = by_vertex
        .into_iter()
        .map(|(vertex, patches)| (vertex, served_corner(result, &patches)))
        .collect();
    // Corners that fell back to a mitre carry no patch; the warning
    // names them.
    for (_, warning) in &result.edge_warnings {
        if let BlendWarning::VertexModeFallback {
            vertex: Some(vertex),
            ..
        } = warning
            && !corners.iter().any(|(v, _)| v == vertex)
        {
            corners.push((
                *vertex,
                CornerEvidence {
                    sides: 0,
                    construction: "classic:mitre".to_string(),
                    achieved: ContinuityLevel::G0,
                    patches: 0,
                    centroid: None,
                },
            ));
        }
    }
    corners.sort_by_key(|(vertex, _)| *vertex);
    corners.into_iter().map(|(_, corner)| corner).collect()
}

fn served_corner(result: &BlendResult, patches: &[&CornerPatchInfo]) -> CornerEvidence {
    let fitted = patches.iter().any(|p| {
        p.sub_patch_provenance
            .as_ref()
            .is_some_and(|prov| prov.family == CornerSubPatchFamily::VariationalFit)
    });
    let construction = if fitted {
        "fitted".to_string()
    } else if patches.iter().any(|p| {
        p.method == popkern_blend::CornerMethod::RollingBallSphere && p.cap_sphere.is_some()
    }) {
        "rolling_ball".to_string()
    } else if patches.iter().any(|p| p.setback_cap) {
        "setback".to_string()
    } else if patches.iter().any(|p| p.exact_boundary_cap) {
        "classic:exact_boundary".to_string()
    } else {
        format!("classic:{}", method_spelling(patches[0].method))
    };
    let achieved = patches
        .iter()
        .map(|p| p.achieved_continuity)
        .min_by_key(|level| continuity_rank(*level))
        .unwrap_or(ContinuityLevel::G0);
    let points: Vec<[f64; 3]> = patches
        .iter()
        .flat_map(|p| result.brep.face_outer_loop_positions(p.face))
        .map(|p| [p.x, p.y, p.z])
        .collect();
    let centroid = (!points.is_empty()).then(|| {
        let n = points.len() as f64;
        let sum = points.iter().fold([0.0; 3], |acc, p| {
            [acc[0] + p[0], acc[1] + p[1], acc[2] + p[2]]
        });
        [sum[0] / n, sum[1] / n, sum[2] / n]
    });
    CornerEvidence {
        sides: patches.iter().map(|p| p.n_sides).max().unwrap_or(0),
        construction,
        achieved,
        patches: patches.len(),
        centroid,
    }
}

fn continuity_rank(level: ContinuityLevel) -> u8 {
    match level {
        ContinuityLevel::G0 => 0,
        ContinuityLevel::G1 => 1,
        ContinuityLevel::G2 => 2,
    }
}

fn method_spelling(method: CornerMethod) -> &'static str {
    match method {
        CornerMethod::RollingBallSphere => "rolling_ball_sphere",
        CornerMethod::PolarNurbs => "polar_nurbs",
        CornerMethod::Mitred => "mitre",
        CornerMethod::Coons => "coons",
        CornerMethod::Gregory => "gregory",
        CornerMethod::PolyhedralNetSpline => "polyhedral_net_spline",
        CornerMethod::TangentPlane => "tangent_plane",
        CornerMethod::PlateMin => "plate_min",
    }
}
