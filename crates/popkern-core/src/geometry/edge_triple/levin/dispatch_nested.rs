use super::dispatch::{
    classify_branch_kind, derive_linear_affine_b, intersect_plane_with_target, provenance_kind,
};
use crate::geometry::Point3;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::EdgeAnchors;
use crate::geometry::edge_triple::types::{
    BranchSelector, EdgeTriple, EdgeTripleProvenance, EdgeTripleQuality, ParameterMap,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// Nested-radical fallback for Levin's algorithm. Builds the rank-2
/// fibre symbolically in `Q`, extracts plane normals through the
/// nested surd field, then converts the selected plane to f64 at the
/// boundary for the existing plane/surface dispatchers.
pub(crate) fn try_levin_pair_nested(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let a_rat = super::rational::surface_to_quadric_rat(surf_a)?;
    let b_rat = super::rational::surface_to_quadric_rat(surf_b)?;
    const MAX_INPUT_BITS: u64 = 32;
    if !super::qsqrtd_nested::matrix_entries_within_bit_budget(&a_rat, MAX_INPUT_BITS)
        || !super::qsqrtd_nested::matrix_entries_within_bit_budget(&b_rat, MAX_INPUT_BITS)
    {
        return None;
    }
    let fiber = super::qsqrtd_nested::find_small_rational_fiber(&a_rat, &b_rat)?;
    let planes_qsqrtd = super::qsqrtd_nested::extract_planes_qsqrtd_nested(&fiber)?;
    let p_a = super::qsqrtd_nested::qsqrtd_sqrt_l_array_to_f64(&planes_qsqrtd[0])?;
    let p_b = super::qsqrtd_nested::qsqrtd_sqrt_l_array_to_f64(&planes_qsqrtd[1])?;

    let p_a = normalise_plane(p_a)?;
    let p_b = normalise_plane(p_b)?;
    let planes: Vec<[f64; 4]> = vec![p_a, p_b];

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let anchor_tol = 100.0 * tol.linear();
    let plane = *planes
        .iter()
        .find(|p| plane_distance(p, start) <= anchor_tol && plane_distance(p, end) <= anchor_tol)?;
    let plane_surf = super::decompose::plane_to_surface(plane);

    let triple_a = intersect_plane_with_target(&plane_surf, surf_a, anchors, tol)?.ok()?;
    let triple_b = intersect_plane_with_target(&plane_surf, surf_b, anchors, tol)?.ok()?;

    let total_branches = planes.len() as u8;
    let mut branch_kinds = Vec::with_capacity(planes.len());
    for p in &planes {
        if let Some(kind) = classify_branch_kind(p, surf_a, anchors, tol) {
            branch_kinds.push(kind);
        }
    }
    let cardinality = Some(crate::geometry::edge_triple::types::BranchCardinality {
        total: total_branches,
        kinds: branch_kinds,
        tangent: false,
    });

    let (a_b_coef, b_b_coef) = derive_linear_affine_b(&triple_a, &triple_b)?;
    let triple = EdgeTriple {
        curve_3d: triple_a.curve_3d.clone(),
        pcurve_a: triple_a.pcurve_a.clone(),
        pcurve_b: triple_b.pcurve_a.clone(),
        t_start: triple_a.t_start,
        t_end: triple_a.t_end,
        topology: triple_a.topology,
        parameter_map: ParameterMap::LinearAffine {
            a_a: 0.0,
            b_a: 1.0,
            a_b: a_b_coef,
            b_b: b_b_coef,
        },
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: provenance_kind(surf_a, surf_b),
            branch: BranchSelector::OnlyBranch,
            cardinality,
        },
        quality: EdgeTripleQuality::default(),
    };
    Some(super::super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

fn normalise_plane(p: [f64; 4]) -> Option<[f64; 4]> {
    let n_len2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
    if n_len2 < 1e-30 {
        return None;
    }
    let inv = 1.0 / n_len2.sqrt();
    Some([p[0] * inv, p[1] * inv, p[2] * inv, p[3] * inv])
}

fn plane_distance(p: &[f64; 4], q: Point3) -> f64 {
    (p[0] * q.x + p[1] * q.y + p[2] * q.z + p[3]).abs()
}
