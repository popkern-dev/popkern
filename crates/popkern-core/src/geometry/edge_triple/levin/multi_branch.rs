//! Stage 2b.50: multi-branch enumeration for Levin quadric pairs.
//!
//! When the rank-2 fibre of the pencil decomposes into **two** real
//! planes (the canonical Levin case for cyl×cyl skew, cyl×cone
//! non-coaxial, cone×cone non-coaxial), the intersection has two
//! branches — one curve in each plane. The single-branch
//! `try_levin_pair` picks whichever plane contains the user-supplied
//! anchors and discards the other. This module enumerates **all**
//! branches by synthesising anchors at the **branch points** for the
//! plane that doesn't match the user's anchors.
//!
//! See `branch_points.rs` for the geometry of branch-point recovery.

#![allow(dead_code)]

use crate::geometry::Point3;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::EdgeAnchors;
use crate::geometry::edge_triple::types::{
    BranchKind, BranchSelector, EdgeTriple, EdgeTripleProvenance, EdgeTripleQuality, ParameterMap,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

use super::branch_points::branch_points_on_kernel_line;
use super::decompose::{DegenerateShape, classify_degenerate, plane_to_surface};
use super::dispatch::{
    classify_branch_kind, derive_linear_affine_b, intersect_plane_with_target, matrix_max_abs,
    provenance_kind,
};
use super::pencil::find_degenerate_fiber;
use super::quadric::surface_to_quadric;

/// Multi-branch counterpart to `try_levin_pair`. Returns one
/// `EdgeTriple` per real plane fibre — typically 2 for rank-2
/// mixed-sign fibres, 1 for rank-1 (single doubled plane).
///
/// For each plane:
/// - If `primary_anchors` lies on the plane within `100 · tol.linear()`,
///   reuse those anchors (the user's intended branch).
/// - Otherwise, synthesise anchors via `branch_points_on_kernel_line`
///   — the two points where the kernel line of the rank-2 fibre meets
///   `surf_a`. These points lie on both planes by construction, so
///   they're valid anchors for the curve in either plane.
///
/// Returns `None` if Levin doesn't apply (non-quadric inputs, no
/// degenerate fibre, rank ≥ 3 fibre, branch-point synthesis fails).
/// On success the returned `Vec` is non-empty.
pub(crate) fn try_levin_pair_all_branches(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    primary_anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Vec<EdgeTriple>> {
    let a_q = surface_to_quadric(surf_a)?;
    let b_q = surface_to_quadric(surf_b)?;
    let fiber = find_degenerate_fiber(&a_q, &b_q)?;
    let scale = matrix_max_abs(&a_q.0).max(matrix_max_abs(&b_q.0)).max(1.0);
    let shape = classify_degenerate(&fiber.quadric, scale);
    let planes: Vec<[f64; 4]> = match shape {
        DegenerateShape::SinglePlane(p) => vec![p],
        DegenerateShape::TwoPlanes(p1, p2) => vec![p1, p2],
        _ => return None,
    };

    let (user_start, user_end) = match primary_anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let plane_distance =
        |p: &[f64; 4], q: Point3| (p[0] * q.x + p[1] * q.y + p[2] * q.z + p[3]).abs();
    let anchor_tol = 100.0 * tol.linear();

    // Precompute branch_kinds for the cardinality report (shared
    // across all triples in the pair).
    let mut branch_kinds = Vec::with_capacity(planes.len());
    for p in &planes {
        if let Some(kind) = classify_branch_kind(p, surf_a, primary_anchors, tol) {
            branch_kinds.push(kind);
        } else {
            branch_kinds.push(BranchKind::Nurbs);
        }
    }
    let total_branches = planes.len() as u8;

    let mut triples = Vec::with_capacity(planes.len());
    for (i, plane) in planes.iter().enumerate() {
        let on_user_plane = plane_distance(plane, user_start) <= anchor_tol
            && plane_distance(plane, user_end) <= anchor_tol;
        let plane_anchors = if on_user_plane {
            primary_anchors.clone()
        } else if planes.len() >= 2 {
            // Synthesise via branch points. `surf_a` is fine for the
            // synthesis source (the kernel line meets surf_a's
            // intersection set, and any branch point on surf_a is
            // also on surf_b by construction since both quadrics
            // contain both branches).
            let other = &planes[1 - i];
            let (p_start, p_end) = branch_points_on_kernel_line(plane, other, surf_a)?;
            EdgeAnchors::Open {
                start: p_start,
                end: p_end,
            }
        } else {
            // Single-plane case but anchors don't match — caller
            // probably passed bad anchors. Skip this branch.
            continue;
        };

        let plane_surf = plane_to_surface(*plane);
        let triple_a =
            intersect_plane_with_target(&plane_surf, surf_a, &plane_anchors, tol)?.ok()?;
        let triple_b =
            intersect_plane_with_target(&plane_surf, surf_b, &plane_anchors, tol)?.ok()?;
        let stitched = stitch_levin_branch(
            triple_a,
            triple_b,
            surf_a,
            surf_b,
            tol,
            total_branches,
            branch_kinds.clone(),
            i as u32,
        )?
        .ok()?;
        triples.push(stitched);
    }

    if triples.is_empty() {
        return None;
    }
    Some(triples)
}

/// Build the per-branch `EdgeTriple` from the plane×surf intersection
/// pair (`triple_a`, `triple_b`) — same stitching `try_levin_pair`
/// uses, factored out so multi-branch enumeration can call it per
/// plane.
#[allow(clippy::too_many_arguments)]
fn stitch_levin_branch(
    triple_a: EdgeTriple,
    triple_b: EdgeTriple,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tol: &Tolerance,
    total_branches: u8,
    branch_kinds: Vec<BranchKind>,
    branch_index: u32,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let cardinality = Some(crate::geometry::edge_triple::types::BranchCardinality {
        total: total_branches,
        kinds: branch_kinds,
        // Multi-branch Levin produces transverse intersections (two
        // distinct real planes). Tangent (rank-1) is routed via
        // `try_tangent_pair` upstream.
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
            branch: if total_branches > 1 {
                BranchSelector::Indexed(branch_index)
            } else {
                BranchSelector::OnlyBranch
            },
            cardinality,
        },
        quality: EdgeTripleQuality::default(),
    };
    Some(super::super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}
