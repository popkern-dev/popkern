//! Live outer-seam refresh for materialized exact-boundary sectors.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::super::sector::ExactBoundaryCapSectorInput;
use super::{ExactBoundaryTopologyDeclineReason, ExactBoundaryTopologyDiagnostic};

pub(super) fn refresh_sector_outer_seams(
    brep: &mut BRep,
    sector: &mut ExactBoundaryCapSectorInput,
) -> Result<(), ExactBoundaryTopologyDiagnostic> {
    let fillet = live_outer_seam(
        brep,
        sector,
        sector.fillet_end_edge,
        sector.fillet_support_face,
    )?;
    sector.fillet_support_face = fillet.1;
    sector.fillet_seam.cap_half_edge = fillet.0;
    sector.fillet_seam.support_face = fillet.1;
    sector.fillet_seam.support_half_edge = fillet.2;
    sector.fillet_seam.cap_pcurve_present = brep
        .half_edges
        .get(fillet.0)
        .is_some_and(|half_edge| half_edge.pcurve.is_some());
    sector.fillet_seam.support_pcurve_present = brep
        .half_edges
        .get(fillet.2)
        .is_some_and(|half_edge| half_edge.pcurve.is_some());
    // Welded sectors (shared-corner contract) have no bridge seam.
    let (Some(face_edge), Some(bridge_support_face)) =
        (sector.face_edge, sector.bridge_support_face)
    else {
        return Ok(());
    };
    let face = live_outer_seam(brep, sector, face_edge, bridge_support_face)?;
    sector.bridge_support_face = Some(face.1);
    if let Some(bridge_seam) = sector.bridge_seam.as_mut() {
        bridge_seam.cap_half_edge = face.0;
        bridge_seam.support_face = face.1;
        bridge_seam.support_half_edge = face.2;
        bridge_seam.cap_pcurve_present = brep
            .half_edges
            .get(face.0)
            .is_some_and(|half_edge| half_edge.pcurve.is_some());
        bridge_seam.support_pcurve_present = brep
            .half_edges
            .get(face.2)
            .is_some_and(|half_edge| half_edge.pcurve.is_some());
    }
    Ok(())
}

fn live_outer_seam(
    brep: &mut BRep,
    sector: &ExactBoundaryCapSectorInput,
    edge: EdgeHandle,
    expected_support_face: FaceHandle,
) -> Result<(HalfEdgeHandle, FaceHandle, HalfEdgeHandle), ExactBoundaryTopologyDiagnostic> {
    rehydrate_unsurfaced_support_face(brep, edge, sector.cap_face, expected_support_face);
    let edge_data = brep.edges.get(edge).ok_or_else(|| {
        ExactBoundaryTopologyDiagnostic::new(
            ExactBoundaryTopologyDeclineReason::TopologyWiringFailed,
            Some(sector.sector_name.clone()),
            None,
        )
    })?;
    let primary = edge_data.half_edge;
    let twin = brep.half_edges.get(primary).map(|half_edge| half_edge.twin);
    for half_edge in [Some(primary), twin] {
        let Some(half_edge) = half_edge else {
            continue;
        };
        if face_for_half_edge(brep, half_edge) != Some(sector.cap_face) {
            continue;
        }
        let support_half_edge = brep.half_edges.get(half_edge).map(|data| data.twin);
        let support_face = support_half_edge.and_then(|support_half_edge| {
            face_for_half_edge(brep, support_half_edge).map(|face| (support_half_edge, face))
        });
        let Some((support_half_edge, support_face)) = support_face else {
            if let Some((semantic_support_half_edge, semantic_support_face)) =
                geometrically_valid_support_side(
                    brep,
                    edge,
                    sector.cap_face,
                    half_edge,
                    expected_support_face,
                )
            {
                return Ok((half_edge, semantic_support_face, semantic_support_half_edge));
            }
            if let Some(semantic_support_half_edge) =
                same_oriented_side_on_face(brep, half_edge, expected_support_face)
            {
                return Ok((half_edge, expected_support_face, semantic_support_half_edge));
            }
            let expected_support_is_cylinder = brep
                .faces
                .get(expected_support_face)
                .and_then(|face| face.surface.as_ref())
                .is_some_and(|surface| {
                    matches!(surface.geometry, SurfaceGeometry::Cylinder { .. })
                });
            if expected_support_is_cylinder
                && let Some((semantic_support_half_edge, semantic_support_face)) =
                    same_oriented_cylinder_support_side(brep, sector.cap_face, half_edge)
            {
                return Ok((half_edge, semantic_support_face, semantic_support_half_edge));
            }
            return Err(ExactBoundaryTopologyDiagnostic::new(
                ExactBoundaryTopologyDeclineReason::TopologyWiringFailed,
                Some(sector.sector_name.clone()),
                None,
            ));
        };
        if support_face != sector.cap_face {
            let tolerance = brep.tolerance.linear().max(1.0e-6);
            let support_is_geometric = support_edge_residual(brep, edge, support_face)
                .is_some_and(|residual| residual <= tolerance);
            let expected_is_geometric = support_edge_residual(brep, edge, expected_support_face)
                .is_some_and(|residual| residual <= tolerance);
            if (!support_is_geometric || support_face != expected_support_face)
                && let Some((semantic_support_half_edge, semantic_support_face)) =
                    geometrically_valid_support_side(
                        brep,
                        edge,
                        sector.cap_face,
                        half_edge,
                        expected_support_face,
                    )
                && (!support_is_geometric
                    || expected_is_geometric
                    || semantic_support_face == expected_support_face)
            {
                return Ok((half_edge, semantic_support_face, semantic_support_half_edge));
            }
            if support_face != expected_support_face
                && let Some(semantic_support_half_edge) =
                    same_oriented_side_on_face(brep, half_edge, expected_support_face)
            {
                return Ok((half_edge, expected_support_face, semantic_support_half_edge));
            }
            let expected_support_is_cylinder = brep
                .faces
                .get(expected_support_face)
                .and_then(|face| face.surface.as_ref())
                .is_some_and(|surface| {
                    matches!(surface.geometry, SurfaceGeometry::Cylinder { .. })
                });
            if expected_support_is_cylinder
                && brep
                    .faces
                    .get(support_face)
                    .and_then(|face| face.surface.as_ref())
                    .is_none()
                && let Some((semantic_support_half_edge, semantic_support_face)) =
                    same_oriented_cylinder_support_side(brep, sector.cap_face, half_edge)
            {
                return Ok((half_edge, semantic_support_face, semantic_support_half_edge));
            }
            return Ok((half_edge, support_face, support_half_edge));
        }
    }
    Err(ExactBoundaryTopologyDiagnostic::new(
        ExactBoundaryTopologyDeclineReason::TopologyWiringFailed,
        Some(sector.sector_name.clone()),
        None,
    ))
}

fn rehydrate_unsurfaced_support_face(
    brep: &mut BRep,
    edge: EdgeHandle,
    cap_face: FaceHandle,
    support_face: FaceHandle,
) {
    let Some(face_data) = brep.faces.get(support_face) else {
        return;
    };
    if face_data.surface.is_some() {
        return;
    }
    let Some(loop_sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
        return;
    };
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    let Some(surface) =
        best_surface_for_loop(brep, edge, cap_face, support_face, &loop_sides, tolerance)
    else {
        return;
    };
    if let Some(face_data) = brep.faces.get_mut(support_face) {
        face_data.surface = Some(surface);
    }
}

fn best_surface_for_loop(
    brep: &BRep,
    edge: EdgeHandle,
    cap_face: FaceHandle,
    support_face: FaceHandle,
    loop_sides: &[HalfEdgeHandle],
    tolerance: f64,
) -> Option<SurfaceOnFace> {
    brep.faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cap_face || face == support_face {
                return None;
            }
            let surface = face_data.surface.as_ref()?;
            if matches!(
                surface.geometry,
                SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_)
            ) {
                return None;
            }
            let seam_residual = support_edge_residual(brep, edge, face)?;
            if seam_residual > tolerance {
                return None;
            }
            let loop_residual = loop_surface_residual(brep, loop_sides, &surface.geometry)?;
            (loop_residual <= tolerance).then_some((surface.clone(), loop_residual))
        })
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(surface, _)| surface)
}

fn loop_surface_residual(
    brep: &BRep,
    loop_sides: &[HalfEdgeHandle],
    surface: &SurfaceGeometry,
) -> Option<f64> {
    loop_sides
        .iter()
        .filter_map(|half_edge| {
            let edge = brep.edges.get(brep.half_edges.get(*half_edge)?.edge)?;
            let curve = edge.curve.as_ref()?;
            Some(
                edge_surface_deviation(&curve.geometry, (curve.t_start, curve.t_end), surface, 9)
                    .max_deviation,
            )
        })
        .reduce(f64::max)
}

fn geometrically_valid_support_side(
    brep: &BRep,
    edge: EdgeHandle,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
    expected_support_face: FaceHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    if let Some(side) = same_oriented_side_on_face(brep, cap_half_edge, expected_support_face)
        && support_edge_residual(brep, edge, expected_support_face).is_some_and(|r| r <= tolerance)
    {
        return Some((side, expected_support_face));
    }
    let mut best = None;
    for (face, face_data) in brep.faces.iter() {
        if face == cap_face {
            continue;
        }
        let Some(residual) = support_edge_residual(brep, edge, face) else {
            continue;
        };
        if residual > tolerance {
            continue;
        }
        let Some(sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in sides {
            if coincident_side(brep, cap_half_edge, side) {
                let replace = best
                    .as_ref()
                    .is_none_or(|(_, _, best_residual)| residual < *best_residual);
                if replace {
                    best = Some((side, face, residual));
                }
            }
        }
    }
    best.map(|(side, face, _)| (side, face))
}

fn support_edge_residual(brep: &BRep, edge: EdgeHandle, face: FaceHandle) -> Option<f64> {
    let edge_data = brep.edges.get(edge)?;
    let curve = edge_data.curve.as_ref()?;
    let surface = brep.faces.get(face)?.surface.as_ref()?;
    let nurbs = curve
        .geometry
        .to_nurbs_in_range((curve.t_start, curve.t_end));
    Some(
        edge_surface_deviation(&nurbs, nurbs.parameter_range(), &surface.geometry, 17)
            .max_deviation,
    )
}

fn same_oriented_side_on_face(
    brep: &BRep,
    cap_half_edge: HalfEdgeHandle,
    expected_support_face: FaceHandle,
) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(expected_support_face)?;
    let sides = collect_loop_half_edges(brep, face_data.outer_loop)?;
    sides
        .into_iter()
        .find(|side| same_oriented_side(brep, cap_half_edge, *side))
}

fn same_oriented_cylinder_support_side(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    for (face, face_data) in brep.faces.iter() {
        if face == cap_face {
            continue;
        }
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if !matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        let Some(sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in sides {
            if same_oriented_side(brep, cap_half_edge, side) {
                return Some((side, face));
            }
        }
    }
    None
}

fn same_oriented_side(brep: &BRep, cap_side: HalfEdgeHandle, support_side: HalfEdgeHandle) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(support_start) = brep.half_edges.get(support_side).map(|data| data.origin) else {
        return false;
    };
    if cap_start != support_start {
        return false;
    }
    let Some(cap_end) = super::half_edge_destination_vertex(brep, cap_side) else {
        return false;
    };
    let Some(support_end) = super::half_edge_destination_vertex(brep, support_side) else {
        return false;
    };
    cap_end == support_end
}

fn coincident_side(brep: &BRep, cap_side: HalfEdgeHandle, support_side: HalfEdgeHandle) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(support_start) = brep.half_edges.get(support_side).map(|data| data.origin) else {
        return false;
    };
    let Some(cap_end) = super::half_edge_destination_vertex(brep, cap_side) else {
        return false;
    };
    let Some(support_end) = super::half_edge_destination_vertex(brep, support_side) else {
        return false;
    };
    (cap_start == support_start && cap_end == support_end)
        || (cap_start == support_end && cap_end == support_start)
}

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let loop_data = brep.loops.get(loop_handle)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut cycle = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        cycle.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(cycle);
        }
    }
    None
}

fn face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<FaceHandle> {
    brep.half_edges
        .get(half_edge)
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_data| loop_data.face)
}
