//! Mixed and boundary corner surgery (PK-326): the two-sector peel
//! with transversal / rim end-circle connectors, plus the shared
//! post-surgery curve re-trim.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::arc_join::offset_corner_position;
use crate::arc_join_corner::CornerFanEntry;
use crate::arc_join_transversal::transversal_arc;
use crate::join_classification::EdgeJoinClass;
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Mixed-corner surgery (PK-326 E1): partition the fan into exactly
/// two sectors — cuts at strip wedges and at middle faces between
/// consecutive classified interior edges — peel the non-retained
/// sector with one `split_vertex`, author the connector as the
/// transversal intersection of the two cut-wedge surfaces, and re-trim
/// every incident curve to the final vertex positions.
#[allow(clippy::too_many_arguments)]
pub(crate) fn construct_mixed_corner(
    skin: &mut BRep,
    v: VertexHandle,
    fan: &[CornerFanEntry],
    strip_faces: &BTreeSet<FaceHandle>,
    classes: &BTreeMap<EdgeHandle, EdgeJoinClass>,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
    reverse_edge: &HashMap<EdgeHandle, EdgeHandle>,
    outer: &BRep,
    p_orig: Point3,
    tolerance: &Tolerance,
    retrim: &mut BTreeSet<EdgeHandle>,
    rim_splits: &mut Vec<(VertexHandle, VertexHandle)>,
) -> Result<(), String> {
    let m = fan.len();
    let classified = |eh: EdgeHandle| {
        matches!(
            classes.get(&eh),
            Some(EdgeJoinClass::Resolvable | EdgeJoinClass::Smooth)
        )
    };
    // Cut wedges: wedge i sits between fan edge i and fan edge i+1.
    // Cuts fall at strip wedges, at sheet-boundary wedges (rim
    // corners), at retained removed-face wedges (midplane openings —
    // their face keeps its ORIGINAL surface and the connector lands in
    // it), and at middle faces between two classified edges.
    let unclassified_wedge = |wedge: Option<FaceHandle>| match wedge {
        None => true,
        Some(f) if strip_faces.contains(&f) => false,
        Some(f) => reverse_face
            .get(&f)
            .is_none_or(|of| !offset_surfaces.contains_key(of)),
    };
    let cuts: Vec<usize> = (0..m)
        .filter(|&i| {
            unclassified_wedge(fan[i].wedge)
                || fan[i].wedge.is_some_and(|f| strip_faces.contains(&f))
                || (classified(fan[i].edge) && classified(fan[(i + 1) % m].edge))
        })
        .collect();
    if cuts.len() != 2 {
        return Err(format!(
            "mixed corner fan has {} cut wedges (2-sector peel needs exactly 2)",
            cuts.len()
        ));
    }

    // Two sectors: entries after cut wedge i up to and including the
    // entry whose after-wedge is the next cut.
    let sector = |from_cut: usize, to_cut: usize| -> Vec<usize> {
        let mut ids = Vec::new();
        let mut i = (from_cut + 1) % m;
        loop {
            ids.push(i);
            if i == to_cut {
                break;
            }
            i = (i + 1) % m;
        }
        ids
    };
    let sectors = [sector(cuts[0], cuts[1]), sector(cuts[1], cuts[0])];

    // Per-sector target: nearest point to the original vertex on the
    // sector's constraining entity.
    let mut targets: Vec<Point3> = Vec::with_capacity(2);
    for ids in &sectors {
        let interior: Vec<usize> = ids
            .iter()
            .copied()
            .filter(|&i| classified(fan[i].edge))
            .collect();
        let target = match interior.len() {
            0 => {
                // Single-face sector: wedge faces inside the sector
                // must agree on one non-strip face.
                let mut faces: BTreeSet<FaceHandle> = BTreeSet::new();
                for &i in &ids[..ids.len() - 1] {
                    if let Some(f) = fan[i].wedge {
                        faces.insert(f);
                    }
                }
                if faces.len() != 1 {
                    return Err("mixed-corner sector spans multiple faces without a \
                                classified edge"
                        .into());
                }
                let face = *faces.iter().next().ok_or("empty mixed-corner sector")?;
                let outer_face = *reverse_face
                    .get(&face)
                    .ok_or("mixed-corner sector face missing from the skin map")?;
                offset_corner_position(offset_surfaces, outer_face, &p_orig)
                    .map(|(q, _)| q)
                    .ok_or("mixed-corner sector projection failed")?
            }
            1 => {
                let p = interior[0];
                let (Some(f1), Some(f2)) = (fan[(p + m - 1) % m].wedge, fan[p].wedge) else {
                    return Err("classified edge beside a sheet boundary".into());
                };
                let of1 = *reverse_face
                    .get(&f1)
                    .ok_or("mixed-corner pair face missing from the skin map")?;
                let of2 = *reverse_face
                    .get(&f2)
                    .ok_or("mixed-corner pair face missing from the skin map")?;
                let (Some(s1), Some(s2)) = (offset_surfaces.get(&of1), offset_surfaces.get(&of2))
                else {
                    return Err("mixed-corner pair faces have no offset surfaces".into());
                };
                crate::curved_recompute::solve_vertex_on_surfaces(&[s1, s2], p_orig, tolerance)
                    .ok_or("mixed-corner pair intersection solve did not converge")?
            }
            n => {
                return Err(format!(
                    "mixed-corner sector holds {n} classified edges (expected at most 1)"
                ));
            }
        };
        targets.push(target);
    }

    // Retained sector: the one whose target the vertex recompute
    // already placed the vertex at.
    let v_pos = skin.vertices[v].position;
    let retained_idx = if targets[0].distance_to(&v_pos) <= targets[1].distance_to(&v_pos) {
        0
    } else {
        1
    };
    let moved_idx = 1 - retained_idx;
    let moved = &sectors[moved_idx];
    let first = fan[moved[0]].out;
    let last = fan[moved[moved.len() - 1]].out;

    let res = skin
        .split_vertex(first, last, targets[moved_idx])
        .map_err(|e| format!("mixed-corner split_vertex failed: {e}"))?;
    skin.vertices[v].position = targets[retained_idx];

    // Connector geometry. Interior corners: the transversal
    // intersection of the two cut-wedge surfaces. Boundary (rim)
    // corners: one cut wedge is the sheet boundary, and the connector
    // is the canal's end cross-section circle — a great circle of the
    // vertex sphere through both targets (the strip's natural end).
    let cut_before = fan[(moved[0] + m - 1) % m].wedge;
    let cut_after = fan[moved[moved.len() - 1]].wedge;
    let arc = match (cut_before, cut_after) {
        (Some(f1), Some(f2)) => {
            let s1 = skin
                .faces
                .get(f1)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.clone())
                .ok_or("mixed-corner cut wedge face has no surface")?;
            let s2 = skin
                .faces
                .get(f2)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.clone())
                .ok_or("mixed-corner cut wedge face has no surface")?;
            transversal_arc(
                &s1,
                &s2,
                targets[retained_idx],
                targets[moved_idx],
                tolerance,
            )
            .ok_or(
                "transversal connector between the cut-wedge surfaces failed \
                    (non-conic pair and the marched fit missed its residual gate)",
            )?
        }
        _ => {
            rim_splits.push((res.vertex, v));
            // Opening rims (removed-face shells): the flush contract
            // pins the connector INTO the removed face's ORIGINAL
            // surface — author it as the strip × removed-surface
            // transversal. Free rims (sheets, no face ever existed on
            // the far side): the canal's end cross-section circle, a
            // great circle of the vertex sphere.
            let removed_surface =
                removed_surface_at_boundary(skin, fan, reverse_face, reverse_edge, outer);
            let strip_geom = [cut_before, cut_after]
                .into_iter()
                .flatten()
                .find(|f| strip_faces.contains(f))
                .and_then(|f| skin.faces.get(f))
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.clone());
            match (removed_surface, strip_geom) {
                (Some(removed), Some(strip)) => transversal_arc(
                    &strip,
                    &removed,
                    targets[retained_idx],
                    targets[moved_idx],
                    tolerance,
                )
                .ok_or(
                    "opening-rim connector (strip × removed original surface) failed                      — the flush contract cannot hold at this corner",
                )?,
                _ => {
                    let d = 0.5
                        * (targets[0].distance_to(&p_orig) + targets[1].distance_to(&p_orig));
                    crate::arc_join_geometry::great_arc(
                        p_orig,
                        d,
                        targets[retained_idx],
                        targets[moved_idx],
                        tolerance,
                    )
                    .ok_or("rim end-circle connector is antipodal or degenerate")?
                }
            }
        }
    };
    skin.edges[res.edge].curve = Some(arc);

    // Schedule every incident curve (rails, resolvable edges) at both
    // corner vertices for re-trim: mixed targets are not the rails'
    // natural similarity images, so parameter ranges must follow the
    // vertices. The re-trim itself runs AFTER all corner surgeries —
    // an edge's far endpoint may belong to a corner not yet peeled.
    for vertex in [v, res.vertex] {
        let Some(anchor) = skin.vertices.get(vertex).map(|vd| vd.outgoing) else {
            continue;
        };
        let mut cursor = anchor;
        for _ in 0..skin.half_edges.len() {
            let Some(he) = skin.half_edges.get(cursor) else {
                break;
            };
            if he.origin != vertex {
                break; // walked off the vertex (boundary fan inconsistency)
            }
            if he.edge != res.edge {
                retrim.insert(he.edge);
            }
            let Some(twin) = skin.half_edges.get(he.twin) else {
                break;
            };
            let next = twin.next;
            if next == anchor {
                break;
            }
            cursor = next;
        }
    }
    Ok(())
}

/// Re-trim an edge's curve parameter range so its endpoints coincide
/// with its vertices' positions (closed curves and curveless edges are
/// left alone). Declines when a vertex does not lie on the curve.
pub(crate) fn retrim_edge_to_vertices(
    skin: &mut BRep,
    eh: EdgeHandle,
    tolerance: &Tolerance,
) -> Result<(), String> {
    use popkern_core::geometry::curve::CurveEval;
    let Some(edge) = skin.edges.get(eh) else {
        return Ok(());
    };
    let Some(curve) = edge.curve.clone() else {
        return Ok(());
    };
    let he = edge.half_edge;
    let Some(he_data) = skin.half_edges.get(he) else {
        return Ok(());
    };
    let Some(twin) = skin.half_edges.get(he_data.twin) else {
        return Ok(());
    };
    if he_data.origin == twin.origin {
        return Ok(()); // closed edge: nothing to anchor
    }
    let a = skin.vertices[he_data.origin].position;
    let b = skin.vertices[twin.origin].position;
    let t_a = curve
        .geometry
        .closest_point(&a)
        .ok_or("re-trim: closest-point solve failed")?;
    let t_b = curve
        .geometry
        .closest_point(&b)
        .ok_or("re-trim: closest-point solve failed")?;
    let gate = (tolerance.linear() * 100.0).max(a.distance_to(&b) * 1.0e-6);
    let d_a = curve.geometry.point_at(t_a).distance_to(&a);
    let d_b = curve.geometry.point_at(t_b).distance_to(&b);
    if d_a > gate || d_b > gate {
        return Err(format!(
            "re-trim: edge {eh:?} vertex leaves its curve (origin off by {d_a:.3e}, \
             dest by {d_b:.3e}, gate {gate:.3e}, origin {a:?}, dest {b:?})"
        ));
    }
    if (t_a - t_b).abs() < 1.0e-12 {
        return Err("re-trim collapsed the edge to a point".into());
    }
    let (t_start, t_end) = if t_a <= t_b { (t_a, t_b) } else { (t_b, t_a) };
    if let Some(edge) = skin.edges.get_mut(eh)
        && let Some(c) = edge.curve.as_mut()
    {
        c.t_start = t_start;
        c.t_end = t_end;
    }
    Ok(())
}

/// The removed face's ORIGINAL surface at a boundary corner, when one
/// exists (opening rims): a rim fan edge maps back to an outer-body
/// edge whose far side face has no image in the skin. Free sheet rims
/// return `None` — no face ever existed there.
fn removed_surface_at_boundary(
    skin: &BRep,
    fan: &[CornerFanEntry],
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
    reverse_edge: &HashMap<EdgeHandle, EdgeHandle>,
    outer: &BRep,
) -> Option<SurfaceGeometry> {
    let skinned_outer: BTreeSet<FaceHandle> = reverse_face.values().copied().collect();
    for entry in fan {
        // Rim edges have a loopless side in the skin.
        let edge = skin.edges.get(entry.edge)?;
        let he = skin.half_edges.get(edge.half_edge)?;
        let twin = skin.half_edges.get(he.twin)?;
        let rim = !skin.loops.contains_key(he.loop_) || !skin.loops.contains_key(twin.loop_);
        if !rim {
            continue;
        }
        let outer_eh = reverse_edge.get(&entry.edge)?;
        let o_edge = outer.edges.get(*outer_eh)?;
        let o_he = outer.half_edges.get(o_edge.half_edge)?;
        let o_twin = outer.half_edges.get(o_he.twin)?;
        for lh in [o_he.loop_, o_twin.loop_] {
            if let Some(l) = outer.loops.get(lh)
                && !skinned_outer.contains(&l.face)
            {
                return outer
                    .faces
                    .get(l.face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| s.geometry.clone());
            }
        }
    }
    None
}
