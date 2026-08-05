//! Rail-spur recut for interior saddle-ring junctions (PK-234).
//!
//! Strip insertion leaves the saddle corner as a mitre complex: on
//! the reflex side a strip's support boundary is a straight CHORD
//! from the far rail point to the mitre tip (measured on the
//! notched-cube saddle: f10's cube-face boundary runs
//! (1,1,0.52)→(1,0.675,0.65), 2.3e-2 off the strip's own cylinder
//! mid-span), so the pushed boundary junction lands INTERIOR to the
//! chord-bounded trim and there is no rail edge to split.
//!
//! The recut restores the true rail as topology:
//!
//! 1. `mev` a rail spur from the chord's on-rail far endpoint to the
//!    junction foot (exact rail-line segment);
//! 2. the caller `mef`s the arc from the foot — the strip splits
//!    into the main strip and a corner sliver (overreach band +
//!    between-station piece);
//! 3. `kef` the chord — the sliver returns to the SUPPORT face, the
//!    spur becomes the restored strip↔support rail, and the junction
//!    sits on the support face's loop for the bridge `mef`.
//!
//! Everything inside the future ring (the sliver's corner portion,
//! the mitre wedges) is absorbed into the blend face by the region
//! merge afterwards; the transient sliver-in-support ownership never
//! survives to a served result.

use super::*;
use popkern_core::geometry::Vec3;

/// How far a candidate chord endpoint may sit off the rail line and
/// still anchor the spur (the far endpoint is where the carve was
/// still clean, so genuine anchors measure ~1e-12; the mitre-tip
/// endpoint sits ~1e-1 off).
const SPUR_ANCHOR_ON_RAIL_BAND: f64 = 1.0e-6;

/// Result of a rail-spur recut: the junction vertex (spur tip on the
/// restored rail) plus the chord edge to `kef` after the arc `mef`
/// (the bridge phase re-derives the support face topologically).
pub(super) struct RailSpurRecut {
    pub vertex: VertexHandle,
    pub chord_edge: EdgeHandle,
}

/// Attempt the recut on `strip_face` for a junction `point` whose
/// rail edge is missing. The rail direction comes from the strip
/// surface's iso direction (cylinder axis for the analytic strip
/// class); the chord is identified geometrically: the unique
/// strip-loop edge departing the rail line (one endpoint on it, one
/// off) whose far-side face's surface contains the junction.
pub(super) fn rail_spur_recut(
    brep: &mut BRep,
    strip_face: FaceHandle,
    strip_surface: &popkern_core::geometry::SurfaceGeometry,
    point: Point3,
) -> Result<RailSpurRecut, String> {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    let rail_dir = rail_direction(strip_surface)
        .ok_or_else(|| "strip surface has no rail direction".to_string())?;

    // The chord: departs the rail line, and its far-side face is the
    // support the junction/bridge live on.
    let mut chord: Option<(EdgeHandle, VertexHandle, f64, FaceHandle)> = None;
    for he in face_loop_half_edges(brep, strip_face)? {
        let hd = &brep.half_edges[he];
        let edge = hd.edge;
        let he0 = brep.edges[edge].half_edge;
        let v_a = brep.half_edges[he0].origin;
        let v_b = brep.half_edges[brep.half_edges[he0].twin].origin;
        let d_a = distance_to_line(&brep.vertices[v_a].position, &point, &rail_dir);
        let d_b = distance_to_line(&brep.vertices[v_b].position, &point, &rail_dir);
        let (v_on, d_on, d_off) = if d_a <= d_b {
            (v_a, d_a, d_b)
        } else {
            (v_b, d_b, d_a)
        };
        if d_on > SPUR_ANCHOR_ON_RAIL_BAND || d_off <= SPUR_ANCHOR_ON_RAIL_BAND {
            continue;
        }
        let other = brep.loops[brep.half_edges[hd.twin].loop_].face;
        let Some(surface) = brep.faces.get(other).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let contains = project_point_to_uv(&surface.geometry, &point)
            .map(|(u, v)| {
                (&surface.geometry as &dyn SurfaceEval)
                    .point_at(u, v)
                    .distance_to(&point)
            })
            .is_some_and(|d| d <= 1.0e-6);
        if !contains {
            continue;
        }
        if let Some((existing, _, _, _)) = chord {
            return Err(format!(
                "ambiguous rail-departing chords on {strip_face:?} ({existing:?} vs {edge:?})"
            ));
        }
        chord = Some((edge, v_on, d_on, other));
    }
    let (chord_edge, v_far, _d_far, _support_face) = chord.ok_or_else(|| {
        format!(
            "no rail-departing chord with a junction-bearing far side on {strip_face:?} \
             — the mitre complex here is not the chord class"
        )
    })?;

    // Junction foot: exact projection onto the rail line through the
    // anchor (the pushed junction carries the trim fit's ~1e-3
    // wobble; the restored rail is the line through the clean far
    // anchor).
    let p_far = brep.vertices[v_far].position;
    let foot = project_to_line(&point, &p_far, &rail_dir);

    let he_far = half_edge_at_vertex_on_face(brep, strip_face, v_far)?;
    let mev = brep
        .mev(he_far, foot)
        .map_err(|e| format!("mev rail spur: {e}"))?;
    if let Some(edge) = brep.edges.get_mut(mev.edge) {
        edge.curve = popkern_core::geometry::curve::CurveOnEdge::line_segment(p_far, foot, 1.0e-12);
    }
    Ok(RailSpurRecut {
        vertex: mev.vertex,
        chord_edge,
    })
}

/// The strip's rail direction: the iso direction along the sweep
/// (cylinder/cone axis). General swept strips decline — their recut
/// class needs a measured fixture first.
fn rail_direction(surface: &popkern_core::geometry::SurfaceGeometry) -> Option<Vec3> {
    use popkern_core::geometry::SurfaceGeometry;
    match surface {
        SurfaceGeometry::Cylinder { axis, .. } | SurfaceGeometry::Cone { axis, .. } => {
            axis.normalized()
        }
        _ => None,
    }
}

fn distance_to_line(p: &Point3, on_line: &Point3, dir: &Vec3) -> f64 {
    let w = *p - *on_line;
    let along = w.dot(dir);
    (w - *dir * along).length()
}

fn project_to_line(p: &Point3, on_line: &Point3, dir: &Vec3) -> Point3 {
    let w = *p - *on_line;
    let along = w.dot(dir);
    *on_line + *dir * along
}
