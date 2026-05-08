//! Diagnostic dumps for topology repair.

use super::topology_repair_geometry::off_cyl_distance;
use super::*;

pub(super) fn dump_cylinder_face(brep: &BRep, fh: FaceHandle) {
    let Some(face) = brep.faces.get(fh) else {
        return;
    };
    let cyl = face.surface.as_ref().map(|sof| &sof.geometry);
    let Some(SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    }) = cyl
    else {
        return;
    };
    let outer_loop = face.outer_loop;
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return;
    };
    let start_he = loop_data.half_edge;
    eprintln!(
        "[topo_dump] cyl face {fh:?} cyl(origin=({:.3},{:.3},{:.3}) axis=({:.3},{:.3},{:.3}) r={:.4}):",
        origin.x, origin.y, origin.z, axis.x, axis.y, axis.z, radius
    );
    let mut cursor = start_he;
    let mut safety = 64;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        let v = he.origin;
        let p = brep.vertices[v].position;
        let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
        let twin_surf = brep
            .faces
            .get(twin_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| match &s.geometry {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cyl",
                SurfaceGeometry::Sphere { .. } => "Sph",
                SurfaceGeometry::Cone { .. } => "Cone",
                SurfaceGeometry::Torus { .. } => "Torus",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                SurfaceGeometry::Gregory(_) => "Gregory",
            })
            .unwrap_or("?");
        let off = off_cyl_distance(&p, *origin, *axis, *radius);
        let on = off.abs() <= radius * 0.05;
        eprintln!(
            "[topo_dump]   v{:?} p=({:+.4},{:+.4},{:+.4}) off_cyl={:+.5} {} twin_face={twin_face:?}({twin_surf})",
            v,
            p.x,
            p.y,
            p.z,
            off,
            if on { "ON" } else { "OFF" }
        );
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }
}
