//! Trace helpers for three-quad cap sub-face diagnostics.

use popkern_core::topology::{BRep, FaceHandle};

pub(super) fn trace_cap_subface_neighbors(brep: &BRep, label: &str, faces: &[FaceHandle]) {
    let kind = |face: FaceHandle| -> &'static str {
        use popkern_core::geometry::surface::SurfaceGeometry;
        match brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        {
            Some(SurfaceGeometry::Plane { .. }) => "Plane",
            Some(SurfaceGeometry::Cylinder { .. }) => "Cyl",
            Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
            Some(SurfaceGeometry::Cone { .. }) => "Cone",
            Some(SurfaceGeometry::Torus { .. }) => "Torus",
            Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
            Some(SurfaceGeometry::Gregory(_)) => "Gregory",
            None => "None",
        }
    };
    for (idx, &face) in faces.iter().enumerate() {
        let Some(face_data) = brep.faces.get(face) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut entries = Vec::new();
        for _ in 0..(brep.half_edges.len() + 1) {
            let Some(he) = brep.half_edges.get(cursor) else {
                break;
            };
            let twin_face = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .map(|l| l.face);
            let twin_desc = twin_face
                .map(|f| format!("{f:?}:{}", kind(f)))
                .unwrap_or_else(|| "None".to_string());
            entries.push(format!("{:?}->{twin_desc}", he.edge));
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
        eprintln!(
            "[direct-apex-trace] K-P poke {label} subface[{idx}] {face:?}: {}",
            entries.join(", ")
        );
    }
}
