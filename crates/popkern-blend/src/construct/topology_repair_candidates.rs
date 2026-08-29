//! Chord-pair candidate discovery for topology repair.

use super::*;

pub(super) fn is_cylinder_face(brep: &BRep, fh: FaceHandle) -> bool {
    brep.faces
        .get(fh)
        .and_then(|f| f.surface.as_ref())
        .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }))
}

#[derive(Debug, Clone, Copy)]
pub(super) struct ChordPairCandidate {
    /// The off-cylinder vertex (the "elbow" of the chord pair).
    pub(super) off_cyl_vertex: VertexHandle,
    /// Half-edge in the cylinder face's outer loop ending AT the
    /// off-cylinder vertex (incoming).
    pub(super) he_in: HalfEdgeHandle,
    /// Half-edge in the cylinder face's outer loop starting FROM
    /// the off-cylinder vertex (outgoing).
    pub(super) he_out: HalfEdgeHandle,
    /// Vertices flanking the off-cylinder vertex (anchor vertices
    /// on the cylinder surface).
    pub(super) anchor_in: VertexHandle,
    pub(super) anchor_out: VertexHandle,
}

pub(super) fn find_chord_pair_candidates(brep: &BRep, fh: FaceHandle) -> Vec<ChordPairCandidate> {
    let Some(face) = brep.faces.get(fh) else {
        return Vec::new();
    };
    let outer_loop = face.outer_loop;
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return Vec::new();
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let cylinder = face.surface.as_ref().map(|sof| sof.geometry.clone());
    let Some(SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    }) = cylinder
    else {
        return Vec::new();
    };

    // Walk the outer loop and collect each vertex's incoming half-
    // edge along with the next half-edge (outgoing).
    let mut hes: Vec<HalfEdgeHandle> = Vec::new();
    let mut cursor = start_he;
    let mut safety = 64;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cursor);
        let next = brep.half_edges[cursor].next;
        cursor = next;
        if cursor == start_he {
            break;
        }
    }
    if hes.len() < 3 {
        return Vec::new();
    }

    // For each vertex (the origin of each half-edge): if the vertex
    // is OFF the cylinder, it's a candidate elbow. The incoming
    // half-edge for that vertex is hes[(i + n - 1) % n] (the prev),
    // the outgoing is hes[i].
    let n = hes.len();
    let mut out = Vec::new();
    for i in 0..n {
        let he_out = hes[i];
        let he_in = hes[(i + n - 1) % n];
        let Some(out_data) = brep.half_edges.get(he_out) else {
            continue;
        };
        let Some(in_data) = brep.half_edges.get(he_in) else {
            continue;
        };
        let v = out_data.origin;
        let Some(v_data) = brep.vertices.get(v) else {
            continue;
        };
        let p = v_data.position;
        if vertex_on_cylinder(&p, origin, axis, radius) {
            continue;
        }
        let anchor_in = in_data.origin;
        let twin_out = out_data.twin;
        let Some(twin_data) = brep.half_edges.get(twin_out) else {
            continue;
        };
        let anchor_out = twin_data.origin;
        out.push(ChordPairCandidate {
            off_cyl_vertex: v,
            he_in,
            he_out,
            anchor_in,
            anchor_out,
        });
    }
    out
}

pub(super) fn vertex_on_cylinder(
    p: &Point3,
    origin: Point3,
    axis: popkern_core::geometry::Vec3,
    radius: f64,
) -> bool {
    let rel = popkern_core::geometry::Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    let axial = rel.dot(&axis) / axis.length().max(1e-12);
    let axis_unit = if axis.length() > 1e-12 {
        popkern_core::geometry::Vec3::new(
            axis.x / axis.length(),
            axis.y / axis.length(),
            axis.z / axis.length(),
        )
    } else {
        return false;
    };
    let radial = popkern_core::geometry::Vec3::new(
        rel.x - axial * axis_unit.x,
        rel.y - axial * axis_unit.y,
        rel.z - axial * axis_unit.z,
    );
    let r = radial.length();
    (r - radius).abs() <= radius * 0.05
}
