//! Vertex-vertex coincidence detection between two BRep bodies.
//!
//! Two vertices coincide when their positions are within the snap tolerance.
//! This is the lowest level of the intersection hierarchy and is computed
//! first because VV coincidences constrain all higher-level intersections.

use popkern_core::topology::BRep;

use super::result::VertexVertexHit;

/// Find all vertex pairs between bodies A and B that are within `snap_tol`.
///
/// Uses a simple brute-force approach for now (O(n*m) where n, m are vertex
/// counts). For bodies with < 1000 vertices this is fast enough. BVH
/// acceleration will be added in a later pass for large bodies.
pub(crate) fn detect_vertex_vertex(a: &BRep, b: &BRep, snap_tol: f64) -> Vec<VertexVertexHit> {
    let mut hits = Vec::new();
    let snap_tol_sq = snap_tol * snap_tol;

    for (va_h, va) in &a.vertices {
        let pa = va.position;
        for (vb_h, vb) in &b.vertices {
            let pb = vb.position;
            let dx = pa.x - pb.x;
            let dy = pa.y - pb.y;
            let dz = pa.z - pb.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;

            if dist_sq < snap_tol_sq {
                hits.push(VertexVertexHit {
                    vertex_a: va_h,
                    vertex_b: vb_h,
                    distance: dist_sq.sqrt(),
                });
            }
        }
    }

    hits
}
