//! The phantom flap region of a blend×blend pair corner (PK-583).
//!
//! The pair-corner settlement retires the planar flap the two stub
//! chains enclose. The flap is usually one face, but a mitred
//! neighbour can leave it split into coplanar pieces. On the step
//! fixture (two saddles joined by a filleted reflex edge), the mitre
//! chord from (1, 0.94) to (0.94, 1) split the C1 flap into two
//! triangles. The one-face recognizer then skipped the corner, and its
//! stale chords shipped 2.16e-2 off the strips (check 9 declined the
//! whole fillet).
//!
//! The region is recognized when every face bordering the chains is a
//! plane coplanar with the first (normals parallel, vertices within
//! the stale tolerance), has no inner loops, and together they are
//! bounded by exactly the chain edges. Their shared edges must form a
//! tree (one fewer than the faces, so no interior vertex). The
//! settlement then merges the region with `kef` into one flap before
//! its surgery. That runs inside the same snapshot and rollback, and
//! a rejected `kef` declines as `EulerRejected`.

use std::collections::HashMap;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Largest normal sine two flap pieces may differ by and still count
/// as one plane (the pieces are split copies of one support plane,
/// which agree to rounding).
const COPLANAR_SINE: f64 = 1.0e-9;

/// The flap region enclosed by a pair corner's stub chains.
pub(super) struct FlapRegion {
    /// The region's first face (survives a merge).
    pub(super) flap: FaceHandle,
    /// Edges shared by two region faces, to dissolve with `kef`
    /// (empty for a single-face flap).
    pub(super) interior: Vec<EdgeHandle>,
}

/// The flap region bordered by `chains` (half-edges on the strip
/// side), or `None` when it is not a coplanar planar region bounded by
/// exactly those edges.
pub(super) fn flap_region(
    brep: &BRep,
    chains: &[HalfEdgeHandle],
    strips: (FaceHandle, FaceHandle),
    stale_tol: f64,
) -> Option<FlapRegion> {
    let mut faces: Vec<FaceHandle> = Vec::new();
    let mut chain_edges: Vec<EdgeHandle> = Vec::new();
    for &he in chains {
        let h = brep.half_edges.get(he)?;
        let face = brep.loops.get(brep.half_edges.get(h.twin)?.loop_)?.face;
        if face == strips.0 || face == strips.1 {
            return None;
        }
        if !faces.contains(&face) {
            faces.push(face);
        }
        chain_edges.push(h.edge);
    }
    let (origin, normal) = plane_of(brep, *faces.first()?)?;
    // Edge multiplicity over the region's outer loops.
    let mut count: HashMap<EdgeHandle, usize> = HashMap::new();
    for &face in &faces {
        let (o, n) = plane_of(brep, face)?;
        if n.cross(&normal).length() > COPLANAR_SINE || (o - origin).dot(&normal).abs() > stale_tol
        {
            return None;
        }
        let f = brep.faces.get(face)?;
        if !f.inner_loops.is_empty() {
            return None;
        }
        let start = brep.loops.get(f.outer_loop)?.half_edge;
        let mut he = start;
        for _ in 0..brep.half_edges.len() {
            let h = brep.half_edges.get(he)?;
            let p = brep.vertices.get(h.origin)?.position;
            if (p - origin).dot(&normal).abs() > stale_tol {
                return None;
            }
            *count.entry(h.edge).or_default() += 1;
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    let mut boundary: Vec<EdgeHandle> = count
        .iter()
        .filter(|(_, c)| **c == 1)
        .map(|(e, _)| *e)
        .collect();
    let mut interior: Vec<EdgeHandle> = count
        .iter()
        .filter(|(_, c)| **c == 2)
        .map(|(e, _)| *e)
        .collect();
    boundary.sort();
    interior.sort();
    chain_edges.sort();
    chain_edges.dedup();
    if boundary != chain_edges
        || boundary.len() + interior.len() != count.len()
        || interior.len() + 1 != faces.len()
    {
        return None;
    }
    Some(FlapRegion {
        flap: faces[0],
        interior,
    })
}

/// Merge the region's pieces into one flap face; returns the survivor.
pub(super) fn merge_flap_region(
    brep: &mut BRep,
    region: &FlapRegion,
) -> Result<FaceHandle, String> {
    let mut flap = region.flap;
    for &edge in &region.interior {
        let he = brep
            .edges
            .get(edge)
            .map(|e| e.half_edge)
            .ok_or_else(|| format!("flap edge {edge:?} vanished"))?;
        let kef = brep
            .kef(he)
            .map_err(|e| format!("kef flap edge {edge:?}: {e}"))?;
        flap = brep.loops[kef.surviving_loop].face;
    }
    Ok(flap)
}

fn plane_of(brep: &BRep, face: FaceHandle) -> Option<(Point3, Vec3)> {
    match brep.faces.get(face)?.surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Plane { origin, normal, .. }) => Some((*origin, *normal)),
        _ => None,
    }
}
