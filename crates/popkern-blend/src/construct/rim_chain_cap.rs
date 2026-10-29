//! Rim chains whose bevels consume their cap (PK-1002): the inset
//! polygon turned inside out.

use super::{RimChainPlan, cap_wall_trims, line_intersection, trim_line};
use crate::construct::EdgeGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Why a chain's bevels consume its cap, if they do: the inset polygon
/// (adjacent cap-side trims meeting) must keep the cap's shape, every
/// inset edge running the way its rim edge does with positive length.
/// Bevels wider than the cap turn it inside out (a 0.3-tall face
/// ringed at d = 0.2): that is a consumed face, which neither the ring
/// surgery nor the per-edge path builds (the latter served a
/// self-intersecting body); it grows from a smaller ring (PK-1001).
pub(crate) fn chain_consumes_cap(
    brep: &BRep,
    plan: &RimChainPlan,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Option<String> {
    let n = plan.entries.len();
    let join_tol = (tolerance.linear() * 10.0).max(1.0e-9);
    let mut inset = Vec::with_capacity(n);
    for i in 0..n {
        let prev = &edge_geometries[plan.entries[(i + n - 1) % n].geometry_index];
        let cur = &edge_geometries[plan.entries[i].geometry_index];
        let a = trim_line(cap_wall_trims(prev, plan.cap_face).0, join_tol)?;
        let b = trim_line(cap_wall_trims(cur, plan.cap_face).0, join_tol)?;
        inset.push(line_intersection(a, b, join_tol)?);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        let rim = brep.vertices.get(plan.entries[j].corner_vertex)?.position
            - brep.vertices.get(plan.entries[i].corner_vertex)?.position;
        if (inset[j] - inset[i]).dot(&rim.normalized()?) <= join_tol {
            return Some(format!(
                "the bevels consume the cap around rim edge {i} (inset edge reversed or vanished)"
            ));
        }
    }
    None
}
