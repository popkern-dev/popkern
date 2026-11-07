//! Exporting a general body (PK-1108).
//!
//! Solids that touch along edges or at vertices are one general body:
//! their faces share the contact edges and vertices. STEP writes such a
//! body as OCC does, one `CLOSED_SHELL` listing every solid's faces with
//! the shared `EDGE_CURVE`s and `VERTEX_POINT`s used by each, in one
//! `MANIFOLD_SOLID_BREP`; it reads back as manifold pieces meeting at
//! recorded contacts, or glued into the general body again
//! (`StepNonManifold::General`). A face between two solids is listed
//! once, as OCC writes two solids sharing a face (it reads back as both
//! solids). Forms with no such representation yet decline with what
//! stands in the way: a lamina or sheet face, a cavity, wire edges or
//! acorns.

use popkern_core::topology::{BRep, RegionKind};

use crate::error::IoError;
use crate::step::export::{ExportContext, list, ref_, str_param};
use crate::step::parser::types::EntityId;

/// Emit `brep` (a general body) as one manifold solid over one closed
/// shell of every face, or decline.
pub(super) fn emit_general(
    ctx: &mut ExportContext,
    brep: &BRep,
    faces: Vec<Vec<EntityId>>,
) -> Result<Vec<EntityId>, IoError> {
    let decline = |why: &str| {
        Err(IoError::Export(format!(
            "general body has no STEP form yet: {why}"
        )))
    };
    if brep
        .shells
        .values()
        .any(|s| !s.wire_edges.is_empty() || s.acorn.is_some())
    {
        return decline("wire edges or acorns");
    }
    let kind = |s| brep.shell_region(s).and_then(|r| brep.region_kind(r));
    for face in brep.faces.values() {
        match (kind(face.back), kind(face.front)) {
            (Some(RegionKind::Solid), Some(RegionKind::Void | RegionKind::Solid)) => {}
            _ => return decline("a lamina or sheet face"),
        }
    }
    if brep
        .solid_regions()
        .iter()
        .any(|&r| brep.region_shells(r).len() > 1)
    {
        return decline("a cavity");
    }
    let refs: Vec<_> = faces.into_iter().flatten().map(ref_).collect();
    let shell = ctx.emit("CLOSED_SHELL", vec![str_param(""), list(refs)]);
    Ok(vec![ctx.emit(
        "MANIFOLD_SOLID_BREP",
        vec![str_param(""), ref_(shell)],
    )])
}
