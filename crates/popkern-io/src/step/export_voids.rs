//! Solid items for exported shells, cavities included (PK-880).
//!
//! A popkern solid owns its outer shell (`Solid::shell`) and any number
//! of cavity shells (shells whose `solid` is it): a hollow subtraction or
//! a cellular cavity partition. One `MANIFOLD_SOLID_BREP` per shell wrote
//! each cavity as a second, positive solid (OpenCASCADE read a dome with
//! a cylindrical void as two solids, 8.378 and 0.006). A solid with
//! cavities is a `BREP_WITH_VOIDS`: its outer `CLOSED_SHELL` and one
//! `ORIENTED_CLOSED_SHELL` per cavity. The schema requires each void's
//! orientation be `.F.` over a shell oriented as the void region's own
//! outside; popkern's cavity faces face into the void (out of the
//! material), so they are WRITTEN reversed (sense and bounds negated,
//! `ExportContext::reverse_face_uses`) and the two reversals cancel.
//! Wrapping them in `ORIENTED_FACE(.F.)` instead is equivalent by the
//! schema, but OpenCASCADE's reader drops such faces from a void shell
//! (measured: the void imported with no faces, the volume of the outer
//! alone).

use popkern_core::topology::{BRep, ShellHandle};

use crate::step::export::{ExportContext, bool_, list, ref_, str_param};
use crate::step::parser::types::{EntityId, Parameter};

/// Whether `shell` is a cavity: a shell of a solid other than its outer.
pub(super) fn is_cavity(brep: &BRep, shell: ShellHandle) -> bool {
    brep.shells
        .get(shell)
        .and_then(|s| brep.solids.get(s.solid))
        .is_some_and(|solid| solid.shell != shell)
}

/// Emit one solid item per popkern solid (cavity shells as voids of
/// their solid), then one `MANIFOLD_SOLID_BREP` per remaining shell.
/// `shells[i]` is the popkern shell `faces[i]` came from (`None` for
/// orphan faces).
pub(super) fn emit_solids(
    ctx: &mut ExportContext,
    brep: &BRep,
    shells: &[Option<ShellHandle>],
    faces: Vec<Vec<EntityId>>,
) -> Vec<EntityId> {
    let mut items = Vec::new();
    let mut used = vec![false; faces.len()];
    let index_of = |sh: ShellHandle| shells.iter().position(|s| *s == Some(sh));
    let mut solid_handles: Vec<_> = brep.solids.keys().collect();
    solid_handles.sort();
    let closed = |ctx: &mut ExportContext, ids: &[EntityId]| {
        let refs: Vec<_> = ids.iter().map(|&id| ref_(id)).collect();
        ctx.emit("CLOSED_SHELL", vec![str_param(""), list(refs)])
    };
    for so in solid_handles {
        let Some(outer) = brep.solids.get(so).and_then(|s| index_of(s.shell)) else {
            continue;
        };
        let voids: Vec<usize> = shells
            .iter()
            .enumerate()
            .filter(|&(i, sh)| {
                i != outer
                    && sh
                        .and_then(|h| brep.shells.get(h))
                        .is_some_and(|s| s.solid == so)
            })
            .map(|(i, _)| i)
            .collect();
        if used[outer] {
            continue;
        }
        let outer_id = closed(ctx, &faces[outer]);
        used[outer] = true;
        if voids.is_empty() {
            items.push(ctx.emit("MANIFOLD_SOLID_BREP", vec![str_param(""), ref_(outer_id)]));
            continue;
        }
        let mut void_refs = Vec::with_capacity(voids.len());
        for i in voids {
            used[i] = true;
            let shell = closed(ctx, &faces[i]);
            void_refs.push(ref_(ctx.emit(
                "ORIENTED_CLOSED_SHELL",
                vec![str_param(""), Parameter::Derived, ref_(shell), bool_(false)],
            )));
        }
        items.push(ctx.emit(
            "BREP_WITH_VOIDS",
            vec![str_param(""), ref_(outer_id), list(void_refs)],
        ));
    }
    for (i, ids) in faces.iter().enumerate() {
        if !used[i] {
            let id = closed(ctx, ids);
            items.push(ctx.emit("MANIFOLD_SOLID_BREP", vec![str_param(""), ref_(id)]));
        }
    }
    items
}
