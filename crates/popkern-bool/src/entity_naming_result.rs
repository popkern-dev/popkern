//! Names and lineage of a result built from a target body and tool
//! bodies: one tool for a two-operand boolean, several for the
//! partition's operations (split, N-ary boolean, imprint). The target
//! is operand `a`: its pool seeds the result's and its unchanged
//! entities clear of every tool keep their names.

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::{entity_lineage, entity_lineage_carry, entity_naming, entity_naming_edges};

/// Name `result`'s faces, edges and vertices after `a` and `tools`, and
/// record each one's lineage, under `ctx`'s naming scheme.
pub(crate) fn name_result(
    result: &mut BRep,
    a: &BRep,
    tools: &[&BRep],
    ctx: &mut TrackingContext,
    tol: &Tolerance,
) {
    let mut input_a_faces = entity_naming::collect_face_data(a);
    // Each tool's faces, edges and vertices beside the tool they came
    // from (their names come from that tool's pool).
    let mut input_tool_faces: Vec<(usize, _)> = tools
        .iter()
        .enumerate()
        .flat_map(|(k, b)| {
            entity_naming::collect_face_data(b)
                .into_iter()
                .map(move |f| (k, f))
        })
        .collect();
    let mut input_edges: Vec<(Option<usize>, _)> = entity_naming_edges::collect_edge_data(a)
        .into_iter()
        .map(|e| (None, e))
        .collect();
    let mut input_vertices: Vec<(Option<usize>, _)> = entity_naming_edges::collect_vertex_data(a)
        .into_iter()
        .map(|v| (None, v))
        .collect();
    for (k, b) in tools.iter().enumerate() {
        input_edges.extend(
            entity_naming_edges::collect_edge_data(b)
                .into_iter()
                .map(|e| (Some(k), e)),
        );
        input_vertices.extend(
            entity_naming_edges::collect_vertex_data(b)
                .into_iter()
                .map(|v| (Some(k), v)),
        );
    }
    let source = |k: Option<usize>| k.map_or(a, |k| tools[k]);

    // Naming scheme 3 (PK-731): the result pool starts as the target
    // operand's (a shared clone), as every other operation's result
    // carries its input's pool: the target's names keep their ids, so
    // importing them is the identity. Scheme 2 built the pool by
    // importing every operand name, ancestors first, which re-interned
    // the target's whole history on every boolean (49k names at the
    // 64th slot of a disc).
    let seeded = ctx.naming_scheme >= 3;
    if seeded {
        result.name_pool = a.name_pool.clone();
    }
    if ctx.naming_scheme >= 2 {
        // Naming scheme 2 (PK-677): each operand name comes into the
        // result pool by value, ancestors first, so a tool-side name's
        // ancestor ids are this pool's (scheme 1 interned the tool's ids
        // as they were).
        let pool = &mut result.name_pool;
        for (name, ..) in input_a_faces.iter_mut().filter(|_| !seeded) {
            entity_naming::import_name(name, a, pool);
        }
        for (k, (name, ..)) in &mut input_tool_faces {
            entity_naming::import_name(name, tools[*k], pool);
        }
        for (k, (name, ..)) in &mut input_edges {
            if k.is_some() || !seeded {
                entity_naming::import_name(name, source(*k), pool);
            }
        }
        for (k, (name, _)) in &mut input_vertices {
            if k.is_some() || !seeded {
                entity_naming::import_name(name, source(*k), pool);
            }
        }
    }
    let input_tool_faces: Vec<_> = input_tool_faces.into_iter().map(|(_, f)| f).collect();
    let input_edges: Vec<_> = input_edges.into_iter().map(|(_, e)| e).collect();
    let input_vertices: Vec<_> = input_vertices.into_iter().map(|(_, v)| v).collect();
    // Naming scheme 4 (PK-803): entities carried unchanged and clear of
    // every tool keep `a`'s names and lineage.
    let carry = (seeded && ctx.naming_scheme >= 4)
        .then(|| entity_lineage_carry::CarryIndex::new(a, tools, 10.0 * tol.linear()));
    let carry_a = carry.as_ref().map(|c| (c, a));
    entity_naming::assign_boolean_names(
        result,
        &input_a_faces,
        &input_tool_faces,
        ctx,
        tol,
        carry_a,
    );
    // Lineage beside the names (PK-633): which input faces each output
    // face is a piece of, so labels follow trimmed faces.
    entity_lineage::record_face_lineage(result, a, tools, tol, carry.as_ref());
    entity_naming_edges::assign_boolean_edge_vertex_names(
        result,
        &input_edges,
        &input_vertices,
        ctx,
        tol,
        carry_a,
    );
    entity_lineage::record_edge_lineage(result, a, tools, tol, carry.as_ref());
    popkern_core::naming_lineage::prune_unreadable_lineage(result);
}
