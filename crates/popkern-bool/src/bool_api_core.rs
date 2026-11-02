//! The body of [`crate::boolean`] and [`crate::boolean_with_options`]:
//! input checks, the pipeline under a non-manifold policy, and the
//! naming and lineage of the result.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolError;
use crate::{
    BoolOperation, TrackingContext, entity_lineage, entity_lineage_carry, entity_naming,
    entity_naming_edges, operation, options, sheet,
};

/// [`crate::boolean`] under `policy`.
pub(crate) fn boolean_under(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
    policy: options::NonManifoldPolicy,
) -> Result<BRep, BoolError> {
    // §2.3: reject sheet body inputs with a structured error rather
    // than letting volume-based classification produce undefined
    // results. Future work will add boundary-relative semantics; for
    // now this is the only safe response.
    let a_is_sheet = sheet::is_sheet_body(&a);
    let b_is_sheet = sheet::is_sheet_body(&b);
    if a_is_sheet || b_is_sheet {
        return Err(BoolError::OpenShellNotSupported {
            a: a_is_sheet,
            b: b_is_sheet,
        });
    }

    // Collect input face data before the boolean consumes the inputs.
    let mut input_a_faces = if tracking.is_some() {
        entity_naming::collect_face_data(&a)
    } else {
        Vec::new()
    };
    let mut input_b_faces = if tracking.is_some() {
        entity_naming::collect_face_data(&b)
    } else {
        Vec::new()
    };
    // Edges and vertices too (PK-386): a dress-up after a boolean is
    // recorded by stable name, so the survivors must keep theirs.
    // Where operand b's entries start in the combined edge / vertex
    // lists (their names come from b's pool).
    let (mut input_edges, mut input_vertices, b_edges_from, b_vertices_from) = if tracking.is_some()
    {
        let mut edges = entity_naming_edges::collect_edge_data(&a);
        let b_edges_from = edges.len();
        edges.extend(entity_naming_edges::collect_edge_data(&b));
        let mut vertices = entity_naming_edges::collect_vertex_data(&a);
        let b_vertices_from = vertices.len();
        vertices.extend(entity_naming_edges::collect_vertex_data(&b));
        (edges, vertices, b_edges_from, b_vertices_from)
    } else {
        (Vec::new(), Vec::new(), 0, 0)
    };

    // No perturbation retry: plan Stage 7d retired it after it fired on
    // no boolean over three whole-suite runs (`boolean_with_perturbation`
    // remains an explicit opt-in).
    let mut result = operation::execute_boolean_with_policy(op, &a, &b, &tol, policy)?;
    if let Some(ctx) = tracking {
        // Naming scheme 3 (PK-731): the result pool starts as the
        // target operand's (a shared clone), as every other operation's
        // result carries its input's pool: the target's names keep
        // their ids, so importing them is the identity. Scheme 2 built
        // the pool by importing every operand name, ancestors first,
        // which re-interned the target's whole history on every
        // boolean (49k names at the 64th slot of a disc).
        let seeded = ctx.naming_scheme >= 3;
        if seeded {
            result.name_pool = a.name_pool.clone();
        }
        if ctx.naming_scheme >= 2 {
            // Naming scheme 2 (PK-677): each operand name comes into the
            // result pool by value, ancestors first, so a tool-side
            // name's ancestor ids are this pool's (scheme 1 interned the
            // tool's ids as they were).
            let pool = &mut result.name_pool;
            for (name, ..) in input_a_faces.iter_mut().filter(|_| !seeded) {
                entity_naming::import_name(name, &a, pool);
            }
            for (name, ..) in &mut input_b_faces {
                entity_naming::import_name(name, &b, pool);
            }
            for (i, (name, ..)) in input_edges.iter_mut().enumerate() {
                if i >= b_edges_from || !seeded {
                    entity_naming::import_name(name, if i < b_edges_from { &a } else { &b }, pool);
                }
            }
            for (i, (name, _)) in input_vertices.iter_mut().enumerate() {
                if i >= b_vertices_from || !seeded {
                    let source = if i < b_vertices_from { &a } else { &b };
                    entity_naming::import_name(name, source, pool);
                }
            }
        }
        // Naming scheme 4 (PK-803): entities carried unchanged and clear
        // of the tool keep `a`'s names and lineage.
        let carry = (seeded && ctx.naming_scheme >= 4)
            .then(|| entity_lineage_carry::CarryIndex::new(&a, &b, 10.0 * tol.linear()));
        let carry_a = carry.as_ref().map(|c| (c, &a));
        entity_naming::assign_boolean_names(
            &mut result,
            &input_a_faces,
            &input_b_faces,
            ctx,
            &tol,
            carry_a,
        );
        // Lineage beside the names (PK-633): which input faces each
        // output face is a piece of, so labels follow trimmed faces.
        entity_lineage::record_face_lineage(&mut result, &a, &b, &tol, carry.as_ref());
        entity_naming_edges::assign_boolean_edge_vertex_names(
            &mut result,
            &input_edges,
            &input_vertices,
            ctx,
            &tol,
            carry_a,
        );
        entity_lineage::record_edge_lineage(&mut result, &a, &b, &tol, carry.as_ref());
        popkern_core::naming_lineage::prune_unreadable_lineage(&mut result);
    }
    Ok(result)
}
