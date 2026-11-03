//! Bodies from a partition's regions (design 6).

use popkern_core::topology::BRep;

use super::{BodySet, Partition, SetExpr};
use crate::error::{BoolError, ContactKind, NonManifoldContact};
use crate::options::{BoolOptions, NonManifoldPolicy};

impl Partition {
    /// The regularized union of the regions `expr` keeps: a face between
    /// two kept regions is internal and merged away. A non-manifold
    /// result follows `options.non_manifold` as a two-operand boolean
    /// does (declined with its contacts by default, served as manifold
    /// pieces meeting at recorded contacts under `SplitToManifold`). No
    /// region kept is an empty body.
    ///
    /// # Errors
    ///
    /// [`BoolError::DegenerateInput`] when `expr` names a body past
    /// [`Partition::body_count`]; [`BoolError::NonManifoldResult`];
    /// the assembler's and acceptance gates' errors.
    pub fn select(&self, expr: &SetExpr, options: &BoolOptions) -> Result<BRep, BoolError> {
        if let Some(i) = expr.max_body()
            && i >= self.body_count
        {
            return Err(BoolError::DegenerateInput(format!(
                "the expression names body {i}, but the partition holds {}",
                self.body_count
            )));
        }
        let n = self.body_count;
        let regions = &self.built.complex;
        let keep = |m: &crate::complex::operand_set::OperandSet| {
            let set: Vec<usize> = m.iter().collect();
            expr.holds(&BodySet::of(&set), n)
        };
        let selection = crate::complex::select::select_fragments_by(
            &regions.lift,
            keep,
            self.built.fragments.len(),
        );
        let nothing_kept = !self.built.regions.iter().any(|r| expr.holds(&r.bodies, n));
        self.extract(&selection, options.non_manifold, nothing_kept)
    }

    /// Assemble `selection` into a validated body: a non-manifold
    /// selection declines with its contacts or, under `SplitToManifold`,
    /// is served as manifold pieces meeting at recorded contacts.
    /// `nothing_kept` says an empty body is the answer.
    pub(super) fn extract(
        &self,
        selection: &crate::complex::select::ComplexSelection,
        policy: NonManifoldPolicy,
        nothing_kept: bool,
    ) -> Result<BRep, BoolError> {
        let regions = &self.built.complex;
        let c = &regions.lift.complex;
        let mut plan = None;
        if !selection.nonmanifold_edges.is_empty() || !selection.nonmanifold_vertices.is_empty() {
            let split = (policy == NonManifoldPolicy::SplitToManifold).then(|| {
                crate::complex::split_plan::split_plan(
                    c,
                    &regions.sides,
                    &selection.kept_face,
                    &selection.selected_shell,
                )
            });
            match split {
                Some(Ok(p)) => plan = Some(p),
                _ => {
                    let mut contacts: Vec<NonManifoldContact> = selection
                        .nonmanifold_edges
                        .iter()
                        .map(|e| NonManifoldContact {
                            kind: ContactKind::Edge,
                            location: c.vertices[c.edges[e.index()].start.index()].position,
                        })
                        .collect();
                    contacts.extend(selection.nonmanifold_vertices.iter().map(|v| {
                        NonManifoldContact {
                            kind: ContactKind::Vertex,
                            location: c.vertices[v.index()].position,
                        }
                    }));
                    return Err(BoolError::NonManifoldResult { contacts });
                }
            }
        }
        let operands: Vec<&BRep> = self.bodies.iter().collect();
        let mut carried = crate::serve_region::Carried::default();
        let body = crate::cellular::assemble::cellular_assemble(
            &self.built.graph,
            &self.built.fragments,
            &selection.selections,
            &operands,
            self.tolerance,
            &mut carried,
            Default::default(),
            plan.as_ref(),
        )
        .map_err(|d| d.into_error())?;
        let certify_empty = || {
            if nothing_kept {
                Ok(())
            } else {
                Err(BoolError::MergeFailed {
                    stage: "Merge".to_string(),
                    reason: "partition: an empty result where the selection keeps a region"
                        .to_string(),
                })
            }
        };
        let body = crate::operation::serve::finish(
            body,
            (&carried, false),
            &self.tolerance,
            &certify_empty,
            None,
            &mut Vec::new(),
            &mut None,
        )?;
        if plan.is_some() {
            return crate::postprocess::contact_record::accept_split(body, &self.tolerance);
        }
        Ok(body)
    }
}
