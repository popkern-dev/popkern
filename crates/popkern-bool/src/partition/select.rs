//! Bodies from a partition's regions (design 6).

use popkern_core::topology::BRep;

use super::{BodySet, Partition, SetExpr};
use crate::error::{BoolError, ContactKind, NonManifoldContact};
use crate::options::{BoolOptions, NonManifoldPolicy};

impl Partition {
    /// The regularized union of the regions `expr` keeps: a face between
    /// two kept regions is internal and merged away. A non-manifold
    /// result follows `options.non_manifold` as a two-operand boolean
    /// does (declined with its contacts by default, one general body
    /// under `General`, separate bodies under `SplitToManifold`, which
    /// [`Partition::select_bodies`] returns). No region kept is an empty
    /// body.
    ///
    /// # Errors
    ///
    /// [`BoolError::DegenerateInput`] when `expr` names a body past
    /// [`Partition::body_count`]; [`BoolError::NonManifoldResult`];
    /// [`BoolError::SplitIntoBodies`]; [`BoolError::GeneralBodyDeclined`];
    /// the assembler's and acceptance gates' errors.
    pub fn select(&self, expr: &SetExpr, options: &BoolOptions) -> Result<BRep, BoolError> {
        let pieces = self.select_pieces(expr, options);
        crate::split_serve::serve_one(pieces, options.non_manifold, &self.tolerance)
    }

    /// [`Partition::select`] as a list of bodies: one body, except under
    /// `SplitToManifold` a selection whose manifold pieces touch, each
    /// piece its own valid body.
    ///
    /// # Errors
    ///
    /// As [`Partition::select`], never [`BoolError::SplitIntoBodies`].
    pub fn select_bodies(
        &self,
        expr: &SetExpr,
        options: &BoolOptions,
    ) -> Result<Vec<BRep>, BoolError> {
        let pieces = self.select_pieces(expr, options);
        crate::split_serve::serve_bodies(pieces, options.non_manifold, &self.tolerance)
    }

    /// [`Partition::select`] as the pipeline answers it: a non-manifold
    /// selection as manifold pieces in one body under the split and
    /// general policies.
    pub(crate) fn select_pieces(
        &self,
        expr: &SetExpr,
        options: &BoolOptions,
    ) -> Result<BRep, BoolError> {
        if let Some(i) = expr.max_body()
            && i >= self.body_count
        {
            return Err(BoolError::DegenerateInput(format!(
                "the expression names body {i}, but the partition holds {}",
                self.body_count
            )));
        }
        let n = self.body_count;
        self.select_where_pieces(|m| expr.holds(m, n), options)
    }

    /// [`Partition::select`] by a predicate on each region's bodies (the
    /// code-only form; a feature tree records a [`SetExpr`]).
    ///
    /// # Errors
    ///
    /// As [`Partition::select`].
    pub fn select_where(
        &self,
        keep: impl Fn(&BodySet) -> bool,
        options: &BoolOptions,
    ) -> Result<BRep, BoolError> {
        let pieces = self.select_where_pieces(keep, options);
        crate::split_serve::serve_one(pieces, options.non_manifold, &self.tolerance)
    }

    /// [`Partition::select_where`] as a list of bodies (as
    /// [`Partition::select_bodies`]).
    ///
    /// # Errors
    ///
    /// As [`Partition::select`], never [`BoolError::SplitIntoBodies`].
    pub fn select_where_bodies(
        &self,
        keep: impl Fn(&BodySet) -> bool,
        options: &BoolOptions,
    ) -> Result<Vec<BRep>, BoolError> {
        let pieces = self.select_where_pieces(keep, options);
        crate::split_serve::serve_bodies(pieces, options.non_manifold, &self.tolerance)
    }

    /// [`Partition::select_where`] as the pipeline answers it (see
    /// [`Partition::select_pieces`]).
    pub(crate) fn select_where_pieces(
        &self,
        keep: impl Fn(&BodySet) -> bool,
        options: &BoolOptions,
    ) -> Result<BRep, BoolError> {
        let regions = &self.built.complex;
        let accepts = |m: &crate::complex::operand_set::OperandSet| {
            let set: Vec<usize> = m.iter().collect();
            keep(&BodySet::of(&set))
        };
        let selection = crate::complex::select::select_fragments_by(
            &regions.lift,
            accepts,
            self.built.fragments.len(),
        );
        let nothing_kept = !self.built.regions.iter().any(|r| keep(&r.bodies));
        self.extract(&selection, options.non_manifold, nothing_kept)
    }

    /// Assemble `selection` into a body: a non-manifold selection
    /// declines with its contacts or, under `SplitToManifold` (and
    /// `General`, which runs it), is served as manifold pieces in one
    /// body, each valid alone (the entry points serve them on).
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
            let split = policy.serves().then(|| {
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
            return crate::postprocess::split_pieces::accept_split(body, &self.tolerance, policy);
        }
        Ok(body)
    }
}

impl Partition {
    /// Every region inside body `body` as one body with the faces as the
    /// complex splits them: no output merge of co-surfacial faces (an
    /// imprint's whole point), accepted by `validate` with its pcurves
    /// completed.
    pub(crate) fn imprinted(&self, body: usize) -> Result<BRep, BoolError> {
        let regions = &self.built.complex;
        let selection = crate::complex::select::select_fragments_by(
            &regions.lift,
            |m| m.contains(body),
            self.built.fragments.len(),
        );
        let operands: Vec<&BRep> = self.bodies.iter().collect();
        let mut carried = crate::serve_region::Carried::default();
        let mut out = crate::cellular::assemble::cellular_assemble(
            &self.built.graph,
            &self.built.fragments,
            &selection.selections,
            &operands,
            self.tolerance,
            &mut carried,
            Default::default(),
            None,
        )
        .map_err(|d| d.into_error())?;
        out.validate()
            .map_err(|e| BoolError::ValidationFailed(format!("imprint: {e:?}")))?;
        let _ = out.recompute_pcurves_where_missing();
        if let Some(face) = crate::operation::first_missing_pcurve_face(&out) {
            return Err(BoolError::ValidationFailed(format!(
                "imprint: face {face:?} has no complete trim"
            )));
        }
        Ok(out)
    }
}
