//! The part of a sheet the regions it lies in select (PK-978 stage 2).
//!
//! A lamina contributes no membership bit, so both sides of each of its
//! complex faces lie in one membership class: the solids' region the
//! piece of sheet passes through. Sheet ∩ solid keeps the lamina faces
//! inside the solid, sheet − solid those outside, each once, in the
//! sheet's own sense. The assembler leaves a kept face's edge that no
//! other kept face shares as a half-edge twinned to itself (a cut), or
//! with its operand's null twin (the sheet's own open boundary carried
//! through); the result records both as the core's open boundary (a null twin, its edge in
//! `boundary_edges`), the form `validate` exempts from pairing and
//! junction checks and the partition reads back as a sheet. A sheet
//! face lying on another body's face (a coincident patch) is not
//! decided here: it declines, and the caller falls back.

use popkern_core::topology::BRep;

use super::{BodySet, Partition};
use crate::error::BoolError;
use crate::merge::fragment_classify::FragmentSelection;

fn declined(reason: String) -> BoolError {
    BoolError::MergeFailed {
        stage: "Merge".to_string(),
        reason: format!("partition: {reason}"),
    }
}

/// Make boundary half-edge `h` its edge's canonical use, its curve run
/// `h`'s way (the core's boundary convention, as face retirement keeps
/// it): an open edge by which end `h` leaves, a closed edge by the
/// direction `h`'s pcurve leaves the vertex in.
fn run_curve_along(body: &mut BRep, h: popkern_core::topology::HalfEdgeHandle) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let he = &body.half_edges[h];
    let Some(edge) = body.edges.get(he.edge) else {
        return;
    };
    let Some(curve) = edge.curve.as_ref() else {
        return;
    };
    let at = |t: f64| curve.geometry.point_at(t);
    let (a, b) = (at(curve.t_start), at(curve.t_end));
    let origin = body.vertices[he.origin].position;
    let reversed = if a.distance_to(&b) > body.tolerance.linear() {
        origin.distance_to(&b) < origin.distance_to(&a)
    } else {
        // Closed: compare the curve's first step with the half-edge's.
        let surface = body
            .loops
            .get(he.loop_)
            .and_then(|l| body.faces.get(l.face))
            .and_then(|f| f.surface.as_ref());
        let (Some(pcurve), Some(surface)) = (he.pcurve.as_ref(), surface) else {
            return;
        };
        let step = 1e-3;
        let image = |t: f64| {
            let q = pcurve.point_at(t);
            surface.geometry.point_at(q.x, q.y)
        };
        let along =
            image(pcurve.t_start + step * (pcurve.t_end - pcurve.t_start)) - image(pcurve.t_start);
        let ahead = at(curve.t_start + step * (curve.t_end - curve.t_start)) - a;
        along.dot(&ahead) < 0.0
    };
    let e = he.edge;
    body.edges[e].half_edge = h;
    if reversed && let Some(curve) = body.edges[e].curve.as_mut() {
        std::mem::swap(&mut curve.t_start, &mut curve.t_end);
    }
}

/// `body`, assembled from kept faces some of whose edges no other kept
/// face shares, in the core's open-boundary form (module doc), with its
/// pcurves completed and `validate_all` clean.
pub(super) fn open_sheet(mut body: BRep) -> Result<BRep, BoolError> {
    let open: Vec<_> = body
        .half_edges
        .iter()
        .filter(|(h, he)| {
            he.mate() == Some(*h) || !he.mate().is_some_and(|m| body.half_edges.contains_key(m))
        })
        .map(|(h, he)| (h, he.edge))
        .collect();
    for &(h, _) in &open {
        body.half_edges[h].clear_mate();
    }
    body.recompute_pcurves_where_missing()
        .map_err(|e| BoolError::ValidationFailed(format!("sheet pcurves: {e}")))?;
    for &(h, _) in &open {
        run_curve_along(&mut body, h);
    }
    // A sheet bounds no volume: its faces lie in the void, both sides.
    body.settle_open_shells();
    let findings = body.validate_all();
    if !findings.is_empty() {
        return Err(BoolError::ValidationFailed(format!(
            "partition sheet result: {findings:?}"
        )));
    }
    Ok(body)
}

impl Partition {
    /// The faces of sheet operand `sheet` lying in regions whose
    /// bodies `keep` accepts, as one sheet body (empty when none are).
    ///
    /// # Errors
    ///
    /// [`BoolError::DegenerateInput`] when `sheet` names no sheet
    /// operand; [`BoolError::MergeFailed`] when a sheet face lies on
    /// another body's face or its two sides disagree on membership;
    /// [`BoolError::ValidationFailed`] when the assembled sheet does
    /// not validate.
    pub(crate) fn lamina_part(
        &self,
        sheet: usize,
        keep: impl Fn(&BodySet) -> bool,
    ) -> Result<BRep, BoolError> {
        let lift = &self.built.complex.lift;
        let c = &lift.complex;
        if sheet >= self.body_count || !c.lamina_operands.contains(sheet) {
            return Err(BoolError::DegenerateInput(format!(
                "body {sheet} is not a sheet of the partition"
            )));
        }
        let membership = |shell: Option<crate::complex::ids::ShellId>| {
            shell
                .and_then(|s| c.shells[s.index()].membership.as_ref())
                .map(|m| BodySet::of(&m.iter().collect::<Vec<_>>()))
        };
        let mut selections = vec![
            FragmentSelection {
                keep: false,
                flip: false,
            };
            self.built.fragments.len()
        ];
        let mut kept = 0usize;
        for (fi, face) in c.faces.iter().enumerate() {
            let Some(primary) = face.owners.first() else {
                continue;
            };
            if primary.operand != sheet {
                if let Some(other) = face.owners.iter().find(|o| o.operand == sheet) {
                    return Err(declined(format!(
                        "sheet face {:?} lies on a face of body {}",
                        other.face, primary.operand
                    )));
                }
                continue;
            }
            let (ahead, behind) = (
                membership(face.sides[0].shell),
                membership(face.sides[1].shell),
            );
            let Some(bodies) = ahead.filter(|a| Some(*a) == behind) else {
                return Err(declined(format!(
                    "sheet face {fi} has sides in regions of different bodies ({ahead:?}, {behind:?})"
                )));
            };
            if !keep(&bodies) {
                continue;
            }
            if let Some(sel) = face
                .source_fragments
                .first()
                .and_then(|&f| selections.get_mut(f))
            {
                *sel = FragmentSelection {
                    keep: true,
                    flip: false,
                };
                kept += 1;
            }
        }
        if kept == 0 {
            return Ok(BRep::new(self.tolerance));
        }
        let operands: Vec<&BRep> = self.bodies.iter().collect();
        let mut carried = crate::serve_region::Carried::default();
        let body = crate::cellular::assemble::cellular_assemble(
            &self.built.graph,
            &self.built.fragments,
            &selections,
            &operands,
            self.tolerance,
            &mut carried,
            Default::default(),
            None,
        )
        .map_err(|d| d.into_error())?;
        open_sheet(body)
    }
}

#[cfg(test)]
#[path = "lamina_part_tests.rs"]
mod tests;
