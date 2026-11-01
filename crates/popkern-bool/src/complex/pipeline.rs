//! Classifying a boolean's fragments through the cell complex.
//!
//! lift → radial order → shells → membership → selection, with the
//! structural validator between stages. Every decline is typed and
//! carries the complex artifact, so a failure can be inspected
//! without rerunning the boolean.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::ids::ShellId;
use super::lift::{Lift, LiftDecline, lift_fragments};
use super::model::{CellComplex, SideIndex};
use super::radial::{FaceGeometry, RadialDecline, RadialSides, order_radially};
use super::regions::{RegionDecline, assemble_shells, assign_membership};
use super::select::{ComplexSelection, select_fragments};
use super::validate::{ComplexFinding, validate_complex};
use crate::error::{BoolOperation, PointClassification};
use crate::merge::graph::IntersectionGraph;
use crate::merge::walk::Fragment;

/// Why classification through the complex declined.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum ComplexDecline {
    /// The lift declined.
    Lift(LiftDecline),
    /// The structural validator found violations.
    Invalid(Vec<ComplexFinding>),
    /// Radial ordering declined.
    Radial(RadialDecline),
    /// Shell or membership assembly declined.
    Region(RegionDecline),
    /// Fragments lying in a coincident patch that grouped with no
    /// partner copy (see `same_domain_check`).
    SameDomainMismatch(Vec<super::same_domain_check::SameDomainMismatch>),
    /// Shared records that split only one of their two faces (the walk
    /// is inconsistent there; see `record_use`).
    HalfRegisteredRecords(Vec<super::ids::EdgeId>),
}

/// A complex verdict: the lifted complex with shells and membership,
/// the fragment selections, and how many point classifications
/// membership needed.
#[derive(Debug, Clone)]
pub(crate) struct ComplexVerdict {
    /// The complex and its fragment correspondence.
    pub lift: Lift,
    /// Selections and non-manifold findings.
    pub selection: ComplexSelection,
    /// Point classifications run (one per undetermined shell class).
    pub classifications: usize,
    /// Which side of each fin's face faces the next wedge about its
    /// edge (the radial stage's output; the split plan pairs by it).
    pub sides: RadialSides,
}

/// Inputs to complex classification.
pub(crate) struct ComplexInput<'a> {
    /// The intersection graph.
    pub graph: &'a IntersectionGraph,
    /// The chart walk's fragments.
    pub fragments: &'a [Fragment],
    /// Operands A and B.
    pub operands: [&'a BRep; 2],
    /// Face boxes of A and B (for point classification).
    pub boxes: [&'a [(FaceHandle, Aabb)]; 2],
    /// Coincident face pairs (A face, B face).
    pub overlaps: &'a [(FaceHandle, FaceHandle)],
    /// Model tolerance.
    pub tolerance: &'a Tolerance,
    /// Graph edge keys on a scoped operand's frontier (between an
    /// active face and a carried far face); empty for unscoped runs.
    pub frontier: &'a std::collections::BTreeSet<usize>,
}

fn checked(c: &CellComplex) -> Result<(), ComplexDecline> {
    let findings = validate_complex(c);
    if findings.is_empty() {
        Ok(())
    } else {
        Err(ComplexDecline::Invalid(findings))
    }
}

/// With `POPKERN_COMPLEX_ARTIFACT=<dir>`, write the declined complex
/// and its decline to `<dir>` as one JSON file, so the failure can be
/// inspected without rerunning the boolean, and both operands beside it
/// (`….operand0.json`, `….operand1.json`, serde BReps), so it can be
/// replayed alone: a decline inside a script, an assembly or a worker
/// pool reproduces without the code that built its operands.
fn emit_decline_artifact(
    c: &CellComplex,
    decline: &ComplexDecline,
    op: BoolOperation,
    operands: [&BRep; 2],
) {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static SEQUENCE: AtomicUsize = AtomicUsize::new(0);
    let Some(dir) = std::env::var_os("POPKERN_COMPLEX_ARTIFACT") else {
        return;
    };
    let n = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let path =
        std::path::Path::new(&dir).join(format!("complex-decline-{}-{n}.json", std::process::id()));
    let body = format!(
        "{{\"decline\":{:?},\"operation\":\"{op:?}\",\"complex\":{}}}",
        format!("{decline:?}"),
        super::artifact::complex_artifact_json(c)
    );
    if std::fs::write(&path, body).is_ok() {
        eprintln!("[complex] decline artifact: {}", path.display());
    }
    for (k, operand) in operands.iter().enumerate() {
        if let Ok(json) = serde_json::to_vec(operand) {
            let _ = std::fs::write(path.with_extension(format!("operand{k}.json")), json);
        }
    }
}

/// A point just off a face side of `shell`, classified against
/// `operand`: one ray for a whole class of shells no face of the
/// operand touches.
fn classify_shell(
    c: &CellComplex,
    input: &ComplexInput<'_>,
    shell: ShellId,
    operand: usize,
) -> Option<bool> {
    let linear = input.tolerance.linear();
    for &(face, side) in c.shells.get(shell.index())?.sides.iter().take(8) {
        let f = &c.faces[face.index()];
        let owner = f.owners.first()?;
        let surf = input.operands[owner.operand]
            .faces
            .get(owner.face)?
            .surface
            .as_ref()?;
        let at = super::same_domain_check::inner_point(c, input.operands, linear, face)?;
        let (base, (u2, v2), len) = (at.point, at.uv, at.fin_length);
        // Membership asks only about operands that bound none of the
        // class's faces, and a face of another operand has the same
        // membership on both sides: classify the point on the face. An
        // offset probe crosses any gap thinner than the offset
        // (measured: a pin 0.002 clear of its bore, probed 0.03 out,
        // read inside the plate and the intersection kept the pin).
        let probe = if owner.operand != operand {
            base
        } else {
            let n2 = surf.folded_normal(u2, v2)?;
            let outward = if f.folded_is_outward { n2 } else { n2 * -1.0 };
            let toward = if side == SideIndex::Ahead {
                outward
            } else {
                outward * -1.0
            };
            base + toward * (len * 1e-3).max(linear * 50.0)
        };
        match crate::classify::point_in_solid_final(
            &probe,
            input.operands[operand],
            input.boxes[operand],
            input.tolerance,
        ) {
            Ok(PointClassification::Inside) => return Some(true),
            Ok(PointClassification::Outside) => return Some(false),
            _ => continue,
        }
    }
    None
}

/// Classify the fragments of one boolean through the complex.
pub(crate) fn classify_with_complex(
    input: &ComplexInput<'_>,
    op: BoolOperation,
) -> Result<ComplexVerdict, ComplexDecline> {
    let linear = input.tolerance.linear();
    let mut lift = lift_fragments(
        input.graph,
        input.fragments,
        input.operands,
        input.boxes,
        input.overlaps,
        input.tolerance,
    )
    .map_err(ComplexDecline::Lift)?;
    for edge in &mut lift.complex.edges {
        edge.frontier = edge
            .source_key
            .is_some_and(|key| input.frontier.contains(&key));
    }
    let staged = (|| {
        checked(&lift.complex)?;
        let half = super::record_use::half_registered_records(
            &lift.complex,
            input.graph,
            input.operands,
            linear,
        );
        if !half.is_empty() {
            return Err(ComplexDecline::HalfRegisteredRecords(half));
        }
        let unmatched = super::same_domain_check::same_domain_mismatches(
            &lift.complex,
            input.operands,
            input.overlaps,
            linear,
        );
        if !unmatched.is_empty() {
            return Err(ComplexDecline::SameDomainMismatch(unmatched));
        }
        let geo = FaceGeometry {
            operands: input.operands,
            graph: Some(input.graph),
        };
        let sides = order_radially(&mut lift.complex, &geo, input.tolerance.angular(), linear)
            .map_err(ComplexDecline::Radial)?;
        assemble_shells(&mut lift.complex, &sides).map_err(ComplexDecline::Region)?;
        let classifications = assign_membership(&mut lift.complex, |c, shell, k| {
            classify_shell(c, input, shell, k)
        })
        .map_err(ComplexDecline::Region)?;
        checked(&lift.complex)?;
        Ok((classifications, sides))
    })();
    let (classifications, sides) = match staged {
        Ok(staged) => staged,
        Err(decline) => {
            emit_decline_artifact(&lift.complex, &decline, op, input.operands);
            return Err(decline);
        }
    };
    let selection = select_fragments(&lift, op, input.fragments.len());
    Ok(ComplexVerdict {
        lift,
        selection,
        classifications,
        sides,
    })
}

#[cfg(test)]
#[path = "pipeline_tests.rs"]
mod tests;
