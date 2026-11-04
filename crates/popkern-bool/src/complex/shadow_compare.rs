//! Comparing the complex verdict with the cellular selections of the
//! same boolean (shadow mode; see `crate::complex_shadow`).

use super::pipeline::{ComplexInput, classify_with_complex};
use crate::complex_shadow::{
    ComplexShadowEvent, ComplexShadowOutcome, ShadowDifference, record_complex_shadow,
};
use crate::error::BoolOperation;
use crate::merge::fragment_classify::FragmentSelection;

/// Classify through the complex and record how it compares with the
/// cellular `selections`. Never changes what is served.
pub(crate) fn shadow_compare(
    input: &ComplexInput<'_>,
    selections: &[FragmentSelection],
    op: BoolOperation,
) {
    let fragments = input.fragments.len();
    let event = match classify_with_complex(input, op) {
        Ok(verdict) => {
            let ours = &verdict.selection.selections;
            let differing: Vec<usize> = ours
                .iter()
                .zip(selections)
                .enumerate()
                .filter(|(_, (c, s))| c.keep != s.keep || (c.keep && c.flip != s.flip))
                .map(|(i, _)| i)
                .collect();
            if popkern_core::trace::trace_enabled("complex") {
                trace_differing(input, &verdict, selections, &differing);
            }
            let faces = &verdict.lift.complex.faces;
            let differing: Vec<(usize, ShadowDifference)> = differing
                .into_iter()
                .map(|i| {
                    let sources = faces
                        .iter()
                        .find(|f| f.source_fragments.contains(&i))
                        .map(|f| f.source_fragments.as_slice());
                    (i, difference(i, sources, ours, selections))
                })
                .collect();
            ComplexShadowEvent {
                operation: op,
                outcome: if differing.is_empty() {
                    ComplexShadowOutcome::Agree
                } else {
                    ComplexShadowOutcome::Disagree {
                        fragments: differing,
                    }
                },
                fragments,
                classifications: verdict.classifications,
                nonmanifold: (
                    verdict.selection.nonmanifold_edges.len(),
                    verdict.selection.nonmanifold_vertices.len(),
                ),
            }
        }
        Err(decline) => ComplexShadowEvent {
            operation: op,
            outcome: ComplexShadowOutcome::Declined(format!("{decline:?}")),
            fragments,
            classifications: 0,
            nonmanifold: (0, 0),
        },
    };
    record_complex_shadow(event);
}

/// How the verdicts on differing fragment `i` differ, given the source
/// fragments of the complex face lifting it (`None`: no face does).
fn difference(
    i: usize,
    sources: Option<&[usize]>,
    ours: &[FragmentSelection],
    theirs: &[FragmentSelection],
) -> ShadowDifference {
    let kept = |sel: &[FragmentSelection], k: usize| sel.get(k).is_some_and(|s| s.keep);
    if kept(ours, i) && kept(theirs, i) {
        return ShadowDifference::Orientation;
    }
    match sources {
        None => ShadowDifference::HoleFiller,
        // One copy of the patch kept by each, or by neither.
        Some(copies)
            if copies.len() > 1
                && copies.iter().any(|&k| kept(ours, k))
                    == copies.iter().any(|&k| kept(theirs, k)) =>
        {
            ShadowDifference::CoincidentCopy
        }
        Some(_) => ShadowDifference::Membership,
    }
}

/// One line per differing fragment: its face, both verdicts, and the
/// membership of the complex face's two sides.
fn trace_differing(
    input: &ComplexInput<'_>,
    verdict: &super::pipeline::ComplexVerdict,
    selections: &[FragmentSelection],
    differing: &[usize],
) {
    let c = &verdict.lift.complex;
    for &i in differing {
        let face = c.faces.iter().position(|f| f.source_fragments.contains(&i));
        let sides = face.map(|fi| {
            let f = &c.faces[fi];
            let member = |s: Option<super::ids::ShellId>| {
                s.and_then(|s| c.shells[s.index()].membership.as_ref())
                    .map(|m| (m.contains(0), m.contains(1)))
            };
            (
                fi,
                f.owners.len(),
                f.folded_is_outward,
                verdict.lift.output_folded_outward,
                member(f.sides[0].shell),
                member(f.sides[1].shell),
            )
        });
        let ours = &verdict.selection.selections[i];
        let theirs = &selections[i];
        eprintln!(
            "[complex] differs: fragment {i} {:?} cellular keep={} flip={} complex keep={} flip={} face/owners/folded_out/output_out/ahead(a,b)/behind(a,b) {sides:?}",
            input.fragments[i].face, theirs.keep, theirs.flip, ours.keep, ours.flip
        );
    }
}

#[cfg(test)]
#[path = "shadow_compare_tests.rs"]
mod tests;
