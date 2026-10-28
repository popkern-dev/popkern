//! Classifying a boolean's fragments through the cell complex.
//!
//! lift → radial order → shells → membership → selection, with the
//! structural validator between stages. Every decline is typed and
//! carries the complex artifact, so a failure can be inspected
//! without rerunning the boolean.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::ids::ShellId;
use super::lift::{Lift, LiftDecline, lift_fragments};
use super::model::{CellComplex, SideIndex};
use super::radial::{FaceGeometry, RadialDecline, order_radially};
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
/// inspected without rerunning the boolean.
fn emit_decline_artifact(c: &CellComplex, decline: &ComplexDecline) {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static SEQUENCE: AtomicUsize = AtomicUsize::new(0);
    let Some(dir) = std::env::var_os("POPKERN_COMPLEX_ARTIFACT") else {
        return;
    };
    let n = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let path =
        std::path::Path::new(&dir).join(format!("complex-decline-{}-{n}.json", std::process::id()));
    let body = format!(
        "{{\"decline\":{:?},\"complex\":{}}}",
        format!("{decline:?}"),
        super::artifact::complex_artifact_json(c)
    );
    if std::fs::write(&path, body).is_ok() {
        eprintln!("[complex] decline artifact: {}", path.display());
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
        let fin = c.loops.get(f.loops.first()?.index())?.first;
        let edge = &c.edges[c.fins[fin.index()].edge.index()];
        let (pa, pb) = (
            c.vertices[edge.start.index()].position,
            c.vertices[edge.end.index()].position,
        );
        // A point on the face a short way in from its first fin.
        let mid = match &edge.curve {
            Some(cv) => {
                use popkern_core::geometry::curve::CurveEval;
                cv.geometry.point_at(0.5 * (cv.t_start + cv.t_end))
            }
            None => pa + (pb - pa) * 0.5,
        };
        let (u, v) = surf.geometry.closest_point(&mid)?;
        let n = surf.folded_normal(u, v)?;
        let t = (pb - pa).normalized().or_else(|| {
            edge.curve.as_ref().and_then(|cv| {
                use popkern_core::geometry::curve::CurveEval;
                cv.geometry
                    .derivative_at(0.5 * (cv.t_start + cv.t_end), 1)
                    .normalized()
            })
        })?;
        let tf = if c.fins[fin.index()].forward {
            t
        } else {
            t * -1.0
        };
        let len = pa.distance_to(&pb).max(linear * 1e3);
        let into = n.cross(&tf).normalized()?;
        let on_face = mid + into * (len * 1e-2).max(linear * 100.0);
        let (u2, v2) = surf.geometry.closest_point(&on_face)?;
        let base = surf.geometry.point_at(u2, v2);
        let n2 = surf.folded_normal(u2, v2)?;
        let outward = if f.folded_is_outward { n2 } else { n2 * -1.0 };
        let toward = if side == SideIndex::Ahead {
            outward
        } else {
            outward * -1.0
        };
        let probe = base + toward * (len * 1e-3).max(linear * 50.0);
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
        let half = super::record_use::half_registered_records(&lift.complex, input.graph);
        if !half.is_empty() {
            return Err(ComplexDecline::HalfRegisteredRecords(half));
        }
        let geo = FaceGeometry {
            operands: input.operands,
        };
        let sides = order_radially(&mut lift.complex, &geo, input.tolerance.angular(), linear)
            .map_err(ComplexDecline::Radial)?;
        assemble_shells(&mut lift.complex, &sides).map_err(ComplexDecline::Region)?;
        let classifications = assign_membership(&mut lift.complex, |c, shell, k| {
            classify_shell(c, input, shell, k)
        })
        .map_err(ComplexDecline::Region)?;
        checked(&lift.complex)?;
        Ok(classifications)
    })();
    let classifications = match staged {
        Ok(n) => n,
        Err(decline) => {
            emit_decline_artifact(&lift.complex, &decline);
            return Err(decline);
        }
    };
    let selection = select_fragments(&lift, op, input.fragments.len());
    Ok(ComplexVerdict {
        lift,
        selection,
        classifications,
    })
}
