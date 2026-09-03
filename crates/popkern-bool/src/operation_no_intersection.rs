//! No-intersection dispatch for boolean operations — split from
//! `operation_helpers.rs` (500-line cap): spatial-relationship
//! classification (disjoint / containment / coincident) when the
//! intersection stage finds no curves, and the per-relationship
//! result arms.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::classify::point_in_solid;
use crate::error::{BoolError, BoolOperation, PointClassification};
use crate::operation::helpers::{combine_breps, copy_brep_into};

pub(super) fn handle_no_intersection(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    aabbs_a: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    aabbs_b: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    // Classify sample points of A against B. A SINGLE vertex probe is
    // a knife edge: one degenerate ray (a sphere pole, a seam vertex)
    // answers OnBoundary for a clearly-separated body and the whole
    // no-intersection case was declared CoincidentFaces with bogus
    // face handles (measured: a bearing's merged races+balls union
    // with the NEXT ball — every operand pair disjoint by 0.01+).
    // Probe several stations (vertices, then face-interior samples)
    // and take the first decisive Inside/Outside; only a unanimously
    // boundary-touching body is the ambiguous coincident class.
    let classify_probes = |inner: &BRep,
                           outer: &BRep,
                           outer_aabbs: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )]|
     -> Result<PointClassification, BoolError> {
        let mut probes: Vec<popkern_core::geometry::Point3> = inner
            .vertices
            .values()
            .take(4)
            .map(|v| v.position)
            .collect();
        for face in inner.faces.keys().take(4) {
            probes.extend(inner.face_outer_loop_sample_positions(face, 3));
        }
        if probes.is_empty() {
            return Err(BoolError::DegenerateInput(
                "solid has no vertices".to_string(),
            ));
        }
        let mut verdict = PointClassification::OnBoundary;
        for p in &probes {
            match point_in_solid(p, outer, outer_aabbs, tolerance)? {
                PointClassification::OnBoundary => {}
                decisive => {
                    verdict = decisive;
                    break;
                }
            }
        }
        Ok(verdict)
    };
    let classification = classify_probes(&a, &b, aabbs_b)?;

    let result = match classification {
        PointClassification::Outside => {
            // A's samples are outside B. Check reverse: is B inside A?
            let reverse = classify_probes(&b, &a, aabbs_a)?;
            match reverse {
                PointClassification::Inside => {
                    // B is inside A — reverse containment.
                    handle_reverse_containment(op, a, b, tolerance)
                }
                _ => {
                    // Truly disjoint.
                    handle_disjoint(op, a, b, tolerance)
                }
            }
        }
        PointClassification::Inside => {
            // Containment: A is entirely inside B.
            handle_containment(op, a, b, tolerance)
        }
        PointClassification::OnBoundary => {
            // Coincident boundary — ambiguous case.
            Err(BoolError::CoincidentFaces {
                face_a: a.faces.keys().next().unwrap_or_default(),
                face_b: b.faces.keys().next().unwrap_or_default(),
            })
        }
    };

    // Log validation issues on the no-intersection result.
    if let Ok(ref brep) = result {
        super::log_validation_issues(brep, op);
    }

    result
}

/// Handle disjoint solids (A is entirely outside B) with no intersections.
fn handle_disjoint(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    match op {
        BoolOperation::Union => {
            // Combine both BReps into one.
            combine_breps(a, b)
        }
        BoolOperation::Intersection => {
            // Empty result.
            Ok(BRep::new(*tolerance))
        }
        BoolOperation::Subtraction => {
            // A minus nothing = A.
            Ok(a)
        }
    }
}

/// Handle containment (A is entirely inside B) with no intersections.
fn handle_containment(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    match op {
        BoolOperation::Union => {
            // B fully contains A, so union is just B.
            Ok(b)
        }
        BoolOperation::Intersection => {
            // A is inside B, so intersection is A.
            Ok(a)
        }
        BoolOperation::Subtraction => {
            // A − B when A is inside B: A is entirely consumed by B, result is empty.
            Ok(BRep::new(*tolerance))
        }
    }
}

/// Handle reverse containment (B is entirely inside A) with no intersections.
fn handle_reverse_containment(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    match op {
        BoolOperation::Union => {
            // A fully contains B, so union is just A.
            Ok(a)
        }
        BoolOperation::Intersection => {
            // B is inside A, so intersection is B.
            Ok(b)
        }
        BoolOperation::Subtraction => {
            // A − B when B is inside A: hollow solid (A with B-shaped void).
            build_hollow_solid(a, b, tolerance)
        }
    }
}

/// Build a hollow solid: A with a B-shaped void inside.
///
/// Copies all of B's topology into A with flipped surface normals, creating
/// an inner void. The result has A's original outer faces plus B's faces
/// with inverted orientation.
fn build_hollow_solid(a: BRep, b: BRep, _tolerance: &Tolerance) -> Result<BRep, BoolError> {
    copy_brep_into(a, b, true)
}
