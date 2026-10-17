//! No-intersection dispatch for boolean operations — split from
//! `operation_helpers.rs` (500-line cap): spatial-relationship
//! classification (disjoint / containment / coincident) when the
//! intersection stage finds no curves, and the per-relationship
//! result arms.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::classify::point_in_solid_final;
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
            // Boolean operands are FINAL bodies: a ray hit inside one of
            // a face's holes is a void and must not count. Counting it
            // classified an output pin passing through a cycloid disc's
            // window as inside the disc, and the intersection served
            // the whole pin (PK-722).
            match point_in_solid_final(p, outer, outer_aabbs, tolerance)? {
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
/// Copies B's skin into A with its material side mirrored, and joins its
/// shells to the A solid enclosing them as cavity shells — the model the
/// cellular shell partition serves (a cavity shell belongs to its
/// encloser's solid). B's own solids are dropped: copied as solids, the
/// void read as a second, positive body (PK-879).
fn build_hollow_solid(a: BRep, b: BRep, _tolerance: &Tolerance) -> Result<BRep, BoolError> {
    let host = enclosing_solid(&a, &b);
    let a_solids: Vec<_> = a.solids.keys().collect();
    let mut result = copy_brep_into(a, b, true)?;
    let Some(host) = host else {
        return Ok(result);
    };
    let void_solids: Vec<_> = result
        .solids
        .keys()
        .filter(|s| !a_solids.contains(s))
        .collect();
    for shell in result.shells.values_mut() {
        if void_solids.contains(&shell.solid) {
            shell.solid = host;
        }
    }
    for s in void_solids {
        result.solids.remove(s);
    }
    Ok(result)
}

/// The A solid whose face boxes' union holds B's — the smallest when
/// several do (nested bodies); `None` when A has no solid. The caller
/// established that B lies inside A, so a box test picks among A's
/// disjoint solids. Face boxes, not vertex boxes: a revolved body's
/// vertices sit on its seam alone.
fn enclosing_solid(a: &BRep, b: &BRep) -> Option<popkern_core::topology::SolidHandle> {
    use popkern_core::geometry::intersection::surface_subdivision::Aabb;
    let union = |boxes: &mut dyn Iterator<Item = Aabb>| {
        boxes.reduce(|u, x| Aabb {
            min: popkern_core::geometry::Point3::new(
                u.min.x.min(x.min.x),
                u.min.y.min(x.min.y),
                u.min.z.min(x.min.z),
            ),
            max: popkern_core::geometry::Point3::new(
                u.max.x.max(x.max.x),
                u.max.y.max(x.max.y),
                u.max.z.max(x.max.z),
            ),
        })
    };
    let a_boxes = crate::intersect::compute_face_aabbs(a);
    let tool = union(
        &mut crate::intersect::compute_face_aabbs(b)
            .into_iter()
            .map(|(_, x)| x),
    )?;
    let mut best: Option<(f64, popkern_core::topology::SolidHandle)> = None;
    for (sh, _) in &a.solids {
        let Some(u) = union(&mut a_boxes.iter().filter_map(|(f, x)| {
            let face = a.faces.get(*f)?;
            (a.shells.get(face.shell)?.solid == sh).then_some(x.clone())
        })) else {
            continue;
        };
        let holds = u.min.x <= tool.min.x
            && u.min.y <= tool.min.y
            && u.min.z <= tool.min.z
            && u.max.x >= tool.max.x
            && u.max.y >= tool.max.y
            && u.max.z >= tool.max.z;
        let size = (u.max.x - u.min.x) * (u.max.y - u.min.y) * (u.max.z - u.min.z);
        if holds && best.is_none_or(|(s, _)| size < s) {
            best = Some((size, sh));
        }
    }
    best.map(|(_, s)| s).or_else(|| a.solids.keys().next())
}
