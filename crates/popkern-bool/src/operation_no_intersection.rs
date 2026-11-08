//! No-intersection dispatch for boolean operations — split from
//! `operation_helpers.rs` (500-line cap): spatial-relationship
//! classification (disjoint / containment / coincident) when the
//! intersection stage finds no curves, and the per-relationship
//! result arms.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::{BoolError, BoolOperation, PointClassification};
use crate::operation::helpers::{combine_breps, copy_brep_into};
use crate::operation::relation_probes;

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
    // Classify probe points of each body against the other. A SINGLE
    // vertex probe is a knife edge: one degenerate ray (a sphere pole,
    // a seam vertex) answers OnBoundary for a clearly-separated body
    // and the whole no-intersection case was declared CoincidentFaces
    // (measured: a bearing's merged races+balls union with the NEXT
    // ball). Probes are face-interior points first, the decisive ones
    // must agree, and an undecided relation declines (PK-1124): an
    // OnBoundary reverse read as "apart" served an empty intersection
    // of a cube and its inscribed ball, whose poles and seam touch the
    // cube's faces.
    let forward = relation_probes::probe_relation(&a, &b, aabbs_b, tolerance)?;
    let result = match forward.verdict() {
        None => Err(forward.undecided(0)),
        Some(PointClassification::Inside) => {
            // Containment: A is entirely inside B.
            handle_containment(op, a, b, tolerance)
        }
        Some(PointClassification::OnBoundary) => {
            // Every probe of A on B's boundary: coincident boundaries.
            Err(BoolError::CoincidentFaces {
                face_a: a.faces.keys().next().unwrap_or_default(),
                face_b: b.faces.keys().next().unwrap_or_default(),
            })
        }
        Some(PointClassification::Outside) => {
            // A's probes are outside B. Is B inside A, or are they apart?
            let reverse = relation_probes::probe_relation(&b, &a, aabbs_a, tolerance)?;
            match reverse.verdict() {
                Some(PointClassification::Inside) => {
                    handle_reverse_containment(op, a, b, tolerance)
                }
                Some(PointClassification::Outside) => handle_disjoint(op, a, b, tolerance),
                // B's probes all on A's boundary, or disagreeing.
                Some(PointClassification::OnBoundary) | None => Err(reverse.undecided(1)),
            }
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
/// Copies B's skin into A with its material side mirrored, and turns B's
/// regions inside out into the A solid enclosing them
/// ([`BRep::invert_regions_into`]): B's skin becomes that solid's cavity
/// shell and B's material a void region — the model the cellular shell
/// partition serves (a cavity shell belongs to its encloser's solid).
/// Copied as solids, the void read as a second, positive body (PK-879).
fn build_hollow_solid(a: BRep, b: BRep, _tolerance: &Tolerance) -> Result<BRep, BoolError> {
    let host = enclosing_solid(&a, &b);
    let a_regions: Vec<_> = a.regions.keys().collect();
    let mut result = copy_brep_into(a, b, true)?;
    let Some(host) = host else {
        return Ok(result);
    };
    // B's bounded regions, its material and any cavity of it, turn
    // inside out into the host: B's material is the host's cavity.
    let b_regions: Vec<_> = result
        .regions
        .keys()
        .filter(|r| !a_regions.contains(r))
        .collect();
    result.invert_regions_into(&b_regions, host);
    Ok(result)
}

/// The A solid whose face boxes' union holds B's — the smallest when
/// several do (nested bodies); `None` when A has no solid. The caller
/// established that B lies inside A, so a box test picks among A's
/// disjoint solids. Face boxes, not vertex boxes: a revolved body's
/// vertices sit on its seam alone.
fn enclosing_solid(a: &BRep, b: &BRep) -> Option<popkern_core::topology::RegionHandle> {
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
    let mut best: Option<(f64, popkern_core::topology::RegionHandle)> = None;
    for sh in a.solid_regions() {
        let Some(u) = union(&mut a_boxes.iter().filter_map(|(f, x)| {
            (a.shell_region(a.face_back_shell(*f)?)? == sh).then_some(x.clone())
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
    best.map(|(_, s)| s)
        .or_else(|| a.solid_regions().first().copied())
}
