//! Certificates for an EMPTY boolean result (PK-465, PK-966).
//!
//! An empty result is the right answer only in configurations the
//! operands themselves can show: A inside B for a subtraction, disjoint
//! interiors for an intersection, never for a union. A selection that
//! keeps nothing cannot certify itself — a walk that classifies every
//! fragment away looks exactly like disjointness — so the certificate
//! probes the ORIGINAL operands directly.

use popkern_core::topology::BRep;

/// Verify a claimed whole-body containment directly against the
/// ORIGINAL operands: every vertex, sampled loop boundary point and
/// face interior point of `inner` must lie on `outer`'s boundary or
/// classify Inside it. The empty-selection shortcut's claim
/// comes from seed-propagated face classifications, and a silently
/// under-split operand poisons every seed (measured: a rim-straddling
/// key subtract left the disc's top/bottom faces unsplit, one seed
/// wiped each side, and the "contained" flange came back as an EMPTY
/// body that passed every downstream gate vacuously — PK-173).
fn verify_whole_body_containment(inner: &BRep, outer: &BRep) -> bool {
    use crate::classify::point_in_solid_final;
    use crate::error::PointClassification;
    use crate::intersect::compute_face_aabbs;
    let aabbs = compute_face_aabbs(outer);
    let tolerance = outer.tolerance;
    // Vertices alone under-cover analytic bodies (a plain cylinder
    // has only its two seam vertices — both can sit inside a partial
    // tool while most of the body is outside), so every face's
    // boundary is curve-sampled as well.
    // An interior point per face as well: a face whose boundary lies
    // on `outer`'s can still bulge out of it.
    let mut probes: Vec<popkern_core::geometry::Point3> =
        inner.vertices.values().map(|v| v.position).collect();
    for face in inner.faces.keys() {
        probes.extend(inner.face_outer_loop_sample_positions(face, 8));
        probes.extend(super::contact::face_interior_point(inner, face));
    }
    if probes.is_empty() {
        return false;
    }
    for probe in &probes {
        // The ray vote has no on-boundary verdict of its own (a corner
        // of a coincident body votes either way), so boundary points
        // are recognised geometrically first.
        if crate::postprocess::zero_thickness::on_body_boundary(probe, outer, &aabbs, &tolerance) {
            continue;
        }
        match point_in_solid_final(probe, outer, &aabbs, &tolerance) {
            Ok(PointClassification::Inside) | Ok(PointClassification::OnBoundary) => {}
            _ => return false,
        }
    }
    true
}

/// Verify a claimed EMPTY intersection directly against the original
/// operands: their interiors are disjoint, so no probe of either —
/// vertices, sampled loop boundaries, one interior point per face —
/// lies strictly inside the other, and at least one lies strictly
/// outside (operands whose probes all sit on each other's boundary
/// can be coincident, and their intersection is the whole body). The
/// selection's own verdict cannot certify itself: a walk that
/// classifies every fragment away looks exactly like disjointness.
fn verify_disjoint_interiors(a: &BRep, b: &BRep) -> bool {
    use crate::classify::point_in_solid_final;
    use crate::error::PointClassification;
    use crate::intersect::compute_face_aabbs;
    let mut outside = false;
    for (inner, outer) in [(a, b), (b, a)] {
        let aabbs = compute_face_aabbs(outer);
        let tolerance = outer.tolerance;
        let mut probes: Vec<popkern_core::geometry::Point3> =
            inner.vertices.values().map(|v| v.position).collect();
        for face in inner.faces.keys() {
            probes.extend(inner.face_outer_loop_sample_positions(face, 8));
            probes.extend(super::contact::face_interior_point(inner, face));
        }
        for probe in &probes {
            // The ray vote has no on-boundary verdict of its own (a
            // point on a contact line votes either way), so contact
            // points are recognised geometrically first.
            if crate::postprocess::zero_thickness::on_body_boundary(
                probe, outer, &aabbs, &tolerance,
            ) {
                continue;
            }
            match point_in_solid_final(probe, outer, &aabbs, &tolerance) {
                Ok(PointClassification::Outside) => outside = true,
                Ok(PointClassification::OnBoundary) => {}
                _ => return false,
            }
        }
    }
    outside
}

/// EMPTY-SERVE LIE GATE (PK-465): an empty subtraction result is
/// legitimate ONLY when A is wholly inside B, an empty intersection
/// only when the operands' interiors are disjoint, and an empty union
/// never. A garbage walk that classifies every fragment away otherwise
/// SERVES an empty body (measured: the coarse 1e-2 sphere−cube
/// envelope pin walked the sphere chart into a fused mega-orbit and
/// served volume 0.0 against truth 0.202).
pub(super) fn certify_empty_result(
    op: crate::error::BoolOperation,
    original_a: &BRep,
    original_b: &BRep,
) -> Result<(), crate::error::BoolError> {
    let legitimate = match op {
        crate::error::BoolOperation::Subtraction => {
            verify_whole_body_containment(original_a, original_b)
        }
        crate::error::BoolOperation::Intersection => {
            verify_disjoint_interiors(original_a, original_b)
        }
        crate::error::BoolOperation::Union => false,
    };
    if legitimate {
        return Ok(());
    }
    Err(crate::error::BoolError::MergeFailed {
        stage: "Merge".to_string(),
        reason: "empty result from operands that are not in the containment \
                 configuration an empty result requires; declining rather than \
                 serving a vanished body"
            .to_string(),
    })
}

/// The contacts of a failed union whose operands only touch: their
/// interiors are disjoint and, set side by side, their shells meet.
/// That union is non-manifold wherever they meet (or, where they share
/// area, would have merged — but the pipeline already failed, so only
/// the decline's reason changes). Never a body: the caller declines
/// with these contacts in place of the pipeline's failure.
///
/// Operands that share area (a face lying on a face of the other) are
/// not merely touching: their union merges there into one manifold
/// body, so a pipeline failure on them is that failure, not a contact,
/// and `None` keeps it (PK-1124).
pub(super) fn touching_union_contacts(
    a: &BRep,
    b: &BRep,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Option<Vec<crate::error::NonManifoldContact>> {
    if !verify_disjoint_interiors(a, b) || shares_area(a, b, tolerance) {
        return None;
    }
    let side_by_side = super::helpers::combine_breps(a.clone(), b.clone()).ok()?;
    let contacts = crate::postprocess::shell_contact::shell_contacts(&side_by_side, tolerance);
    (!contacts.is_empty()).then_some(contacts)
}

/// Whether a face of `a` lies on `b`'s boundary over an area: at least
/// two of its three face-interior probes on `b`'s boundary. One is not
/// enough, since a point contact can sit exactly at a face's centre
/// (a cube's inscribed ball touches each face there).
pub(super) fn shares_area(
    a: &BRep,
    b: &BRep,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    let aabbs_b = crate::intersect::compute_face_aabbs(b);
    a.faces.keys().any(|f| {
        super::relation_probes::face_interior_points(a, f, 3)
            .iter()
            .filter(|p| {
                crate::postprocess::zero_thickness::on_body_boundary(p, b, &aabbs_b, tolerance)
            })
            .count()
            >= 2
    })
}

#[cfg(test)]
#[path = "operation_empty_tests.rs"]
mod tests;
