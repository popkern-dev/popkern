//! The verified coincident-body containment shortcut (split from
//! `operation.rs` under the 500-line file budget).

use crate::error::BoolOperation;
use crate::intersect::IntersectionResult;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Which operand the containment shortcut verified inside the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Contained {
    /// A lies inside B.
    AInB,
    /// B lies inside A.
    BInA,
}

/// Coincident-body shortcut. The old gate was a FACE-COUNT
/// heuristic — "overlaps ≥ half the smaller body's faces ⇒ full
/// containment" — and it silently swallowed a flush-sided gusset
/// whose 2 coplanar caps tripped the count on a 5-face tool: the
/// union returned the bare arm, watertight and χ-clean, minus
/// 44 cm³ of material. Containment is now VERIFIED geometrically:
/// every face-interior probe of the candidate body must classify
/// Inside the other. Bodies that merely SHARE some coplanar faces
/// fall through to the full pipeline. Returns which body the probes
/// verified inside the other, the one answer the caller dispatches on:
/// re-probing the relation reads every probe of coincident boundaries
/// On and declined `CoincidentFaces` (a cone and itself turned about its
/// axis).
///
/// Containment alone decides the result only where the shared
/// faces cannot matter: A ∪ B = A and A ∩ B = B when B ⊂ A, and
/// A − B = ∅ when A ⊂ B. A − B with B ⊂ A is NOT A with a cavity
/// once B shares a face with A — the shared face opens the void to
/// the outside (a corner notch flush on three faces was served as
/// an outer shell plus a cavity coincident with it, validator-clean;
/// PK-965). That case takes the full pipeline.
pub(crate) fn coincident_containment(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    intersection_result: &IntersectionResult,
    tolerance: &Tolerance,
) -> Option<Contained> {
    if intersection_result.coincidence.face_overlaps.is_empty() {
        return None;
    }
    let n_a = a.faces.len();
    let n_b = b.faces.len();
    let n_overlaps = intersection_result.coincidence.face_overlaps.len();
    if n_overlaps < n_a.min(n_b) / 2 {
        return None;
    }
    // Probe FACE INTERIOR samples, not vertices: corner
    // vertices sit on several faces of the other body at
    // once and the ray votes degrade (two identical cubes
    // misclassified their own corners as Outside). A face
    // interior point of a genuinely contained body is
    // Inside-or-On the other; a poking face has interior
    // points strictly Outside.
    // Each probe is offset inward from its face sample along the
    // surface normal (`inward_probes`): an on-boundary sample of a
    // coincident face has ill-defined ray parity, but its inward
    // offset is strictly interior to the candidate — Inside the other
    // body iff containment holds. Any face that cannot be probed
    // declines the shortcut (conservative: the full pipeline runs).
    let contained = |inner: &BRep,
                     inner_aabbs: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
                     outer: &BRep,
                     outer_aabbs: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )]|
     -> bool {
        let eps = tolerance.linear() * 50.0;
        for fh in inner.faces.keys() {
            let probes = inward_probes(inner, fh, inner_aabbs, eps, tolerance);
            if probes.is_empty() {
                return false;
            }
            for probe in probes {
                let verdict =
                    crate::classify::point_in_solid_final(&probe, outer, outer_aabbs, tolerance);
                if popkern_core::trace::trace_enabled("containment") {
                    eprintln!(
                        "[contain] {fh:?} probe ({:.4},{:.4},{:.4}) -> {verdict:?}",
                        probe.x, probe.y, probe.z
                    );
                }
                if !matches!(verdict, Ok(crate::error::PointClassification::Inside)) {
                    return false;
                }
            }
        }
        true
    };
    let a_in_b = n_a <= n_b
        && contained(
            a,
            &intersection_result.aabbs_a,
            b,
            &intersection_result.aabbs_b,
        );
    // Every overlap counted here is a coincident face pair, so a
    // B ⊂ A subtraction always has B touching A's boundary.
    let b_in_a = !a_in_b
        && !matches!(op, BoolOperation::Subtraction)
        && n_b <= n_a
        && contained(
            b,
            &intersection_result.aabbs_b,
            a,
            &intersection_result.aabbs_a,
        );
    if a_in_b {
        return Some(Contained::AInB);
    }
    if b_in_a {
        return Some(Contained::BInA);
    }
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "[boolean] coincident shortcut declined: {n_overlaps} overlaps but \
             containment not verified — full pipeline"
        );
    }
    None
}

/// Points `eps` inside `body` off samples spread over its face `face`
/// (`face_grid_samples`; one picked sample where the grid keeps none),
/// each along the surface normal there, the side chosen by the body
/// itself (the probe classifying Inside it). The face's outer-loop vertex
/// polygon cannot give that side on a curved face: a cone's lateral loop
/// holds its rim vertex and its apex only, its Newell vector vanished, and
/// a cone united with itself turned about its axis never took the shortcut
/// (32 scenes of the generated corpus's coincident cell). One sample per
/// face decides only near itself: a flush cylinder's lateral face sampled
/// once inside the other read the pair as one inside the other. Samples
/// where neither side classifies Inside are skipped; none left declines.
fn inward_probes(
    body: &BRep,
    face: popkern_core::topology::FaceHandle,
    aabbs: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    eps: f64,
    tolerance: &Tolerance,
) -> Vec<popkern_core::geometry::Point3> {
    use popkern_core::geometry::surface::SurfaceEval as _;
    let Some(surface) = body.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return Vec::new();
    };
    let surface = &surface.geometry;
    let mut uvs = crate::classify::face_grid_samples(body, face, 4);
    if uvs.is_empty() {
        let picked = crate::classify::pick_sample_point(body, face);
        uvs.extend(popkern_core::geometry::uv_projection::project_point_to_uv(
            surface, &picked,
        ));
    }
    uvs.into_iter()
        .filter_map(|(u, v)| {
            let at = surface.point_at(u, v);
            let n = surface.normal_at(u, v).normalized()?;
            [1.0, -1.0]
                .into_iter()
                .map(|k| at + n * (k * eps))
                .find(|q| {
                    matches!(
                        crate::classify::point_in_solid_final(q, body, aabbs, tolerance),
                        Ok(crate::error::PointClassification::Inside)
                    )
                })
        })
        .collect()
}

#[cfg(test)]
#[path = "operation_containment_tests.rs"]
mod tests;
