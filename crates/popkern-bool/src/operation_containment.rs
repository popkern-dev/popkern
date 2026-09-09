//! The verified coincident-body containment shortcut (split from
//! `operation.rs` under the 500-line file budget).

use crate::intersect::IntersectionResult;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Coincident-body shortcut. The old gate was a FACE-COUNT
/// heuristic — "overlaps ≥ half the smaller body's faces ⇒ full
/// containment" — and it silently swallowed a flush-sided gusset
/// whose 2 coplanar caps tripped the count on a 5-face tool: the
/// union returned the bare arm, watertight and χ-clean, minus
/// 44 cm³ of material. Containment is now VERIFIED geometrically:
/// every face-interior probe of the candidate body must classify
/// Inside the other. Bodies that merely SHARE some coplanar faces
/// fall through to the full pipeline. Returns `true` when the
/// no-intersection containment path applies.
pub(crate) fn coincident_containment_applies(
    a: &BRep,
    b: &BRep,
    intersection_result: &IntersectionResult,
    tolerance: &Tolerance,
) -> bool {
    if intersection_result.coincidence.face_overlaps.is_empty() {
        return false;
    }
    let n_a = a.faces.len();
    let n_b = b.faces.len();
    let n_overlaps = intersection_result.coincidence.face_overlaps.len();
    if n_overlaps < n_a.min(n_b) / 2 {
        return false;
    }
    // Probe FACE INTERIOR samples, not vertices: corner
    // vertices sit on several faces of the other body at
    // once and the ray votes degrade (two identical cubes
    // misclassified their own corners as Outside). A face
    // interior point of a genuinely contained body is
    // Inside-or-On the other; a poking face has interior
    // points strictly Outside.
    // Each probe is offset INWARD from its face sample (the
    // Newell winding of the outer loop gives the outward
    // side): an on-boundary sample of a coincident face has
    // ill-defined ray parity, but its inward offset is
    // strictly interior to the candidate — Inside the other
    // body iff containment holds. Any face that cannot be
    // probed declines the shortcut (conservative: the full
    // pipeline runs).
    let contained = |inner: &BRep,
                     outer: &BRep,
                     outer_aabbs: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )]|
     -> bool {
        let eps = tolerance.linear() * 50.0;
        for (fh, face) in inner.faces.iter() {
            let sample = crate::classify::pick_sample_point(inner, fh);
            // Outward normal from the outer loop's Newell
            // vector (loops wind CCW seen from outside).
            let Some(loop_data) = inner.loops.get(face.outer_loop) else {
                return false;
            };
            let start = loop_data.half_edge;
            if !inner.half_edges.contains_key(start) {
                return false;
            }
            let mut pts = Vec::new();
            let mut he = start;
            loop {
                pts.push(inner.vertices[inner.half_edges[he].origin].position);
                he = inner.half_edges[he].next;
                if he == start {
                    break;
                }
            }
            if pts.len() < 3 {
                return false;
            }
            let mut n = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
            for i in 0..pts.len() {
                let p = pts[i];
                let q = pts[(i + 1) % pts.len()];
                n.x += (p.y - q.y) * (p.z + q.z);
                n.y += (p.z - q.z) * (p.x + q.x);
                n.z += (p.x - q.x) * (p.y + q.y);
            }
            let len = n.length();
            if len < 1.0e-12 {
                return false;
            }
            // Measured convention: this kernel's loop Newell
            // vector points INWARD (a unit cube's bottom
            // face yielded +z), so the interior probe adds.
            let probe = sample + n * (eps / len);
            let verdict =
                crate::classify::point_in_solid_final(&probe, outer, outer_aabbs, tolerance);
            if popkern_core::trace::trace_enabled("containment") {
                eprintln!(
                    "[contain] {fh:?} probe ({:.4},{:.4},{:.4}) -> {verdict:?}",
                    probe.x, probe.y, probe.z
                );
            }
            match verdict {
                Ok(crate::error::PointClassification::Inside) => {}
                Ok(_) => return false,
                Err(_) => return false,
            }
        }
        true
    };
    let a_in_b = n_a <= n_b && contained(a, b, &intersection_result.aabbs_b);
    let b_in_a = !a_in_b && n_b <= n_a && contained(b, a, &intersection_result.aabbs_a);
    if a_in_b || b_in_a {
        true
    } else {
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!(
                "[boolean] coincident shortcut declined: {n_overlaps} overlaps but \
                 containment not verified — full pipeline"
            );
        }
        false
    }
}
