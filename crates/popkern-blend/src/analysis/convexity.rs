//! Winding-free edge convexity classification.
//!
//! Extracted from `analysis::mod` (file-size cap). Samples the
//! into-face/outward-normal invariant along the edge; see
//! [`classify_convexity`] for the sign convention and the trust
//! model around stored orientation flags, shell orientation, and
//! the pcurve-free fallback for surgery output.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, SurfaceEval, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::{Convexity, shell_outward_sign, surface_normal_at_point_natural};

const CONVEXITY_SAMPLES: usize = 7;

/// Classify the convexity of an edge between two surfaces.
///
/// Samples normals at [`CONVEXITY_SAMPLES`] points along the edge and
/// uses the sign of (N1 x T) . N2 at each sample:
///   - positive -> convex
///   - negative -> concave
///   - near zero -> indeterminate (skipped)
///
/// If all definite samples agree, returns [`Convexity::Convex`] or
/// [`Convexity::Concave`]. If samples disagree, returns [`Convexity::Mixed`].
pub(super) fn classify_convexity(
    brep: &BRep,
    edge_handle: EdgeHandle,
    surface_1: &SurfaceGeometry,
    surface_2: &SurfaceGeometry,
    face_1: FaceHandle,
    face_2: FaceHandle,
) -> Result<Convexity, String> {
    let edge = brep
        .edges
        .get(edge_handle)
        .ok_or_else(|| format!("edge {edge_handle:?} not found"))?;

    let curve = edge
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_handle:?} has no curve geometry"))?;

    const MIXED_THRESHOLD: f64 = 1e-6;

    let mut convex_count: usize = 0;
    let mut concave_count: usize = 0;

    // Winding-free convexity: OUTWARD normals (natural normal
    // corrected by the winding-derived effective orientation) plus
    // INTO-FACE tangent directions (toward each face's uv-interior).
    // For a convex edge the material wedge is < π and each face's
    // interior direction opposes the OTHER face's outward normal;
    // for a concave (reentrant) edge it agrees:
    //     dot(d1, n2) + dot(d2, n1)  < 0  ⇔ convex, > 0 ⇔ concave.
    // (Box top edge: d1=+y into the top face, n2=−y outward front
    // ⇒ −1. Pocket floor/wall: d1=+y, n2=+y ⇒ +1.) The previous
    // formula used the RAW curve tangent — untied to either face's
    // winding — and stored orientation flags, so its sign was
    // producer-dependent: every convex plate edge on boolean output
    // classified Concave, sending the rolling ball to the exterior
    // weld-bead mirror.
    // Stored orientation flags are the manifold truth: primitives
    // encode outwardness in the flag alone (their loop windings are
    // unenforced), and boolean/sweep output is winding-flag
    // consistent since the PK-099 source fixes (validated: zero
    // mismatches across the flagship union). The former shoelace
    // override compensated for a since-fixed producer defect.
    let stored = |face: FaceHandle| -> bool {
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|sof| sof.orientation)
            .unwrap_or(true)
    };
    // GLOBAL shell orientation: primitive constructors build
    // INWARD-oriented shells (extrude keeps the base face's
    // winding), sweep/boolean output builds OUTWARD ones — the flag
    // is only honest relative to its shell. One signed-volume
    // measurement disambiguates.
    let shell = shell_outward_sign(brep);

    // Winding-free side witness (PK-264): nudge the sample point a
    // step along ±d in 3D, project into the face chart, and keep the
    // sign whose nudged point classifies INSIDE the face's trimmed
    // region (even-odd containment — loop direction never enters, and
    // holes classify correctly). Decisive only when one side is
    // Interior and the other is not; tangency, projection failure, or
    // boundary-band ambiguity keeps the winding-derived answer.
    //
    // This witness is what makes the interior-left rule below honest
    // on primitive output: cone/cylinder constructors encode
    // outwardness in the orientation flag alone and wind their
    // single-vertex base-circle loops AGAINST the flag, so the raw
    // interior-left direction mirrors there (measured: the cone base
    // rim classified Concave — both d1 and d2 inverted — sending the
    // walker to the exterior mirror ball; the cylinder base's Convex
    // verdict rested on a −0.000 knife-edge).
    let trim_side_witness =
        |face: FaceHandle, point: &popkern_core::geometry::Point3, dir: &Vec3| -> Option<Vec3> {
            let face_data = brep.faces.get(face)?;
            let surface = &face_data.surface.as_ref()?.geometry;
            let step = 1.0e-3;
            let classify = |p: &popkern_core::geometry::Point3| {
                let (u, v) =
                    popkern_core::geometry::uv_projection::project_point_to_uv(surface, p)?;
                brep.face_uv_trim_containment(face, surface, u, v)
            };
            use popkern_core::validation::FaceUvContainment as C;
            let plus = classify(&(*point + *dir * step));
            let minus = classify(&(*point - *dir * step));
            match (plus, minus) {
                (Some(C::Interior), Some(C::Outside | C::Boundary)) => Some(*dir),
                (Some(C::Outside | C::Boundary), Some(C::Interior)) => Some(*dir * -1.0),
                _ => None,
            }
        };

    // Into-face direction at a boundary sample: walking a face's own
    // half-edge with the face's OUTWARD normal up, the interior is on
    // the LEFT — d = n_outward × t̂, with t̂ the half-edge traversal
    // direction from its pcurve. Holds for outer AND inner loops
    // (that is why inner loops wind opposite), unlike any
    // uv-box-center heuristic (the plate top's box center sits in
    // its own hole). The trim-side witness above overrides the sign
    // whenever it is decisive; this rule remains the tie-breaker.
    let into_face_dir = |face: FaceHandle,
                         point: &popkern_core::geometry::Point3,
                         outward: &Vec3,
                         t_sample: f64|
     -> Option<Vec3> {
        let face_data = brep.faces.get(face)?;
        let surface = &face_data.surface.as_ref()?.geometry;
        // The half-edge of THIS edge whose loop belongs to `face`.
        let edge_data = brep.edges.get(edge_handle)?;
        let he_a = edge_data.half_edge;
        let he_b = brep.half_edges.get(he_a)?.twin;
        let mut he_face = None;
        for cand in [he_a, he_b] {
            let d = brep.half_edges.get(cand)?;
            if brep.loops.get(d.loop_)?.face == face {
                he_face = Some(cand);
                break;
            }
        }
        let he_face = he_face?;
        let d = brep.half_edges.get(he_face)?;
        let Some(pcurve) = d.pcurve.as_ref() else {
            // Pcurve-free fallback: the traversal tangent comes from
            // the EDGE curve (oriented along the primary half-edge by
            // convention), flipped when this face uses the twin.
            // Surgery output (unfillet-restored edges) can lack
            // pcurves entirely; without this fallback every sample
            // is indeterminate and the edge misclassifies Mixed —
            // which sent the recreated rolling ball to the exterior
            // mirror (measured: strip centre (−0.2,−0.2) on a unit
            // cube re-fillet after unfillet).
            let mut tangent = curve.geometry.derivative_at(t_sample, 1).normalized()?;
            // Orient along THIS half-edge by vertex positions — no
            // curve-direction convention is assumed (surgery output
            // violates the primary-half-edge convention).
            let origin = brep.vertices.get(d.origin)?.position;
            let dest = brep
                .vertices
                .get(brep.half_edges.get(d.twin)?.origin)?
                .position;
            let along = dest - origin;
            if along.dot(&tangent) < 0.0 {
                tangent = tangent * -1.0;
            }
            let mut dir = (outward.cross(&tangent)).normalized()?;
            // Winding-free sign witness: surgery output (the only
            // producer of pcurve-less half-edges) does not honor the
            // winding-flag convention the interior-left rule rests
            // on, so verify the side against the face's outer-loop
            // centroid projected into the tangent plane. Ambiguous
            // witnesses (centroid on the edge) keep the winding
            // answer.
            let loop_points = brep.face_outer_loop_positions(face);
            if loop_points.len() >= 3 {
                let m = loop_points.len() as f64;
                let centroid = popkern_core::geometry::Point3::new(
                    loop_points.iter().map(|p| p.x).sum::<f64>() / m,
                    loop_points.iter().map(|p| p.y).sum::<f64>() / m,
                    loop_points.iter().map(|p| p.z).sum::<f64>() / m,
                );
                let w = centroid - *point;
                let w_t = w - *outward * outward.dot(&w);
                if w_t.length() > 1e-9 && w_t.dot(&dir) < 0.0 {
                    dir = dir * -1.0;
                }
            }
            return Some(dir);
        };
        let (t0, t1) = pcurve.parameter_range();
        // Pcurve parameter nearest the sample point (in 3D).
        let mut best = (f64::INFINITY, t0);
        for k in 0..=16 {
            let t = t0 + (t1 - t0) * k as f64 / 16.0;
            let uv = pcurve.point_at(t);
            let q = surface.point_at(uv.x, uv.y);
            let dist = q.distance_to(point);
            if dist < best.0 {
                best = (dist, t);
            }
        }
        let t = best.1;
        let h = (t1 - t0) * 1.0e-4;
        let (ta, tb) = if t + h <= t1 { (t, t + h) } else { (t - h, t) };
        let uv_a = pcurve.point_at(ta);
        let uv_b = pcurve.point_at(tb);
        let pa = surface.point_at(uv_a.x, uv_a.y);
        let pb = surface.point_at(uv_b.x, uv_b.y);
        let mut tangent = (pb - pa).normalized()?;
        // Direction verification (no repo-wide pcurve direction
        // contract; see effective_outward_orientation).
        let origin = brep.vertices.get(d.origin)?.position;
        let dest = brep
            .vertices
            .get(brep.half_edges.get(d.twin)?.origin)?
            .position;
        let uv0 = pcurve.point_at(t0);
        let q0 = surface.point_at(uv0.x, uv0.y);
        if q0.distance_to(&dest) + 1.0e-12 < q0.distance_to(&origin) {
            tangent = tangent * -1.0;
        }
        (outward.cross(&tangent)).normalized()
    };

    for i in 0..CONVEXITY_SAMPLES {
        let frac = (i as f64) / ((CONVEXITY_SAMPLES - 1).max(1) as f64);
        let t = curve.t_start + (curve.t_end - curve.t_start) * frac;
        let point = curve.geometry.point_at(t);

        // Outward normals fold in the shell sign; the interior-left
        // rule below is a property of the face's OWN manifold
        // normal (flag-corrected natural, shell-independent) — on an
        // inward-oriented shell the two differ and mixing them
        // mirrors the into-face direction.
        let nat_1 = surface_normal_at_point_natural(surface_1, &point);
        let nat_2 = surface_normal_at_point_natural(surface_2, &point);
        let flag_sign = |f: bool| if f { 1.0 } else { -1.0 };
        let m1 = nat_1 * flag_sign(stored(face_1));
        let m2 = nat_2 * flag_sign(stored(face_2));
        let n1 = m1 * shell;
        let n2 = m2 * shell;

        let Some(d1) = into_face_dir(face_1, &point, &m1, t) else {
            continue;
        };
        let Some(d2) = into_face_dir(face_2, &point, &m2, t) else {
            continue;
        };
        let d1 = trim_side_witness(face_1, &point, &d1).unwrap_or(d1);
        let d2 = trim_side_witness(face_2, &point, &d2).unwrap_or(d2);
        let dot = d1.dot(&n2) + d2.dot(&n1);
        if popkern_core::trace::trace_enabled("convexity") && i == 0 {
            eprintln!(
                "  [cvx] {edge_handle:?} p({:.2},{:.2},{:.2}) n1({:.2},{:.2},{:.2}) d1({:.2},{:.2},{:.2}) n2({:.2},{:.2},{:.2}) d2({:.2},{:.2},{:.2}) dot {dot:.3}",
                point.x,
                point.y,
                point.z,
                n1.x,
                n1.y,
                n1.z,
                d1.x,
                d1.y,
                d1.z,
                n2.x,
                n2.y,
                n2.z,
                d2.x,
                d2.y,
                d2.z
            );
        }

        if dot < -MIXED_THRESHOLD {
            convex_count += 1;
        } else if dot > MIXED_THRESHOLD {
            concave_count += 1;
        }
        // Near-zero samples are indeterminate — skipped.
    }

    if convex_count > 0 && concave_count > 0 {
        Ok(Convexity::Mixed)
    } else if convex_count > 0 {
        Ok(Convexity::Convex)
    } else if concave_count > 0 {
        Ok(Convexity::Concave)
    } else {
        // All samples indeterminate.
        Ok(Convexity::Mixed)
    }
}
