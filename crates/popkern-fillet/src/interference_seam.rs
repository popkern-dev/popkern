//! PK-159 mid-face strip-interference completion (blend×blend
//! intersection trimming, stage 3 of
//! `docs/architecture/blend-blend-intersection-trimming.md`).
//!
//! Two fillets on vertex-disjoint edges of one support face whose
//! bands collide mid-face used to decline transactionally
//! (`FilletInterference`, PK-140). For PARALLEL cylinder bands the
//! true mutual boundary is analytic: in the cross-section the two
//! circles intersect, and the seam is the straight line through the
//! material-side intersection point, parallel to both axes. Each
//! fillet's trim on the shared face is REPLACED by that seam line
//! before the apply phase: the two fillets then consume
//! complementary regions of the shared face and share the seam edge,
//! with each strip trimmed to the seam (the seam lies ON both
//! cylinders by construction).
//!
//! **Fallback:** any pair that is not two parallel cylinders with a
//! transversal circle intersection keeps the structured decline —
//! this pass only upgrades the configurations whose seam is exact.
//! **Tolerance:** the caller's blend tolerance for the tangency
//! degeneracy gate; the seam itself is closed-form.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::EdgeGeometry;
use crate::overlap_strip_interference::StripInterference;

/// Attempt to upgrade every detected interference pair to a seam
/// completion by patching the shared-face trims in place. Returns
/// `Err` with the first decline reason when ANY hit cannot be
/// patched (the caller declines or drops per its policy, leaving
/// the geometries untouched semantically — a declined operation
/// never reaches apply).
pub(crate) fn try_complete_interference(
    edge_geometries: &mut [EdgeGeometry],
    brep: &BRep,
    hits: &[StripInterference],
    tolerance: &Tolerance,
) -> Result<(), String> {
    let trace = popkern_core::trace::trace_enabled("interference");
    for hit in hits {
        let (ia, ib) = (hit.fillet_a, hit.fillet_b);
        if ia >= edge_geometries.len() || ib >= edge_geometries.len() {
            return Err("interference pair index out of range".into());
        }
        let (cyl_a, cyl_b) = match (
            cylinder_of(&edge_geometries[ia]),
            cylinder_of(&edge_geometries[ib]),
        ) {
            (Some(a), Some(b)) => (a, b),
            _ => {
                if trace {
                    eprintln!("[iseam] pair ({ia},{ib}): non-cylinder strip — decline stands");
                }
                return Err(
                    "non-cylinder strip — the analytic seam completion only covers \
                     parallel-cylinder pairs"
                        .into(),
                );
            }
        };
        let axis = match cyl_a.1.normalized() {
            Some(a) => a,
            None => return Err("degenerate strip axis".into()),
        };
        if axis.cross(&cyl_b.1).length() > 1e-6 * cyl_b.1.length() {
            if trace {
                eprintln!("[iseam] pair ({ia},{ib}): axes not parallel — decline stands");
            }
            return Err(
                "strip axes not parallel — the analytic seam completion only covers \
                 parallel-cylinder pairs"
                    .into(),
            );
        }
        // Shared support face and which trim side each fillet has on it.
        let Some((shared, side_a, side_b)) =
            shared_face_sides(&edge_geometries[ia], &edge_geometries[ib])
        else {
            return Err("no shared support face resolved for the pair".into());
        };
        let Some(face_surface) = brep
            .faces
            .get(shared)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
        else {
            return Err("shared support face has no surface".into());
        };
        // Cross-section circle intersection.
        let (c_a, _, r_a) = cyl_a;
        let (c_b, _, r_b) = cyl_b;
        let sep = c_b - c_a;
        let w = sep - axis * sep.dot(&axis);
        let d = w.length();
        let gate = tolerance.linear().max(1e-9) * 10.0;
        if d <= gate || d >= r_a + r_b - gate || d <= (r_a - r_b).abs() + gate {
            if trace {
                eprintln!(
                    "[iseam] pair ({ia},{ib}): degenerate circle pair d={d:.4} r=({r_a:.4},{r_b:.4}) — decline stands"
                );
            }
            return Err(format!(
                "degenerate cross-section circle pair (d = {d:.4}, r = ({r_a:.4}, {r_b:.4}))"
            ));
        }
        let w_hat = match w.normalized() {
            Some(v) => v,
            None => return Err("degenerate axis separation".into()),
        };
        let u_hat = axis.cross(&w_hat);
        let alpha = (d * d + r_a * r_a - r_b * r_b) / (2.0 * d);
        let h = (r_a * r_a - alpha * alpha).sqrt();
        let base = c_a + w_hat * alpha;
        let p_plus = base + u_hat * h;
        let p_minus = base - u_hat * h;
        // Material side: the candidate nearer the shared face plane.
        let plane_dist = |p: Point3| -> f64 {
            match &face_surface {
                SurfaceGeometry::Plane { origin, normal, .. } => {
                    let n = normal.normalized().unwrap_or(*normal);
                    ((p - *origin).dot(&n)).abs()
                }
                _ => f64::INFINITY,
            }
        };
        let seam_pt = if plane_dist(p_plus) <= plane_dist(p_minus) {
            p_plus
        } else {
            p_minus
        };
        if plane_dist(seam_pt).is_infinite() {
            if trace {
                eprintln!("[iseam] pair ({ia},{ib}): shared face not planar — decline stands");
            }
            return Err("shared support face is not planar".into());
        }
        if trace {
            eprintln!(
                "[iseam] pair ({ia},{ib}) shared {shared:?}: seam through ({:.4},{:.4},{:.4}) depth {:.4e}",
                seam_pt.x,
                seam_pt.y,
                seam_pt.z,
                plane_dist(seam_pt)
            );
        }
        // Patch both fillets' shared-side trims to the seam line,
        // spanning each old trim's own axial extent.
        for (idx, side) in [(ia, side_a), (ib, side_b)] {
            if !patch_trim_to_seam(
                &mut edge_geometries[idx],
                side,
                seam_pt,
                axis,
                &face_surface,
            ) {
                return Err("seam trim patch failed (UV projection)".into());
            }
        }
    }
    Ok(())
}

/// The strip cylinder of one computed fillet geometry.
fn cylinder_of(eg: &EdgeGeometry) -> Option<(Point3, Vec3, f64)> {
    match &eg.surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some((*origin, *axis, *radius)),
        _ => None,
    }
}

/// The face both fillets support on, with each fillet's trim side
/// (1 or 2) on it.
fn shared_face_sides(a: &EdgeGeometry, b: &EdgeGeometry) -> Option<(FaceHandle, u8, u8)> {
    for (fa, sa) in [(a.analysis.face_1, 1u8), (a.analysis.face_2, 2u8)] {
        for (fb, sb) in [(b.analysis.face_1, 1u8), (b.analysis.face_2, 2u8)] {
            if fa == fb {
                return Some((fa, sa, sb));
            }
        }
    }
    None
}

/// Replace one fillet's shared-face trim (3-D curve + pcurve + exact
/// circle) with the seam line, spanning the old trim's axial extent.
fn patch_trim_to_seam(
    eg: &mut EdgeGeometry,
    side: u8,
    seam_pt: Point3,
    axis: Vec3,
    face_surface: &SurfaceGeometry,
) -> bool {
    use popkern_core::geometry::CurveEval;
    let old = if side == 1 {
        &eg.trim_3d_1
    } else {
        &eg.trim_3d_2
    };
    let (t0, t1) = old.parameter_range();
    let (e0, e1) = (old.point_at(t0), old.point_at(t1));
    let project = |e: Point3| -> Point3 { seam_pt + axis * (e - seam_pt).dot(&axis) };
    let (q0, q1) = (project(e0), project(e1));
    let line = line_nurbs(q0, q1);
    let (Some(uv0), Some(uv1)) = (
        project_point_to_uv(face_surface, &q0),
        project_point_to_uv(face_surface, &q1),
    ) else {
        return false;
    };
    let pline = line_nurbs(
        Point3::new(uv0.0, uv0.1, 0.0),
        Point3::new(uv1.0, uv1.1, 0.0),
    );
    if side == 1 {
        eg.trim_3d_1 = line;
        eg.pcurve_1 = pline;
        eg.trim_circle_1 = None;
    } else {
        eg.trim_3d_2 = line;
        eg.pcurve_2 = pline;
        eg.trim_circle_2 = None;
    }
    true
}

fn line_nurbs(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

/// Detect interference, upgrade the analytic parallel class to a
/// seam completion (recorded as an overlap warning), and resolve
/// the rest per the caller's [`InterferencePolicy`]: decline with
/// the structured `FilletInterference` error (default), or drop
/// both members of every interfering pair and serve the rest,
/// reporting each drop (PK-185).
pub(crate) fn complete_or_decline_interference(
    edge_geometries: &mut Vec<EdgeGeometry>,
    brep: &BRep,
    tolerance: &Tolerance,
    policy: crate::types::InterferencePolicy,
    overlap_warnings: &mut Vec<crate::error::BlendWarning>,
) -> Result<(), crate::error::BlendError> {
    let hits = crate::overlap_strip_interference::detect_strip_interference(
        edge_geometries,
        brep,
        tolerance,
    );
    if hits.is_empty() {
        return Ok(());
    }
    match try_complete_interference(edge_geometries, brep, &hits, tolerance) {
        Ok(()) => {
            for hit in &hits {
                overlap_warnings.push(crate::error::BlendWarning::FilletOverlap {
                    fillet_a: hit.fillet_a,
                    fillet_b: hit.fillet_b,
                });
            }
            Ok(())
        }
        Err(reason) => match policy {
            crate::types::InterferencePolicy::Decline => {
                Err(crate::error::BlendError::FilletInterference {
                    fillet_a: hits[0].fillet_a,
                    fillet_b: hits[0].fillet_b,
                })
            }
            crate::types::InterferencePolicy::ServeNonInterfering => {
                // Drop BOTH members of every hit pair; indices in
                // the warnings refer to the ORIGINAL spec order
                // (edge_geometries is 1:1 with specs here — the
                // transactional compute gate ran first).
                let mut dropped = std::collections::BTreeSet::new();
                for hit in &hits {
                    dropped.insert(hit.fillet_a);
                    dropped.insert(hit.fillet_b);
                    overlap_warnings.push(crate::error::BlendWarning::FilletPairDropped {
                        fillet_a: hit.fillet_a,
                        fillet_b: hit.fillet_b,
                        reason: reason.clone(),
                    });
                }
                let mut idx = 0usize;
                edge_geometries.retain(|_| {
                    let keep = !dropped.contains(&idx);
                    idx += 1;
                    keep
                });
                Ok(())
            }
        },
    }
}
