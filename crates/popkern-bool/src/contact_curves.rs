//! Tangent records as contact curves (W3c, tangency design §4.8): a
//! record along which the operands' faces are tangent at the model's
//! resolution is where they touch, not where one crosses into the
//! other. Where every record of a boolean is one, the operands meet in
//! no crossing: the records leave the arrangement and go to the contact
//! reading with the touches, which imprints each on a reading copy
//! (`BRep::mark_contact_curve`) and reads the line contact from its
//! topology; the general glue imprints it in the body itself.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, PortableMath};
use popkern_core::topology::BRep;

use crate::intersect::IntersectionRecord;

/// What construction found the operands touch at, beside its body: the
/// contact reading's input (`operation::answer_contacts`).
#[derive(Debug, Clone, Default)]
pub(crate) struct Touches {
    /// Isolated touches.
    pub points: Vec<Point3>,
    /// Curves the faces are tangent along.
    pub curves: Vec<CurveOnEdge>,
}

/// Stations along a record where its tangency is read.
const STATIONS: [f64; 5] = [0.1, 0.3, 0.5, 0.7, 0.9];

/// Whether `record`'s faces are tangent along it at the model's
/// resolution: at every station the normals are parallel within
/// `angular`, or the angle `θ` between them and their cross-section
/// curvatures across the record, `Δκ` apart, enclose a lens no thicker
/// than the linear tolerance (`θ² ≤ 2·linear·Δκ`, the radial order's
/// rule): two faces overlapping by 5e-7 along a ruling meet in a record
/// crossing at 1.4e-3 rad, tangent within the tolerance.
pub(crate) fn is_contact_curve(
    record: &IntersectionRecord,
    (a, b): (&BRep, &BRep),
    (linear, angular): (f64, f64),
) -> bool {
    let (Some(sa), Some(sb)) = (surface_of(a, record.face_a), surface_of(b, record.face_b)) else {
        return false;
    };
    let c = &record.curve;
    let (t0, t1) = c.curve_3d.parameter_range();
    STATIONS.iter().all(|&f| {
        let s = t0 + (t1 - t0) * f;
        let (q1, q2) = (c.curve_uv1.point_at(s), c.curve_uv2.point_at(s));
        let (Some(n1), Some(n2), Some(t)) = (
            sa.normal_at(q1.x, q1.y).normalized(),
            sb.normal_at(q2.x, q2.y).normalized(),
            c.curve_3d.tangent_at(s).normalized(),
        ) else {
            return false;
        };
        let theta = n1.cross(&n2).length().min(1.0).pasin();
        if theta <= angular {
            return true;
        }
        let Some(d) = n1.cross(&t).normalized() else {
            return false;
        };
        let k =
            |sf, uv| popkern_core::topology::radial_order::cross_section_curvature(sf, uv, d, t);
        let (Some(k1), Some(k2)) = (k(sa, (q1.x, q1.y)), k(sb, (q2.x, q2.y))) else {
            return false;
        };
        theta * theta <= 2.0 * linear * (k1 - k2).abs()
    })
}

fn surface_of(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> Option<&popkern_core::geometry::surface::SurfaceGeometry> {
    Some(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}

/// `record` as a contact curve: its 3-D curve over its range.
pub(crate) fn as_contact_curve(record: &IntersectionRecord) -> CurveOnEdge {
    let curve = &record.curve.curve_3d;
    let (t_start, t_end) = curve.parameter_range();
    CurveOnEdge {
        geometry: popkern_core::geometry::curve::CurveGeometry::Nurbs(curve.clone()),
        t_start,
        t_end,
    }
}

/// `curves` less those whose middle lies on one of `records` within the
/// linear reach: a record there carries the contact.
pub(crate) fn not_on_records(
    curves: Vec<CurveOnEdge>,
    records: &[IntersectionRecord],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Vec<CurveOnEdge> {
    let reach = tolerance.linear_reach();
    curves
        .into_iter()
        .filter(|c| {
            let mid = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
            !records.iter().any(|r| {
                let curve = &r.curve.curve_3d;
                crate::curve_projection::project_onto_curve_dense(curve, &mid)
                    .is_some_and(|t| curve.point_at(t).distance_to(&mid) <= reach)
            })
        })
        .collect()
}

#[cfg(test)]
#[path = "contact_curves_tests.rs"]
mod tests;
