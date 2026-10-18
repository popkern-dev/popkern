//! Clip gaps versus marching slop between two pieces of one face
//! pair's intersection curve (PK-680).
//!
//! Consolidation welds the facing ends of a face pair's pieces when
//! they lie within the marching slop class (3.5e-4..1.2e-3 measured).
//! Pieces of one curve also end facing each other where the CLIP cut
//! the curve at two distinct junctions — and on fine features those
//! junctions can sit inside the slop class too (measured: an M2.5
//! rod's slot, whose top edge passes 2.0e-3 under the point where the
//! crest helix crosses the slot wall; the two pieces of the major
//! cylinder × wall line were averaged onto the midpoint, the crest
//! edge was split there, and its end left the flank by 8.7e-4). The
//! distance cannot tell the two apart; the trims can: between two
//! traces of one junction the curve stays inside both faces, while a
//! clip gap is exactly where it leaves one of them.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceEval as _;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::IntersectionRecord;

/// True when the stretch of the intersection curve between the
/// facing ends of `ra` (`a_end`: its end, else its start) and `rb`
/// leaves the trim of either face: the ends are two clip junctions,
/// not two traces of one, and must keep their own positions.
pub(super) fn gap_leaves_a_face(
    a: &BRep,
    b: &BRep,
    ra: &IntersectionRecord,
    a_end: bool,
    rb: &IntersectionRecord,
    b_end: bool,
    tolerance: &Tolerance,
) -> bool {
    let outside = |brep: &BRep, face: FaceHandle, ca: &NurbsCurve, cb: &NurbsCurve| {
        let Some(mid) = chart_midpoint(brep, face, ca, a_end, cb, b_end) else {
            return false;
        };
        !crate::split_inner_loop::point_inside_face_region_uv_with_density(
            brep, face, mid, tolerance, 64,
        )
    };
    outside(a, ra.face_a, &ra.curve.curve_uv1, &rb.curve.curve_uv1)
        || outside(b, ra.face_b, &ra.curve.curve_uv2, &rb.curve.curve_uv2)
}

/// The chart midpoint of two pieces' facing ends on `face`, with the
/// second end brought onto the first's period branch. `None` when the
/// face has no surface.
fn chart_midpoint(
    brep: &BRep,
    face: FaceHandle,
    ca: &NurbsCurve,
    a_end: bool,
    cb: &NurbsCurve,
    b_end: bool,
) -> Option<Point3> {
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let end_of = |c: &NurbsCurve, end: bool| {
        let (t0, t1) = c.parameter_range();
        c.point_at(if end { t1 } else { t0 })
    };
    let (pa, mut pb) = (end_of(ca, a_end), end_of(cb, b_end));
    let (u_period, v_period) = surface.parameter_periodicity();
    if let Some(p) = u_period {
        pb.x -= p * ((pb.x - pa.x) / p).round();
    }
    if let Some(p) = v_period {
        pb.y -= p * ((pb.y - pa.y) / p).round();
    }
    Some(Point3::new(0.5 * (pa.x + pb.x), 0.5 * (pa.y + pb.y), 0.0))
}

/// True when `weld` lies off either face pair surface by more than the
/// linear tolerance (or than the ends it replaces already did): moving
/// the facing ends there would leave the records off their own
/// intersection (PK-788). Two traces of one junction weld onto their
/// midpoint within the curve's sagitta of both surfaces; two DISTINCT
/// exact junctions inside the slop gate do not — the cap of a boss
/// 1e-5 inside a ring's extreme meets the ring wall in two exact lines
/// 0.068 apart, and their averaged start sat 1e-5 off the ring wall
/// at a point no junction occupies.
pub(super) fn weld_leaves_a_surface(
    a: &BRep,
    b: &BRep,
    record: &IntersectionRecord,
    ends: [Point3; 2],
    weld: Point3,
    tolerance: &Tolerance,
) -> bool {
    let off = |brep: &BRep, face: FaceHandle, p: &Point3| -> Option<f64> {
        let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
        let (u, v) = surface.closest_point(p)?;
        Some(surface.point_at(u, v).distance_to(p))
    };
    [(a, record.face_a), (b, record.face_b)]
        .into_iter()
        .any(|(brep, face)| {
            let Some(d) = off(brep, face, &weld) else {
                return false;
            };
            // The solvers' own tangency band: a lens exactly one
            // tolerance deep has a sagitta of exactly one tolerance, so
            // its welded midpoint and its secant crossings must be read
            // against the same open bound.
            let scale = [ends[0], ends[1], weld]
                .iter()
                .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
                .fold(0.0, f64::max);
            let allowed = ends
                .iter()
                .filter_map(|p| off(brep, face, p))
                .fold(tolerance.tangency_band(scale), f64::max);
            d >= allowed
        })
}
