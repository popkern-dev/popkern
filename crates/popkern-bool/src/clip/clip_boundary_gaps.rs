//! Gaps in a record running along a face's boundary line.
//!
//! A record collinear with a face's boundary edges can pass through a
//! notch in that boundary: the edges before and after the notch lie on
//! the record, the notch's mouth does not, and the record touches the
//! boundary only at the notch corners in between. The crossing finder
//! sees those touches as too few crossings for the per-interval
//! classification, and the tuned crossing hull keeps the whole span,
//! mouth included (measured: a boss cap flush with a web end cuts the
//! web's bottom face along its x = -335 boundary line; the bottom face
//! is notched over z 137..188, the clip kept 132.5..192.5, and the
//! record served as a slit across the notch on that face alone).
//!
//! The kept window is cut at every boundary vertex the record passes
//! through, and a piece the face's exact trim puts outside at all of
//! its samples, its middle clear of every boundary curve by ten
//! tolerances, is dropped: a gap is outside however it was reached,
//! while a piece riding a boundary edge is on it and stays.
//!
//! Periodic charts too (PK-1027: an arm's groove cylinder, v 0..8, kept
//! the record along its boundary ruling over v -150..60): the corners
//! are tried at their period translates, the record's uv branch need
//! not be the trim's, and the outside verdict classifies the piece's
//! 3D points by the face's multi-sheet projection.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

use super::clip_uv::UvClipResult;

/// Samples per piece for the exact verdict.
const PIECE_SAMPLES: usize = 5;

/// Whether the piece `(a, b)` of `curve_uv` is a gap of `face`. On a
/// periodic chart the record's uv need not sit on the trim's branch,
/// so the verdict classifies the piece's 3D points instead, by the
/// face's multi-sheet projection.
fn is_gap(brep: &BRep, face: FaceHandle, curve_uv: &NurbsCurve, (a, b): (f64, f64)) -> bool {
    let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let periodic = surface.geometry.periodicity() != (None, None);
    let outside = (1..=PIECE_SAMPLES).all(|k| {
        let p = curve_uv.point_at(a + (b - a) * k as f64 / (PIECE_SAMPLES + 1) as f64);
        let verdict = if periodic {
            let q = surface.geometry.point_at(p.x, p.y);
            brep.face_point_trim_containment(face, &q)
                .filter(|c| c.distance <= brep.tolerance.linear())
                .map(|c| c.containment)
        } else {
            brep.face_uv_trim_containment(face, &surface.geometry, p.x, p.y)
        };
        verdict == Some(FaceUvContainment::Outside)
    });
    let mid = curve_uv.point_at(0.5 * (a + b));
    outside
        && !super::clip_windows_clearance::on_any_boundary_exact(
            brep,
            face,
            (mid.x, mid.y),
            10.0 * brep.tolerance.linear(),
        )
}

/// `corner` and its translates by one period either way along each
/// periodic direction of `periods`.
fn translates(corner: Point3, periods: (Option<f64>, Option<f64>)) -> Vec<Point3> {
    let shifts = |p: Option<f64>| p.map_or(vec![0.0], |p| vec![-p, 0.0, p]);
    let mut out = Vec::new();
    for du in shifts(periods.0) {
        for dv in shifts(periods.1) {
            out.push(Point3::new(corner.x + du, corner.y + dv, corner.z));
        }
    }
    out
}

/// `result` with every boundary gap cut out of an open record's
/// windows (see the module docs). Closed records keep their windows:
/// their circular per-interval classification already sees gaps.
pub(super) fn cut_boundary_gaps(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    closed: bool,
    result: UvClipResult,
) -> UvClipResult {
    let UvClipResult::Windows(windows) = result else {
        return result;
    };
    let Some(periods) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.periodicity())
    else {
        return UvClipResult::Windows(windows);
    };
    if closed {
        return UvClipResult::Windows(windows);
    }
    let Some(boundary) = crate::clip::face_uv_boundary(brep, face) else {
        return UvClipResult::Windows(windows);
    };
    let band = super::clip_uv::uv_boundary_band(brep, curve_uv);
    // A periodic chart's corners may sit a period off the record's
    // branch: each corner is tried at its translates.
    let corners: Vec<_> = boundary
        .iter()
        .flat_map(|(_, pc)| translates(pc.point_at(pc.parameter_range().0), periods))
        .collect();
    let mut out = Vec::new();
    for (a, b) in windows {
        let mut cuts: Vec<f64> = corners
            .iter()
            .filter_map(|q| {
                let t = curve_uv.closest_point(q)?;
                let p = curve_uv.point_at(t);
                let d = ((p.x - q.x).powi(2) + (p.y - q.y).powi(2)).sqrt();
                (d < band && t > a && t < b).then_some(t)
            })
            .collect();
        if cuts.is_empty() {
            out.push((a, b));
            continue;
        }
        cuts.sort_by(f64::total_cmp);
        let mut start = a;
        let mut run: Option<(f64, f64)> = None;
        for end in cuts.into_iter().chain(std::iter::once(b)) {
            if end - start > (b - a) * 1e-9 && is_gap(brep, face, curve_uv, (start, end)) {
                if let Some(w) = run.take() {
                    out.push(w);
                }
            } else {
                run = Some(run.map_or((start, end), |(s, _)| (s, end)));
            }
            start = end;
        }
        if let Some(w) = run {
            out.push(w);
        }
    }
    if out.is_empty() {
        UvClipResult::Outside
    } else {
        UvClipResult::Windows(out)
    }
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::transform::Transform3;
    use popkern_core::geometry::{Point3, Vec3};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    use super::*;

    /// A record riding a periodic face's boundary ruling past both its
    /// ends keeps only the stretch on the face: the rest is a gap of the
    /// trim however the tuned window reached it (PK-1027: the arm's
    /// groove cylinder, v 0..8, kept its boundary-line record over v
    /// -150..60, and the slit beyond the groove's end declined the
    /// complex).
    #[test]
    fn a_record_along_a_cylinder_boundary_ruling_keeps_only_the_face() {
        let tol = Tolerance::default();
        // An r = 4.25 groove 8 long in a block's top, its centre 1.502
        // above the top: the groove meets z = 10 at x = ±3.976.
        let cutter = BRep::cylinder(4.25, 8.0, tol)
            .expect("cutter")
            .transform(&Transform3::rotation_x(-std::f64::consts::FRAC_PI_2))
            .expect("cutter along y")
            .transform(&Transform3::translation(0.0, -5.0, 11.502))
            .expect("placed cutter");
        let block = BRep::prism(
            &[
                Point3::new(-10.0, -20.0, 0.0),
                Point3::new(10.0, -20.0, 0.0),
                Point3::new(10.0, 20.0, 0.0),
                Point3::new(-10.0, 20.0, 0.0),
            ],
            Vec3::new(0.0, 0.0, 1.0),
            10.0,
            tol,
        )
        .expect("block");
        let grooved = crate::boolean(crate::BoolOperation::Subtraction, block, cutter, tol, None)
            .expect("grooved block");
        let groove = grooved
            .faces
            .iter()
            .find(|(_, f)| {
                f.surface
                    .as_ref()
                    .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
            })
            .map(|(h, _)| h)
            .expect("groove face");
        // The face's chart: its boundary ruling u0 and its v span.
        let boundary = crate::clip::face_uv_boundary(&grooved, groove).expect("boundary");
        let corners: Vec<Point3> = boundary
            .iter()
            .map(|(_, pc)| pc.point_at(pc.parameter_range().0))
            .collect();
        let u0 = corners.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let (v_lo, v_hi) = corners
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                (lo.min(p.y), hi.max(p.y))
            });
        assert!(
            (v_hi - v_lo - 8.0).abs() < 1e-6,
            "groove v span {v_lo}..{v_hi}"
        );
        // The record along u0, 10 beyond the face either way.
        let (a, b) = (v_lo - 10.0, v_hi + 10.0);
        let record = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(u0, a, 0.0), Point3::new(u0, b, 0.0)],
            weights: vec![1.0; 2],
        };
        let UvClipResult::Windows(kept) = cut_boundary_gaps(
            &grooved,
            groove,
            &record,
            false,
            UvClipResult::Windows(vec![(0.0, 1.0)]),
        ) else {
            panic!("the record keeps a window");
        };
        let on_face = ((v_lo - a) / (b - a), (v_hi - a) / (b - a));
        assert_eq!(kept.len(), 1, "{kept:?}");
        assert!(
            (kept[0].0 - on_face.0).abs() < 1e-6 && (kept[0].1 - on_face.1).abs() < 1e-6,
            "kept {kept:?}, the face spans {on_face:?}"
        );
    }
}
