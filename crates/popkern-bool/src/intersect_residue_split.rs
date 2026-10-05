//! Verdict-uniform trim residues (PK-553).
//!
//! A trim residue (PK-465) is judged ONCE, at its midpoint, against
//! the other operand: its two flanks classify identically at every
//! point, so one verdict decides whether the slit partitions a kept
//! region (informationless — dropped) or a discarded one. That holds
//! only while the residue stays on one side of the other operand's
//! boundary. A residue can run THROUGH the other solid's wall: the
//! keyed bore's key-wall plane crosses the boss top from the key slot,
//! through the bore (discarded) and on to the boss rim (kept
//! material). Judged at its midpoint in the bore, the whole residue
//! survived and slit the kept boss top from the bore to the rim —
//! hidden only while the coarse boundary conversion (PK-553) left its
//! rim end 5e-5 off the rim, so the walk filtered it as an interior
//! graze.
//!
//! Before the keep/drop judgement, each residue is therefore cut
//! where it crosses the other operand's boundary AND its verdict
//! changes there. Every such crossing on the residue's owner face F
//! lies on an ordinary record F × G (G the crossed face), so the
//! crossings are the residue's intersections, in F's uv chart, with
//! the uv curves of the ordinary records on F (tried one period either
//! way on a periodic chart, pruned by uv and 3D hull overlap — a full
//! curve × face-surface scan put a torus × cylinder union from 5.8 s
//! to 33 s). A crossing becomes a cut only when the pieces on its two
//! sides take definite, different verdicts (a graze or a crossing of
//! an inner wall leaves the verdict alone and is not cut). The pieces
//! are then judged individually by the existing rule. A residue with
//! no such cut is left untouched — no behaviour change for the
//! uniform classes the residue gates were tuned on.
//!
//! Why exact crossings, not sampled verdicts: every residue lies ON
//! its pair face's surface and starts at that face's trim end, where
//! ray parity is knife-edge (measured on the keyed bore: a
//! bore-interior point 4.5e-3 from the key wall's trim read Outside),
//! and a cut must land on the other operand's wall to meet the
//! ordinary record there.

use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::PointClassification;
use crate::intersect::IntersectionRecord;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

/// Crossings closer than this fraction of the residue's parameter span
/// to an end, or to each other, are merged into it (no sliver pieces).
/// The residue's junction walls touch it exactly at its ends, so the
/// wall prefilter looks only at the interior this margin leaves.
const MIN_PIECE_FRACTION: f64 = 1e-3;

/// The other operand of a residue owned by `owner_is_a`.
pub(crate) struct OtherOperand<'a> {
    pub(crate) brep_a: &'a BRep,
    pub(crate) brep_b: &'a BRep,
    pub(crate) aabbs_a: &'a [(FaceHandle, Aabb)],
    pub(crate) aabbs_b: &'a [(FaceHandle, Aabb)],
}

impl OtherOperand<'_> {
    /// Inside (`true`) / outside (`false`) the other operand; `None` on
    /// its boundary (within `100 x` the linear tolerance of one of its
    /// faces, the residue filter's boundary gate) or inconclusive.
    fn verdict(
        &self,
        owner_is_a: bool,
        point: &popkern_core::geometry::Point3,
        tolerance: &Tolerance,
    ) -> Option<bool> {
        let (brep, aabbs) = if owner_is_a {
            (self.brep_b, self.aabbs_b)
        } else {
            (self.brep_a, self.aabbs_a)
        };
        if on_operand_boundary(point, brep, aabbs, 100.0 * tolerance.linear(), tolerance) {
            return None;
        }
        // A FINAL operand: final-body gates (see drop_kept_flank_residues).
        match crate::classify::point_in_solid_final(point, brep, aabbs, tolerance) {
            Ok(PointClassification::Inside) => Some(true),
            Ok(PointClassification::Outside) => Some(false),
            _ => None,
        }
    }
}

/// Cut every residue in `records` at the parameters where its
/// classification against the other operand changes (see the module
/// docs). Ordinary records and verdict-uniform residues are kept as
/// they are; a split residue is replaced in place by its pieces, in
/// parameter order.
pub(crate) fn split_residues_at_verdict_changes(
    records: &mut Vec<IntersectionRecord>,
    other: &OtherOperand<'_>,
    tolerance: &Tolerance,
) {
    let trace = popkern_core::trace::trace_enabled("boolean");
    // A residue on owner face F can cross the other operand's boundary
    // only where a face G of that boundary meets F — ON an ordinary
    // record F × G. So G is a candidate only when one of its F × G
    // records' control hulls overlaps the residue's (scanning every
    // AABB-overlapping face with the general solver put a torus ×
    // cylinder union from 5.8 s to 33 s).
    let gate = 100.0 * tolerance.linear();
    let mut walls: std::collections::HashMap<(bool, FaceHandle), Vec<(Hull, NurbsCurve)>> =
        std::collections::HashMap::new();
    for r in records.iter().filter(|r| r.residue_owner_a.is_none()) {
        let hull = Hull::of(&r.curve.curve_3d, gate);
        walls
            .entry((true, r.face_a))
            .or_default()
            .push((hull, r.curve.curve_uv1.clone()));
        walls
            .entry((false, r.face_b))
            .or_default()
            .push((hull, r.curve.curve_uv2.clone()));
    }
    let mut out = Vec::with_capacity(records.len());
    for record in records.drain(..) {
        let Some(owner_is_a) = record.residue_owner_a else {
            out.push(record);
            continue;
        };
        let owner = if owner_is_a {
            record.face_a
        } else {
            record.face_b
        };
        let reach = Hull::of(&record.curve.curve_3d, gate);
        let candidates: Vec<&NurbsCurve> = walls
            .get(&(owner_is_a, owner))
            .into_iter()
            .flatten()
            .filter(|(hull, _)| hull.overlaps(&reach))
            .map(|(_, uv)| uv)
            .collect();
        let cuts = verdict_transitions(&record, owner_is_a, &candidates, other, tolerance);
        if cuts.is_empty() {
            out.push(record);
            continue;
        }
        let (t0, t1) = record.curve.curve_3d.parameter_range();
        let bounds: Vec<f64> = std::iter::once(t0)
            .chain(cuts.iter().copied())
            .chain(std::iter::once(t1))
            .collect();
        let pieces: Vec<IntersectionRecord> = bounds
            .windows(2)
            .filter_map(|w| {
                let curve = restrict(&record.curve, w[0], w[1])?;
                Some(IntersectionRecord {
                    curve,
                    ..record.clone()
                })
            })
            .collect();
        if pieces.len() != bounds.len() - 1 {
            // A piece failed to extract: keep the residue whole rather
            // than drop part of a partition.
            out.push(record);
            continue;
        }
        if trace {
            let p = |t: f64| record.curve.curve_3d.point_at(t);
            let cut_points: Vec<_> = cuts
                .iter()
                .map(|&t| {
                    let q = p(t);
                    (q.x, q.y, q.z)
                })
                .collect();
            eprintln!(
                "  [residue-split] A{:?} x B{:?} owner={} cut at {cut_points:.4?}",
                record.face_a,
                record.face_b,
                if owner_is_a { "A" } else { "B" },
            );
        }
        out.extend(pieces);
    }
    *records = out;
}

/// A padded axis-aligned box around a NURBS curve's control points
/// (which contain the curve).
#[derive(Clone, Copy)]
struct Hull {
    lo: popkern_core::geometry::Point3,
    hi: popkern_core::geometry::Point3,
}

impl Hull {
    fn of(curve: &NurbsCurve, pad: f64) -> Self {
        use popkern_core::geometry::Point3;
        let (lo, hi) = curve.control_points.iter().fold(
            (
                Point3::new(f64::MAX, f64::MAX, f64::MAX),
                Point3::new(f64::MIN, f64::MIN, f64::MIN),
            ),
            |(lo, hi), p| {
                (
                    Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z)),
                    Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z)),
                )
            },
        );
        Self {
            lo: Point3::new(lo.x - pad, lo.y - pad, lo.z - pad),
            hi: Point3::new(hi.x + pad, hi.y + pad, hi.z + pad),
        }
    }

    fn overlaps(&self, other: &Hull) -> bool {
        self.lo.x <= other.hi.x
            && other.lo.x <= self.hi.x
            && self.lo.y <= other.hi.y
            && other.lo.y <= self.hi.y
            && self.lo.z <= other.hi.z
            && other.lo.z <= self.hi.z
    }
}

/// Parameters, ascending, where `record`'s residue crosses the other
/// operand's boundary and its verdict changes.
fn verdict_transitions(
    record: &IntersectionRecord,
    owner_is_a: bool,
    walls: &[&NurbsCurve],
    other: &OtherOperand<'_>,
    tolerance: &Tolerance,
) -> Vec<f64> {
    let curve = &record.curve;
    let (t0, t1) = curve.curve_3d.parameter_range();
    let span = t1 - t0;
    if span <= 0.0 {
        return Vec::new();
    }
    let mut crossings = boundary_crossings(record, owner_is_a, walls, other, tolerance);
    crossings.retain(|&t| t - t0 > MIN_PIECE_FRACTION * span && t1 - t > MIN_PIECE_FRACTION * span);
    crossings.sort_by(f64::total_cmp);
    crossings.dedup_by(|b, a| *b - *a <= MIN_PIECE_FRACTION * span);
    if crossings.is_empty() {
        return Vec::new();
    }
    let bounds: Vec<f64> = std::iter::once(t0)
        .chain(crossings.iter().copied())
        .chain(std::iter::once(t1))
        .collect();
    let verdicts: Vec<Option<bool>> = bounds
        .windows(2)
        .map(|w| {
            let mid = curve.curve_3d.point_at(0.5 * (w[0] + w[1]));
            other.verdict(owner_is_a, &mid, tolerance)
        })
        .collect();
    // Cut between definite, different verdicts only; carry the last
    // definite verdict across undecided pieces.
    let mut cuts = Vec::new();
    let mut held: Option<bool> = None;
    for (k, verdict) in verdicts.iter().enumerate() {
        let Some(v) = *verdict else {
            continue;
        };
        if held.is_some_and(|h| h != v) {
            cuts.push(bounds[k]);
        }
        held = Some(v);
    }
    cuts
}

/// Parameters where the residue crosses one of the `walls` — the uv
/// curves, on the residue's owner face, of the ordinary records there
/// (every crossing of the other operand's boundary lies on one). On a
/// periodic owner chart each wall is also tried one period either way.
fn boundary_crossings(
    record: &IntersectionRecord,
    owner_is_a: bool,
    walls: &[&NurbsCurve],
    other: &OtherOperand<'_>,
    tolerance: &Tolerance,
) -> Vec<f64> {
    use popkern_core::geometry::curve::CurveGeometry;
    use popkern_core::geometry::intersection::{CurveIntersection, intersect_curves};
    use popkern_core::geometry::surface::SurfaceEval as _;
    let (owner_brep, owner_face, uv) = if owner_is_a {
        (other.brep_a, record.face_a, &record.curve.curve_uv1)
    } else {
        (other.brep_b, record.face_b, &record.curve.curve_uv2)
    };
    let periods = owner_brep
        .faces
        .get(owner_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.parameter_periodicity())
        .unwrap_or((None, None));
    let shifts = |period: Option<f64>| match period {
        Some(p) => vec![0.0, p, -p],
        None => vec![0.0],
    };
    let residue = CurveGeometry::Nurbs(uv.clone());
    let range = uv.parameter_range();
    // uv hulls (control points contain the curve): a shifted wall whose
    // hull misses the residue's cannot cross it.
    let pad = 1e-9;
    let span = range.1 - range.0;
    let Some(interior) = crate::clip::clip_curve_to_crossings(
        uv,
        range.0 + MIN_PIECE_FRACTION * span,
        range.1 - MIN_PIECE_FRACTION * span,
    ) else {
        return Vec::new();
    };
    let reach = Hull::of(&interior, pad);
    let mut hits = Vec::new();
    for wall in walls {
        for du in shifts(periods.0) {
            for dv in shifts(periods.1) {
                let mut shifted = (*wall).clone();
                for p in &mut shifted.control_points {
                    p.x += du;
                    p.y += dv;
                }
                if !Hull::of(&shifted, pad).overlaps(&reach) {
                    continue;
                }
                let wall_range = shifted.parameter_range();
                if wall_range.1 <= wall_range.0 {
                    continue;
                }
                let Ok(result) = intersect_curves(
                    &residue,
                    range,
                    &CurveGeometry::Nurbs(shifted),
                    wall_range,
                    tolerance,
                ) else {
                    continue;
                };
                hits.extend(result.intersections.iter().filter_map(|i| match i {
                    CurveIntersection::Point { t1, .. } => Some(*t1),
                    CurveIntersection::Overlap { .. } => None,
                }));
            }
        }
    }
    hits
}

/// `curve` restricted to `[t_enter, t_exit]` on all three carriers.
fn restrict(curve: &CurveOnSurface, t_enter: f64, t_exit: f64) -> Option<CurveOnSurface> {
    use crate::clip::clip_curve_to_crossings as clip;
    Some(CurveOnSurface {
        curve_3d: clip(&curve.curve_3d, t_enter, t_exit)?,
        curve_uv1: clip(&curve.curve_uv1, t_enter, t_exit)?,
        curve_uv2: clip(&curve.curve_uv2, t_enter, t_exit)?,
        is_closed: false,
    })
}

/// Whether `point` lies within `gate` of a face of `brep`, inside that
/// face's trim region.
pub(crate) fn on_operand_boundary(
    point: &popkern_core::geometry::Point3,
    brep: &BRep,
    aabbs: &[(FaceHandle, Aabb)],
    gate: f64,
    tolerance: &Tolerance,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval as _;
    aabbs.iter().any(|(fh, aabb)| {
        let inside_box = point.x >= aabb.min.x - gate
            && point.x <= aabb.max.x + gate
            && point.y >= aabb.min.y - gate
            && point.y <= aabb.max.y + gate
            && point.z >= aabb.min.z - gate
            && point.z <= aabb.max.z + gate;
        if !inside_box {
            return false;
        }
        let Some(sof) = brep.faces.get(*fh).and_then(|f| f.surface.as_ref()) else {
            return false;
        };
        let Some((u, v)) = sof.geometry.closest_point(point) else {
            return false;
        };
        sof.geometry.point_at(u, v).distance_to(point) < gate
            && crate::split_inner_loop::point_inside_face_region_uv(
                brep,
                *fh,
                popkern_core::geometry::Point3::new(u, v, 0.0),
                tolerance,
            )
    })
}

#[cfg(test)]
#[path = "intersect_residue_split_tests.rs"]
mod tests;
