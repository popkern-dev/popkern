//! Consequences of a resolution for the records (plan §4.3a): a record
//! whose two ends resolved into one vertex is a micro edge when its
//! whole extent stays within that vertex's tolerance (it collapses into
//! the vertex), a spur when it runs out and back along one line (it
//! bounds nothing, so it carries no boundary either), and a loop
//! through the vertex otherwise. A closed record's ends meet by
//! construction and are none of these.
//!
//! Measured (suite 19, table on): the only records whose ends met
//! were spurs, 1.8e-4 out and back along a cylinder ruling where a
//! wall plane grazes a rim fillet; their extent alone called them
//! loops.
//!
//! **Stage.** Measured in shadow beside the graph's point-contact guard
//! (`insert_intersection_curves_at`), which skips a record whose ends
//! share a node and whose extent stays within 500 tolerances (capped by
//! the model size): the table decides by the vertex's own tolerance
//! instead, and the counts show which records the two disagree on
//! before either decides the graph.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use super::candidates::{CandidateVertex, Incidence, SourceKind};
use super::resolve::Resolution;
use crate::intersect::IntersectionRecord;

/// Samples along a record for its extent and width.
const SHAPE_SAMPLES: usize = 32;

/// Each record's (start, end) candidate, `None` where it has none.
pub(crate) fn record_end_candidates(
    cands: &[CandidateVertex],
    records: usize,
) -> Vec<[Option<usize>; 2]> {
    let mut ends: Vec<[Option<usize>; 2]> = vec![[None; 2]; records];
    for (i, c) in cands.iter().enumerate() {
        if c.source != SourceKind::RecordEndpoint {
            continue;
        }
        if let Some(&Incidence::Record {
            record,
            end: Some(end),
            ..
        }) = c.incidence.first()
            && let Some(slot) = ends.get_mut(record)
        {
            slot[usize::from(end)] = Some(i);
        }
    }
    ends
}

/// A record whose ends resolved into one vertex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CollapsedEnds {
    /// The record.
    pub record: usize,
    /// The farthest the record strays from the vertex.
    pub extent: f64,
    /// The farthest the record strays from the line through the vertex
    /// and the record's farthest point: zero for an out-and-back spur,
    /// on either side of the vertex.
    pub width: f64,
    /// The vertex's tolerance.
    pub tolerance: f64,
}

/// What a record whose ends met at one vertex becomes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Collapse {
    /// The whole record lies within the vertex: it collapses into it.
    MicroEdge,
    /// The record runs out and back along one line through the
    /// vertex: it bounds nothing, so it carries no boundary.
    Spur,
    /// The record encloses width: a loop through the vertex.
    Loop,
}

impl CollapsedEnds {
    /// The record's fate, by the vertex's own tolerance.
    pub(crate) fn kind(&self) -> Collapse {
        if self.extent <= self.tolerance {
            Collapse::MicroEdge
        } else if self.width <= self.tolerance {
            Collapse::Spur
        } else {
            Collapse::Loop
        }
    }
}

/// The farthest `curve` strays from `anchor`, and from the line through
/// `anchor` and that farthest point.
fn shape(curve: &NurbsCurve, anchor: Point3) -> (f64, f64) {
    let (t0, t1) = curve.parameter_range();
    let samples: Vec<Point3> = (0..=SHAPE_SAMPLES)
        .map(|k| curve.point_at(t0 + (t1 - t0) * k as f64 / SHAPE_SAMPLES as f64))
        .collect();
    let (extent, far) =
        samples
            .iter()
            .map(|q| (q.distance_to(&anchor), *q))
            .fold(
                (0.0_f64, anchor),
                |best, c| if c.0 > best.0 { c } else { best },
            );
    let width = samples
        .iter()
        .map(|q| line_distance(q, &anchor, &far))
        .fold(0.0_f64, f64::max);
    (extent, width)
}

/// What a record on `curve` whose ends met at the vertex `anchor`
/// becomes, by the vertex's tolerance.
pub(crate) fn collapse_of(curve: &NurbsCurve, anchor: Point3, tolerance: f64) -> Collapse {
    let (extent, width) = shape(curve, anchor);
    CollapsedEnds {
        record: 0,
        extent,
        width,
        tolerance,
    }
    .kind()
}

/// Distance from `p` to the line through `a` and `b`.
fn line_distance(p: &Point3, a: &Point3, b: &Point3) -> f64 {
    let ab = *b - *a;
    let len2 = ab.length_squared();
    if len2 <= f64::MIN_POSITIVE {
        return p.distance_to(a);
    }
    let s = (*p - *a).dot(&ab) / len2;
    p.distance_to(&(*a + ab * s))
}

/// Every open record of `records` whose two ends resolved into one
/// vertex.
pub(crate) fn collapsed_ends(
    cands: &[CandidateVertex],
    res: &Resolution,
    records: &[IntersectionRecord],
) -> Vec<CollapsedEnds> {
    let mut out = Vec::new();
    for (ri, [s, e]) in record_end_candidates(cands, records.len())
        .into_iter()
        .enumerate()
    {
        let (Some(s), Some(e)) = (s, e) else {
            continue;
        };
        if records[ri].curve.is_closed {
            continue;
        }
        let group = res.group_of[s];
        if group != res.group_of[e] {
            continue;
        }
        let vertex = &res.groups[group];
        let anchor = cands[vertex.anchor].position;
        let (extent, width) = shape(&records[ri].curve.curve_3d, anchor);
        out.push(CollapsedEnds {
            record: ri,
            extent,
            width,
            tolerance: vertex.tolerance,
        });
    }
    out
}
