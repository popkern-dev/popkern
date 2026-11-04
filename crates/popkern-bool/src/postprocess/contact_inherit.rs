//! Contacts a boolean result inherits from its operands (PK-986).
//!
//! A body served as manifold pieces meeting at recorded contacts
//! (`BRep::topology_contacts`) touches itself on purpose. A later
//! boolean on it rebuilds the topology, so the records' handles die,
//! and the result-side zero-thickness check would read the same touch
//! as a new non-manifold contact and decline. A contact of the result
//! that lies on a contact an operand recorded is that operand's
//! deliberate touch, carried: it is recorded again on the result and
//! excused. Any other contact is new and follows the boolean's policy.
//!
//! "Lies on" is geometric, within a band of ten linear tolerances
//! beyond the sampled locus's own chord error: every test point of the
//! result's contact (its location, its vertices, its edges' ends and
//! midpoints) lies on the locus of an operand contact (that contact's
//! vertices and sampled edges, or every face it names).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, ContactEntity, EdgeHandle, FaceHandle, TopologyContact};

use crate::error::{ContactKind, NonManifoldContact};

/// Samples per contact edge (each span's chord error joins the band).
const EDGE_SAMPLES: usize = 32;
/// The band around an operand's contact locus, in linear tolerances.
const BAND_TOLERANCES: f64 = 10.0;

/// One operand contact's locus: points, sampled edge spans with their
/// chord error, and the faces it names.
struct Locus<'a> {
    body: &'a BRep,
    points: Vec<Point3>,
    spans: Vec<(Point3, Point3, f64)>,
    faces: Vec<FaceHandle>,
}

fn edge_point(body: &BRep, e: EdgeHandle, s: f64) -> Option<Point3> {
    let edge = body.edges.get(e)?;
    if let Some(c) = &edge.curve {
        return Some(c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * s));
    }
    let he = body.half_edges.get(edge.half_edge)?;
    let a = body.vertices.get(he.origin)?.position;
    let b = body
        .vertices
        .get(body.half_edges.get(he.twin)?.origin)?
        .position;
    Some(a + (b - a) * s)
}

fn locus<'a>(body: &'a BRep, contact: &TopologyContact) -> Locus<'a> {
    let mut l = Locus {
        body,
        points: vec![contact.location],
        spans: Vec::new(),
        faces: Vec::new(),
    };
    for entity in &contact.entities {
        match *entity {
            ContactEntity::Vertex(v) => l.points.extend(body.vertices.get(v).map(|v| v.position)),
            ContactEntity::Edge(e) => {
                let at = |k: usize| edge_point(body, e, k as f64 / EDGE_SAMPLES as f64);
                for k in 0..EDGE_SAMPLES {
                    let (Some(a), Some(b)) = (at(k), at(k + 1)) else {
                        continue;
                    };
                    let mid = edge_point(body, e, (k as f64 + 0.5) / EDGE_SAMPLES as f64);
                    let sag = mid.map_or(0.0, |m| m.distance_to(&(a + (b - a) * 0.5)));
                    l.spans.push((a, b, sag));
                }
            }
            ContactEntity::Face(f) => l.faces.push(f),
        }
    }
    l
}

fn segment_distance(p: &Point3, a: &Point3, b: &Point3) -> f64 {
    let ab = *b - *a;
    let len2 = ab.dot(&ab);
    let t = if len2 > 0.0 {
        ((*p - *a).dot(&ab) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    p.distance_to(&(*a + ab * t))
}

impl Locus<'_> {
    fn holds(&self, p: &Point3, band: f64, boxes: &[(FaceHandle, Aabb)], tol: &Tolerance) -> bool {
        if self.points.iter().any(|q| q.distance_to(p) <= band)
            || self
                .spans
                .iter()
                .any(|(a, b, sag)| segment_distance(p, a, b) <= band + sag)
        {
            return true;
        }
        // A touch of faces (shells meeting) is known by its faces: the
        // point lies on every one of them.
        !self.faces.is_empty()
            && self.faces.iter().all(|f| {
                boxes.iter().any(|(h, aabb)| {
                    h == f && super::zero_thickness::on_face(p, self.body, *f, aabb, tol)
                })
            })
    }
}

/// The test points of a result contact: its location, its vertices, and
/// its edges' ends and midpoints.
fn test_points(result: &BRep, contact: &TopologyContact) -> Vec<Point3> {
    let mut points = vec![contact.location];
    for entity in &contact.entities {
        match *entity {
            ContactEntity::Vertex(v) => points.extend(result.vertices.get(v).map(|v| v.position)),
            ContactEntity::Edge(e) => {
                points.extend(
                    [0.0, 0.5, 1.0]
                        .iter()
                        .filter_map(|&s| edge_point(result, e, s)),
                );
            }
            ContactEntity::Face(_) => {}
        }
    }
    points
}

/// Whether some operand recorded contacts (the result may inherit them).
pub(crate) fn operands_record_contacts(operands: &[&BRep]) -> bool {
    operands.iter().any(|b| !b.topology_contacts.is_empty())
}

/// The policy to run an operation under: a declining operation over
/// operands with recorded contacts runs split, so the result records
/// its contacts and [`inherited_or_declined`] can sort them; `None`
/// when nothing is inherited (run as asked).
pub(crate) fn run_policy(
    policy: crate::options::NonManifoldPolicy,
    operands: &[&BRep],
) -> Option<crate::options::NonManifoldPolicy> {
    (policy == crate::options::NonManifoldPolicy::Decline && operands_record_contacts(operands))
        .then_some(crate::options::NonManifoldPolicy::SplitToManifold)
}

/// `result` (run under [`run_policy`]) when every contact it records was
/// recorded by an operand; else the decline naming the new ones.
pub(crate) fn inherited_or_declined(
    result: BRep,
    operands: &[&BRep],
    tolerance: &Tolerance,
) -> Result<BRep, crate::error::BoolError> {
    let contacts = uninherited(&result, operands, tolerance);
    if contacts.is_empty() {
        Ok(result)
    } else {
        Err(crate::error::BoolError::NonManifoldResult { contacts })
    }
}

/// The contacts recorded on `result` that no operand recorded: each one
/// with a test point on no operand contact's locus.
pub(crate) fn uninherited(
    result: &BRep,
    operands: &[&BRep],
    tolerance: &Tolerance,
) -> Vec<NonManifoldContact> {
    let band = BAND_TOLERANCES * tolerance.linear();
    let boxes: Vec<Vec<(FaceHandle, Aabb)>> = operands
        .iter()
        .map(|b| crate::adaptive_precision::body_face_aabbs(b))
        .collect();
    let loci: Vec<(usize, Locus<'_>)> = operands
        .iter()
        .enumerate()
        .flat_map(|(i, b)| b.topology_contacts.iter().map(move |c| (i, locus(b, c))))
        .collect();
    result
        .topology_contacts
        .iter()
        .filter(|contact| {
            // A result contact may join several operand records along
            // one line: each point lies on some recorded locus.
            !test_points(result, contact).iter().all(|p| {
                loci.iter()
                    .any(|(i, l)| l.holds(p, band, &boxes[*i], tolerance))
            })
        })
        .map(|c| NonManifoldContact {
            kind: match c.kind {
                popkern_core::topology::TopologyContactKind::Edge => ContactKind::Edge,
                popkern_core::topology::TopologyContactKind::Vertex => ContactKind::Vertex,
            },
            location: c.location,
        })
        .collect()
}
