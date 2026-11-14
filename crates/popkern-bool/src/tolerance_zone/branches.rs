//! Records of one face pair (plan §4.3a, §5.4).
//!
//! Two records of one face pair are two branches of that pair's
//! intersection, or two pieces of one branch cut at a hit, a seam or a
//! band lattice. A record end is not given an interior incidence on
//! another record of its own pair (`incidence`): where two branches of
//! one pair come within the tolerance of each other, the face
//! arrangement's meeting is the vertex (`meetings`), one point on both.
//!
//! The rule that kept such a pass as two vertices (two branches the SSI
//! marched through) is retired (PK-1130): within the tolerance two
//! points are one, and no served boolean depended on it (measured with
//! the vertex table on: no record's interior met another of its pair
//! outside the rule's own fixture).

use crate::intersect::IntersectionRecord;

/// Whether records `ra` and `rb` are distinct records of one face pair.
pub(crate) fn same_pair(records: &[IntersectionRecord], ra: usize, rb: usize) -> bool {
    match (records.get(ra), records.get(rb)) {
        (Some(a), Some(b)) => ra != rb && a.face_a == b.face_a && a.face_b == b.face_b,
        _ => false,
    }
}
