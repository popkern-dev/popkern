//! Static k-d tree over 3-D points answering axis-aligned box queries.
//!
//! Built once over a point set and queried with many boxes. The
//! orientation votes use it to find, for each B-rep edge, the boundary
//! segments of a face that start near the edge, instead of testing
//! every segment against every edge, which was quadratic on a thread
//! rod (PK-794).

use popkern_core::geometry::Point3;

/// Points stored as an implicit k-d tree: each range's median on the
/// split axis (x, y, z cycling with depth) sits at the range's middle.
pub(crate) struct PointBoxIndex {
    nodes: Vec<(Point3, usize)>,
    /// Points with a non-finite coordinate. No box contains them, so
    /// every query returns them and the caller's exact test decides.
    non_finite: Vec<usize>,
}

impl PointBoxIndex {
    /// Index `points`; a query returns positions in this sequence.
    pub(crate) fn build(points: impl IntoIterator<Item = Point3>) -> PointBoxIndex {
        let mut nodes = Vec::new();
        let mut non_finite = Vec::new();
        for (i, p) in points.into_iter().enumerate() {
            if p.x.is_finite() && p.y.is_finite() && p.z.is_finite() {
                nodes.push((p, i));
            } else {
                non_finite.push(i);
            }
        }
        split(&mut nodes, 0);
        PointBoxIndex { nodes, non_finite }
    }

    /// Append to `out` the index of every point inside the closed box
    /// `[lo, hi]`, plus every non-finite point.
    pub(crate) fn query(&self, lo: Point3, hi: Point3, out: &mut Vec<usize>) {
        out.extend_from_slice(&self.non_finite);
        query(&self.nodes, 0, lo, hi, out);
    }
}

fn coord(p: &Point3, axis: usize) -> f64 {
    match axis {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    }
}

fn split(nodes: &mut [(Point3, usize)], depth: usize) {
    if nodes.len() <= 1 {
        return;
    }
    let axis = depth % 3;
    let mid = nodes.len() / 2;
    nodes.select_nth_unstable_by(mid, |a, b| coord(&a.0, axis).total_cmp(&coord(&b.0, axis)));
    let (left, right) = nodes.split_at_mut(mid);
    split(left, depth + 1);
    split(&mut right[1..], depth + 1);
}

fn query(nodes: &[(Point3, usize)], depth: usize, lo: Point3, hi: Point3, out: &mut Vec<usize>) {
    if nodes.is_empty() {
        return;
    }
    let axis = depth % 3;
    let mid = nodes.len() / 2;
    let (p, i) = nodes[mid];
    if (lo.x..=hi.x).contains(&p.x) && (lo.y..=hi.y).contains(&p.y) && (lo.z..=hi.z).contains(&p.z)
    {
        out.push(i);
    }
    let c = coord(&p, axis);
    if coord(&lo, axis) <= c {
        query(&nodes[..mid], depth + 1, lo, hi, out);
    }
    if c <= coord(&hi, axis) {
        query(&nodes[mid + 1..], depth + 1, lo, hi, out);
    }
}

#[cfg(test)]
mod tests {
    use super::PointBoxIndex;
    use popkern_core::geometry::Point3;

    /// Every query returns exactly the points a scan finds in the box.
    #[test]
    fn box_queries_match_a_scan() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut unit = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        // Repeated coordinates exercise ties at the split medians.
        let points: Vec<Point3> = (0..2000)
            .map(|_| {
                Point3::new(
                    (unit() * 20.0).round() * 0.5,
                    unit() * 3.0,
                    (unit() * 4.0).round(),
                )
            })
            .collect();
        let index = PointBoxIndex::build(points.iter().copied());
        for _ in 0..300 {
            let (a, b) = (
                Point3::new(unit() * 10.0, unit() * 3.0, unit() * 4.0),
                Point3::new(unit() * 10.0, unit() * 3.0, unit() * 4.0),
            );
            let lo = Point3::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z));
            let hi = Point3::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z));
            let mut got = Vec::new();
            index.query(lo, hi, &mut got);
            got.sort_unstable();
            let want: Vec<usize> = (0..points.len())
                .filter(|&i| {
                    let p = points[i];
                    (lo.x..=hi.x).contains(&p.x)
                        && (lo.y..=hi.y).contains(&p.y)
                        && (lo.z..=hi.z).contains(&p.z)
                })
                .collect();
            assert_eq!(got, want);
        }
    }

    /// A non-finite point is in every answer.
    #[test]
    fn non_finite_points_are_always_returned() {
        let index = PointBoxIndex::build([
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(f64::NAN, 0.0, 0.0),
            Point3::new(5.0, 5.0, 5.0),
        ]);
        let mut got = Vec::new();
        index.query(
            Point3::new(4.0, 4.0, 4.0),
            Point3::new(6.0, 6.0, 6.0),
            &mut got,
        );
        got.sort_unstable();
        assert_eq!(got, vec![1, 2]);
    }
}
