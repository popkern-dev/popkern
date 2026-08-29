//! Exact polynomial restriction of degree-9 boundary rows.

use popkern_core::geometry::Point3;

/// Exact restriction of a degree-9 Bezier row to `[start, end]`
/// (oriented start → end), via two de Casteljau splits.
/// Polynomial subdivision is exact, so the returned row traces
/// the same curve segment to floating-point roundoff.
pub(super) fn bezier_row_segment(row: &[Point3; 10], start: f64, end: f64) -> [Point3; 10] {
    let lo = start.min(end);
    let hi = start.max(end);
    let (_, tail) = de_casteljau_split_row(row, lo);
    let local = if (1.0 - lo).abs() < 1.0e-12 {
        0.0
    } else {
        (hi - lo) / (1.0 - lo)
    };
    let (mut segment, _) = de_casteljau_split_row(&tail, local);
    if start > end {
        segment.reverse();
    }
    segment
}

fn de_casteljau_split_row(row: &[Point3; 10], t: f64) -> ([Point3; 10], [Point3; 10]) {
    let mut left = [Point3::new(0.0, 0.0, 0.0); 10];
    let mut right = [Point3::new(0.0, 0.0, 0.0); 10];
    let mut working = *row;
    left[0] = working[0];
    right[9] = working[9];
    for level in 1..10 {
        for index in 0..(10 - level) {
            working[index] = working[index] + (working[index + 1] - working[index]) * t;
        }
        left[level] = working[0];
        right[9 - level] = working[9 - level];
    }
    (left, right)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_row() -> [Point3; 10] {
        let mut row = [Point3::new(0.0, 0.0, 0.0); 10];
        for (index, point) in row.iter_mut().enumerate() {
            let t = index as f64;
            *point = Point3::new(t, t * t * 0.1, (t - 4.0) * (t - 7.0) * 0.05);
        }
        row
    }

    fn eval_row(row: &[Point3; 10], t: f64) -> Point3 {
        let mut working = *row;
        for level in 1..10 {
            for index in 0..(10 - level) {
                working[index] = working[index] + (working[index + 1] - working[index]) * t;
            }
        }
        working[0]
    }

    #[test]
    fn segment_matches_source_curve() {
        let row = sample_row();
        let (start, end) = (0.25, 0.75);
        let segment = bezier_row_segment(&row, start, end);
        for sample in 0..=8 {
            let local = sample as f64 / 8.0;
            let global = start + (end - start) * local;
            let expected = eval_row(&row, global);
            let actual = eval_row(&segment, local);
            assert!(
                expected.distance_to(&actual) < 1.0e-12,
                "segment diverged at local {local}: {expected:?} vs {actual:?}"
            );
        }
    }

    #[test]
    fn reversed_segment_walks_backwards() {
        let row = sample_row();
        let segment = bezier_row_segment(&row, 0.5, 0.0);
        let expected_start = eval_row(&row, 0.5);
        let actual_start = eval_row(&segment, 0.0);
        assert!(expected_start.distance_to(&actual_start) < 1.0e-12);
        let expected_end = eval_row(&row, 0.0);
        let actual_end = eval_row(&segment, 1.0);
        assert!(expected_end.distance_to(&actual_end) < 1.0e-12);
    }
}
