//! Surface subdivision infrastructure.

use crate::geometry::Point3;
use serde::{Deserialize, Serialize};

/// Axis-aligned bounding box.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Aabb {
    /// Minimum corner of the bounding box.
    pub min: Point3,
    /// Maximum corner of the bounding box.
    pub max: Point3,
}

impl Aabb {
    /// Compute AABB from a set of points.
    pub fn from_points(points: &[Point3]) -> Option<Self> {
        if points.is_empty() {
            return None;
        }
        let mut min = points[0];
        let mut max = points[0];
        for p in points.iter().skip(1) {
            min = Point3::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
            max = Point3::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
        }
        Some(Self { min, max })
    }

    /// Check if two AABBs overlap with tolerance padding.
    pub fn overlaps(&self, other: &Aabb, padding: f64) -> bool {
        self.min.x - padding <= other.max.x
            && self.max.x + padding >= other.min.x
            && self.min.y - padding <= other.max.y
            && self.max.y + padding >= other.min.y
            && self.min.z - padding <= other.max.z
            && self.max.z + padding >= other.min.z
    }

    /// Diagonal length.
    pub fn diagonal(&self) -> f64 {
        (self.max - self.min).length()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_points_empty_returns_none() {
        assert!(Aabb::from_points(&[]).is_none());
    }

    #[test]
    fn from_points_single_point() {
        let p = Point3::new(1.0, 2.0, 3.0);
        let aabb = Aabb::from_points(&[p]).unwrap();
        assert_eq!(aabb.min.x, 1.0);
        assert_eq!(aabb.max.z, 3.0);
    }

    #[test]
    fn from_points_multiple_points() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 3.0),
            Point3::new(-1.0, 5.0, 1.0),
        ];
        let aabb = Aabb::from_points(&points).unwrap();
        assert_eq!(aabb.min.x, -1.0);
        assert_eq!(aabb.min.y, 0.0);
        assert_eq!(aabb.min.z, 0.0);
        assert_eq!(aabb.max.x, 1.0);
        assert_eq!(aabb.max.y, 5.0);
        assert_eq!(aabb.max.z, 3.0);
    }

    #[test]
    fn overlaps_touching_boxes_no_padding() {
        let a = Aabb {
            min: Point3::new(0.0, 0.0, 0.0),
            max: Point3::new(1.0, 1.0, 1.0),
        };
        let b = Aabb {
            min: Point3::new(1.0, 0.0, 0.0),
            max: Point3::new(2.0, 1.0, 1.0),
        };
        // Touching at x=1 — counts as overlapping (<=/>= comparisons)
        assert!(a.overlaps(&b, 0.0));
    }

    #[test]
    fn overlaps_separated_boxes_no_padding() {
        let a = Aabb {
            min: Point3::new(0.0, 0.0, 0.0),
            max: Point3::new(1.0, 1.0, 1.0),
        };
        let b = Aabb {
            min: Point3::new(2.0, 0.0, 0.0),
            max: Point3::new(3.0, 1.0, 1.0),
        };
        assert!(!a.overlaps(&b, 0.0));
    }

    #[test]
    fn overlaps_separated_with_padding_bridges_gap() {
        let a = Aabb {
            min: Point3::new(0.0, 0.0, 0.0),
            max: Point3::new(1.0, 1.0, 1.0),
        };
        let b = Aabb {
            min: Point3::new(1.5, 0.0, 0.0),
            max: Point3::new(2.5, 1.0, 1.0),
        };
        // Gap of 0.5 in x; padding of 0.6 should bridge it
        assert!(a.overlaps(&b, 0.6));
        // Padding of 0.4 should not
        assert!(!a.overlaps(&b, 0.4));
    }

    #[test]
    fn diagonal_unit_cube() {
        let aabb = Aabb {
            min: Point3::new(0.0, 0.0, 0.0),
            max: Point3::new(1.0, 1.0, 1.0),
        };
        let expected = 3_f64.sqrt();
        assert!((aabb.diagonal() - expected).abs() < 1e-12);
    }

    #[test]
    fn diagonal_zero_for_point_aabb() {
        let p = Point3::new(5.0, 5.0, 5.0);
        let aabb = Aabb::from_points(&[p]).unwrap();
        assert!(aabb.diagonal() < f64::EPSILON);
    }
}
