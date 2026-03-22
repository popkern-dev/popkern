//! Surface subdivision infrastructure.

use crate::geometry::nurbs_surface::NurbsSurface;
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

/// A Bezier surface patch with its bounding box.
#[derive(Debug, Clone)]
pub struct SurfacePatch {
    /// The sub-surface for this patch.
    pub surface: NurbsSurface,
    /// Parameter range in u covered by this patch.
    pub u_range: (f64, f64),
    /// Parameter range in v covered by this patch.
    pub v_range: (f64, f64),
    /// Axis-aligned bounding box of the patch's control points.
    pub bbox: Aabb,
}

/// Compute AABB from all control points of a NurbsSurface.
///
/// Returns a degenerate unit box at the origin if the surface has no control points.
fn aabb_from_surface(surface: &NurbsSurface) -> Aabb {
    Aabb::from_points(&surface.control_points).unwrap_or(Aabb {
        min: Point3::new(0.0, 0.0, 0.0),
        max: Point3::new(0.0, 0.0, 0.0),
    })
}

/// Decompose a NURBS surface into sub-patches by splitting at all interior knots.
///
/// Interior knots within `u_range` and `v_range` define the split locations.
/// Each resulting patch has its own AABB computed from its control points.
pub(crate) fn decompose_to_patches(
    surface: &NurbsSurface,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> Vec<SurfacePatch> {
    // Collect distinct interior knots strictly inside u_range.
    let mut u_knots: Vec<f64> = surface
        .knots_u
        .iter()
        .copied()
        .filter(|&k| k > u_range.0 && k < u_range.1)
        .collect();
    u_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    u_knots.dedup();

    // Collect distinct interior knots strictly inside v_range.
    let mut v_knots: Vec<f64> = surface
        .knots_v
        .iter()
        .copied()
        .filter(|&k| k > v_range.0 && k < v_range.1)
        .collect();
    v_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v_knots.dedup();

    // Split at each u knot to get u-strips.
    // Each entry is (sub_surface, u_start, u_end).
    let mut u_strips: Vec<(NurbsSurface, f64, f64)> = Vec::new();
    let mut current = surface.clone();
    let mut current_u_start = u_range.0;

    for &uk in &u_knots {
        match current.split_at_u(uk) {
            Ok((left, right)) => {
                u_strips.push((left, current_u_start, uk));
                current = right;
                current_u_start = uk;
            }
            Err(_) => {
                // Skip knots that fail (e.g. already at boundary).
            }
        }
    }
    u_strips.push((current, current_u_start, u_range.1));

    // For each u-strip, split at each v knot.
    let mut patches = Vec::new();
    for (u_surf, u_start, u_end) in u_strips {
        let mut current_v = u_surf.clone();
        let mut current_v_start = v_range.0;

        for &vk in &v_knots {
            match current_v.split_at_v(vk) {
                Ok((bottom, top)) => {
                    let bbox = aabb_from_surface(&bottom);
                    patches.push(SurfacePatch {
                        surface: bottom,
                        u_range: (u_start, u_end),
                        v_range: (current_v_start, vk),
                        bbox,
                    });
                    current_v = top;
                    current_v_start = vk;
                }
                Err(_) => {
                    // Skip knots that fail.
                }
            }
        }
        let bbox = aabb_from_surface(&current_v);
        patches.push(SurfacePatch {
            surface: current_v,
            u_range: (u_start, u_end),
            v_range: (current_v_start, v_range.1),
            bbox,
        });
    }

    patches
}

/// Split a patch at the midpoint of both u and v, returning 4 sub-patches.
pub(crate) fn refine_patch(patch: &SurfacePatch) -> Vec<SurfacePatch> {
    let u_mid = (patch.u_range.0 + patch.u_range.1) / 2.0;
    let v_mid = (patch.v_range.0 + patch.v_range.1) / 2.0;

    // Split at u_mid → left covers [u_start, u_mid], right covers [u_mid, u_end].
    let (left, right) = match patch.surface.split_at_u(u_mid) {
        Ok(pair) => pair,
        Err(_) => {
            // Fallback: return the original patch as a single item (degenerate input).
            return vec![patch.clone()];
        }
    };

    let mut result = Vec::with_capacity(4);

    for (u_surf, u_start, u_end) in [
        (left, patch.u_range.0, u_mid),
        (right, u_mid, patch.u_range.1),
    ] {
        match u_surf.split_at_v(v_mid) {
            Ok((bottom, top)) => {
                let bbox_b = aabb_from_surface(&bottom);
                result.push(SurfacePatch {
                    surface: bottom,
                    u_range: (u_start, u_end),
                    v_range: (patch.v_range.0, v_mid),
                    bbox: bbox_b,
                });
                let bbox_t = aabb_from_surface(&top);
                result.push(SurfacePatch {
                    surface: top,
                    u_range: (u_start, u_end),
                    v_range: (v_mid, patch.v_range.1),
                    bbox: bbox_t,
                });
            }
            Err(_) => {
                // Fallback: push unsplit patch.
                let bbox = aabb_from_surface(&u_surf);
                result.push(SurfacePatch {
                    surface: u_surf,
                    u_range: (u_start, u_end),
                    v_range: patch.v_range,
                    bbox,
                });
            }
        }
    }

    result
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

    // ---- SurfacePatch / decompose_to_patches / refine_patch tests ----

    fn bilinear_surface() -> NurbsSurface {
        NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0, 1.0, 1.0, 1.0],
        }
    }

    #[test]
    fn decompose_bilinear_one_patch() {
        let s = bilinear_surface();
        let patches = decompose_to_patches(&s, (0.0, 1.0), (0.0, 1.0));
        assert_eq!(patches.len(), 1);
    }

    #[test]
    fn aabb_contains_surface_points() {
        use crate::geometry::surface::SurfaceEval;
        let s = bilinear_surface();
        let patches = decompose_to_patches(&s, (0.0, 1.0), (0.0, 1.0));
        for patch in &patches {
            for &u in &[0.0f64, 0.5, 0.99] {
                for &v in &[0.0f64, 0.5, 0.99] {
                    let tu = patch.u_range.0 + u * (patch.u_range.1 - patch.u_range.0);
                    let tv = patch.v_range.0 + v * (patch.v_range.1 - patch.v_range.0);
                    let p = patch.surface.point_at(tu, tv);
                    assert!(
                        p.x >= patch.bbox.min.x - 1e-10 && p.x <= patch.bbox.max.x + 1e-10,
                        "x={} not in [{}, {}]",
                        p.x,
                        patch.bbox.min.x,
                        patch.bbox.max.x
                    );
                    assert!(
                        p.y >= patch.bbox.min.y - 1e-10 && p.y <= patch.bbox.max.y + 1e-10,
                        "y={} not in [{}, {}]",
                        p.y,
                        patch.bbox.min.y,
                        patch.bbox.max.y
                    );
                }
            }
        }
    }

    #[test]
    fn refine_produces_four_subpatches() {
        let s = bilinear_surface();
        let patches = decompose_to_patches(&s, (0.0, 1.0), (0.0, 1.0));
        assert_eq!(patches.len(), 1);
        let refined = refine_patch(&patches[0]);
        assert_eq!(refined.len(), 4);

        // The union of sub-patch parameter ranges should cover the original range.
        let u0 = patches[0].u_range.0;
        let u1 = patches[0].u_range.1;
        let v0 = patches[0].v_range.0;
        let v1 = patches[0].v_range.1;
        let u_mid = (u0 + u1) / 2.0;
        let v_mid = (v0 + v1) / 2.0;

        // Each quadrant should appear exactly once.
        let ranges: Vec<_> = refined.iter().map(|p| (p.u_range, p.v_range)).collect();
        assert!(ranges.contains(&((u0, u_mid), (v0, v_mid))));
        assert!(ranges.contains(&((u0, u_mid), (v_mid, v1))));
        assert!(ranges.contains(&((u_mid, u1), (v0, v_mid))));
        assert!(ranges.contains(&((u_mid, u1), (v_mid, v1))));
    }
}
