//! Engraving: pockets cut into planar polyhedron faces from 2D
//! polygon sets (text glyphs or arbitrary profiles).
//!
//! A pocket is exact topology: an inner loop on the face, straight
//! walls, a flat floor — and for polygons with holes, islands rising
//! back to the surface. Everything stays planes and lines.

use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

use crate::error::PartError;
use crate::poly::PolyFace;

/// A simple polygon with optional holes, in 2D layout coordinates
/// (x right, y up), all contours CCW. Text glyphs produce these;
/// any profile works.
#[derive(Debug, Clone)]
pub struct HoledPolygon {
    /// Outer contour, CCW.
    pub outer: Vec<(f64, f64)>,
    /// Hole contours, CCW.
    pub holes: Vec<Vec<(f64, f64)>>,
}

/// Bounding box `(min_x, max_x, min_y, max_y)` of a polygon set.
pub fn polygons_bbox(polys: &[HoledPolygon]) -> (f64, f64, f64, f64) {
    let (mut min_x, mut max_x) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut min_y, mut max_y) = (f64::INFINITY, f64::NEG_INFINITY);
    for g in polys {
        for &(x, y) in g.outer.iter().chain(g.holes.iter().flatten()) {
            min_x = min_x.min(x);
            max_x = max_x.max(x);
            min_y = min_y.min(y);
            max_y = max_y.max(y);
        }
    }
    (min_x, max_x, min_y, max_y)
}

/// The in-plane frame of a planar face: centroid, outward normal
/// (Newell over the outer cycle), and the inradius (nearest outer
/// edge to the centroid).
#[derive(Debug, Clone, Copy)]
pub struct FaceFrame {
    /// Outer-cycle centroid.
    pub centroid: Point3,
    /// Outward unit normal.
    pub normal: Vec3,
    /// Distance from the centroid to the nearest outer edge.
    pub inradius: f64,
}

/// Compute a face's engraving frame.
pub fn face_frame(points: &[Point3], face: &PolyFace) -> FaceFrame {
    let outline: Vec<Point3> = face.outer.iter().map(|&i| points[i]).collect();
    let n = outline.len();
    let zero = Point3::new(0.0, 0.0, 0.0);
    let centroid = zero
        + outline
            .iter()
            .fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| acc + (*p - zero))
            * (1.0 / n as f64);
    let mut normal = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let p = outline[i];
        let q = outline[(i + 1) % n];
        normal.x += (p.y - q.y) * (p.z + q.z);
        normal.y += (p.z - q.z) * (p.x + q.x);
        normal.z += (p.x - q.x) * (p.y + q.y);
    }
    let normal = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let mut inradius = f64::INFINITY;
    for i in 0..n {
        let a = outline[i];
        let b = outline[(i + 1) % n];
        let ab = b - a;
        let len2 = ab.dot(&ab);
        let t = ((centroid - a).dot(&ab) / len2).clamp(0.0, 1.0);
        let foot = a + ab * t;
        inradius = inradius.min(centroid.distance_to(&foot));
    }
    FaceFrame {
        centroid,
        normal,
        inradius,
    }
}

/// Placement of a polygon set on a face.
#[derive(Debug, Clone, Copy)]
pub struct EngravePlacement {
    /// Pocket centre in 3D (must lie on the face plane).
    pub center: Point3,
    /// In-plane text-right direction.
    pub right: Vec3,
    /// In-plane text-up direction.
    pub up: Vec3,
    /// Layout-to-model scale factor.
    pub scale: f64,
}

impl EngravePlacement {
    /// Centre the polygons on the face, upright (text-up = the
    /// projection of world-up onto the face; +Y for horizontal
    /// faces), scaled so the layout's half-diagonal fits
    /// `fit × inradius`.
    pub fn centered(frame: &FaceFrame, polys: &[HoledPolygon], fit: f64) -> EngravePlacement {
        let up_hint = if frame.normal.z.abs() > 0.9 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let up = (up_hint - frame.normal * up_hint.dot(&frame.normal))
            .normalized()
            .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
        let right = up.cross(&frame.normal);
        let (min_x, max_x, min_y, max_y) = polygons_bbox(polys);
        let half_diag = (((max_x - min_x) / 2.0).powi(2) + ((max_y - min_y) / 2.0).powi(2)).sqrt();
        EngravePlacement {
            center: frame.centroid,
            right,
            up,
            scale: fit * frame.inradius / half_diag.max(1e-12),
        }
    }
}

/// Engrave `polys` into `faces[face_idx]` as pockets of `depth`,
/// appending wall/floor/island faces and inner loops in place.
///
/// # Errors
///
/// [`PartError::NotApplicable`] when the scaled layout does not fit
/// inside the face's outer boundary (overflowing pockets would
/// self-intersect — the check carries the fitting metrics).
pub fn engrave_polygons(
    points: &mut Vec<Point3>,
    faces: &mut Vec<PolyFace>,
    face_idx: usize,
    polys: &[HoledPolygon],
    placement: &EngravePlacement,
    depth: f64,
) -> Result<(), PartError> {
    let frame = face_frame(points, &faces[face_idx]);
    let (min_x, max_x, min_y, max_y) = polygons_bbox(polys);
    let layout_center = ((min_x + max_x) / 2.0, (min_y + max_y) / 2.0);
    // Overflow guard: every laid-out point must fall inside the face
    // outer boundary with a small margin (overflowing pockets would
    // self-intersect the face). Exact 2D point-in-polygon in the
    // face plane, so off-centre placements (dice corner numbers) are
    // judged against the real outline, not the inradius.
    {
        let outline2d: Vec<(f64, f64)> = faces[face_idx]
            .outer
            .iter()
            .map(|&i| {
                let d = points[i] - frame.centroid;
                (d.dot(&placement.right), d.dot(&placement.up))
            })
            .collect();
        let margin = frame.inradius * 0.02;
        let to_2d = |x: f64, y: f64| -> (f64, f64) {
            let p = placement.center
                + placement.right * ((x - layout_center.0) * placement.scale)
                + placement.up * ((y - layout_center.1) * placement.scale);
            let d = p - frame.centroid;
            (d.dot(&placement.right), d.dot(&placement.up))
        };
        let inside = |pt: (f64, f64)| -> bool {
            let n = outline2d.len();
            let mut in_poly = false;
            let mut min_edge_dist = f64::INFINITY;
            for i in 0..n {
                let (a, b) = (outline2d[i], outline2d[(i + 1) % n]);
                if (a.1 > pt.1) != (b.1 > pt.1) {
                    let x = a.0 + (pt.1 - a.1) / (b.1 - a.1) * (b.0 - a.0);
                    if x > pt.0 {
                        in_poly = !in_poly;
                    }
                }
                let (dx, dy) = (b.0 - a.0, b.1 - a.1);
                let l2 = (dx * dx + dy * dy).max(1e-300);
                let t = (((pt.0 - a.0) * dx + (pt.1 - a.1) * dy) / l2).clamp(0.0, 1.0);
                let (qx, qy) = (a.0 + dx * t, a.1 + dy * t);
                min_edge_dist = min_edge_dist.min((pt.0 - qx).phypot(pt.1 - qy));
            }
            in_poly && min_edge_dist >= margin
        };
        for glyph in polys {
            for &(x, y) in glyph.outer.iter().chain(glyph.holes.iter().flatten()) {
                let pt = to_2d(x, y);
                if !inside(pt) {
                    return Err(PartError::NotApplicable(format!(
                        "engraving does not fit: a pocket point lands outside the face                          (or within {margin:.4} of its boundary); scale {:.4}, face inradius {:.4}",
                        placement.scale, frame.inradius
                    )));
                }
            }
        }
    }

    let ring = |poly: &[(f64, f64)], points: &mut Vec<Point3>| -> (Vec<usize>, Vec<usize>) {
        let top: Vec<usize> = poly
            .iter()
            .map(|&(x, y)| {
                let p = placement.center
                    + placement.right * ((x - layout_center.0) * placement.scale)
                    + placement.up * ((y - layout_center.1) * placement.scale);
                points.push(p);
                points.len() - 1
            })
            .collect();
        let bottom: Vec<usize> = top
            .iter()
            .map(|&i| {
                let p = points[i] + frame.normal * (-depth);
                points.push(p);
                points.len() - 1
            })
            .collect();
        (top, bottom)
    };
    for glyph in polys {
        let (top, bottom) = ring(&glyph.outer, points);
        let m = top.len();
        let mut inner = top.clone();
        inner.reverse();
        faces[face_idx].inners.push(inner);
        for i in 0..m {
            let (a, b) = (top[i], top[(i + 1) % m]);
            let (af, bf) = (bottom[i], bottom[(i + 1) % m]);
            faces.push(PolyFace::simple(vec![a, b, bf, af]));
        }
        let mut floor = PolyFace::simple(bottom);
        for hole in &glyph.holes {
            let (h_top, h_bottom) = ring(hole, points);
            let k = h_top.len();
            faces.push(PolyFace::simple(h_top.clone()));
            for i in 0..k {
                let (a, b) = (h_top[i], h_top[(i + 1) % k]);
                let (af, bf) = (h_bottom[i], h_bottom[(i + 1) % k]);
                faces.push(PolyFace::simple(vec![b, a, af, bf]));
            }
            let mut floor_inner = h_bottom;
            floor_inner.reverse();
            floor.inners.push(floor_inner);
        }
        faces.push(floor);
    }
    Ok(())
}
