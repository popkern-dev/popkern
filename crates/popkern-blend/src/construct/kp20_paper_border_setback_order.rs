//! The angular order of a setback corner's rounds about its vertex
//! (best-fit frame from the trim-curve midpoints), shared by the setback
//! extractor and the user-setback range query (PK-623).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

use crate::fillet_surface::FilletSurfaceResult;

/// The rounds' indices in angular order about the vertex.
pub(super) fn angular_order(fillet_surfaces: &[FilletSurfaceResult]) -> Result<Vec<usize>, String> {
    let n = fillet_surfaces.len();
    let rep_point = |fi: usize| -> Point3 {
        let c = &fillet_surfaces[fi].trim_curve_1;
        let (t0, t1) = c.parameter_range();
        c.point_at(0.5 * (t0 + t1))
    };
    let centroid = {
        let mut acc = [0.0_f64; 3];
        for fi in 0..n {
            let p = rep_point(fi);
            acc[0] += p.x;
            acc[1] += p.y;
            acc[2] += p.z;
        }
        Point3::new(acc[0] / n as f64, acc[1] / n as f64, acc[2] / n as f64)
    };
    let mut axis = Vec3::new(0.0, 0.0, 0.0);
    for fi in 0..n {
        let a = rep_point(fi) - centroid;
        let b = rep_point((fi + 1) % n) - centroid;
        axis = axis + a.cross(&b);
    }
    let axis = axis
        .normalized()
        .ok_or_else(|| "degenerate fillet arrangement".to_string())?;
    let seed = rep_point(0) - centroid;
    let u0 = (seed - axis * seed.dot(&axis))
        .normalized()
        .ok_or_else(|| "degenerate angular frame".to_string())?;
    let v0 = axis.cross(&u0);
    let angle_of = |p: Point3| -> f64 {
        let d = p - centroid;
        d.dot(&v0).patan2(d.dot(&u0))
    };
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        angle_of(rep_point(a))
            .partial_cmp(&angle_of(rep_point(b)))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(order)
}
