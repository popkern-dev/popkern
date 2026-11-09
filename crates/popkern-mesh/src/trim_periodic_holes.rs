//! Holes of a whole periodic chart that cross its seam (PK-1136).
//!
//! A face bounded only by its seam (a whole sphere with holes cut in
//! it) is meshed over its surface's natural rectangle, one period of the
//! chart. A hole's loop is unwrapped continuously, so a hole straddling
//! the seam runs past the rectangle's edge (the lobes of Viviani's
//! figure eight cut from a sphere: u from -pi/2 to pi/2), and the part
//! beyond it trimmed nothing: the face kept that stretch of the hole as
//! material (measured: 1142.41 meshed of a face of 1028.32, one lobe).
//! Each such hole is also cut shifted by one period, so the rectangle
//! sees the whole of it; a copy's part beyond the rectangle cuts no
//! grid cell.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use crate::trim::TrimLoop;

/// Keep `loops`' holes, dropping its outer loops (the caller supplies
/// the rectangle), and add each hole shifted by a period wherever it
/// runs past the surface's natural range in that direction.
pub(crate) fn retain_holes_wrapped(loops: &mut Vec<TrimLoop>, surface: &SurfaceGeometry) {
    loops.retain(|l| !l.is_outer);
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let (pu, pv) = surface.parameter_periodicity();
    let mut copies = Vec::new();
    for l in loops.iter() {
        let (mut lo, mut hi) = (
            (f64::INFINITY, f64::INFINITY),
            (f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for &(u, v) in &l.points {
            lo = (lo.0.min(u), lo.1.min(v));
            hi = (hi.0.max(u), hi.1.max(v));
        }
        let mut shift = |du: f64, dv: f64| {
            copies.push(TrimLoop {
                points: l.points.iter().map(|&(u, v)| (u + du, v + dv)).collect(),
                points_3d: l.points_3d.clone(),
                is_outer: false,
                from_pcurves: l.from_pcurves,
            });
        };
        if let Some(p) = pu.filter(|p| *p > 0.0) {
            if lo.0 < u0 {
                shift(p, 0.0);
            }
            if hi.0 > u1 {
                shift(-p, 0.0);
            }
        }
        if let Some(p) = pv.filter(|p| *p > 0.0) {
            if lo.1 < v0 {
                shift(0.0, p);
            }
            if hi.1 > v1 {
                shift(0.0, -p);
            }
        }
    }
    loops.extend(copies);
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::{Point3, Vec3};

    fn hole(us: &[f64]) -> TrimLoop {
        TrimLoop {
            points: us.iter().map(|&u| (u, 1.5)).chain([(0.0, 1.0)]).collect(),
            points_3d: vec![Point3::new(0.0, 0.0, 0.0); us.len() + 1],
            is_outer: false,
            from_pcurves: true,
        }
    }

    /// A sphere's hole running to u = -0.5, past its chart's start, is
    /// also cut a period on, at u up to 2 pi; one inside the chart, and
    /// the outer loop, are left as they are and dropped respectively.
    #[test]
    fn a_hole_across_the_seam_is_cut_a_period_on() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
            axis: Vec3::new(0.0, 0.0, 1.0),
        };
        let mut outer = hole(&[0.0, 6.0]);
        outer.is_outer = true;
        let mut loops = vec![outer, hole(&[-0.5, 0.5]), hole(&[1.0, 2.0])];
        retain_holes_wrapped(&mut loops, &sphere);
        assert_eq!(loops.len(), 3, "two holes and one wrapped copy");
        assert!(loops.iter().all(|l| !l.is_outer));
        let tau = std::f64::consts::TAU;
        assert!(
            loops
                .iter()
                .any(|l| (l.points[0].0 - (tau - 0.5)).abs() < 1e-12)
        );
    }
}
