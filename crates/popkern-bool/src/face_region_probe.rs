//! Reusable face-region containment for repeated queries on one face
//! (PK-814).
//!
//! [`crate::split_inner_loop::point_inside_face_region_uv_with_density`]
//! rebuilds the face's outer boundary polygon and every hole's box on
//! each call. The clip's zero-crossing vote asks it five times per
//! record, and on a disc carrying 63 slot holes every late slot cut
//! sampled the rim 65 times and boxed 252 hole edges per query
//! (measured: ~2 ms of an 11 ms cut, all in the clip). The probe builds
//! the planar polygon and the hole boxes once and each hole polygon only
//! when a query lands in its box. Answers are identical by construction:
//! it is the same sampling, the same box gate and the same even-odd test,
//! and every chart the planar path does not serve (periodic charts, full
//! chart covers, chart cells) defers to the original function per query.

use std::cell::OnceCell;

use popkern_core::geometry::Point3;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::split_inner_loop::ray_cast_contains;

/// Samples per hole pcurve, as `hole_contour_contains` samples them.
const HOLE_SAMPLES_PER_EDGE: usize = 32;

/// One hole of a planar face: its conservative box and, once a query
/// lands in the box, its sampled polygon.
struct Hole<'a> {
    contour: Vec<(HalfEdgeHandle, &'a Pcurve)>,
    uv_box: Option<crate::pcurve_box::UvBox>,
    polygon: OnceCell<Vec<Point3>>,
}

/// Containment of chart points in one face's region (outer loop minus
/// holes), built once and queried many times.
pub(crate) struct FaceRegionProbe<'a> {
    brep: &'a BRep,
    face: FaceHandle,
    tolerance: &'a Tolerance,
    density: usize,
    planar: Option<(Vec<Point3>, Vec<Hole<'a>>)>,
}

impl<'a> FaceRegionProbe<'a> {
    /// A probe answering exactly as
    /// `point_inside_face_region_uv_with_density(brep, face, _, tolerance, density)`.
    pub(crate) fn new(
        brep: &'a BRep,
        face: FaceHandle,
        tolerance: &'a Tolerance,
        density: usize,
    ) -> Self {
        let planar = Self::planar_parts(brep, face, tolerance, density);
        Self {
            brep,
            face,
            tolerance,
            density,
            planar,
        }
    }

    /// The outer polygon and holes when the face takes the planar
    /// even-odd path; `None` defers every query to the original.
    fn planar_parts(
        brep: &'a BRep,
        face: FaceHandle,
        tolerance: &Tolerance,
        density: usize,
    ) -> Option<(Vec<Point3>, Vec<Hole<'a>>)> {
        let surface = brep.faces.get(face)?.surface.as_ref()?;
        let (u_period, v_period) = crate::clip::classification_chart_periods(&surface.geometry);
        if u_period.is_some()
            || v_period.is_some()
            || crate::clip::face_covers_full_chart(brep, face)
            || crate::clip::chart_cell::face_chart_cell(brep, face, tolerance).is_some()
        {
            return None;
        }
        let boundary = crate::clip::face_uv_boundary(brep, face)?;
        let mut outer = Vec::new();
        for (_he, pcurve) in &boundary {
            let (t0, t1) = pcurve.parameter_range();
            let range = t1 - t0;
            if range.abs() < tolerance.linear() {
                continue;
            }
            for i in 0..density {
                let q = pcurve.point_at(t0 + range * i as f64 / density as f64);
                outer.push(Point3::new(q.x, q.y, 0.0));
            }
        }
        let holes = crate::clip::face_uv_inner_boundary_refs(brep, face)
            .into_iter()
            .map(|contour| Hole {
                uv_box: crate::pcurve_box::contour_uv_box(contour.iter().map(|(_, pc)| *pc)),
                contour,
                polygon: OnceCell::new(),
            })
            .collect();
        Some((outer, holes))
    }

    /// Whether `uv_point` lies in the face region.
    pub(crate) fn contains(&self, uv_point: Point3) -> bool {
        let Some((outer, holes)) = &self.planar else {
            return crate::split_inner_loop::point_inside_face_region_uv_with_density(
                self.brep,
                self.face,
                uv_point,
                self.tolerance,
                self.density,
            );
        };
        if !ray_cast_contains(outer, uv_point.x, uv_point.y) {
            return false;
        }
        let mut inside = true;
        for hole in holes {
            if hole
                .uv_box
                .as_ref()
                .is_some_and(|bx| crate::pcurve_box::outside(bx, uv_point.x, uv_point.y, 0.0))
            {
                continue;
            }
            let polygon = hole.polygon.get_or_init(|| {
                let mut polygon = Vec::new();
                for (_he, pcurve) in &hole.contour {
                    let (t0, t1) = pcurve.parameter_range();
                    let range = t1 - t0;
                    if range.abs() < self.tolerance.linear() {
                        continue;
                    }
                    for i in 0..HOLE_SAMPLES_PER_EDGE {
                        let t = t0 + range * i as f64 / HOLE_SAMPLES_PER_EDGE as f64;
                        let q = pcurve.point_at(t);
                        polygon.push(Point3::new(q.x, q.y, 0.0));
                    }
                }
                polygon
            });
            if ray_cast_contains(polygon, uv_point.x, uv_point.y) {
                inside = !inside;
            }
        }
        inside
    }
}

#[cfg(test)]
#[path = "face_region_probe_tests.rs"]
mod tests;
