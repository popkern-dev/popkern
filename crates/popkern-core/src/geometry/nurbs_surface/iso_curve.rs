//! Exact iso-parameter curve extraction.
//!
//! The iso curve of a NURBS surface at a fixed parameter is itself a
//! NURBS curve whose homogeneous control points are the de Boor
//! combination of the grid's rows (or columns) at that parameter —
//! no sampling, no fitting. Callers that need an exact on-surface
//! cross-section (fillet station arcs, seam sections, trim rebuilds)
//! use this instead of a Hermite fit, whose ~1e-3 curvature-class
//! error fails on-surface postconditions in the 1e-5 band.

use super::NurbsSurface;
use crate::geometry::NurbsCurve;
use crate::geometry::nurbs_ops::{Homogeneous, de_boor_1d};

impl NurbsSurface {
    /// The exact iso curve at fixed `u`, parameterized by `v` (same
    /// degree and knot vector as the surface's v direction). Every
    /// point of the returned curve satisfies
    /// `curve.point_at(v) == surface.point_at(u, v)` to roundoff.
    pub fn iso_curve_at_u(&self, u: f64) -> NurbsCurve {
        let mut control_points = Vec::with_capacity(self.num_v);
        let mut weights = Vec::with_capacity(self.num_v);
        let mut column: Vec<Homogeneous> = Vec::with_capacity(self.num_u);
        for i_v in 0..self.num_v {
            column.clear();
            for i_u in 0..self.num_u {
                column.push(Homogeneous::from_point_weight(
                    self.cp(i_u, i_v),
                    self.wt(i_u, i_v),
                ));
            }
            let q = de_boor_1d(self.degree_u, &self.knots_u, &column, u);
            control_points.push(q.to_point3());
            weights.push(q.w);
        }
        NurbsCurve {
            degree: self.degree_v,
            knots: self.knots_v.clone(),
            control_points,
            weights,
        }
    }

    /// The exact iso curve at fixed `v`, parameterized by `u` (same
    /// degree and knot vector as the surface's u direction). Every
    /// point of the returned curve satisfies
    /// `curve.point_at(u) == surface.point_at(u, v)` to roundoff.
    pub fn iso_curve_at_v(&self, v: f64) -> NurbsCurve {
        let mut control_points = Vec::with_capacity(self.num_u);
        let mut weights = Vec::with_capacity(self.num_u);
        let mut row: Vec<Homogeneous> = Vec::with_capacity(self.num_v);
        for i_u in 0..self.num_u {
            row.clear();
            for i_v in 0..self.num_v {
                row.push(Homogeneous::from_point_weight(
                    self.cp(i_u, i_v),
                    self.wt(i_u, i_v),
                ));
            }
            let q = de_boor_1d(self.degree_v, &self.knots_v, &row, v);
            control_points.push(q.to_point3());
            weights.push(q.w);
        }
        NurbsCurve {
            degree: self.degree_u,
            knots: self.knots_u.clone(),
            control_points,
            weights,
        }
    }
}
