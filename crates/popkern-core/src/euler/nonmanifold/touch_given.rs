//! The glue at given contacts (split from `touch.rs` under the 500-line
//! cap): a construction's touches glued as they are, or completed by the
//! detection where none of them lies.

use super::EulerError;
use super::touch::TouchGlue;
use super::touch_detect::TouchPlan;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::topology::BRep;

impl BRep {
    /// [`Self::glue_touching`] at the contacts a construction found,
    /// detecting none: each point a vertex every touching piece uses,
    /// each curve (its ends among the points; a closed one in two halves)
    /// an edge every touching piece uses, glued. A boolean's general
    /// policy glues at its construction's touches (tangency design §4.8),
    /// where detecting them again from the geometry placed a contact
    /// curve's end 3.5e-6 off the vertex it ends at.
    ///
    /// # Errors
    ///
    /// As [`Self::glue_touching`].
    pub fn glue_touching_at(
        &mut self,
        points: &[Point3],
        curves: &[CurveOnEdge],
    ) -> Result<TouchGlue, EulerError> {
        let plan = self.given_plan(points, curves);
        self.glue_plan(plan)
    }

    /// [`Self::glue_touching_at`] completed by detection: the given
    /// contacts first (exact, a construction's), then whatever the
    /// geometric detection finds that none of them covers (a point
    /// within the tolerance of a given one, a curve with the same ends
    /// and middle, is the given one). A boolean's construction does not
    /// yet trace a tangent curve its analytic arm is silent on (a rod
    /// lying on a face with its seam off the line, W3c), which the
    /// detection finds in closed form.
    ///
    /// # Errors
    ///
    /// As [`Self::glue_touching`].
    pub fn glue_touching_with(
        &mut self,
        points: &[Point3],
        curves: &[CurveOnEdge],
    ) -> Result<TouchGlue, EulerError> {
        let tol = self.tolerance.linear_reach();
        let mut plan = self.given_plan(points, curves);
        let found = self.touch_plan();
        let ends = |c: &CurveOnEdge| [c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end)];
        // On a given curve: within the tolerance of its closest point in
        // its run. A detected curve whose middle lies on a given one is
        // that curve (detection's ends are only as exact as its run
        // search), and a detected point on one is on it.
        let given: Vec<CurveOnEdge> = plan.curves.clone();
        let on_given = |p: &Point3| {
            given.iter().any(|g| {
                let (lo, hi) = (g.t_start.min(g.t_end), g.t_start.max(g.t_end));
                g.geometry.closest_point(p).is_some_and(|t| {
                    t >= lo - 1e-9 && t <= hi + 1e-9 && g.geometry.point_at(t).distance_to(p) <= tol
                })
            })
        };
        for c in found.curves {
            let mid = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
            if !on_given(&mid) {
                plan.add_curve(c, tol);
            }
        }
        for p in found.points {
            if !on_given(&p) {
                plan.add_point(p, tol);
            }
        }
        let curve_ends: Vec<Point3> = plan.curves.iter().flat_map(ends).collect();
        for p in curve_ends {
            plan.add_point(p, tol);
        }
        self.glue_plan(plan)
    }

    /// The plan of the given contacts: each point, each curve with its
    /// ends among the points, a closed curve in two halves.
    fn given_plan(&self, points: &[Point3], curves: &[CurveOnEdge]) -> TouchPlan {
        let tol = self.tolerance.linear_reach();
        let mut plan = TouchPlan {
            points: points.to_vec(),
            curves: Vec::new(),
        };
        let add_point = |plan: &mut TouchPlan, p: Point3| {
            if plan.points.iter().all(|q| q.distance_to(&p) > tol) {
                plan.points.push(p);
            }
        };
        for c in curves {
            let at = |t: f64| c.geometry.point_at(t);
            let (p0, p1) = (at(c.t_start), at(c.t_end));
            if p0.distance_to(&p1) <= tol {
                let tm = 0.5 * (c.t_start + c.t_end);
                add_point(&mut plan, p0);
                add_point(&mut plan, at(tm));
                plan.curves.push(CurveOnEdge {
                    t_end: tm,
                    ..c.clone()
                });
                plan.curves.push(CurveOnEdge {
                    t_start: tm,
                    ..c.clone()
                });
            } else {
                add_point(&mut plan, p0);
                add_point(&mut plan, p1);
                plan.curves.push(c.clone());
            }
        }
        plan
    }
}
