use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::curve_inversion::invert_curve_at_point;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, SplitEdgeResult};

impl BRep {
    /// Split an edge at the parameter that places `point` on the edge's curve.
    ///
    /// Wraps [`crate::geometry::curve_inversion::invert_curve_at_point`]
    /// + [`BRep::split_edge`]. Per ADR §14.2 (constructive apex builder).
    pub fn split_edge_at_point(
        &mut self,
        edge: EdgeHandle,
        point: Point3,
        tolerance: &Tolerance,
    ) -> Result<SplitEdgeResult, CoreError> {
        let edge_data = self
            .edges
            .get(edge)
            .ok_or_else(|| CoreError::InvalidHandle("edge handle is invalid".into()))?;
        let curve_on_edge = edge_data.curve.as_ref().ok_or_else(|| {
            CoreError::PreconditionFailed("split_edge_at_point: edge has no curve geometry".into())
        })?;

        let he_ab = edge_data.half_edge;
        let he_ab_data = self
            .half_edges
            .get(he_ab)
            .ok_or_else(|| CoreError::InvalidHandle("half_edge on edge is invalid".into()))?;
        let vertex_a = he_ab_data.origin;
        let pos_a = self.vertices[vertex_a].position;

        let t_curve = invert_curve_at_point(curve_on_edge, point, tolerance).map_err(|e| {
            CoreError::PreconditionFailed(format!(
                "split_edge_at_point: point inversion failed: {e:?}"
            ))
        })?;

        let t_range = curve_on_edge.t_end - curve_on_edge.t_start;
        if t_range.abs() < 1e-30 {
            return Err(CoreError::PreconditionFailed(
                "split_edge_at_point: edge has zero parameter range".into(),
            ));
        }
        let t_curve_frac = (t_curve - curve_on_edge.t_start) / t_range;

        let p_t_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
        let p_t_end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);
        let curve_aligned = pos_a.distance_to(&p_t_start) <= pos_a.distance_to(&p_t_end);
        let t_he = if curve_aligned {
            t_curve_frac
        } else {
            1.0 - t_curve_frac
        };

        let endpoint_eps = 1e-9;
        if t_he <= endpoint_eps || t_he >= 1.0 - endpoint_eps {
            return Err(CoreError::PreconditionFailed(format!(
                "split_edge_at_point: point too close to existing edge endpoint \
                 (t_he={t_he}); use the existing vertex instead"
            )));
        }

        self.split_edge(edge, t_he)
    }
}
