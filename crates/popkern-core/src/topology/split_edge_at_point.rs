use crate::error::CoreError;
use crate::geometry::Point3;
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

        // The curve's direction along the half-edge as `split_edge` reads
        // it (a closed edge's ends are one point, so its pcurves witness
        // the direction): comparing the ends' distances to the start
        // vertex broke the tie by rounding on a closed circle and split
        // at the point's mirror image (PK-1189).
        let pos_b = self
            .half_edges
            .get(he_ab_data.next)
            .and_then(|next| self.vertices.get(next.origin))
            .map(|v| v.position)
            .unwrap_or(pos_a);
        let he_ba = he_ab_data.mate().unwrap_or(he_ab);
        let curve_aligned = super::split_position::curve_alignment_for_split(
            self,
            he_ab,
            he_ba,
            curve_on_edge,
            pos_a,
            pos_b,
        );
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

#[cfg(test)]
#[path = "split_edge_at_point_tests.rs"]
mod tests;
