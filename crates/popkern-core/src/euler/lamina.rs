//! A closed lamina loop as one transaction (PK-709): the `mvfs`, the
//! `mev` chain and the closing `mef` a sketch profile compiles to,
//! verified once at the boundary instead of after every step.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::contract::{EulerDelta, EulerTopologyCounts};

/// What [`BRep::closed_lamina`] built.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LaminaResult {
    /// One vertex per point, in order.
    pub vertices: Vec<VertexHandle>,
    /// Edge `k` runs `vertices[k]` → `vertices[k + 1]`; the last closes
    /// the loop back to `vertices[0]`.
    pub edges: Vec<EdgeHandle>,
    /// The face the closing `mef` enclosed.
    pub mef_face: FaceHandle,
    /// The `mvfs` face (the lamina's other side).
    pub mvfs_face: FaceHandle,
}

impl BRep {
    /// Build a new lamina (a two-faced sheet) whose loop runs through
    /// `points` in order and closes back to the first: `mvfs`, `mev` per
    /// further point, one closing `mef`. Geometry is the caller's: the
    /// edges carry no curves yet.
    ///
    /// One transaction: the steps run on a working copy that must show
    /// exactly their combined Euler delta and pass the structural checks
    /// (half-edge pairing, cycles, vertex anchors, manifold edges) before
    /// it joins the body. Per-step whole-model checks made a 1000-point
    /// profile quadratic (PK-709).
    ///
    /// # Errors
    ///
    /// [`CoreError::PreconditionFailed`] for fewer than two points or a
    /// failed postcondition (the body is then unchanged).
    pub fn closed_lamina(&mut self, points: &[Point3]) -> Result<LaminaResult, CoreError> {
        if points.len() < 2 {
            return Err(CoreError::PreconditionFailed(
                "a closed lamina needs at least two points".into(),
            ));
        }
        let before = EulerTopologyCounts::from_brep(self);
        let mut work = self.clone();
        let mvfs = work.mvfs(points[0])?;
        let first = work.apply_mev_at_vertex(mvfs.vertex, mvfs.loop_, points[1])?;
        let mut vertices = vec![mvfs.vertex, first.vertex];
        let mut edges = vec![first.edge];
        let mut last_twin = first.half_edge_twin;
        for &p in &points[2..] {
            let mev = work.apply_mev(last_twin, p)?;
            last_twin = mev.half_edge_twin;
            vertices.push(mev.vertex);
            edges.push(mev.edge);
        }
        let mef = work.apply_mef(last_twin, first.half_edge)?;
        edges.push(mef.edge);
        let steps = points.len() as isize - 1;
        EulerDelta::verify(
            "closed_lamina",
            before,
            EulerTopologyCounts::from_brep(&work),
            EulerDelta::mvfs()
                .then(EulerDelta::mev().times(steps))
                .then(EulerDelta::mef()),
        )?;
        work.validate_half_edge_pairing().map_err(CoreError::from)?;
        work.validate_half_edge_cycles().map_err(CoreError::from)?;
        work.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        work.validate_manifold_edges().map_err(CoreError::from)?;
        *self = work;
        Ok(LaminaResult {
            vertices,
            edges,
            mef_face: mef.face,
            mvfs_face: mvfs.face,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    fn square() -> Vec<Point3> {
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ]
    }

    #[test]
    fn a_square_lamina_has_two_faces_and_matches_the_stepwise_build() {
        let mut brep = BRep::new(Tolerance::default());
        let lamina = brep.closed_lamina(&square()).expect("lamina");
        assert_eq!(lamina.vertices.len(), 4);
        assert_eq!(lamina.edges.len(), 4);
        assert_eq!(brep.faces.len(), 2);
        assert_eq!(brep.edges.len(), 4);
        // Edge k runs vertex k -> k+1.
        for (k, &e) in lamina.edges.iter().enumerate() {
            let he = brep.edges[e].half_edge;
            let a = brep.half_edges[he].origin;
            let b = brep.half_edges[brep.half_edges[he].twin].origin;
            let (v0, v1) = (lamina.vertices[k], lamina.vertices[(k + 1) % 4]);
            assert!((a, b) == (v0, v1) || (a, b) == (v1, v0), "edge {k}");
        }
        // The same steps through the verified operators give the same
        // topology counts.
        let mut stepwise = BRep::new(Tolerance::default());
        let pts = square();
        let mvfs = stepwise.mvfs(pts[0]).expect("mvfs");
        let first = stepwise
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, pts[1])
            .expect("mev");
        let mut last = first.half_edge_twin;
        for &p in &pts[2..] {
            last = stepwise.mev(last, p).expect("mev").half_edge_twin;
        }
        stepwise.mef(last, first.half_edge).expect("mef");
        assert_eq!(stepwise.faces.len(), brep.faces.len());
        assert_eq!(stepwise.half_edges.len(), brep.half_edges.len());
    }

    #[test]
    fn a_bigon_closes_and_one_point_is_refused() {
        let mut brep = BRep::new(Tolerance::default());
        let pts = [Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let lamina = brep.closed_lamina(&pts).expect("bigon");
        assert_eq!(lamina.edges.len(), 2);
        let before = brep.faces.len();
        assert!(brep.closed_lamina(&pts[..1]).is_err());
        assert_eq!(
            brep.faces.len(),
            before,
            "a refused lamina leaves the body unchanged"
        );
    }
}
