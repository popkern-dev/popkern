//! Sheet body constructors — planar open shells.
//!
//! A **sheet body** in popkern-core is a single-face or multi-face
//! BRep that describes an open shell rather than a closed solid.
//! Sheet bodies are the natural input for §5.7 2D fillets, §8.3
//! free-boundary cap tests, and any sheet-metal / thin-wall
//! workflow where the geometry is not volumetric.
//!
//! ## Manifold representation
//!
//! popkern's BRep is strictly manifold: every edge has exactly two
//! half-edges, each belonging to exactly one face. To represent a
//! sheet body inside this model, every "free-boundary" edge is
//! backed by a **backside face** on the opposite orientation of the
//! sheet. The resulting topology has twice as many faces as a
//! visually-single-sided sheet, and the [`Face::surface`] on the
//! backside face is left as `None` — downstream consumers that
//! distinguish front and back can filter on that.
//!
//! For a single planar quad this produces a 2-face, 4-edge,
//! 4-vertex body that satisfies `validate_manifold_edges`,
//! `validate_loop_closure`, and the Euler-Poincaré identity for a
//! sphere-topology shell: `V − E + F = 2` (4 − 4 + 2 = 2). The
//! "solid" wrapper is the same single-shell container popkern uses
//! for closed solids; sheet bodies simply happen to have all their
//! faces coplanar.
//!
//! ## Constructors
//!
//! - [`BRep::sheet_quad`] — axis-aligned planar quad
//! - [`BRep::sheet_triangle`] — planar triangle
//! - [`BRep::sheet_polygon`] — planar polygon with n ≥ 3 vertices
//!
//! Each constructor attaches a [`SurfaceGeometry::Plane`] to the
//! front face (the one holding the walked-in-order outer loop)
//! and leaves the backside face's surface as `None`. Edges carry
//! straight-line [`CurveGeometry::Line`] geometry and each
//! half-edge on the front face gets a 2D pcurve via
//! `attach_planar_pcurves`.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

impl BRep {
    /// Construct an axis-aligned planar quad sheet body in the
    /// `z = z_plane` plane, spanning
    /// `x ∈ [x_min, x_max]`, `y ∈ [y_min, y_max]`.
    ///
    /// Walks the 4 corners counter-clockwise when viewed from `+z`:
    /// `(x_min, y_min) → (x_max, y_min) → (x_max, y_max) → (x_min, y_max)`.
    /// The outward normal of the front face points in `+z`.
    ///
    /// # Errors
    ///
    /// - [`CoreError::PreconditionFailed`] if `x_min >= x_max` or
    ///   `y_min >= y_max`.
    /// - [`CoreError::Validation`] if any downstream Euler operator
    ///   or validation pass fails.
    pub fn sheet_quad(
        x_min: f64,
        x_max: f64,
        y_min: f64,
        y_max: f64,
        z_plane: f64,
        tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        if x_min >= x_max || y_min >= y_max {
            return Err(CoreError::PreconditionFailed(format!(
                "sheet_quad: invalid extents x=[{x_min},{x_max}] y=[{y_min},{y_max}]"
            )));
        }
        let corners = [
            Point3::new(x_min, y_min, z_plane),
            Point3::new(x_max, y_min, z_plane),
            Point3::new(x_max, y_max, z_plane),
            Point3::new(x_min, y_max, z_plane),
        ];
        Self::sheet_polygon(&corners, tolerance)
    }

    /// Construct a planar triangle sheet body with the three given
    /// vertex positions. The three points must not be colinear.
    ///
    /// The outward normal of the front face is `(p1 − p0) × (p2 − p0)`
    /// normalised. The 3 points are walked in the order given; callers
    /// wanting a specific outward-facing orientation should pass them
    /// in counter-clockwise order as viewed from the intended outside.
    ///
    /// # Errors
    ///
    /// - [`CoreError::PreconditionFailed`] if the three points are
    ///   colinear (cross product magnitude below the tolerance's
    ///   linear threshold).
    /// - [`CoreError::Validation`] if any downstream Euler operator
    ///   or validation pass fails.
    pub fn sheet_triangle(
        p0: Point3,
        p1: Point3,
        p2: Point3,
        tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        let e1 = p1 - p0;
        let e2 = p2 - p0;
        let cross = e1.cross(&e2);
        if cross.length() < tolerance.linear() {
            return Err(CoreError::PreconditionFailed(
                "sheet_triangle: three points are colinear".into(),
            ));
        }
        Self::sheet_polygon(&[p0, p1, p2], tolerance)
    }

    /// Construct a planar polygon sheet body from `n ≥ 3` vertex
    /// positions. Points must be coplanar; colinear / degenerate
    /// configurations are rejected by the downstream attach and
    /// validate passes.
    ///
    /// The front face's outward normal is derived from the first
    /// three non-colinear vertices via the right-hand rule.
    ///
    /// # Errors
    ///
    /// - [`CoreError::PreconditionFailed`] if `vertices.len() < 3`.
    /// - [`CoreError::Validation`] if any downstream Euler operator
    ///   or validation pass fails.
    pub fn sheet_polygon(vertices: &[Point3], tolerance: Tolerance) -> Result<Self, CoreError> {
        if vertices.len() < 3 {
            return Err(CoreError::PreconditionFailed(format!(
                "sheet_polygon: need ≥ 3 vertices, got {}",
                vertices.len()
            )));
        }

        let mut brep = BRep::new(tolerance);

        // 1. mvfs at the first vertex.
        let mvfs = brep.mvfs(vertices[0])?;

        // 2. mev chain through vertices[1..n-1] (leaving the last
        //    vertex to be closed by the final mef). After the
        //    chain, the outer loop holds 2(n-1) half-edges: n-1
        //    forward from vertices[0] → vertices[n-1] plus n-1
        //    back, all on the same face.
        let first_mev = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, vertices[1])?;
        let mut last_twin = first_mev.half_edge_twin;
        for v in vertices.iter().skip(2) {
            let mev = brep.mev(last_twin, *v)?;
            last_twin = mev.half_edge_twin;
        }

        // 3. mef to close the loop, splitting the face into a
        //    "front" face (keeping mvfs.face) and a "backside" face
        //    (the new face returned by mef). The mef result's
        //    `face` is the backside face because mef keeps the old
        //    face on the first argument's side of the split.
        brep.mef(last_twin, first_mev.half_edge)?;

        // Attach planar geometry to the front face and straight-line
        // curves to every edge. `attach_geometry` walks the BRep
        // and inserts the standard analytic geometry for faces
        // produced by Euler operators on planar loops.
        brep.attach_geometry()?;
        // Pcurves on the planar face's half-edges.
        brep.attach_planar_pcurves()?;

        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}

#[cfg(test)]
#[path = "sheet_tests.rs"]
mod tests;
