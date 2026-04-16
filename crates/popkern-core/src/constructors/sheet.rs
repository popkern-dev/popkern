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

use crate::constructors::sheet_helpers::rewrite_edge_as_arc;
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

    /// Construct an "arched triangle" sheet body — a half-disk
    /// represented as a 3-vertex triangle where two edges are
    /// quarter-circle arcs and one edge is the straight diameter.
    ///
    /// Vertex layout for the upper-half disk in the `z = z_plane`
    /// plane:
    ///
    /// - `v0 = (−radius, 0, z_plane)` — left end of the diameter
    /// - `v1 = (0, radius, z_plane)` — top of the arc
    /// - `v2 = (radius, 0, z_plane)` — right end of the diameter
    ///
    /// Edges:
    /// - `v0 → v1` — quarter-circle arc (left upper quadrant)
    /// - `v1 → v2` — quarter-circle arc (right upper quadrant)
    /// - `v2 → v0` — straight diameter
    ///
    /// **Corner angles** (useful for testing 2D corner ops with
    /// mixed line-arc support):
    /// - `v0` and `v2`: sharp **90°** corners between the diameter
    ///   and the tangent of the arc (which is vertical at the
    ///   diameter endpoints)
    /// - `v1`: **180°** straight-through — the two adjacent arcs
    ///   share a common tangent at the top. 2D corner operations
    ///   correctly reject `v1` as "near-colinear," so it is not
    ///   a usable corner target.
    ///
    /// **Test-vehicle purpose.** This shape is the smallest sheet
    /// body in popkern with genuine line-arc sharp corners. It
    /// exists specifically so future §5.7 extensions (curved-edge
    /// chamfer and fillet) have a testable input without
    /// requiring parallel-edge topology between the same two
    /// vertices.
    ///
    /// # Errors
    ///
    /// - [`CoreError::PreconditionFailed`] if `radius ≤ 0`.
    /// - [`CoreError::Validation`] if any downstream Euler
    ///   operator or validation pass fails.
    pub fn sheet_arched_triangle(
        radius: f64,
        z_plane: f64,
        tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        if radius <= 0.0 || !radius.is_finite() {
            return Err(CoreError::PreconditionFailed(format!(
                "sheet_arched_triangle: radius must be positive and finite, got {radius}"
            )));
        }

        // Build the triangular topology with straight-line curves
        // on every edge first — `sheet_polygon` handles all the
        // Euler operators and planar surface attachment for us.
        let v0 = Point3::new(-radius, 0.0, z_plane);
        let v1 = Point3::new(0.0, radius, z_plane);
        let v2 = Point3::new(radius, 0.0, z_plane);
        let mut brep = Self::sheet_polygon(&[v0, v1, v2], tolerance)?;

        // Now rewrite two of the edges' curves from Line to Arc.
        // We need to locate the edges `v0 → v1` and `v1 → v2`.
        // The third edge `v2 → v0` stays as a Line — that's the
        // diameter.
        //
        // The circle the arcs lie on: centre at the origin of the
        // diameter midpoint (0, 0, z_plane), normal (0, 0, -1)
        // (matches the plane normal from `plane_from_points`,
        // which for our counter-clockwise walk of v0, v1, v2
        // computes `(v1 - v0) × (v2 - v0) = (0, 0, −2r²) → normal
        // (0, 0, −1)`), and our chosen `radius`.
        //
        // With `normal = (0, 0, −1)` the `circle_basis` helper
        // picks `u = (0, −1, 0)` and `v = (−1, 0, 0)` (see
        // `popkern_core::geometry::analytic_curves::circle_basis`).
        // Angles of the three vertices in that frame:
        // - v0 (−r, 0, 0) → `u·(−r,0,0)=0`, `v·(−r,0,0)=r` →
        //   atan2(r, 0) = π/2
        // - v1 (0,  r, 0) → `u·(0,r,0)=−r`, `v·(0,r,0)=0` →
        //   atan2(0, −r) = π
        // - v2 (r,  0, 0) → `u·(r,0,0)=0`, `v·(r,0,0)=−r` →
        //   atan2(−r, 0) = 3π/2 (or −π/2)
        //
        // So the arc `v0 → v1` spans angles `[π/2, π]` (increasing),
        // and the arc `v1 → v2` spans angles `[π, 3π/2]`
        // (also increasing). Beautiful monotonic parameter
        // ranges.
        let arc_center = Point3::new(0.0, 0.0, z_plane);
        let arc_normal = crate::geometry::Vec3::new(0.0, 0.0, -1.0);
        rewrite_edge_as_arc(&mut brep, v0, v1, arc_center, arc_normal, radius)?;
        rewrite_edge_as_arc(&mut brep, v1, v2, arc_center, arc_normal, radius)?;

        // The curves changed; re-run the full validate to catch
        // any watertight / manifold regressions before returning.
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }

    /// Construct a **biased** arched triangle — like
    /// [`Self::sheet_arched_triangle`] but with two *different* arc
    /// radii so the apex is a genuine non-tangent-matched arc-arc
    /// corner.
    ///
    /// Fixed vertex layout:
    /// - `v0 = (−1, 0, z_plane)`
    /// - `v1 = (0, 1.5, z_plane)` (apex)
    /// - `v2 = (1, 0, z_plane)`
    ///
    /// Edges:
    /// - `v0 → v1` — arc of radius `radius_left`, bulging upper-left
    ///   (centre on the chord's right-hand side)
    /// - `v1 → v2` — arc of radius `radius_right`, bulging upper-right
    /// - `v2 → v0` — straight diameter-like line at `y = 0`
    ///
    /// **Corner at `v1`** is the arc-arc test vehicle for §5.7
    /// 2D corner ops on curved supports: when `radius_left ==
    /// radius_right` the two arcs share a tangent at `v1` and the
    /// corner is degenerate; any non-equal pair yields a real
    /// (non-colinear) corner that fillet/chamfer can attack.
    ///
    /// Chord length from either `v0 → v1` or `v1 → v2` is
    /// `sqrt(1 + 1.5²) = sqrt(3.25) ≈ 1.803`, so both radii must
    /// be strictly greater than half that (`≈ 0.901`) for the
    /// geometry to be feasible.
    ///
    /// # Errors
    ///
    /// - [`CoreError::PreconditionFailed`] if either radius is
    ///   non-finite, non-positive, or smaller than half the chord
    ///   length.
    /// - [`CoreError::Validation`] on any downstream Euler /
    ///   validate failure.
    pub fn sheet_biased_arched_triangle(
        radius_left: f64,
        radius_right: f64,
        z_plane: f64,
        tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        const APEX_Y: f64 = 1.5;
        let chord_half = (1.0_f64 + APEX_Y * APEX_Y).sqrt() * 0.5;

        for (name, r) in [("radius_left", radius_left), ("radius_right", radius_right)] {
            if !r.is_finite() || r <= 0.0 {
                return Err(CoreError::PreconditionFailed(format!(
                    "sheet_biased_arched_triangle: {name} must be positive and finite, got {r}"
                )));
            }
            if r <= chord_half + 1.0e-12 {
                return Err(CoreError::PreconditionFailed(format!(
                    "sheet_biased_arched_triangle: {name} {r} must exceed chord/2 = {chord_half}"
                )));
            }
        }

        let v0 = Point3::new(-1.0, 0.0, z_plane);
        let v1 = Point3::new(0.0, APEX_Y, z_plane);
        let v2 = Point3::new(1.0, 0.0, z_plane);
        let mut brep = Self::sheet_polygon(&[v0, v1, v2], tolerance)?;

        // Compute arc centres. The perpendicular bisector of each
        // chord carries the centre at distance
        // `sqrt(r² − (chord/2)²)` from the chord midpoint. We
        // offset toward the **exterior** of the straight-triangle
        // walk so the arcs bulge **outward** away from the third
        // vertex — the shape is still sphere-topology and the
        // bottom "line" stays the diameter.
        let arc_normal = crate::geometry::Vec3::new(0.0, 0.0, -1.0);

        let (center_a, center_b) = biased_arch_centers(APEX_Y, radius_left, radius_right, z_plane);

        rewrite_edge_as_arc(&mut brep, v0, v1, center_a, arc_normal, radius_left)?;
        rewrite_edge_as_arc(&mut brep, v1, v2, center_b, arc_normal, radius_right)?;

        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}

/// Compute the arc centres for
/// [`BRep::sheet_biased_arched_triangle`]. Shared with
/// [`BRep::sheet_biased_arched_triangle`] so the test-file
/// constants can be kept in one place.
fn biased_arch_centers(
    apex_y: f64,
    radius_left: f64,
    radius_right: f64,
    z_plane: f64,
) -> (Point3, Point3) {
    let chord_half_sq = 0.25 * (1.0 + apex_y * apex_y);

    // Arc a: chord (-1, 0) → (0, apex_y). Midpoint (-0.5, apex_y/2).
    // Chord direction normalised (1, apex_y) / chord_len. Perpendicular
    // rotated 90° clockwise (that's the **right** side of the walk
    // from v0 to v1, which coincides with the exterior of the
    // counter-clockwise triangle walk v0 → v1 → v2). We place the
    // centre on that side to make the arc bulge outward on the
    // opposite side.
    let chord_len = (1.0_f64 + apex_y * apex_y).sqrt();
    let mid_a = Point3::new(-0.5, apex_y * 0.5, z_plane);
    let perp_a = crate::geometry::Vec3::new(apex_y / chord_len, -1.0 / chord_len, 0.0);
    let offset_a = (radius_left * radius_left - chord_half_sq).sqrt();
    let center_a = mid_a + perp_a * offset_a;

    let mid_b = Point3::new(0.5, apex_y * 0.5, z_plane);
    // Arc b: chord (0, apex_y) → (1, 0). Chord direction (1, -apex_y)/chord_len.
    // The exterior perpendicular (90° CW of the chord direction)
    // points to `(-apex_y, -1)/chord_len` — down-left.
    let perp_b = crate::geometry::Vec3::new(-apex_y / chord_len, -1.0 / chord_len, 0.0);
    let offset_b = (radius_right * radius_right - chord_half_sq).sqrt();
    let center_b = mid_b + perp_b * offset_b;

    (center_a, center_b)
}

#[cfg(test)]
#[path = "sheet_tests.rs"]
mod tests;
