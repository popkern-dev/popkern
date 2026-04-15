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
}

/// Locate the edge whose canonical half-edge runs from `p_from`
/// to `p_to` (within linear tolerance) and replace its `Line`
/// curve with an `Arc` on the specified circle.
///
/// Preserves the canonical half-edge's direction: the arc
/// parameterisation `(t_start, t_end)` matches the canonical
/// endpoint order.
fn rewrite_edge_as_arc(
    brep: &mut BRep,
    p_from: Point3,
    p_to: Point3,
    center: Point3,
    normal: crate::geometry::Vec3,
    radius: f64,
) -> Result<(), CoreError> {
    let edge_h = find_edge_between_points(brep, p_from, p_to, 1.0e-9).ok_or_else(|| {
        CoreError::PreconditionFailed(format!(
            "arched_triangle: cannot find edge between {p_from:?} and {p_to:?}"
        ))
    })?;

    // The canonical half-edge's origin gives the "first" vertex
    // of the edge in parameter-range order. We use that to
    // compute the angle at t_start.
    let canonical_he = brep.edges[edge_h].half_edge;
    let canonical_origin_vh = brep.half_edges[canonical_he].origin;
    let canonical_twin_h = brep.half_edges[canonical_he].twin;
    let canonical_twin_origin_vh = brep.half_edges[canonical_twin_h].origin;
    let canonical_origin_pos = brep.vertices[canonical_origin_vh].position;
    let canonical_twin_origin_pos = brep.vertices[canonical_twin_origin_vh].position;

    let (u_basis, v_basis) = circle_basis_local(normal);
    let t_start = angle_in_basis(canonical_origin_pos - center, u_basis, v_basis);
    let t_end_raw = angle_in_basis(canonical_twin_origin_pos - center, u_basis, v_basis);
    // Choose the "short" arc by ensuring `t_end` is within `π`
    // of `t_start` (add or subtract `2π` as needed). For our
    // arched-triangle case both endpoints are adjacent on the
    // circle so the natural difference is already `< π`.
    let t_end = wrap_nearest(t_start, t_end_raw);

    let arc = crate::geometry::curve::CurveGeometry::Arc {
        center,
        normal,
        radius,
        start_angle: t_start,
        end_angle: t_end,
    };
    brep.edges[edge_h].curve = Some(crate::geometry::curve::CurveOnEdge {
        geometry: arc,
        t_start,
        t_end,
    });
    Ok(())
}

fn find_edge_between_points(
    brep: &BRep,
    p_from: Point3,
    p_to: Point3,
    tol_eps: f64,
) -> Option<crate::topology::EdgeHandle> {
    for (eh, edge) in &brep.edges {
        let he_h = edge.half_edge;
        let twin_h = brep.half_edges[he_h].twin;
        let a = brep.vertices[brep.half_edges[he_h].origin].position;
        let b = brep.vertices[brep.half_edges[twin_h].origin].position;
        let matches = (a.distance_to(&p_from) < tol_eps && b.distance_to(&p_to) < tol_eps)
            || (a.distance_to(&p_to) < tol_eps && b.distance_to(&p_from) < tol_eps);
        if matches {
            return Some(eh);
        }
    }
    None
}

fn angle_in_basis(
    rel: crate::geometry::Vec3,
    u: crate::geometry::Vec3,
    v: crate::geometry::Vec3,
) -> f64 {
    rel.dot(&v).atan2(rel.dot(&u))
}

/// Wrap `t` so it's within `±π` of `anchor`. Used to pick the
/// short-arc parameter for an otherwise-ambiguous endpoint angle.
fn wrap_nearest(anchor: f64, t: f64) -> f64 {
    use std::f64::consts::{PI, TAU};
    let mut result = t;
    while result - anchor > PI {
        result -= TAU;
    }
    while result - anchor < -PI {
        result += TAU;
    }
    result
}

/// Local replica of the (private) `circle_basis` helper in
/// `analytic_curves`, used here so the arc parameterisation we
/// compute matches the frame `CurveGeometry::Arc::point_at`
/// evaluates in at runtime.
fn circle_basis_local(
    normal: crate::geometry::Vec3,
) -> (crate::geometry::Vec3, crate::geometry::Vec3) {
    use crate::geometry::Vec3;
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let perp = normal.cross(&axis);
    let u = perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

#[cfg(test)]
#[path = "sheet_tests.rs"]
mod tests;
