//! Phase 6c-C.2 — apex-region deletion *planning*.
//!
//! The surgery itself is non-trivial because the existing
//! cube_uniform topology has a structure where each rolling-ball
//! footprint vertex sits at the junction of TWO different fillet
//! cylinders.  A simple `kev` can't merge its two cap-side edges
//! (they bound different cylinder faces).  This module builds a
//! structured `DeletionPlan` that enumerates exactly what to
//! delete and what cap-tangent arcs to insert in their place.
//! Phase 6c-C.2b (separate commit) will execute the plan via
//! Euler ops; this commit only computes it.
//!
//! # Domain model (cube_uniform after the centroid fix)
//!
//! For a uniform-radius cube apex, the post-fillet topology has:
//! * **Cap face** (Sphere primitive) with a 6-edge outer loop
//!   alternating "sphere-fit endpoints" and "rolling-ball
//!   footprints":
//!   * 3 sphere-fit vertices V₁, V₃, V₅ at face-interior points
//!     `(r,0,r)`, `(r,r,0)`, `(0,r,r)` — on the cap sphere.
//!   * 3 rolling-ball-footprint vertices V₀, V₂, V₄ at cube-edge
//!     points `(0,0,r)`, `(r,0,0)`, `(0,r,0)` — *not* on the
//!     cap sphere.
//! * **6 cap-cylinder edges**, all of them straight cylinder
//!   *generators* (not cap-tangent arcs).  Pairs of generators
//!   meet at the rolling-ball-footprint vertices.
//! * **3 fillet cylinder faces**, each with a portion of its
//!   outer loop running through one footprint vertex via two
//!   generators + one cube-face boundary edge segment.
//!
//! # Target topology (post Phase 6c-C)
//!
//! * Cap face's outer loop becomes a triangle: V₁ → V₃ → V₅ → V₁
//!   with each edge a **cap-tangent half-circle** (the curve
//!   where the cap sphere meets one of the three fillet
//!   cylinders).
//! * Each fillet cylinder face's apex-side boundary becomes the
//!   single cap-tangent arc (no generators, no footprint vertex).
//! * Rolling-ball-footprint vertices V₀/V₂/V₄ are removed.
//!
//! # Deletion + insertion sequence
//!
//! 1. **For each footprint Vₖ**: enumerate its 2 cap-side
//!    half-edges and their twins.  The two edges bound the cap
//!    on one side; their twins bound *two different* cylinder
//!    faces.
//! 2. **Insert** 3 new cap-tangent arc edges between the
//!    sphere-fit vertices, one per fillet cylinder.  Curve
//!    geometry: half-circle in the appropriate plane (the
//!    sphere-cylinder intersection).
//! 3. **Replace** each cylinder face's apex-side outer-loop
//!    fragment (2 generators + a cube-face edge) with the new
//!    cap-tangent arc.  This involves `mef` to insert the new
//!    edge across the loop, then `kef` to remove the old
//!    fragment.
//! 4. **`kev`** each rolling-ball-footprint vertex.  At this
//!    point, each footprint is only incident to two edges (the
//!    obsolete cube-face boundary segment between adjacent
//!    cylinder trims) — those become a single merged edge that
//!    no longer borders the cap face.
//!
//! Phase 6c-C.2b will encode this sequence.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Structured plan for converting the existing 6-vertex cap
/// loop into the target 3-vertex sphere-fit topology.  Returned
/// by `plan_apex_region_deletion`.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub(crate) struct DeletionPlan {
    /// The cap face that this plan targets.
    pub(crate) cap_face: Option<FaceHandle>,
    /// Cap-sphere parameters (centre, radius). Drives the
    /// cap-tangent arc geometry for the new edges.
    pub(crate) cap_sphere: Option<(Point3, f64)>,
    /// Sphere-fit endpoint vertices that survive the rebuild.
    /// In cap-loop walking order; length N for the N-fillet
    /// apex (3 for cube uniform, 4 for square pyramid uniform).
    pub(crate) keep_vertices: Vec<VertexHandle>,
    /// Rolling-ball-footprint vertices to remove.  Each is
    /// incident to exactly 2 cap-cylinder generator edges
    /// (which become obsolete) and possibly cube-face-boundary
    /// edges that may need re-routing.
    pub(crate) delete_vertices: Vec<VertexHandle>,
    /// Cap-cylinder generator edges to remove (6 for cube
    /// uniform: two per cylinder).
    pub(crate) delete_cap_cylinder_generators: Vec<EdgeHandle>,
    /// New cap-tangent arc edges to insert.  Each entry is
    /// `(start_vertex, end_vertex, cylinder_surface_for_arc)`
    /// — the arc lies on (cap_sphere ∩ cylinder).
    pub(crate) insert_cap_tangent_arcs: Vec<(VertexHandle, VertexHandle, SurfaceGeometry)>,
}

impl DeletionPlan {
    /// True when the plan is complete enough to feed to
    /// 6c-C.2b's executor.  Empty plans (e.g. from non-uniform
    /// fixtures) return false; callers fall back to today's path.
    #[allow(dead_code)]
    pub(crate) fn is_actionable(&self) -> bool {
        self.cap_face.is_some()
            && self.cap_sphere.is_some()
            && !self.keep_vertices.is_empty()
            && self.delete_vertices.len() == self.delete_cap_cylinder_generators.len() / 2
            && self.insert_cap_tangent_arcs.len() == self.keep_vertices.len()
    }
}

/// Build the deletion plan for `cap_face` post-fillet.  Returns
/// an empty plan when the cap face's structure isn't the
/// canonical 6-edge alternating sphere-fit / rolling-ball-foot
/// topology (e.g. variable-radius corners take a different
/// shape and need a different plan).
#[allow(dead_code)]
pub(crate) fn plan_apex_region_deletion(brep: &BRep, cap_face: FaceHandle) -> DeletionPlan {
    let mut plan = DeletionPlan::default();

    // 1. Cap surface must be a Sphere primitive (centroid fix
    //    in `surface_fix.rs::try_sphere_for_corner` ensures
    //    this for uniform-radius corners).
    let Some(sof) = brep.faces.get(cap_face).and_then(|f| f.surface.as_ref()) else {
        return plan;
    };
    let (cap_centre, cap_radius) = match &sof.geometry {
        SurfaceGeometry::Sphere { center, radius } => (*center, *radius),
        _ => return plan,
    };
    plan.cap_face = Some(cap_face);
    plan.cap_sphere = Some((cap_centre, cap_radius));

    // 2. Walk the cap face's outer loop and classify each
    //    vertex as "sphere-fit" (lies on the cap sphere within
    //    tolerance) or "rolling-ball footprint" (off the
    //    sphere).  Tolerance: 1% of cap_radius.
    let Some(face_data) = brep.faces.get(cap_face) else {
        return plan;
    };
    let outer_loop = face_data.outer_loop;
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return plan;
    };
    let start = loop_data.half_edge;

    let mut loop_vertices: Vec<(VertexHandle, EdgeHandle)> = Vec::new();
    let mut cursor = start;
    for _ in 0..32 {
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        loop_vertices.push((he.origin, he.edge));
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    if loop_vertices.is_empty() {
        return plan;
    }

    let on_sphere_tol = (cap_radius * 0.01).max(1e-6);
    let on_sphere = |v: VertexHandle| -> bool {
        let Some(p) = brep.vertices.get(v).map(|vd| vd.position) else {
            return false;
        };
        (p.distance_to(&cap_centre) - cap_radius).abs() < on_sphere_tol
    };

    // 3. Build the keep / delete vertex sets and the cap-cylinder
    //    generator edges to delete (each footprint contributes
    //    its 2 incident loop edges).
    for (v, e) in &loop_vertices {
        if on_sphere(*v) {
            if !plan.keep_vertices.contains(v) {
                plan.keep_vertices.push(*v);
            }
        } else if !plan.delete_vertices.contains(v) {
            plan.delete_vertices.push(*v);
        }
        if !plan.delete_cap_cylinder_generators.contains(e) {
            plan.delete_cap_cylinder_generators.push(*e);
        }
    }

    // 4. For each consecutive pair of sphere-fit vertices in
    //    the loop walking order, insert one cap-tangent arc
    //    edge.  The cylinder surface for the arc is whichever
    //    fillet cylinder that pair of sphere-fit endpoints sits
    //    on (the sphere-cylinder tangent circle that contains
    //    both endpoints).  We don't pick the cylinder here —
    //    6c-C.2b looks it up by querying which cylinder face's
    //    outer loop currently contains both vertices.
    let n = plan.keep_vertices.len();
    for i in 0..n {
        let a = plan.keep_vertices[i];
        let b = plan.keep_vertices[(i + 1) % n];
        // Placeholder cylinder surface — 6c-C.2b derives the real
        // one from the BRep topology around the two vertices.
        let placeholder = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
            u_axis: popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0),
            v_axis: popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0),
        };
        plan.insert_cap_tangent_arcs.push((a, b, placeholder));
    }
    plan
}
