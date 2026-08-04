//! Blend region boundary identification.
//!
//! At a vertex where multiple fillets meet, this module computes the closed
//! boundary loop around the vertex. The boundary alternates between fillet-end
//! arcs (the transverse edge of each fillet surface near the vertex) and
//! original face edge segments between adjacent fillets.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::boundary_hex::{build_face_edge_bridge_between, reorient_arcs_for_hex_closure};
use super::boundary_saddle_stations::apply_saddle_station_pushback;
use crate::fillet_surface::FilletSurfaceResult;

// Re-export sampling types/functions so existing call sites that use
// `super::boundary::sample_boundary_by_arc_length` /
// `super::boundary::BoundarySample` continue to resolve.  The
// implementation lives in `super::boundary_sampling` to keep this
// file under the 500-line cap.
pub use super::boundary_sampling::{BoundarySample, sample_boundary_by_arc_length};

/// Boundary of a vertex blend region.
#[derive(Debug, Clone)]
pub struct BlendBoundary {
    /// The vertex being blended.
    pub vertex: VertexHandle,
    /// Position of the vertex (saved before vertex removal by blend insertion).
    pub vertex_position: Point3,
    /// Boundary curve segments (alternating fillet-end arcs and face edge segments).
    pub segments: Vec<BoundarySegment>,
    /// §5.3 item #6 — names of support surfaces where the
    /// surface-hugging bridge construction failed and the code
    /// fell back to a straight-line Bezier. Used by
    /// `insert_vertex_blend` to emit
    /// [`BlendWarning::SurfaceBridgeFallback`](crate::error::BlendWarning::SurfaceBridgeFallback).
    /// Empty in the common all-planar case.
    pub bridge_fallbacks: Vec<&'static str>,
    /// Phase 6c-A: `true` when sphere-fit substitution succeeded
    /// for this boundary. Phase 6c-D's edge-wiring helper reads
    /// this flag to decide whether to overwrite `edge.curve` —
    /// without it, the helper would propagate the original
    /// hermite-cubic synthesis (stored in `segments` when
    /// substitution wasn't performed) into BRep edges, *worsening*
    /// the cap-outer-boundary metric on fixtures where Phase 6c-A
    /// declined to run.
    pub sphere_fit_substituted: bool,
    /// Phase 6c-A continuation §4.2: per-fillet `(cap_centre_k,
    /// r_k)` pairs in the same order as `FilletEnd` segments
    /// (i.e., index `k` matches `segments[2*k]`'s
    /// `fillet_index`). Populated when sphere-fit substitution
    /// succeeds; empty otherwise. Downstream cap construction
    /// can use these to assign each sub-face a per-fillet
    /// `Sphere` surface — perfectly G1-tangent to that fillet's
    /// cylinder by construction, sidestepping the L-M solver
    /// stall on variable-radius corners.
    pub per_fillet_caps: Vec<(Point3, f64)>,
    /// Phase 6c-D Stage 1b: per-sub-patch hex-side outer boundary
    /// CPs, derived from the substituted `FilletEnd` (cap-tangent
    /// arc on `sphere ∩ cylinder`) and `FaceEdge` (bridge across
    /// slant face) curves and degree-elevated to 9.
    ///
    /// `hex_boundary_cps_deg9[k] = (row0, col9)` where:
    /// * `row0` — deg-9 Bernstein CPs of the sub-patch k's left
    ///   hex edge (`u=0`, `v ∈ [0, 1]`), sourced from
    ///   `segments[2*k]` (FilletEnd[k]).
    /// * `col9` — deg-9 Bernstein CPs of sub-patch k's top hex
    ///   edge (`v=1`, `u ∈ [0, 1]`), sourced from
    ///   `segments[2*k + 1]` (FaceEdge[k]).
    ///
    /// Length = N (one entry per fillet / sub-patch). Populated
    /// when sphere-fit substitution succeeds; empty otherwise.
    /// Downstream cap construction reads this when promoting to
    /// the K-P bi-nonic ring builder, where it bypasses the
    /// `k_p_boundary` Hermite synthesis that would otherwise
    /// produce a corner-Hermite blend instead of the analytic
    /// `sphere ∩ cylinder` arc.
    ///
    /// See `docs/decisions/2026-04-29-phase-6c-d-stage-1-boundary-injection.md`.
    pub hex_boundary_cps_deg9: Vec<([Point3; 10], [Point3; 10])>,
    /// Phase 6c-D Stage 4.7c: per-fillet analytic arc parameters
    /// for the substituted FilletEnd.  Populated when sphere-fit
    /// substitution succeeds AND the special case applies (cap
    /// centre on cylinder axis, cap_radius = cyl_radius — pyramid
    /// path).  Empty otherwise.
    ///
    /// Used by the wire (`sphere_fit_edge_wire`) to write
    /// `CurveGeometry::Ellipse` (with the cylinder's exact
    /// orthonormal basis) to BRep edges instead of a NURBS
    /// approximation, plus `Pcurve::Linear` with parameter range
    /// matching the analytic angle.  Gives bit-exact agreement
    /// between cap-side and twin (Cyl) surface evaluations on
    /// cap-vs-Cyl edges.
    pub fillet_end_arcs: Vec<FilletEndArcParams>,
}

/// Per-fillet analytic arc parameters captured during sphere-fit
/// substitution.  Stores the centre, cylinder basis, radius, and
/// angular range needed to reconstruct the substituted FilletEnd
/// as `CurveGeometry::Ellipse` (parameterised by angle directly,
/// matching the cylinder's UV trajectory).
#[derive(Debug, Clone)]
pub struct FilletEndArcParams {
    /// Arc centre (= cap centre, on the cylinder axis).
    pub center: Point3,
    /// First basis vector (length = radius), aligned with the
    /// cylinder surface's `perp1` from `orthonormal_basis(axis)`.
    pub u_basis: popkern_core::geometry::Vec3,
    /// Second basis vector (length = radius), aligned with `perp2`.
    pub v_basis: popkern_core::geometry::Vec3,
    /// Start angle of the substituted arc (radians).
    pub start_angle: f64,
    /// End angle of the substituted arc (radians).
    pub end_angle: f64,
    /// Axial coordinate on the cylinder (constant for this arc).
    pub v_const: f64,
}

/// A single segment of a blend boundary.
#[derive(Debug, Clone)]
pub enum BoundarySegment {
    /// End arc of a fillet surface at the vertex.
    FilletEnd {
        /// Index of the fillet in the global fillet list.
        fillet_index: usize,
        /// The boundary curve segment.
        curve: NurbsCurve,
    },
    /// Edge segment on an original face between two adjacent fillets.
    FaceEdge {
        /// The face this segment lies on.
        face: FaceHandle,
        /// The boundary curve segment.
        curve: NurbsCurve,
    },
}

/// Per-call options for [`compute_blend_boundary_with_options`].
///
/// Today's only field, `disable_sphere_fit_substitution`, lets the
/// caller override the global `POPKERN_SPHERE_FIT_BOUNDARY` env-var
/// gate.  Used by §6c-X.A.7's auto-dispatch: at uniform-radius
/// apices the sphere-fit substitution produces degenerate cap-loop
/// edges that the apex vertex-relocation repair can't resolve cleanly
/// (see ADR / journal entries, 2026-05-01), so the caller
/// suppresses substitution and lets the original rolling-ball
/// topology stand.
#[derive(Debug, Default, Clone, Copy)]
pub struct BlendBoundaryOptions {
    /// Force-disable sphere-fit substitution for this call,
    /// regardless of the `POPKERN_SPHERE_FIT_BOUNDARY` env var.
    /// `None` honours the env var (today's default behaviour).
    pub disable_sphere_fit_substitution: Option<bool>,
    /// PK-234: saddle (`OppositeSide`) station pushback. The default
    /// closest-to-vertex arc stations are only valid when every
    /// bridge region exists beyond the strips' tangency rails — true
    /// at convex corners, FALSE at a saddle's rail-crossing
    /// junctions, where each neighbor strip consumes the shared
    /// support band up to the rails' crossing (measured on the
    /// notched-cube saddle: the default station x=0.9875 sits past
    /// the carved rail end at x=0.9, and the emitted bridge runs
    /// through the neighbor strip's band, 3.3e-2 off the surface
    /// that actually owns the region). With this flag the arc
    /// stations are pushed back past every in-strip rail crossing
    /// with a setback margin, so all junctions land on real rails
    /// and all bridges cross real face regions.
    pub saddle_station_pushback: bool,
}

/// Compute the blend boundary at a vertex where multiple fillets meet.
///
/// Collects the end arcs of each fillet surface near the vertex, plus the
/// original face edge segments between fillets. These form a closed loop
/// around the vertex that defines the blend region.
///
/// # Errors
///
/// Returns an error if any fillet index is out of range or if the vertex
/// position cannot be determined.
pub fn compute_blend_boundary(
    vertex: VertexHandle,
    fillet_indices: &[usize],
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    tolerance: &Tolerance,
) -> Result<BlendBoundary, String> {
    compute_blend_boundary_with_options(
        vertex,
        fillet_indices,
        fillet_surfaces,
        brep,
        tolerance,
        BlendBoundaryOptions::default(),
    )
}

/// Variant of [`compute_blend_boundary`] that accepts per-call
/// [`BlendBoundaryOptions`].  Today's only option lets the caller
/// override the global `POPKERN_SPHERE_FIT_BOUNDARY` env-var
/// substitution gate.
pub fn compute_blend_boundary_with_options(
    vertex: VertexHandle,
    fillet_indices: &[usize],
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    tolerance: &Tolerance,
    options: BlendBoundaryOptions,
) -> Result<BlendBoundary, String> {
    if fillet_indices.is_empty() {
        return Err("no fillet indices provided for blend boundary".into());
    }

    // Resolve the vertex position.
    let vertex_data = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found in BRep"))?;
    let vertex_pos = vertex_data.position;

    // Build each fillet-end arc and remember its two 3D endpoints.
    let n = fillet_indices.len();
    let mut fillet_arcs: Vec<NurbsCurve> = Vec::with_capacity(n);
    let mut fillet_endpoints: Vec<[Point3; 2]> = Vec::with_capacity(n);
    for &fi in fillet_indices {
        let fillet = fillet_surfaces
            .get(fi)
            .ok_or_else(|| format!("fillet index {fi} out of range"))?;
        let end_curve = extract_fillet_end_near_vertex(fillet, &vertex_pos, tolerance)?;
        let (t0, t1) = end_curve.parameter_range();
        let ep = [end_curve.point_at(t0), end_curve.point_at(t1)];
        fillet_arcs.push(end_curve);
        fillet_endpoints.push(ep);
    }

    // Soft-fail: an unpushable station keeps the default boundary
    // (today's behavior — the saddle ring carve then declines with
    // the junction mismatch recorded) instead of killing the whole
    // vertex blend into a Mitred fallback.
    if options.saddle_station_pushback
        && let Err(reason) = apply_saddle_station_pushback(
            fillet_indices,
            fillet_surfaces,
            &vertex_pos,
            &mut fillet_arcs,
            &mut fillet_endpoints,
        )
        && std::env::var_os("PK_CARVE_TRACE").is_some()
    {
        eprintln!("  [saddle-stations] pushback declined: {reason}");
    }

    // Collect the supporting surface of every face incident to the
    // vertex. §5.3 item #6 lifts the prior planes-only restriction
    // — the face-edge-bridge builder now handles any C² surface
    // (Plane, Cylinder, Sphere, Cone, Torus, NURBS, Gregory) via
    // [`build_surface_bridge`](super::surface_bridge::build_surface_bridge).
    let faces_at_vertex = faces_at_vertex_set(brep, vertex);
    let face_supports: Vec<(FaceHandle, popkern_core::geometry::SurfaceGeometry)> = faces_at_vertex
        .iter()
        .filter_map(|&fh| {
            let face = brep.faces.get(fh)?;
            let sof = face.surface.as_ref()?;
            Some((fh, sof.geometry.clone()))
        })
        .collect();

    // Emit alternating FilletEnd / FaceEdge segments. For
    // hexagonal closure the bridge between arc[k] and arc[k+1]
    // MUST connect `arc[k].end` to `arc[k+1].start` (not
    // arbitrary endpoint pair) — anything else breaks the
    // closed-loop invariant. Each fillet's arc endpoints come
    // from the trim-curve closest-point construction in
    // `extract_fillet_end_near_vertex`: index 0 is the arc start
    // (UV-side from trim_curve_1's nearest point) and index 1 is
    // the arc end (trim_curve_2's nearest point).
    //
    // `try_reorient_arcs_for_closure` re-orders fillet indexing
    // so consecutive arcs actually share a nearby endpoint via a
    // cube face — the caller's `fillet_indices` order may not be
    // CCW around the vertex.
    let mut segments = Vec::with_capacity(n * 2);
    let mut bridge_fallbacks: Vec<&'static str> = Vec::new();
    // §5.3 item #6 — previously the hex-closure path required
    // `all_planar`. Now it runs whenever every face at the vertex
    // has a defined supporting surface (any variant).
    let has_face_supports =
        !face_supports.is_empty() && face_supports.len() == faces_at_vertex.len();

    let sphere_fit_enabled = match options.disable_sphere_fit_substitution {
        Some(true) => false,
        Some(false) => true,
        None => std::env::var("POPKERN_SPHERE_FIT_BOUNDARY")
            .map(|v| v == "1")
            .unwrap_or(false),
    };

    if has_face_supports && n >= 2 {
        let reorder = reorient_arcs_for_hex_closure(
            fillet_indices,
            &fillet_arcs,
            &fillet_endpoints,
            &face_supports,
            tolerance.linear().max(1e-6),
        );
        if let Some((order, flip, arc_by_k)) = reorder {
            // `arc_by_k[k]` is the arc to emit at position `k`,
            // already re-oriented; `flip[k]` records whether we
            // reversed arc `order[k]` (so the bridge builder can
            // use the correct endpoints).
            for k in 0..n {
                let fi = fillet_indices[order[k]];
                segments.push(BoundarySegment::FilletEnd {
                    fillet_index: fi,
                    curve: arc_by_k[k].clone(),
                });
                let next_k = (k + 1) % n;
                // Bridge connects current arc's end to next arc's
                // start (post-reorientation endpoints).
                let current_endpoints = if flip[k] {
                    [fillet_endpoints[order[k]][1], fillet_endpoints[order[k]][0]]
                } else {
                    fillet_endpoints[order[k]]
                };
                let next_endpoints = if flip[next_k] {
                    [
                        fillet_endpoints[order[next_k]][1],
                        fillet_endpoints[order[next_k]][0],
                    ]
                } else {
                    fillet_endpoints[order[next_k]]
                };
                let arc_end = current_endpoints[1];
                let next_start = next_endpoints[0];
                if let Some((bridge, fallback_kind)) = build_face_edge_bridge_between(
                    &arc_end,
                    &next_start,
                    &face_supports,
                    tolerance.linear().max(1e-6),
                ) {
                    segments.push(bridge);
                    if let Some(kind) = fallback_kind {
                        bridge_fallbacks.push(kind);
                    }
                }
            }
            let mut boundary = BlendBoundary {
                vertex,
                vertex_position: vertex_pos,
                segments,
                bridge_fallbacks,
                sphere_fit_substituted: false,
                per_fillet_caps: Vec::new(),
                hex_boundary_cps_deg9: Vec::new(),
                fillet_end_arcs: Vec::new(),
            };
            if sphere_fit_enabled {
                boundary.sphere_fit_substituted =
                    super::sphere_fit_boundary::apply_sphere_fit_substitution(
                        &mut boundary,
                        fillet_surfaces,
                        brep,
                        tolerance,
                    )
                    .is_ok();
            }
            return Ok(boundary);
        }
        // Fall through to the FilletEnd-only path if reorientation
        // couldn't determine a valid hex order (e.g., some arcs
        // don't share faces with neighbours — non-planar geometry,
        // degenerate corner, etc.). This preserves pre-Round 8
        // behaviour.
    }

    for k in 0..n {
        let fi = fillet_indices[k];
        segments.push(BoundarySegment::FilletEnd {
            fillet_index: fi,
            curve: fillet_arcs[k].clone(),
        });
    }

    Ok(BlendBoundary {
        vertex,
        vertex_position: vertex_pos,
        segments,
        bridge_fallbacks,
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    })
}

/// Collect the set of face handles incident to `vertex` by walking
/// all half-edges originating from it.
fn faces_at_vertex_set(brep: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
    let mut seen: Vec<FaceHandle> = Vec::new();
    for he in brep.half_edges.values() {
        if he.origin != vertex {
            continue;
        }
        let Some(loop_data) = brep.loops.get(he.loop_) else {
            continue;
        };
        let face = loop_data.face;
        if !seen.contains(&face) {
            seen.push(face);
        }
    }
    seen
}

/// Extract the end curve of a fillet surface closest to a vertex.
///
/// Takes the trim curves of the fillet and finds the endpoint nearest to the
/// vertex position, then builds a short transverse arc between the two trim
/// curve endpoints.
fn extract_fillet_end_near_vertex(
    fillet: &FilletSurfaceResult,
    vertex_pos: &Point3,
    _tolerance: &Tolerance,
) -> Result<NurbsCurve, String> {
    // §5.3 Round 7 C4b.1a — find the closest POINT on each trim
    // curve to the vertex (not just one of the two endpoints).
    // The fillet's trim curves typically span the entire filleted
    // edge including some internal extended-rolling-ball overshoot
    // (e.g., from `x = -1` to `x = 2` on a unit cube at `x ∈ [0, 1]`),
    // so the nearest trim-curve points to a vertex at the edge
    // endpoint are in the INTERIOR of the trim curves, not at
    // their parametric endpoints. The legacy `d_start` / `d_end`
    // selection picked whichever pair was "less far" and would
    // land on the far-overshoot plane, producing boundary arcs
    // geometrically disconnected from the vertex.
    let p1 = closest_point_on_nurbs_curve(&fillet.trim_curve_1, vertex_pos, 64);
    let p2 = closest_point_on_nurbs_curve(&fillet.trim_curve_2, vertex_pos, 64);

    // Project both points onto the fillet surface so we can trace
    // an ARC ALONG THE SURFACE between them (the geometrically
    // correct fillet-end arc) instead of a straight line through
    // space.
    let Some((u1, v1)) = project_point_to_uv(&fillet.surface, &p1) else {
        return Err(format!("failed to project p1 = {p1:?} onto fillet surface"));
    };
    let Some((u2, v2)) = project_point_to_uv(&fillet.surface, &p2) else {
        return Err(format!("failed to project p2 = {p2:?} onto fillet surface"));
    };

    // Periodicity handling: for Cylinder / Cone / Torus the `u`
    // parameter is angular (period `TAU`). Walking linearly from
    // `u1` to `u2` may take the long way around the cylinder; pick
    // the shorter arc direction.
    let (u_periodicity, _) = fillet.surface.periodicity();
    let mut du = u2 - u1;
    if let Some(period) = u_periodicity
        && du.abs() > period * 0.5
    {
        du -= period * du.signum();
    }

    // Exact chart-boundary iso first (swept strips' terminal
    // sections — exact for circular AND conic profiles); the Hermite
    // cubic stays the general-surface fallback.
    if let Some(exact) = super::boundary_hex::exact_boundary_iso_between(&fillet.surface, &p1, &p2)
    {
        return Ok(exact);
    }
    super::boundary_hex::fit_hermite_cubic_arc_on_surface(&fillet.surface, u1, v1, u2, v2, du)
}

/// Find the closest point on `curve` to `target` by sampling at
/// `n_samples + 1` uniform parameter values and picking the
/// minimum-distance sample. Chord-length-accurate to about
/// `curve_length / n_samples`; for our N=3-corner fixtures with
/// `n_samples = 64`, accuracy is sub-millimetre on unit-scale
/// models — well under the rolling-ball radius.
///
/// Could be improved with Newton polishing on the best sample,
/// but the extra precision hasn't been needed in practice.
fn closest_point_on_nurbs_curve(curve: &NurbsCurve, target: &Point3, n_samples: usize) -> Point3 {
    let (t0, t1) = curve.parameter_range();
    let mut best_d = f64::INFINITY;
    let mut best_p = curve.point_at(t0);
    let samples = n_samples.max(1);
    for i in 0..=samples {
        let t = t0 + (t1 - t0) * (i as f64) / (samples as f64);
        let p = curve.point_at(t);
        let d = p.distance_to(target);
        if d < best_d {
            best_d = d;
            best_p = p;
        }
    }
    best_p
}
