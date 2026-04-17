//! §5.3 Round 8 C3 — three-quad decomposition of a hexagonal
//! vertex-blend boundary into 3 rectangular Gregory sub-patches.
//!
//! Given a closed hexagonal [`BlendBoundary`] (3 FilletEnd + 3
//! FaceEdge alternating, as emitted by Round 8 C1), this module
//! computes a central point `C` and constructs 3
//! [`RectangularGregory`] patches meeting at `C`. Each sub-patch
//! covers two consecutive hex sides (one arc + one bridge) plus
//! two **spoke** curves from `C` to adjacent hex vertices.
//!
//! ```text
//!                 V1
//!                  /\
//!         arc0    /  \    bridge0
//!                /    \
//!           V0 /------\ V2
//!             /    C    \
//!    bridge2 /           \ arc1
//!           / ----------- \
//!         V5      ...      V3
//!             arc2 / bridge1
//! ```
//!
//! Each quad has 4 CCW boundaries:
//! 1. `segment[2k]` — the arc boundary
//! 2. `segment[2k+1]` — the face-edge bridge
//! 3. Spoke from hex-vertex `V_{2k+2 mod 6}` to `C`
//! 4. Spoke from `C` back to hex-vertex `V_{2k}`
//!
//! Three spokes (V0↔C, V2↔C, V4↔C) are each shared between two
//! sub-patches — a future BRep-surgery commit will promote them
//! to real half-edges so the three Gregory patches land as three
//! faces meeting at a new `C` vertex.
//!
//! **Scope caveat.** This module only produces the decomposition;
//! it does not modify topology. Callers that want to replace the
//! single merged vertex-blend face with three sub-faces need to
//! perform the topology surgery separately (R8 C4+ follow-up).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::RectangularGregory;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::gregory::fit_cubic_bezier_through;
use crate::fillet_surface::FilletSurfaceResult;

/// Three rectangular Gregory sub-patches meeting at a central
/// point, decomposing a hexagonal vertex-blend boundary.
#[allow(dead_code)] // diagnostic / follow-up: wired in R8 C4+
#[derive(Debug)]
pub(crate) struct ThreeQuadDecomposition {
    /// The central point shared by all three sub-patches.
    pub central_point: Point3,
    /// Sub-patches in hex-CCW order starting from arc 0.
    pub sub_patches: [RectangularGregory; 3],
}

/// Build a three-quad decomposition of a hexagonal boundary.
///
/// The boundary must have exactly 6 segments alternating
/// `FilletEnd` (even indices) / `FaceEdge` (odd indices). Returns
/// `None` if the structure doesn't match or any sub-patch's
/// cubic-Bezier boundary fit fails the corner-coincidence check.
///
/// `fillet_surfaces` are the incident fillet strips (one per
/// `FilletEnd` segment), passed through so hex-corner inner CPs
/// can be projected onto the strip tangent plane for G1 parity
/// with the adjacent fillet surfaces (§5.3 C4b.2 twist
/// extraction). Pass an empty slice to skip that refinement and
/// fall back to the pure-parallelogram minimum-twist baseline.
#[allow(dead_code)] // wired in R8 C4+
pub(crate) fn build_three_quad_decomposition(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<ThreeQuadDecomposition> {
    if boundary.segments.len() != 6 {
        return None;
    }
    // Verify alternation.
    for (i, seg) in boundary.segments.iter().enumerate() {
        match (i % 2, seg) {
            (0, BoundarySegment::FilletEnd { .. }) | (1, BoundarySegment::FaceEdge { .. }) => {}
            _ => return None,
        }
    }

    // Extract 6 hex vertices (one per segment start).
    let mut hex_vertices: Vec<Point3> = Vec::with_capacity(6);
    for seg in &boundary.segments {
        let curve = match seg {
            BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => {
                curve
            }
        };
        let (t0, _) = curve.parameter_range();
        hex_vertices.push(curve.point_at(t0));
    }
    if hex_vertices.len() != 6 {
        return None;
    }

    // Central point = centroid of the 6 hex vertices. Inside the
    // hex for convex boundaries (always true for N=3 trihedral
    // corners where the bridge-arc-bridge chains form a ~convex
    // hexagon in the offset-direction plane).
    let c = Point3::new(
        hex_vertices.iter().map(|p| p.x).sum::<f64>() / 6.0,
        hex_vertices.iter().map(|p| p.y).sum::<f64>() / 6.0,
        hex_vertices.iter().map(|p| p.z).sum::<f64>() / 6.0,
    );

    let mut sub_patches: Vec<RectangularGregory> = Vec::with_capacity(3);
    for k in 0..3usize {
        let b0_curve = match &boundary.segments[2 * k] {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            _ => return None,
        };
        let b1_curve = match &boundary.segments[2 * k + 1] {
            BoundarySegment::FaceEdge { curve, .. } => curve,
            _ => return None,
        };
        let v_start = hex_vertices[2 * k];
        let v_end = hex_vertices[(2 * k + 2) % 6];

        let spoke_v_end_to_c = build_spoke(v_end, c)?;
        let spoke_c_to_v_start = build_spoke(c, v_start)?;

        let b0 = fit_cubic_bezier_through(b0_curve);
        let b1 = fit_cubic_bezier_through(b1_curve);
        let b2 = fit_cubic_bezier_through(&spoke_v_end_to_c);
        let b3 = fit_cubic_bezier_through(&spoke_c_to_v_start);

        // C4b.2 / R8 C8c proper — parallelogram minimum-twist
        // inner CPs with fillet-strip tangent-plane projection
        // at hex corners. Replaces the pure-bilinear no-twist
        // baseline. See `compute_parallelogram_inner_cps` for
        // the derivation. This is the G1-at-hex-corners half of
        // the Chiyokura–Kimura construction; the central vertex
        // stays as pure parallelogram until the cross-patch G1
        // solve (C4b.2b) lands — that caveat is documented in
        // `ThreeQuadDecomposition`.
        //
        // Sub-patch k's hex corners (0, 1, 2) sit on fillet
        // strips — specifically:
        //   corner 0 = V_{2k}        on fillet[k-1 mod 3] and fillet[k]
        //   corner 1 = V_{2k+1}      on fillet[k]
        //   corner 2 = V_{2k+2}      on fillet[k] and fillet[k+1 mod 3]
        // Corner 3 = C is interior (no fillet strip).
        let strip_current = fillet_surfaces.get(k);
        let strip_prev = fillet_surfaces.get((k + 2) % 3);
        let strip_next = fillet_surfaces.get((k + 1) % 3);
        let corner_strips: [Option<&FilletSurfaceResult>; 4] =
            [strip_prev, strip_current, strip_next, None];

        let inner_cps = compute_parallelogram_inner_cps([&b0, &b1, &b2, &b3], corner_strips);
        let patch =
            RectangularGregory::from_boundary_beziers([b0, b1, b2, b3], inner_cps, 1e-3).ok()?;
        sub_patches.push(patch);
    }

    if sub_patches.len() != 3 {
        return None;
    }
    let mut it = sub_patches.into_iter();
    let p0 = it.next()?;
    let p1 = it.next()?;
    let p2 = it.next()?;
    Some(ThreeQuadDecomposition {
        central_point: c,
        sub_patches: [p0, p1, p2],
    })
}

/// C4b.2 — compute minimum-twist inner-CP pairs from the 4
/// boundary Beziers using the parallelogram extension formula,
/// then project hex-corner inner CPs onto the incident fillet
/// strip's tangent plane.
///
/// **Parallelogram baseline.** For a bicubic Bezier patch, the
/// interior CP that captures each corner's tangent plane with
/// zero twist is the parallelogram completion:
/// `P_11 = P_01 + P_10 − P_00` (and rotations for the other three
/// corners). This is the classical "no-twist" construction from
/// Farin and Hoschek — equivalent to extending the two incident
/// boundary tangent vectors through the corner.
///
/// **Tangent-plane projection at hex corners.** Each hex corner
/// sits on a fillet strip surface. The boundary tangent *along
/// the arc* at the corner already lies in the strip tangent plane
/// (the arc is on the strip). But the boundary tangent *along the
/// spoke* does not — spokes are interior construction lines. The
/// parallelogram P therefore drifts off the strip tangent plane
/// by the spoke's normal component. Projecting back onto the plane
/// `P' = P − ((P − corner) · n̂) · n̂` restores G1 at the hex
/// corner without changing the boundary shape.
///
/// Sets `P_u = P_v` at every corner — equal candidates give a
/// pure-parallelogram (G1) patch, not a twisted (G2) one. The
/// cross-derivative `∂²S/∂u∂v` still differs per corner because
/// the parallelogram value itself differs; that's what makes this
/// strictly better than the 4-corner-bilinear baseline.
///
/// Corner 3 is the central vertex C — no incident fillet strip;
/// pure parallelogram (caveat documented on the type).
fn compute_parallelogram_inner_cps(
    boundaries: [&[Point3; 4]; 4],
    corner_strips: [Option<&FilletSurfaceResult>; 4],
) -> [(Point3, Point3); 4] {
    let b0 = boundaries[0];
    let b1 = boundaries[1];
    let b2 = boundaries[2];
    let b3 = boundaries[3];

    // Bezier grid slots (see `from_boundary_beziers` layout):
    // corner 0 = grid[0][0] = b0[0] = b3[3]
    // corner 1 = grid[0][3] = b0[3] = b1[0]
    // corner 2 = grid[3][3] = b1[3] = b2[0]
    // corner 3 = grid[3][0] = b2[3] = b3[0]
    let corners = [b0[0], b0[3], b1[3], b2[3]];

    // Parallelogram inner CP = adjacent_boundary_cp_along_u
    //                         + adjacent_boundary_cp_along_v
    //                         − corner
    // At each corner, "along u" and "along v" are the two incident
    // boundary CPs one step away from the corner:
    //
    //   corner 0 (grid[0][0]): grid[1][0] + grid[0][1] − grid[0][0]
    //                         = b3[2]    + b0[1]    − b0[0]
    //   corner 1 (grid[0][3]): grid[0][2] + grid[1][3] − grid[0][3]
    //                         = b0[2]    + b1[1]    − b0[3]
    //   corner 2 (grid[3][3]): grid[3][2] + grid[2][3] − grid[3][3]
    //                         = b2[1]    + b1[2]    − b1[3]
    //   corner 3 (grid[3][0]): grid[2][0] + grid[3][1] − grid[3][0]
    //                         = b3[1]    + b2[2]    − b2[3]
    let raw_parallelogram = [
        parallelogram_add(b3[2], b0[1], b0[0]),
        parallelogram_add(b0[2], b1[1], b0[3]),
        parallelogram_add(b2[1], b1[2], b1[3]),
        parallelogram_add(b3[1], b2[2], b2[3]),
    ];

    let mut inner_cps: [(Point3, Point3); 4] =
        [(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)); 4];
    for (i, slot) in inner_cps.iter_mut().enumerate() {
        let corner = corners[i];
        let p_raw = raw_parallelogram[i];
        let p_refined = if let Some(strip) = corner_strips[i] {
            project_onto_strip_tangent_plane(&strip.surface, corner, p_raw).unwrap_or(p_raw)
        } else {
            // Central vertex (or any corner with no incident
            // strip) keeps the raw parallelogram value. C4b.2b
            // cross-patch G1 solve will refine this later.
            p_raw
        };
        *slot = (p_refined, p_refined);
    }
    inner_cps
}

fn parallelogram_add(p_a: Point3, p_b: Point3, base: Point3) -> Point3 {
    Point3::new(
        p_a.x + p_b.x - base.x,
        p_a.y + p_b.y - base.y,
        p_a.z + p_b.z - base.z,
    )
}

/// Project `p` onto the tangent plane of `surface` at the
/// location nearest to `corner`. The tangent plane is
/// `tangent ⊥ normal` through `corner`; subtracting the normal
/// component of `(p − corner)` lands the result back in the
/// plane.
fn project_onto_strip_tangent_plane(
    surface: &SurfaceGeometry,
    corner: Point3,
    p: Point3,
) -> Option<Point3> {
    let (u, v) = project_point_to_uv(surface, &corner)?;
    let n_hat = surface.normal_at(u, v).normalized()?;
    let delta = Vec3::new(p.x - corner.x, p.y - corner.y, p.z - corner.z);
    let off_plane = delta.dot(&n_hat);
    Some(Point3::new(
        p.x - off_plane * n_hat.x,
        p.y - off_plane * n_hat.y,
        p.z - off_plane * n_hat.z,
    ))
}

/// Build a straight-line spoke as a degree-1 NURBS from `a` to `b`
/// through a midpoint. Returns `None` if the two points coincide.
fn build_spoke(a: Point3, b: Point3) -> Option<NurbsCurve> {
    if a.distance_to(&b) < 1e-15 {
        return None;
    }
    let mid = Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5);
    crate::center_curve::fit_nurbs_through_points(&[a, mid, b]).ok()
}
