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
use popkern_core::geometry::{NurbsCurve, Point3};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::gregory::fit_cubic_bezier_through;

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
#[allow(dead_code)] // wired in R8 C4+
pub(crate) fn build_three_quad_decomposition(
    boundary: &BlendBoundary,
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

        let patch =
            RectangularGregory::from_boundary_beziers_no_twist([b0, b1, b2, b3], 1e-3).ok()?;
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

/// Build a straight-line spoke as a degree-1 NURBS from `a` to `b`
/// through a midpoint. Returns `None` if the two points coincide.
fn build_spoke(a: Point3, b: Point3) -> Option<NurbsCurve> {
    if a.distance_to(&b) < 1e-15 {
        return None;
    }
    let mid = Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5);
    crate::center_curve::fit_nurbs_through_points(&[a, mid, b]).ok()
}
