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
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregory, RectangularGregoryDeg5};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::gregory::{fit_cubic_bezier_through, fit_quintic_bezier_through};
use super::gregory_parallelogram::{
    compute_parallelogram_inner_cps, compute_parallelogram_inner_cps_deg5, sub3,
};
use super::gregory_vertex_enclosure::shirman_sequin_vertex_enclosure_solve;
use crate::fillet_surface::FilletSurfaceResult;
use crate::vertex_classifier::VertexBlendClassification;

/// `N` rectangular Gregory sub-patches meeting at a central point,
/// decomposing a `2N`-segment vertex-blend boundary (alternating
/// FilletEnd / FaceEdge). `N = 3` is the trihedral hex case; `N ≥
/// 4` covers pyramid apex, hub-and-spoke fittings, etc.
///
/// Each sub-patch is wrapped in [`GregoryPatch`] so cubic
/// (bidegree 3) and quintic (bidegree 5) sub-patches can flow
/// through the same dispatch pipeline.
#[allow(dead_code)] // diagnostic / follow-up
#[derive(Debug)]
pub(crate) struct NQuadDecomposition {
    /// The central point shared by all `N` sub-patches.
    pub central_point: Point3,
    /// Sub-patches in CCW order starting from arc 0. Length `N`.
    pub sub_patches: Vec<GregoryPatch>,
}

/// Back-compat alias for the N=3 trihedral case.
#[allow(dead_code)]
pub(crate) type ThreeQuadDecomposition = NQuadDecomposition;

/// Build an N-quad decomposition of a `2N`-segment boundary.
///
/// The boundary must have an even number `2N ≥ 6` of segments
/// alternating `FilletEnd` (even indices) / `FaceEdge` (odd
/// indices). Returns `None` if the structure doesn't match or any
/// sub-patch's cubic-Bezier boundary fit fails the corner-
/// coincidence check.
///
/// `fillet_surfaces` are the incident fillet strips (one per
/// `FilletEnd` segment), passed through so corner inner CPs can be
/// projected onto the strip tangent plane for G1 parity with the
/// adjacent fillet surfaces (§5.3 C4b.2 twist extraction). Pass an
/// empty slice to skip that refinement and fall back to the pure-
/// parallelogram minimum-twist baseline.
#[allow(dead_code)]
pub(crate) fn build_n_quad_decomposition(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<NQuadDecomposition> {
    build_n_quad_decomposition_with_degree(boundary, fillet_surfaces, false)
}

/// Same as [`build_n_quad_decomposition`] but with explicit control
/// over the Bezier degree of the sub-patches. `use_deg5 = true`
/// emits [`RectangularGregoryDeg5`] quintic sub-patches; `false`
/// emits cubic [`RectangularGregory`]. Quintic closes the ~2%
/// cubic-Bezier arc approximation error that caps §2.4.
#[allow(dead_code)]
pub(crate) fn build_n_quad_decomposition_with_degree(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    use_deg5: bool,
) -> Option<NQuadDecomposition> {
    build_n_quad_decomposition_with_options(
        boundary,
        fillet_surfaces,
        use_deg5,
        VertexBlendClassification::SameSideConvex,
    )
}

/// `build_n_quad_decomposition_with_degree` + classification hint so
/// the builder can route saddle (`OppositeSide`) configurations
/// through a saddle-specific central-point computation (midpoint of
/// the N fillet-end arc midpoints rather than the 2N hex-vertex
/// centroid, which sits incorrectly inside the convex side for
/// saddles).
#[allow(dead_code)]
pub(crate) fn build_n_quad_decomposition_with_options(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    use_deg5: bool,
    classification: VertexBlendClassification,
) -> Option<NQuadDecomposition> {
    let two_n = boundary.segments.len();
    if two_n < 6 || !two_n.is_multiple_of(2) {
        return None;
    }
    let n = two_n / 2;
    // Verify alternation.
    for (i, seg) in boundary.segments.iter().enumerate() {
        match (i % 2, seg) {
            (0, BoundarySegment::FilletEnd { .. }) | (1, BoundarySegment::FaceEdge { .. }) => {}
            _ => return None,
        }
    }

    // Extract the `2N` polygon vertices (one per segment start).
    let mut hex_vertices: Vec<Point3> = Vec::with_capacity(two_n);
    for seg in &boundary.segments {
        let curve = match seg {
            BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => {
                curve
            }
        };
        let (t0, _) = curve.parameter_range();
        hex_vertices.push(curve.point_at(t0));
    }
    if hex_vertices.len() != two_n {
        return None;
    }

    // Central point selection — classification-dependent.
    //
    // **SameSideConvex / SameSideConcave:** hex centroid (average
    // of the 2N hex vertices). The hex is convex so the centroid
    // is inside the blend — the sub-patches fan out from here.
    //
    // **OppositeSide (saddle):** §5.3 item #5. A saddle has 1+
    // convex fillet(s) + 1+ concave fillet(s) meeting at the
    // vertex. The convex arc midpoints sit on one side of the
    // saddle ridge; concave arc midpoints on the other. The
    // 2N-hex centroid is pulled off-ridge toward the convex side
    // (more positive-curvature contributions). Averaging only
    // the N **arc midpoints** (the FilletEnd segments — not the
    // bridge corners) is symmetric across the ridge because each
    // arc midpoint lies on the spine of its fillet, equidistant
    // from the saddle ridge regardless of sign.
    let denom_hex = two_n as f64;
    let hex_centroid = Point3::new(
        hex_vertices.iter().map(|p| p.x).sum::<f64>() / denom_hex,
        hex_vertices.iter().map(|p| p.y).sum::<f64>() / denom_hex,
        hex_vertices.iter().map(|p| p.z).sum::<f64>() / denom_hex,
    );
    let c = if matches!(classification, VertexBlendClassification::OppositeSide) {
        let arc_mids: Vec<Point3> = boundary
            .segments
            .iter()
            .filter_map(|s| match s {
                BoundarySegment::FilletEnd { curve, .. } => {
                    let (t0, t1) = curve.parameter_range();
                    Some(curve.point_at(0.5 * (t0 + t1)))
                }
                _ => None,
            })
            .collect();
        if arc_mids.len() == n {
            let denom_arc = n as f64;
            Point3::new(
                arc_mids.iter().map(|p| p.x).sum::<f64>() / denom_arc,
                arc_mids.iter().map(|p| p.y).sum::<f64>() / denom_arc,
                arc_mids.iter().map(|p| p.z).sum::<f64>() / denom_arc,
            )
        } else {
            hex_centroid
        }
    } else {
        hex_centroid
    };

    // Central-vertex normal. For convex configurations, the sum
    // of cross products of consecutive hex-vertex-to-centroid
    // vectors gives a consistent outward normal. For saddles the
    // sum can cancel to near-zero; in that case we fall back to
    // averaging the fillet-surface normals at the arc midpoints
    // (a robust choice that respects the saddle ridge direction).
    let n_c = central_vertex_normal(&hex_vertices, c)
        .or_else(|| saddle_normal_fallback(boundary, fillet_surfaces))
        .unwrap_or(Vec3::new(0.0, 0.0, 0.0));

    // Pass 1: build each sub-patch's boundaries and raw parallelogram
    // inner CPs (with C-K twist applied at corners 0–2). Corner 3
    // (central vertex) stays at its raw parallelogram position; the
    // Shirman-Séquin cross-patch solve in pass 2 corrects it.
    //
    // `use_deg5` selects between cubic (bidegree 3, 4 boundary CPs)
    // and quintic (bidegree 5, 6 boundary CPs) sub-patch construction.
    enum StagedBoundaries {
        Cubic([[Point3; 4]; 4]),
        Deg5([[Point3; 6]; 4]),
    }
    type SubPatchStaged = (StagedBoundaries, [(Point3, Point3); 4]);
    let mut per_patch: Vec<SubPatchStaged> = Vec::with_capacity(n);
    for k in 0..n {
        let b0_curve = match &boundary.segments[2 * k] {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            _ => return None,
        };
        let b1_curve = match &boundary.segments[2 * k + 1] {
            BoundarySegment::FaceEdge { curve, .. } => curve,
            _ => return None,
        };
        let v_start = hex_vertices[2 * k];
        let v_end = hex_vertices[(2 * k + 2) % two_n];

        let spoke_v_end_to_c = build_spoke(v_end, c)?;
        let spoke_c_to_v_start = build_spoke(c, v_start)?;

        // Sub-patch k's outer corners (0, 1, 2) sit on fillet strips:
        //   corner 0 = V_{2k}        on fillet[k-1 mod N] and fillet[k]
        //   corner 1 = V_{2k+1}      on fillet[k]
        //   corner 2 = V_{2k+2}      on fillet[k] and fillet[k+1 mod N]
        // Corner 3 = C is interior (handled by SS in pass 2).
        let strip_current = fillet_surfaces.get(k);
        let strip_prev = fillet_surfaces.get((k + n - 1) % n);
        let strip_next = fillet_surfaces.get((k + 1) % n);
        let corner_strips: [Option<&FilletSurfaceResult>; 4] =
            [strip_prev, strip_current, strip_next, None];

        if use_deg5 {
            let b0 = fit_quintic_bezier_through(b0_curve);
            let b1 = fit_quintic_bezier_through(b1_curve);
            let b2 = fit_quintic_bezier_through(&spoke_v_end_to_c);
            let b3 = fit_quintic_bezier_through(&spoke_c_to_v_start);
            let inner_cps_raw =
                compute_parallelogram_inner_cps_deg5([&b0, &b1, &b2, &b3], corner_strips);
            per_patch.push((StagedBoundaries::Deg5([b0, b1, b2, b3]), inner_cps_raw));
        } else {
            let b0 = fit_cubic_bezier_through(b0_curve);
            let b1 = fit_cubic_bezier_through(b1_curve);
            let b2 = fit_cubic_bezier_through(&spoke_v_end_to_c);
            let b3 = fit_cubic_bezier_through(&spoke_c_to_v_start);
            let inner_cps_raw =
                compute_parallelogram_inner_cps([&b0, &b1, &b2, &b3], corner_strips);
            per_patch.push((StagedBoundaries::Cubic([b0, b1, b2, b3]), inner_cps_raw));
        }
    }

    // Pass 2: Shirman-Séquin spoke-wise G1 solve at the central
    // vertex. Spoke `k` is shared between sub-patches `k` and
    // `(k+1) mod N`; its 3-D direction at `C` points toward hex
    // vertex `V_{2(k+1) mod 2N}` (sub-patch `k`'s `v_end` /
    // sub-patch `(k+1)`'s `v_start`). Returns `None` for
    // degenerate normals or when all spokes project to zero in
    // the tangent plane; caller leaves raw parallelogram values.
    let raw_p21: Vec<Point3> = per_patch.iter().map(|(_, cps)| cps[3].0).collect();
    let spoke_dirs: Vec<Vec3> = (0..n)
        .map(|k| {
            let v_shared = hex_vertices[(2 * k + 2) % two_n];
            sub3(v_shared, c)
        })
        .collect();
    let ss = shirman_sequin_vertex_enclosure_solve(c, n_c, &raw_p21, &spoke_dirs);

    let mut sub_patches: Vec<GregoryPatch> = Vec::with_capacity(n);
    for (k, (boundaries, mut inner_cps)) in per_patch.into_iter().enumerate() {
        if let Some(ref shifted) = ss {
            inner_cps[3] = (shifted[k], shifted[k]);
        }
        match boundaries {
            StagedBoundaries::Cubic([b0, b1, b2, b3]) => {
                let patch =
                    RectangularGregory::from_boundary_beziers([b0, b1, b2, b3], inner_cps, 1e-3)
                        .ok()?;
                sub_patches.push(GregoryPatch::Rectangular(Box::new(patch)));
            }
            StagedBoundaries::Deg5([b0, b1, b2, b3]) => {
                let patch = RectangularGregoryDeg5::from_boundary_beziers(
                    [b0, b1, b2, b3],
                    inner_cps,
                    1e-3,
                )
                .ok()?;
                sub_patches.push(GregoryPatch::RectangularDeg5(Box::new(patch)));
            }
        }
    }

    if sub_patches.len() != n {
        return None;
    }
    Some(NQuadDecomposition {
        central_point: c,
        sub_patches,
    })
}

/// Back-compat helper: `build_three_quad_decomposition` is just the
/// N = 3 case of [`build_n_quad_decomposition`]. Kept as a named
/// entry point so callers that only handle the trihedral case keep
/// compiling while the N ≥ 4 dispatch is wired in.
#[allow(dead_code)]
pub(crate) fn build_three_quad_decomposition(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<NQuadDecomposition> {
    let d = build_n_quad_decomposition(boundary, fillet_surfaces)?;
    if d.sub_patches.len() != 3 {
        return None;
    }
    Some(d)
}

/// Approximate outward normal of the hex plane at the central
/// vertex `C`. Sums cross products of consecutive
/// centroid-to-hex-vertex vectors — for a convex hex in CCW
/// order the sum points along the plane normal. Returns `None`
/// if the 6 cross products sum to a near-zero vector
/// (degenerate hex, e.g., all 6 vertices on a line, or a saddle
/// where the positive + negative curvature contributions cancel).
fn central_vertex_normal(hex_vertices: &[Point3], c: Point3) -> Option<Vec3> {
    if hex_vertices.len() != 6 {
        return None;
    }
    let mut nx = 0.0_f64;
    let mut ny = 0.0_f64;
    let mut nz = 0.0_f64;
    for k in 0..6 {
        let a = sub3(hex_vertices[k], c);
        let b = sub3(hex_vertices[(k + 1) % 6], c);
        let cross = a.cross(&b);
        nx += cross.x;
        ny += cross.y;
        nz += cross.z;
    }
    let n = Vec3::new(nx, ny, nz);
    n.normalized()
}

/// §5.3 item #5 — saddle central-vertex normal fallback.
///
/// When [`central_vertex_normal`] returns `None` (saddle-case
/// cancellation of the cross-product sum), pick the normal
/// direction from the fillet-surface normals at each fillet-end
/// arc midpoint. Averaging the absolute-valued normals (to avoid
/// sign cancellation across the ridge) gives a direction
/// perpendicular to the saddle ridge — the correct common
/// tangent-plane normal at C.
fn saddle_normal_fallback(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<Vec3> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    let mut n = 0usize;
    let mut seed: Option<Vec3> = None;
    for (i, seg) in boundary.segments.iter().enumerate() {
        let BoundarySegment::FilletEnd { curve, .. } = seg else {
            continue;
        };
        // Pair each FilletEnd arc (even index 2k) with the fillet
        // strip at the same index.
        let strip_idx = i / 2;
        let Some(strip) = fillet_surfaces.get(strip_idx) else {
            continue;
        };
        let (t0, t1) = curve.parameter_range();
        let mid = curve.point_at(0.5 * (t0 + t1));
        let Some((u, v)) = strip.surface.closest_point(&mid) else {
            continue;
        };
        let normal = strip.surface.normal_at(u, v);
        let Some(unit) = normal.normalized() else {
            continue;
        };
        // Flip the contribution if it's anti-aligned with the
        // first sampled normal — this ensures the average is in
        // one consistent half-space across the saddle ridge.
        let aligned = match seed {
            None => {
                seed = Some(unit);
                unit
            }
            Some(ref s) => {
                if unit.dot(s) >= 0.0 {
                    unit
                } else {
                    Vec3::new(-unit.x, -unit.y, -unit.z)
                }
            }
        };
        sum = Vec3::new(sum.x + aligned.x, sum.y + aligned.y, sum.z + aligned.z);
        n += 1;
    }
    if n == 0 {
        return None;
    }
    sum.normalized()
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
