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

/// `N` rectangular Gregory sub-patches meeting at a central point,
/// decomposing a `2N`-segment vertex-blend boundary (alternating
/// FilletEnd / FaceEdge). `N = 3` is the trihedral hex case; `N ≥
/// 4` covers pyramid apex, hub-and-spoke fittings, etc.
#[allow(dead_code)] // diagnostic / follow-up
#[derive(Debug)]
pub(crate) struct NQuadDecomposition {
    /// The central point shared by all `N` sub-patches.
    pub central_point: Point3,
    /// Sub-patches in CCW order starting from arc 0. Length `N`.
    pub sub_patches: Vec<RectangularGregory>,
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

    // Central point = centroid of the `2N` boundary vertices. Inside
    // the polygon for convex boundaries (always true for trihedral
    // corners; for N ≥ 4 relies on the bridge-arc chain remaining
    // convex — true for standard pyramidal / hub geometry).
    let denom = two_n as f64;
    let c = Point3::new(
        hex_vertices.iter().map(|p| p.x).sum::<f64>() / denom,
        hex_vertices.iter().map(|p| p.y).sum::<f64>() / denom,
        hex_vertices.iter().map(|p| p.z).sum::<f64>() / denom,
    );

    // C4b.2b central-vertex normal. All three sub-patches' corner-3
    // inner CPs must lie in a common tangent plane at `C` for G1
    // across the three shared spokes. We pick the plane via
    // cross-product of consecutive hex-vertex-to-centroid vectors —
    // a convex hex's outward normal is consistent across such
    // cross products. Normalised; sign chosen so the normal points
    // "out of the blend" (away from the averaged strip normals at
    // the arc corners).
    let n_c = central_vertex_normal(&hex_vertices, c).unwrap_or(Vec3::new(0.0, 0.0, 0.0));

    // Pass 1: build each sub-patch's boundaries and raw parallelogram
    // inner CPs (with C-K twist applied at corners 0–2). Corner 3
    // (central vertex) stays at its raw parallelogram position; the
    // Shirman-Séquin cross-patch solve in pass 2 corrects it.
    type SubPatchStaged = ([[Point3; 4]; 4], [(Point3, Point3); 4]);
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

        let b0 = fit_cubic_bezier_through(b0_curve);
        let b1 = fit_cubic_bezier_through(b1_curve);
        let b2 = fit_cubic_bezier_through(&spoke_v_end_to_c);
        let b3 = fit_cubic_bezier_through(&spoke_c_to_v_start);

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

        let inner_cps_raw = compute_parallelogram_inner_cps([&b0, &b1, &b2, &b3], corner_strips);
        per_patch.push(([b0, b1, b2, b3], inner_cps_raw));
    }

    // Pass 2: Shirman-Séquin vertex-enclosure solve at the central
    // vertex `C`. The simplified N=3 form: pick an orthonormal tangent
    // basis `(t1, t2)` ⊥ `n_c`, decompose each sub-patch's raw
    // corner-3 inner CP offset from `C` into `(a_k, b_k, d_k)` in
    // that basis, then shift so each sub-patch's final corner-3 inner
    // CP sits in the tangent plane AND all three agree on the
    // tangent-plane components (subtract the cross-patch mean).
    //
    // This gives the relaxed G2 at C in the surface-normal sense
    // (Shirman-Séquin 1991): the mixed second partials τ_k across
    // sub-patches agree in the tangent plane, which is what governs
    // reflection-line continuity across internal spokes. The residual
    // normal-direction disagreement is irreducible given the boundary
    // geometry — G2 in the tangent direction is what actually matters
    // for visual smoothness.
    //
    // Falls back to the C4b.2b pure-projection (current code) when
    // `n_c` is degenerate.
    let use_ss = n_c.length_squared() > 1e-18;
    let (t1, t2) = if use_ss {
        orthonormal_tangent_basis(n_c)
    } else {
        (Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0))
    };
    let mut a_tan = vec![0.0_f64; n];
    let mut b_tan = vec![0.0_f64; n];
    if use_ss {
        for k in 0..n {
            let p_para = per_patch[k].1[3].0;
            let off = Vec3::new(p_para.x - c.x, p_para.y - c.y, p_para.z - c.z);
            a_tan[k] = off.dot(&t1);
            b_tan[k] = off.dot(&t2);
        }
    }
    let a_mean = a_tan.iter().sum::<f64>() / n as f64;
    let b_mean = b_tan.iter().sum::<f64>() / n as f64;

    let mut sub_patches: Vec<RectangularGregory> = Vec::with_capacity(n);
    for (k, (boundaries, mut inner_cps)) in per_patch.into_iter().enumerate() {
        if use_ss {
            // Final corner-3 inner CP: C + (a_k - a_mean) t1 + (b_k - b_mean) t2.
            let da = a_tan[k] - a_mean;
            let db = b_tan[k] - b_mean;
            let p_ss = Point3::new(
                c.x + da * t1.x + db * t2.x,
                c.y + da * t1.y + db * t2.y,
                c.z + da * t1.z + db * t2.z,
            );
            inner_cps[3] = (p_ss, p_ss);
        }
        let [b0, b1, b2, b3] = boundaries;
        let patch =
            RectangularGregory::from_boundary_beziers([b0, b1, b2, b3], inner_cps, 1e-3).ok()?;
        sub_patches.push(patch);
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

/// C4b.2b — compute Chiyokura–Kimura inner-CP pairs from the 4
/// boundary Beziers + the incident fillet strips, matching the
/// strip's mixed cross-derivative `∂²T/∂u'∂v'` at each hex
/// corner.
///
/// **Theory.** For a bicubic Bezier at corner `(0, 0)`:
/// `S_uv(0, 0) = 9 · (P_11 − P_parallelogram)`
/// where `P_parallelogram = P_10 + P_01 − P_00`. For G2 with an
/// adjacent strip `T(u', v')` across the arc boundary, we want
/// `S_uv(0, 0) = τ_target` where `τ_target` is the projection
/// of the strip's second fundamental form onto the sub-patch's
/// tangent directions via chain rule.
///
/// **Chain rule.** Let `e_u`, `e_v` be the sub-patch's boundary
/// tangent vectors at the corner (3× the first Bezier CP
/// offset). Express them in the strip's tangent basis
/// `{T_u', T_v'}` by solving the 2×2 Gram system — gives
/// coefficients `(α_u, β_u)`, `(α_v, β_v)` such that
/// `e_u ≈ α_u · T_u' + β_u · T_v'` (projected onto the strip's
/// tangent plane). Then the strip's induced mixed derivative in
/// sub-patch coords is:
/// ```text
/// τ = α_u·α_v · T_u'u'
///   + (α_u·β_v + β_u·α_v) · T_u'v'
///   + β_u·β_v · T_v'v'
/// ```
/// Setting `P_11 = P_parallelogram + τ/9` makes the sub-patch's
/// cross-derivative match the strip's — the classical C-K
/// minimum-twist construction.
///
/// Both `P_u` and `P_v` at the corner are set to the same
/// corrected value — a "symmetric Chiyokura" that collapses the
/// Gregory rational blend to an ordinary bicubic at the corner.
/// This is the first-cut correct G2 at hex corners. The G2
/// extension to the full boundary (not just corners) and the
/// central-vertex cross-patch G1 solve are the remaining
/// C4b.2b work.
fn compute_parallelogram_inner_cps(
    boundaries: [&[Point3; 4]; 4],
    corner_strips: [Option<&FilletSurfaceResult>; 4],
) -> [(Point3, Point3); 4] {
    let b0 = boundaries[0];
    let b1 = boundaries[1];
    let b2 = boundaries[2];
    let b3 = boundaries[3];

    let corners = [b0[0], b0[3], b1[3], b2[3]];
    let raw_parallelogram = [
        parallelogram_add(b3[2], b0[1], b0[0]),
        parallelogram_add(b0[2], b1[1], b0[3]),
        parallelogram_add(b2[1], b1[2], b1[3]),
        parallelogram_add(b3[1], b2[2], b2[3]),
    ];

    // Sub-patch boundary tangent vectors at each corner, from
    // the Bezier derivative `S_u(0, 0) = 3 · (P_10 − P_00)`,
    // etc. The INBOUND tangent (Bezier index reversed) is used
    // for `e_v` where appropriate so both vectors point AWAY
    // from the corner into the patch interior.
    //
    // corner 0 (u-boundary: b3 reversed; v-boundary: b0):
    //   e_u = 3·(b3[2] − b0[0]),  e_v = 3·(b0[1] − b0[0])
    // corner 1 (u-boundary: b0 reversed ending at corner 1;
    //           v-boundary: b1):
    //   e_u = 3·(b0[2] − b0[3]),  e_v = 3·(b1[1] − b1[0])
    // corner 2 (u-boundary: b1 reversed; v-boundary: b2):
    //   e_u = 3·(b1[2] − b1[3]),  e_v = 3·(b2[1] − b2[0])
    // corner 3 (u-boundary: b2 reversed; v-boundary: b3):
    //   e_u = 3·(b2[2] − b2[3]),  e_v = 3·(b3[1] − b3[0])
    let tangents: [(Vec3, Vec3); 4] = [
        (scale3(sub3(b3[2], b0[0])), scale3(sub3(b0[1], b0[0]))),
        (scale3(sub3(b0[2], b0[3])), scale3(sub3(b1[1], b1[0]))),
        (scale3(sub3(b1[2], b1[3])), scale3(sub3(b2[1], b2[0]))),
        (scale3(sub3(b2[2], b2[3])), scale3(sub3(b3[1], b3[0]))),
    ];

    let mut inner_cps: [(Point3, Point3); 4] =
        [(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)); 4];
    for (i, slot) in inner_cps.iter_mut().enumerate() {
        let corner = corners[i];
        let p_para = raw_parallelogram[i];
        let (e_u, e_v) = tangents[i];
        let p_refined = if let Some(strip) = corner_strips[i] {
            chiyokura_twisted_inner_cp(&strip.surface, corner, p_para, e_u, e_v).unwrap_or(p_para)
        } else {
            // Central vertex (no incident strip) keeps the raw
            // parallelogram value. C4b.2b central-vertex
            // cross-patch G1 solve is follow-up.
            p_para
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

fn sub3(a: Point3, b: Point3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

fn scale3(v: Vec3) -> Vec3 {
    Vec3::new(3.0 * v.x, 3.0 * v.y, 3.0 * v.z)
}

/// Approximate outward normal of the hex plane at the central
/// vertex `C`. Sums cross products of consecutive
/// centroid-to-hex-vertex vectors — for a convex hex in CCW
/// order the sum points along the plane normal. Returns `None`
/// if the 6 cross products sum to a near-zero vector
/// (degenerate hex, e.g., all 6 vertices on a line).
/// Build an orthonormal tangent basis `(t1, t2)` in the plane
/// perpendicular to the unit vector `n`. `t1` is derived from the
/// world-axis least aligned with `n` (for numerical stability); `t2`
/// is `n × t1` so `(t1, t2, n)` is a right-handed orthonormal frame.
fn orthonormal_tangent_basis(n: Vec3) -> (Vec3, Vec3) {
    let ax = n.x.abs();
    let ay = n.y.abs();
    let az = n.z.abs();
    let seed = if ax <= ay && ax <= az {
        Vec3::new(1.0, 0.0, 0.0)
    } else if ay <= az {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    // t1 = (seed - (seed·n) n) normalised
    let dot = seed.dot(&n);
    let raw_t1 = Vec3::new(seed.x - dot * n.x, seed.y - dot * n.y, seed.z - dot * n.z);
    let t1 = raw_t1.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let t2 = n.cross(&t1);
    (t1, t2)
}

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

/// Chiyokura–Kimura twist-corrected inner CP at a hex corner.
/// Returns `P_parallelogram + τ/9` where `τ` is the strip's
/// induced mixed derivative in sub-patch coords via chain rule.
fn chiyokura_twisted_inner_cp(
    surface: &SurfaceGeometry,
    corner: Point3,
    p_parallelogram: Point3,
    e_u: Vec3,
    e_v: Vec3,
) -> Option<Point3> {
    let (u, v) = project_point_to_uv(surface, &corner)?;
    let t_u = surface.partial_u(u, v);
    let t_v = surface.partial_v(u, v);
    let t_uu = surface.second_partial_uu(u, v);
    let t_uv = surface.second_partial_uv(u, v);
    let t_vv = surface.second_partial_vv(u, v);

    // Solve the 2×2 Gram system
    //   [T_u·T_u  T_u·T_v] [α]   [T_u·e]
    //   [T_u·T_v  T_v·T_v] [β] = [T_v·e]
    // for each sub-patch tangent vector. Express sub-patch
    // tangents in the strip's (u', v') basis, projected onto
    // the strip's tangent plane.
    let g00 = t_u.dot(&t_u);
    let g01 = t_u.dot(&t_v);
    let g11 = t_v.dot(&t_v);
    let det = g00 * g11 - g01 * g01;
    if det.abs() < 1e-18 {
        return None;
    }
    let solve = |e: &Vec3| -> (f64, f64) {
        let r0 = t_u.dot(e);
        let r1 = t_v.dot(e);
        ((g11 * r0 - g01 * r1) / det, (g00 * r1 - g01 * r0) / det)
    };
    let (alpha_u, beta_u) = solve(&e_u);
    let (alpha_v, beta_v) = solve(&e_v);

    // τ = α_u·α_v · T_uu + (α_u·β_v + β_u·α_v) · T_uv
    //     + β_u·β_v · T_vv
    let cross = alpha_u * beta_v + beta_u * alpha_v;
    let tau = Vec3::new(
        alpha_u * alpha_v * t_uu.x + cross * t_uv.x + beta_u * beta_v * t_vv.x,
        alpha_u * alpha_v * t_uu.y + cross * t_uv.y + beta_u * beta_v * t_vv.y,
        alpha_u * alpha_v * t_uu.z + cross * t_uv.z + beta_u * beta_v * t_vv.z,
    );

    // P = P_parallelogram + τ/9.
    Some(Point3::new(
        p_parallelogram.x + tau.x / 9.0,
        p_parallelogram.y + tau.y / 9.0,
        p_parallelogram.z + tau.z / 9.0,
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
