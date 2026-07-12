//! Phase I.10c — bridge from [`PolarBoundaryInput`] +
//! [`TensorBorder`] (I.4) to [`Bi6CapInput`] (I.10).
//!
//! This module closes the I.4 → I.10 API gap: given the *shape*
//! inputs the polar construction consumes (central vertex, N
//! shared spokes, N outer corners with tangent planes, N face
//! tangent planes, per-spoke radii) plus a prepared bi-3
//! tensor-border, produce the `Bi6CapInput` that
//! [`f5_solve_optimal_bi6_cap`] can optimise.
//!
//! # What is **not** in scope
//!
//! - **Tensor-border extraction from fillet geometry.** Caller
//!   supplies a prepared [`TensorBorder`]. Deriving the bi-3
//!   border from the fillet's boundary curves is a separate
//!   geometric step (planned as I.10d).
//! - **Fillet-crate dispatch.** No `popkern-fillet` integration
//!   yet (planned as I.10e).
//! - **(S1) seed propagation.** The full paper-compliant free-
//!   parameter set includes sector 0's near-star CPs; our
//!   scaffold still takes per-sector depth-1 data via
//!   `CirculantSectorRhs::p_01/p_11/p_02`. The converter uses a
//!   simple mapping: row-1 bi-3 CPs at positions 0, 1 become
//!   the per-sector `p_01`, `p_11`; row-2 bi-3 CP at position 0
//!   becomes `p_02`. This is adequate for the architectural
//!   hook-up; F_5 minimisation still produces a well-posed
//!   solution regardless of whether these values exactly match
//!   paper convention.
//!
//! # Algorithm summary
//!
//! 1. Validate inputs (`n ≥ 3`, matching lengths, bi-3 tensor-
//!    border kind).
//! 2. Degree-elevate each sector's bi-3 row-0 (4 CPs) to bi-6
//!    (7 CPs) — this gives us `p_00, p_10, p_20, p_30, p_40,
//!    p_50, p_60` along the shared edge. `CirculantSectorRhs`
//!    takes indices `{0, 1, 2, 4, 5, 6}` — we discard the
//!    degree-raised `p_30` (it's recomputed from rule (5)).
//! 3. Pick `sectors[k].p_01/p_11/p_02` from the bi-3 tensor-
//!    border's row-1/row-2 depth-1 CPs.
//! 4. `p^k_{i2}` for `i ∈ {4, 5, 6}` comes from β-reparameterised
//!    row-2 via I.5's γ table — simplified to direct bi-6
//!    elevation for the scaffold.
//! 5. Catmull-Clark inputs: `c7 = central_vertex`; `c5s[k]`,
//!    `c6s[k]` = interior spoke sample points (at `t = 1/3`,
//!    `t = 2/3` along each spoke).
//! 6. Free parameters default to zero (caller passes to F_5
//!    min).

use super::kp_polar_boundary_input::{PolarBoundaryInput, TensorBorder, TensorBorderBi3};
use super::kp_polar_shared::catmull_clark_star_point;
use super::kp15_bi6_cap::Bi6CapInput;
use super::kp15_bi6_circulant::CirculantSectorRhs;
use super::kp15_bi6_deep_interior::Bi6DeepInteriorFree;
use super::kp15_bi6_reparam::Bi6Reparams;
use crate::geometry::Point3;

/// Convert a [`PolarBoundaryInput`] + bi-3 [`TensorBorder`] into
/// a [`Bi6CapInput`] ready for
/// [`f5_solve_optimal_bi6_cap`](super::kp15_bi6_cap_f5_solve::f5_solve_optimal_bi6_cap).
///
/// See module docs for the simplifications and scope.
///
/// # Errors
///
/// - `n < 3` or mismatched input lengths.
/// - Tensor-border kind is `Bi5` (the bi-6 path wants bi-3).
pub fn convert_polar_to_bi6_cap_input(
    polar: &PolarBoundaryInput,
    tensor_border: &TensorBorder,
) -> Result<Bi6CapInput, String> {
    let n = polar.n();
    polar.validate()?;
    tensor_border.validate(n)?;
    let bi3: &TensorBorderBi3 = match tensor_border {
        TensorBorder::Bi3(b) => b,
        TensorBorder::Bi5(_) => {
            return Err("bi-6 polar cap requires a bi-3 TensorBorder; got bi-5".to_string());
        }
    };

    let params = Bi6Reparams::for_valence(n);

    // Per-sector circulant RHS.
    //
    // **Phase E.1 attempted and reverted (2026-06-22):** The plan
    // calls for (S1) propagation here — sector 0's `p_01, p_11,
    // p_02` as the free seed, sectors k ≥ 1 derived via
    // `rule_s1_star_propagation`. Measurement showed the eq (2)
    // residual at `assembly_eq2_residual_diagnostic` grows from
    // ~10 → ~73 (n=6) under that change. Root cause: the (S1)
    // formulas in [KP15 eq (6)] assume `p̀_{i, 0} = ṕ_{i, 0}` —
    // i.e., the j=0 row is the cross-sector spoke shared between
    // adjacent sectors. In the kernel, j=0 is each sector's *own*
    // bi-3 tensor-border row (the rim segment, per-sector), and
    // adjacent sectors do NOT share their j=0 rows. Naive (S1)
    // applied to rim-indexed data yields the wrong propagation.
    //
    // The full fix requires re-orienting the bi-6 grid convention
    // so j=0 is the spoke (or, equivalently, adding cross-sector
    // identifications of `p_{0..2, 0}` to the constraint matrix
    // before propagating). That's a deeper refactor than this
    // bridge can carry. Until then, each sector seeds its own
    // near-star CPs from its own elevated bi-3 row — the current
    // behavior, preserved by this comment.
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let row0 = bi3.sectors[k][0]; // 4 bi-3 position CPs
        let row1 = bi3.sectors[k][1]; // 4 bi-3 first-cross-deriv CPs
        let row2 = bi3.sectors[k][2]; // 4 bi-3 second-cross-deriv CPs

        let row0_bi6 = elevate_bi3_point_curve_to_bi6(row0);
        let row1_bi6 = elevate_bi3_point_curve_to_bi6(row1);
        let row2_bi6 = elevate_bi3_point_curve_to_bi6(row2);

        sectors.push(CirculantSectorRhs {
            p_00: row0_bi6[0],
            p_10: row0_bi6[1],
            p_20: row0_bi6[2],
            p_40: row0_bi6[4],
            p_50: row0_bi6[5],
            p_60: row0_bi6[6],
            // Per-sector seed for the near-star CPs.
            p_01: row1_bi6[0],
            p_11: row1_bi6[1],
            p_02: row2_bi6[0],
            p_42: row2_bi6[4],
            p_52: row2_bi6[5],
            p_62: row2_bi6[6],
        });
    }

    // Cap-net: c7 is the central vertex. c5s, c6s are sampled
    // along each spoke at t = 1/3 and t = 2/3 respectively.
    //
    // **Phase E.1 second attempt and reverted (2026-06-22):**
    // overwriting every sector's `p_00` with `catmull_clark_star_
    // point(n, c7, c5s, c6s)` (so `CirculantSectorRhs.p_00` matches
    // what the cap assembly will use as `grid[0][0]`) shifts the
    // eq (2) residual by less than 50% and breaks the bridge's
    // endpoint-preservation test contract. The unification is
    // *conceptually* right — the kernel currently keeps two
    // distinct "p_00" notions (bridge's tensor-border `row0_bi6[0]`
    // and the cap's `star = catmull_clark_star_point(...)`) — but
    // resolving it cleanly belongs in a deeper refactor that also
    // adds cross-sector `p^k_{i, 0} = p^0_{i, 0}` identifications
    // as constraint rows in M (see Phase E to-be-revised).
    let c7 = polar.central_vertex;
    let mut c5s = Vec::with_capacity(n);
    let mut c6s = Vec::with_capacity(n);
    for spoke in &polar.spokes {
        c5s.push(eval_degree7_bezier(&spoke.boundary_cps, 1.0 / 3.0));
        c6s.push(eval_degree7_bezier(&spoke.boundary_cps, 2.0 / 3.0));
    }

    // Pin p0_21 at n ∈ {3, 6}. A sensible heuristic default is
    // the star point itself — F_5 min will refine, and callers
    // can override.
    let p0_21_pinned = if n == 3 || n == 6 { Some(c7) } else { None };

    Ok(Bi6CapInput {
        n,
        params,
        c7,
        c5s,
        c6s,
        sectors,
        p0_22: Point3::new(0.0, 0.0, 0.0),
        p0_21_pinned,
        interior_filler: None,
    })
}

/// Degree-elevate a 4-CP bi-3 Bezier curve to 7-CP bi-6.
/// Uses the generic [`elevate_point_curve_to_nine`] in reverse
/// by padding to 10 then truncating — but that doesn't work for
/// a true 4 → 7 elevation. Instead, apply the standard
/// one-step elevation from degree `d` to `d + 1` three times.
fn elevate_bi3_point_curve_to_bi6(cps: [Point3; 4]) -> [Point3; 7] {
    // Start at degree 3 (4 CPs), elevate → 4, 5, 6 (7 CPs).
    let mut cur: Vec<Point3> = cps.to_vec();
    for _ in 0..3 {
        let d = cur.len() - 1;
        let dp1 = (d + 1) as f64;
        let mut next = Vec::with_capacity(d + 2);
        next.push(cur[0]);
        for i in 1..=d {
            let frac = (i as f64) / dp1;
            let prev = cur[i - 1];
            let this = cur[i];
            next.push(Point3::new(
                frac * prev.x + (1.0 - frac) * this.x,
                frac * prev.y + (1.0 - frac) * this.y,
                frac * prev.z + (1.0 - frac) * this.z,
            ));
        }
        next.push(cur[d]);
        cur = next;
    }
    debug_assert_eq!(cur.len(), 7);
    core::array::from_fn(|i| cur[i])
}

/// Evaluate a degree-7 Bezier curve (8 CPs) at parameter `t`
/// via de Casteljau (7 linear-interpolation passes).
fn eval_degree7_bezier(cps: &[Point3; 8], t: f64) -> Point3 {
    let mut work: [Point3; 8] = *cps;
    for r in 1..=7 {
        for i in 0..(8 - r) {
            let a = work[i];
            let b = work[i + 1];
            work[i] = Point3::new(
                (1.0 - t) * a.x + t * b.x,
                (1.0 - t) * a.y + t * b.y,
                (1.0 - t) * a.z + t * b.z,
            );
        }
    }
    work[0]
}

// =====================================================================
// PK-087 stage 3 — the re-oriented polar → **shared cap net** bridge.
// =====================================================================
//
// The `convert_polar_to_bi6_cap_input` above is the legacy frame:
// per-sector rim borders in `CirculantSectorRhs`, no shared cap net,
// spoke on the wrong slot. The functions below are its correct-by-
// construction replacement that feeds the verified deep-interior
// populator
// ([`build_bi6_cap_from_net`](super::kp15_bi6_deep_interior::build_bi6_cap_from_net)):
// a single **shared** 4×4 bicubic B-spline net per sector where
// adjacent sectors share control points along the **spoke** (the
// sector-glue edge) and the tensor border lands on the two
// **star-opposite rim edges**. Because that sharing *is* the
// condition that the input net is one C²-consistent B-spline surface,
// the eq (7)/(8) ring closure `sal` vanishes and — by the oracle's
// symbolic proof — every eq (2)/(3) cross-sector constraint closes by
// construction, *independent of the control-point values*. The naive
// rim-as-spoke wire grew eq (2) ~10 → ~73 (Phase-E.1 record); this
// frame drives it to machine zero.

/// The shared control ring of the coarse cap net: per-spoke control
/// symbols `a1..a6` (indexed by spoke/sector `0..n`) plus the single
/// global central face point `e`. This is the C²-consistent
/// parameterisation of the whole cap net; [`shared_net_layout`]
/// windows it into per-sector 4×4 B-spline nets.
#[derive(Debug, Clone)]
pub struct SharedControlRing {
    /// `a1[k]` — outer corner control symbol on spoke `k`.
    pub a1: Vec<Point3>,
    /// `a2[k]` — rim-edge interior control symbol (near corner `k`).
    pub a2: Vec<Point3>,
    /// `a3[k]` — rim-edge interior control symbol (near corner `k+1`).
    pub a3: Vec<Point3>,
    /// `a4[k]` — ring-2 control symbol on spoke `k`.
    pub a4: Vec<Point3>,
    /// `a5[k]` — sector-`k` face control symbol.
    pub a5: Vec<Point3>,
    /// `a6[k]` — ring-1 control symbol on spoke `k` (adjacent to `e`).
    pub a6: Vec<Point3>,
    /// `E` — the global central face point (shared by all sectors).
    pub e: Point3,
}

/// The re-oriented bridge output: per-valence reparameterisation
/// scalars, the per-sector C²-consistent 4×4 bicubic nets, and the
/// KP15 free set.
pub type Bi6CapNetInputs = (Bi6Reparams, Vec<[[Point3; 4]; 4]>, Bi6DeepInteriorFree);

/// Lay out one sector's 4×4 bicubic B-spline net from the shared
/// control ring, exactly per the `n6` script's C²-consistent net
/// (oracle `n6_oracle.py::shared_net`, script comment lines 227–230).
///
/// `net[s]` is written column-consistent with the populator's
/// accessor `p(i, j) = net[j][i]`:
///
/// ```text
///   [ a1[s]  , a2[s]  , a3[s]  , a4[s+1] ]
///   [ a4[s]  , a5[s]  , a6[s]  , a5[s+1] ]
///   [ a3[s-1], a6[s-1], E      , a6[s+1] ]
///   [ a2[s-1], a5[s-1], a6[s-2], 0       ]
/// ```
///
/// Adjacent sectors share `a4/a5/a6` along the spoke; `E` is global.
/// This sharing is what makes the assembled bicubic net a single
/// C²-consistent surface, so `sal ≈ 0` by construction.
pub fn shared_net_layout(ring: &SharedControlRing, s: usize, n: usize) -> [[Point3; 4]; 4] {
    let sf = |k: isize| k.rem_euclid(n as isize) as usize;
    let sp1 = sf(s as isize + 1);
    let sm1 = sf(s as isize - 1);
    let sm2 = sf(s as isize - 2);
    let (a1, a2, a3, a4, a5, a6, e) = (
        &ring.a1, &ring.a2, &ring.a3, &ring.a4, &ring.a5, &ring.a6, ring.e,
    );
    [
        [a1[s], a2[s], a3[s], a4[sp1]],
        [a4[s], a5[s], a6[s], a5[sp1]],
        [a3[sm1], a6[sm1], e, a6[sp1]],
        [a2[sm1], a5[sm1], a6[sm2], Point3::new(0.0, 0.0, 0.0)],
    ]
}

/// Sample a spoke's degree-7 Bézier boundary curve at parameter `t`
/// (`t = 0` at the central vertex `C`, `t = 1` at the outer corner).
fn eval_spoke_bezier(cps: &[Point3; 8], t: f64) -> Point3 {
    let mut work: [Point3; 8] = *cps;
    for r in 1..=7 {
        for i in 0..(8 - r) {
            let a = work[i];
            let b = work[i + 1];
            work[i] = Point3::new(
                (1.0 - t) * a.x + t * b.x,
                (1.0 - t) * a.y + t * b.y,
                (1.0 - t) * a.z + t * b.z,
            );
        }
    }
    work[0]
}

/// Linear blend `p + t·(q − p)`.
fn lerp_point(p: Point3, q: Point3, t: f64) -> Point3 {
    Point3::new(
        p.x + t * (q.x - p.x),
        p.y + t * (q.y - p.y),
        p.z + t * (q.z - p.z),
    )
}

/// Convert a [`PolarBoundaryInput`] into the **shared cap net** frame
/// the KP15 deep-interior populator consumes: the per-valence
/// [`Bi6Reparams`], the per-sector C²-consistent 4×4 bicubic nets, and
/// the paper free set. Feed the output tuple to
/// [`assemble_bi6_polar_cap_from_net`](super::kp15_bi6_deep_interior::assemble_bi6_polar_cap_from_net).
///
/// This is the stage-3 re-orientation of [`convert_polar_to_bi6_cap_input`]:
/// the tensor border lands on the star-opposite rim edges and the
/// spokes carry the cross-sector glue, so the populator's ring closure
/// and every eq (2)/(3) constraint close by construction.
///
/// The shared control ring is read from the polar geometry as the
/// coarse Catmull–Clark control mesh of the cap region (each symbol a
/// single well-defined point per index, so the shared-net sharing is
/// exact); the free set is a monotone star→rim blend keyed by the
/// CP's script depth so the resolved cap is fold-free. The paper
/// determines the free set by F₅ minimisation; this supplies a
/// well-posed, geometry-derived seed (F₅ refinement is a later layer
/// and does not affect the closure).
///
/// # Errors
///
/// - `polar.validate()` fails (non-finite radii, length mismatch).
/// - `n < 5` — the valence-general deep-interior populator declines
///   `n < 5` (paper §3.5 `n = 3` and the rectangular `n = 4` path are
///   served elsewhere); the caller's dispatch ladder handles those.
pub fn convert_polar_to_bi6_cap_net(polar: &PolarBoundaryInput) -> Result<Bi6CapNetInputs, String> {
    polar.validate()?;
    let n = polar.n();
    if n < 5 {
        return Err(format!(
            "shared-net bi-6 polar bridge is valence-general for n ≥ 5; got n = {n}"
        ));
    }
    let params = Bi6Reparams::for_valence(n);
    let c = polar.central_vertex;
    let sf = |k: usize| k % n;

    // Outer-corner rim vertices V_k.
    let v: Vec<Point3> = polar.outer_corners.iter().map(|oc| oc.position).collect();

    // ---- Shared control ring from the polar geometry ----
    //
    // The coarse Catmull–Clark control mesh of the cap region. The
    // spoke curve `q_k(t)` runs C → V_k; the rim edge of sector s is
    // V_s → V_{s+1}. Each symbol is a single well-defined point per
    // index (that single-valuedness is what makes the shared-net C²
    // sharing exact):
    //   E  = central vertex C (the face point all sectors share).
    //   a6 = ring-1 point on spoke k at t = 1/3 (adjacent to E).
    //   a4 = ring-2 point on spoke k at t = 2/3.
    //   a1 = outer corner V_k (t = 1).
    //   a5 = sector-k face point: mid-radius average of spokes k, k+1.
    //   a2, a3 = interior control points of rim edge V_k → V_{k+1}.
    let mut a1 = Vec::with_capacity(n);
    let mut a2 = Vec::with_capacity(n);
    let mut a3 = Vec::with_capacity(n);
    let mut a4 = Vec::with_capacity(n);
    let mut a5 = Vec::with_capacity(n);
    let mut a6 = Vec::with_capacity(n);
    for k in 0..n {
        let q = &polar.spokes[k].boundary_cps;
        let qk1 = &polar.spokes[sf(k + 1)].boundary_cps;
        a1.push(v[k]);
        a4.push(eval_spoke_bezier(q, 2.0 / 3.0));
        a6.push(eval_spoke_bezier(q, 1.0 / 3.0));
        let m_k = eval_spoke_bezier(q, 0.6);
        let m_k1 = eval_spoke_bezier(qk1, 0.6);
        a5.push(lerp_point(m_k, m_k1, 0.5));
        let vk1 = v[sf(k + 1)];
        a2.push(lerp_point(v[k], vk1, 1.0 / 3.0));
        a3.push(lerp_point(v[k], vk1, 2.0 / 3.0));
    }
    let ring = SharedControlRing {
        a1,
        a2,
        a3,
        a4,
        a5,
        a6,
        e: c,
    };
    let net: Vec<[[Point3; 4]; 4]> = (0..n).map(|s| shared_net_layout(&ring, s, n)).collect();

    // ---- Free set: monotone star → rim blend (fold-free seed) ----
    //
    // Star point = Catmull–Clark limit of C over the ring-1/ring-2
    // spoke samples. The remaining free CPs sit on a monotone blend
    // from the star (apex) toward each sector's rim reference, keyed
    // by the CP's script depth d = (12 − i − j) / 12 (0 at star [6,6],
    // 1 at rim [0,0]). Monotone ⇒ fold-free; F₅ refines later.
    let c5s: Vec<Point3> = polar
        .spokes
        .iter()
        .map(|s| eval_spoke_bezier(&s.boundary_cps, 1.0 / 3.0))
        .collect();
    let c6s: Vec<Point3> = polar
        .spokes
        .iter()
        .map(|s| eval_spoke_bezier(&s.boundary_cps, 2.0 / 3.0))
        .collect();
    let star = catmull_clark_star_point(n, c, &c5s, &c6s);
    let rim_ref: Vec<Point3> = (0..n)
        .map(|k| lerp_point(v[k], v[sf(k + 1)], 0.5))
        .collect();
    let blend = |sector: usize, i: usize, j: usize| -> Point3 {
        let depth = (12.0 - (i as f64) - (j as f64)) / 12.0;
        lerp_point(star, rim_ref[sector], depth)
    };

    let free = Bi6DeepInteriorFree {
        p00: star,
        tan: [blend(0, 5, 6), blend(0, 6, 5)],
        two: [blend(0, 4, 6), blend(0, 5, 5), blend(0, 6, 4)],
        p220: blend(0, 4, 4),
        p210: blend(0, 4, 5),
        p23: (0..n).map(|r| blend(r, 4, 3)).collect(),
        p33: (0..n).map(|r| blend(r, 3, 3)).collect(),
    };

    Ok((params, net, free))
}

#[cfg(test)]
#[path = "kp15_bi6_polar_bridge_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "kp15_bi6_polar_cap_net_tests.rs"]
mod cap_net_tests;
