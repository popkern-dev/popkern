//! PnS bi-4 cap on the SETBACK border — the G1 sibling of the KP20
//! setback candidate (PK-097 the anchor-split stage).
//!
//! The existing PnS dispatch consumes the mitred blend boundary
//! (apex-end anchors); this entry derives the border from the fillet
//! surfaces' setback arcs — the geometry the cap must actually bound
//! — via the same extractor the KP20 path uses.
//!
//! Design (docs/architecture/g1-fillet-family-setback-cap.md):
//! - Sector k ↔ setback arc k (corner-to-corner span on sweep k).
//! - Position row: degree-2 interpolation of the arc cubic (the PnS
//!   paper wants quadratic borders; exactness is recovered below).
//! - Derivative row: transversal directions from the neighbor
//!   tangent planes — the sweep at the arc midpoint, the lateral
//!   plane at the corners (derived from the two adjacent arc
//!   end-tangents, which span it).
//! - After the paper-faithful assembly, the bi-4 boundary rows are
//!   reconciled onto the neighbor union (sweeps ∪ corner planes) by
//!   the PK-096 Bernstein-interpolation solve at degree 4, with the
//!   derivative row translated in parallel. The paper gates measure
//!   the construction against its own border; the shell gates
//!   measure it against the neighbors. Both are reported.

use popkern_core::geometry::gregory::{
    BetaReparam, Bi2BorderStrip, Bi2TensorBorderRing, CapAssembleInput, assemble_pns_bi4_cap,
    pole_from_boundary, sigma_center_jet_default,
};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, PortableMath, Vec3};

use super::kp20_paper_border_extract::kp20_paper_border_from_fillet_setback;
use crate::fillet_surface::FilletSurfaceResult;

/// Fraction of arc length used as the transversal derivative
/// magnitude (matches the mitred-boundary adapter's convention).
const DERIV_MAG_FRAC: f64 = 0.5;

/// F_κ order for the assembly (same default as the promote path).
const DEFAULT_KAPPA: usize = 1;

/// Measured quality of the candidate.
#[derive(Debug)]
pub struct PnsSetbackReport {
    /// Max distance from each patch's outer border curve to its
    /// setback arc BEFORE reconciliation (the degree-reduction gap).
    pub g0_border_before: f64,
    /// Same, AFTER bi-4 reconciliation onto the neighbor union.
    pub g0_border_after: f64,
    /// Max angle (radians) between the cap normal and the neighbor
    /// surface normal along the border, after reconciliation.
    pub g1_angle_after: f64,
    /// Max positional mismatch across adjacent-sector radial seams.
    pub seam_g0: f64,
    /// Valence.
    pub n: usize,
}

/// The assembled candidate.
pub struct PnsSetbackCandidate {
    /// One bi-4 patch per sector, boundary-reconciled.
    pub patches: Vec<[[Point3; 5]; 5]>,
    /// The same patches as NURBS surfaces (exact conversion).
    pub surfaces: Vec<SurfaceGeometry>,
    /// Per-sector half-arc cubics (border spans) with their fillet
    /// index, in ring order — sector k's border is half_cubics[k].
    pub half_cubics: Vec<([Point3; 4], usize)>,
    /// The cap's central pole point.
    pub pole: Point3,
    /// Quality report.
    pub report: PnsSetbackReport,
}

/// Build the PnS bi-4 cap from the fillet surfaces' setback border.
pub fn try_pns_setback_candidate(
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
) -> Result<PnsSetbackCandidate, String> {
    let (_border, extract_report) =
        kp20_paper_border_from_fillet_setback(fillet_surfaces, vertex_position)?;
    let arcs = &extract_report.setback_arcs;
    let n = arcs.len();
    if n < 3 {
        return Err(format!("setback border has {n} arcs, need >= 3"));
    }

    // ---- ring: two sectors per arc (split at the arc midpoint) ----------
    // One bi-2 strip per whole arc cannot follow the sweep's tangent
    // plane, which rotates by roughly the fillet wrap angle over the
    // arc (measured 0.12 rad mid-arc G1 error). Halving the span
    // quarters the rotation per sector. Valence 2n must stay within
    // PnS's ≤ 12, i.e. n ≤ 6.
    if 2 * n > 12 {
        return Err(format!(
            "valence {n}: split ring would exceed PnS's max valence 12"
        ));
    }
    // Per-arc half cubics, in ring order.
    let half_cubics: Vec<([Point3; 4], usize)> = arcs
        .iter()
        .flat_map(|arc| {
            let (a, b) = split_cubic(&arc.cubic, 0.5);
            [(a, arc.fillet_index), (b, arc.fillet_index)]
        })
        .collect();
    let mut sectors = Vec::with_capacity(2 * n);
    for (k, (cubic, fillet_index)) in half_cubics.iter().enumerate() {
        let at = |t: f64| bezier3_eval(cubic, t);
        let tan = |t: f64| bezier3_tangent(cubic, t);
        let p0 = at(0.0);
        let p_mid = at(0.5);
        let p2 = at(1.0);
        let mid_bb = Point3::new(
            2.0 * p_mid.x - 0.5 * (p0.x + p2.x),
            2.0 * p_mid.y - 0.5 * (p0.y + p2.y),
            2.0 * p_mid.z - 0.5 * (p0.z + p2.z),
        );
        let arc_len = p0.distance_to(&p_mid) + p_mid.distance_to(&p2);
        let mag = arc_len * DERIV_MAG_FRAC;
        let sweep = &fillet_surfaces[*fillet_index].surface;
        let d_mid = surface_transversal(sweep, p_mid, tan(0.5), vertex_position, mag)
            .ok_or_else(|| format!("sector {k}: no transversal at mid"))?;
        // Sector ends: an arc-midpoint join (smooth on the sweep —
        // sample the sweep) or a corner join (lateral plane from the
        // adjacent tangents).
        let prev = &half_cubics[(k + 2 * n - 1) % (2 * n)];
        let next = &half_cubics[(k + 1) % (2 * n)];
        let d0 = if k % 2 == 1 {
            // Start of this half = the arc's interior midpoint m.
            surface_transversal(sweep, p0, tan(0.0), vertex_position, mag)
                .ok_or_else(|| format!("sector {k}: no transversal at m"))?
        } else {
            corner_transversal(
                bezier3_tangent(&prev.0, 1.0),
                tan(0.0),
                p0,
                vertex_position,
                mag,
            )?
        };
        let d2 = if k % 2 == 0 {
            // End of this half = the arc's interior midpoint m.
            surface_transversal(sweep, p2, tan(1.0), vertex_position, mag)
                .ok_or_else(|| format!("sector {k}: no transversal at m"))?
        } else {
            corner_transversal(
                tan(1.0),
                bezier3_tangent(&next.0, 0.0),
                p2,
                vertex_position,
                mag,
            )?
        };
        sectors.push(Bi2BorderStrip {
            pos: [p0, mid_bb, p2],
            deriv: [d0, d_mid, d2],
        });
    }
    let ring = Bi2TensorBorderRing { sectors };

    // ---- paper-faithful assembly ----------------------------------------
    let pole = pole_from_boundary(&ring);
    let sigma = sigma_center_jet_default(n).map_err(|e| format!("PnS σ jet: {e}"))?;
    let beta = BetaReparam::bi4(n);
    let cap = assemble_pns_bi4_cap(&CapAssembleInput {
        boundary: ring,
        beta,
        sigma,
        pole,
        kappa: DEFAULT_KAPPA,
    })
    .map_err(|e| format!("PnS assembly: {e}"))?;
    let mut patches = cap.patches;

    // ---- shell gates: measure, reconcile, re-measure ---------------------
    // Neighbor union for projection: the sweeps plus each corner's
    // tangent plane (origin at the corner, normal from the adjacent
    // arc tangents).
    let corner_planes: Vec<(Point3, Vec3)> = (0..n)
        .filter_map(|k| {
            let prev = &arcs[(k + n - 1) % n];
            let arc = &arcs[k];
            let t_prev = bezier3_tangent(&prev.cubic, 1.0);
            let t_this = bezier3_tangent(&arc.cubic, 0.0);
            let normal = t_prev.cross(&t_this);
            let len = normal.length();
            (len > 1e-12).then(|| (bezier3_eval(&arc.cubic, 0.0), normal * (1.0 / len)))
        })
        .collect();
    let project = |p: Point3| -> Option<Point3> {
        let mut best: Option<(f64, Point3)> = None;
        for fsr in fillet_surfaces {
            if let Some((u, v)) = project_point_to_uv(&fsr.surface, &p) {
                let q = fsr.surface.point_at(u, v);
                let d = q.distance_to(&p);
                if best.map(|(bd, _)| d < bd).unwrap_or(true) {
                    best = Some((d, q));
                }
            }
        }
        for (origin, normal) in &corner_planes {
            let d = (p - *origin).dot(normal);
            let q = Point3::new(p.x - normal.x * d, p.y - normal.y * d, p.z - normal.z * d);
            let dd = d.abs();
            if best.map(|(bd, _)| dd < bd).unwrap_or(true) {
                best = Some((dd, q));
            }
        }
        best.and_then(|(d, q)| (d < 5.0e-3).then_some(q))
    };

    let g0_border_before = border_gap(&patches, &half_cubics);
    // Reconcile: the degree-reduction gap is TANGENTIAL (the bi-2
    // border slides along the surface), so perpendicular projection
    // cannot close it. The bi-4 border CAN represent the cubic arc
    // exactly — replace the border column (j=4) with the arc's
    // degree-4 elevation, translating the j=3 column by the same
    // per-CP deltas (preserving the paper's G1 cross field relative
    // to the border). The border may run against the arc: match by
    // endpoint distance.
    let _ = &project;
    for (k, patch) in patches.iter_mut().enumerate() {
        let mut elevated = elevate_cubic_to_quartic(&half_cubics[k].0);
        let border_start = patch[0][4];
        if elevated[0].distance_to(&border_start) > elevated[4].distance_to(&border_start) {
            elevated.reverse();
        }
        for i in 0..5 {
            let delta = Vec3::new(
                elevated[i].x - patch[i][4].x,
                elevated[i].y - patch[i][4].y,
                elevated[i].z - patch[i][4].z,
            );
            for col in [4usize, 3] {
                let g = patch[i][col];
                patch[i][col] = Point3::new(g.x + delta.x, g.y + delta.y, g.z + delta.z);
            }
        }
    }
    let g0_border_after = border_gap(&patches, &half_cubics);

    // Tangency correction: project the cross-derivative field onto
    // the sweep tangent plane at the three INTERIOR interpolation
    // nodes (end nodes pinned — their j=3 entries lie on the radial
    // seam rows shared with the neighbor patch, and touching them
    // would open the seams). The cross derivative at the border is
    // ∂v ∝ Σ (col4 − col3)·B_i(u); correcting col3 by the Bernstein
    // interpolation of the out-of-plane components zeroes the
    // measured G1 error at the nodes.
    for (k, patch) in patches.iter_mut().enumerate() {
        let sweep = &fillet_surfaces[half_cubics[k].1].surface;
        for _pass in 0..2 {
            let ts: [f64; 5] = core::array::from_fn(|k2| k2 as f64 / 4.0);
            let border: [Point3; 5] = core::array::from_fn(|i| patch[i][4]);
            let col3: [Point3; 5] = core::array::from_fn(|i| patch[i][3]);
            let mut rhs = [[0.0_f64; 3]; 5];
            let mut ok = true;
            for (node, &t) in ts.iter().enumerate() {
                if node == 0 || node == 4 {
                    continue; // pinned
                }
                let p = bezier4_eval(&border, t);
                let Some((u, v)) = project_point_to_uv(sweep, &p) else {
                    ok = false;
                    break;
                };
                let su = sweep.partial_u(u, v);
                let sv = sweep.partial_v(u, v);
                let Some(nrm) = su.cross(&sv).normalized() else {
                    ok = false;
                    break;
                };
                // Cross vector at t: D(t) = Σ (col4−col3)_i B_i(t).
                let bu = bernstein4_all(t);
                let mut d = Vec3::new(0.0, 0.0, 0.0);
                for i in 0..5 {
                    let diff = Vec3::new(
                        border[i].x - col3[i].x,
                        border[i].y - col3[i].y,
                        border[i].z - col3[i].z,
                    );
                    d = d + diff * bu[i];
                }
                // Out-of-plane component to remove from D — i.e. ADD
                // to col3's blend at t.
                let out = nrm * d.dot(&nrm);
                rhs[node] = [out.x, out.y, out.z];
            }
            if !ok {
                break;
            }
            let mut mat = [[0.0_f64; 5]; 5];
            for (node, &t) in ts.iter().enumerate() {
                mat[node] = bernstein4_all(t);
            }
            let Some(deltas) = solve_5(&mat, &rhs) else {
                break;
            };
            for i in 0..5 {
                let g = patch[i][3];
                patch[i][3] =
                    Point3::new(g.x + deltas[i][0], g.y + deltas[i][1], g.z + deltas[i][2]);
            }
        }
    }

    // G1 angle along the border after reconciliation.
    let mut g1_angle_after = 0.0_f64;
    for (k, patch) in patches.iter().enumerate() {
        let sweep = &fillet_surfaces[half_cubics[k].1].surface;
        let border: [Point3; 5] = core::array::from_fn(|i| patch[i][4]);
        for i in 1..8 {
            let t = i as f64 / 8.0;
            let p = bezier4_eval(&border, t);
            let cap_n = patch_normal(patch, t, 1.0);
            let Some((u, v)) = project_point_to_uv(sweep, &p) else {
                continue;
            };
            let su = sweep.partial_u(u, v);
            let sv = sweep.partial_v(u, v);
            let sn = su.cross(&sv);
            if let (Some(a), Some(b)) = (cap_n.normalized(), sn.normalized()) {
                let cosang = a.dot(&b).abs().clamp(-1.0, 1.0);
                g1_angle_after = g1_angle_after.max(cosang.pacos());
            }
        }
    }

    // Radial seam G0 between adjacent sectors.
    let mut seam_g0 = 0.0_f64;
    let np = patches.len();
    for k in 0..np {
        let a = &patches[k];
        let b = &patches[(k + 1) % np];
        // Candidate seam rows/cols of a vs b: compare closest pairing
        // of a's 4 sides against b's 4 sides at 5 samples.
        let sides_a = patch_sides(a);
        let sides_b = patch_sides(b);
        let mut best_pair = f64::INFINITY;
        for sa in &sides_a {
            for sb in &sides_b {
                let mut worst = 0.0_f64;
                for i in 0..5 {
                    let d_fwd = sa[i].distance_to(&sb[i]);
                    let d_rev = sa[i].distance_to(&sb[4 - i]);
                    worst = worst.max(d_fwd.min(d_rev));
                }
                best_pair = best_pair.min(worst);
            }
        }
        seam_g0 = seam_g0.max(best_pair);
    }

    // First-class provenance (PK-200), same pattern as the promote
    // path. The 2026-07-27 revert here was a mis-attribution — the
    // strict-export crash proved pre-existing and unrelated (PK-205)
    // — so the switch was never falsified; it ships now, measured
    // against the cap gates and the corpus. The assembler stores
    // `patch[u_c][v_c]`; the variant's net is `net[j][i]` with `i`
    // the u index — transpose or the shipped chart swaps u and v.
    // Valence is the number of sectors meeting at the pole — this
    // ring splits every arc at its midpoint, so that is 2n, not n.
    let cap_valence = patches.len();
    let surfaces = patches
        .iter()
        .enumerate()
        .map(|(sector, patch)| {
            let net: [[Point3; 5]; 5] =
                core::array::from_fn(|j| core::array::from_fn(|i| patch[i][j]));
            SurfaceGeometry::PolyhedralSpline(Box::new(
                popkern_core::geometry::polyhedral_spline_patch::PolyhedralSplinePatch::new(
                    net,
                    popkern_core::geometry::polyhedral_spline_patch::PolyhedralSplineProvenance {
                        valence: cap_valence,
                        sector,
                    },
                ),
            ))
        })
        .collect();
    Ok(PnsSetbackCandidate {
        patches,
        surfaces,
        half_cubics,
        pole,
        report: PnsSetbackReport {
            g0_border_before,
            g0_border_after,
            g1_angle_after,
            seam_g0,
            n,
        },
    })
}

/// Max distance from each patch's border row curve to its arc cubic.
fn border_gap(patches: &[[[Point3; 5]; 5]], half_cubics: &[([Point3; 4], usize)]) -> f64 {
    let mut worst = 0.0_f64;
    for (k, patch) in patches.iter().enumerate() {
        let cubic = half_cubics[k].0;
        let border: [Point3; 5] = core::array::from_fn(|i| patch[i][4]);
        for i in 0..=16 {
            let t = i as f64 / 16.0;
            let p = bezier4_eval(&border, t);
            let mut best = f64::INFINITY;
            for j in 0..=64 {
                let s = j as f64 / 64.0;
                best = best.min(bezier3_eval(&cubic, s).distance_to(&p));
            }
            worst = worst.max(best);
        }
    }
    worst
}

/// Transversal direction on a surface: in-surface perpendicular to
/// the border tangent, oriented away from the vertex, scaled to
/// `mag`.
fn surface_transversal(
    surface: &SurfaceGeometry,
    p: Point3,
    tangent: Vec3,
    vertex: Point3,
    mag: f64,
) -> Option<Point3> {
    let (u, v) = project_point_to_uv(surface, &p)?;
    let su = surface.partial_u(u, v);
    let sv = surface.partial_v(u, v);
    let normal = su.cross(&sv);
    oriented_perp(normal, tangent, p, vertex, mag)
}

/// Transversal at a corner: the lateral plane's normal is spanned by
/// the two adjacent arc tangents.
fn corner_transversal(
    t_in: Vec3,
    t_out: Vec3,
    p: Point3,
    vertex: Point3,
    mag: f64,
) -> Result<Point3, String> {
    let normal = t_in.cross(&t_out);
    if normal.length() < 1e-12 {
        return Err("degenerate corner tangents".into());
    }
    // Border direction at the corner: mean of the tangents.
    let tangent = t_in + t_out;
    oriented_perp(normal, tangent, p, vertex, mag).ok_or_else(|| "degenerate corner".into())
}

fn oriented_perp(
    normal: Vec3,
    tangent: Vec3,
    p: Point3,
    vertex: Point3,
    mag: f64,
) -> Option<Point3> {
    let perp = normal.cross(&tangent);
    let len = perp.length();
    if len < 1e-12 {
        return None;
    }
    let unit = perp * (1.0 / len);
    let to_vertex = Vec3::new(vertex.x - p.x, vertex.y - p.y, vertex.z - p.z);
    // PnS derivative rows point INTO the surrounding surface (away
    // from the cap), matching the mitred-boundary adapter convention.
    let sign = if unit.dot(&to_vertex) > 0.0 {
        -1.0
    } else {
        1.0
    };
    Some(Point3::new(
        unit.x * mag * sign,
        unit.y * mag * sign,
        unit.z * mag * sign,
    ))
}

/// Normal of a bi-4 patch at (u, v) by Bernstein partials.
fn patch_normal(patch: &[[Point3; 5]; 5], u: f64, v: f64) -> Vec3 {
    // point(u,v) = Σ P[i][j] B_i(u) B_j(v) with i indexing u.
    let bu = bernstein4_all(u);
    let bv = bernstein4_all(v);
    let dbu = bernstein4_deriv_all(u);
    let dbv = bernstein4_deriv_all(v);
    let mut su = Vec3::new(0.0, 0.0, 0.0);
    let mut sv = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..5 {
        for j in 0..5 {
            let p = patch[i][j];
            su = su + Vec3::new(p.x, p.y, p.z) * (dbu[i] * bv[j]);
            sv = sv + Vec3::new(p.x, p.y, p.z) * (bu[i] * dbv[j]);
        }
    }
    su.cross(&sv)
}

fn patch_sides(p: &[[Point3; 5]; 5]) -> [[Point3; 5]; 4] {
    [
        core::array::from_fn(|j| p[0][j]),
        core::array::from_fn(|j| p[4][j]),
        core::array::from_fn(|i| p[i][0]),
        core::array::from_fn(|i| p[i][4]),
    ]
}

#[path = "pns_setback_candidate_bezier.rs"]
mod bezier;
use bezier::*;
