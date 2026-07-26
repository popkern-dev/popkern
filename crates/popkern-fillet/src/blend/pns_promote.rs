//! Promotion: `BlendBoundary` → PnS cap surfaces ready for BRep
//! insertion.
//!
//! End-to-end wrapper over the PnS pipeline:
//!
//! ```text
//! BlendBoundary
//!   → Bi2TensorBorderRing   (pns_adapter)
//!   → pole                  (pns_center eq (5))
//!   → assemble_pns_bi4_cap  (pns_cap_assemble)
//!   → Vec<NurbsSurface>     (pns_bi4_to_nurbs per patch)
//!   → Vec<SurfaceGeometry>  (wrapped)
//! ```
//!
//! The output is `N_fillets · 2 = N_seg` surfaces — one per
//! `BlendBoundary` segment, matching the natural sub-face
//! decomposition of the blend face.
//!
//! Downstream topology insertion (e.g.,
//! `poke_n_gon_face_into_n_quads` + per-sub-face surface
//! assignment) is the caller's responsibility — this module only
//! produces the geometries.

use popkern_core::error::CoreError;
use popkern_core::geometry::gregory::{
    BetaReparam, CapAssembleInput, assemble_pns_bi4_cap, pns_bi4_to_nurbs, pole_from_boundary,
    sigma_center_jet_default,
};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::boundary::BlendBoundary;
use super::pns_adapter::{PnsAdapterError, blend_boundary_to_pns_ring_with_surfaces};
use crate::construct::poke_hex::poke_n_gon_face_into_n_quads;
use crate::fillet_surface::FilletSurfaceResult;

/// Error cases for the PnS promote pipeline.
#[derive(Debug, thiserror::Error)]
pub enum PnsPromoteError {
    /// Adapter step failed (segment count / surface lookup).
    #[error("PnS promote: adapter failed: {0}")]
    Adapter(#[from] PnsAdapterError),
    /// Cap assembly step failed (bad valence, singular F_κ, etc).
    #[error("PnS promote: cap assembly failed: {0}")]
    Assembly(#[from] popkern_core::geometry::gregory::CapAssembleError),
    /// The σ center jet could not be built for this valence.
    #[error("PnS promote: σ center jet failed: {0}")]
    SigmaJet(#[from] popkern_core::geometry::gregory::SigmaJetError),
    /// BRep topology op (poke) failed.
    #[error("PnS promote: BRep operation failed: {0}")]
    Brep(#[from] CoreError),
    /// Mismatched sub-face / surface count after poke.
    #[error(
        "PnS promote: poke produced {n_sub_faces} sub-faces, expected {expected} (adapter surfaces)"
    )]
    SubFaceCountMismatch {
        /// Number of sub-faces after poke.
        n_sub_faces: usize,
        /// Expected count (= number of adapter-produced surfaces).
        expected: usize,
    },
}

/// Result of a full PnS insertion on a blend face.
#[derive(Debug, Clone)]
pub struct PnsInsertResult {
    /// `N` sub-faces after poke, each assigned a PnS NURBS surface.
    pub sub_faces: Vec<FaceHandle>,
    /// Central pole vertex (all N sub-faces meet here).
    pub central_vertex: VertexHandle,
}

/// F_κ order for the cap assembly. 3 is a good default for caps
/// in the supported valence range (n ∈ [5, 10]); could be made
/// configurable if needed.
const DEFAULT_KAPPA: usize = 3;

/// Lowest corner valence the PnS cap is dispatched for.
///
/// Held at 5 for most of this cap's life, for two reasons that have
/// both since been retired:
///
/// - `n = 4` was believed structurally impossible (σ's `1/cos(2π/n)`
///   pole). It is not — the pole is removable and cancels exactly
///   under the paper's Appendix B scalings. See
///   `popkern_core::geometry::gregory::pns_sigma`.
/// - `n = 3` was simply never wired up, even though the adapter has
///   always supported it.
///
/// Both are now measured through the same live-seam gates as the
/// higher valences (`pns_low_valence_gates`), so the floor sits at
/// the construction's true lower bound: a triangle.
///
/// This matters disproportionately because valence tracks corner
/// frequency in reverse — `n = 3` is the cube corner and `n = 4` the
/// pyramid apex, i.e. the two commonest corners in mechanical parts,
/// while `n ≥ 5` is comparatively rare.
pub const PNS_MIN_VALENCE: usize = 3;

/// Produce one `SurfaceGeometry::Nurbs` per `BlendBoundary`
/// segment, ready for assignment to a corresponding sub-face.
///
/// See module docs for the pipeline.
pub fn promote_blend_boundary_to_pns_surfaces(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
) -> Result<Vec<SurfaceGeometry>, PnsPromoteError> {
    let ring = blend_boundary_to_pns_ring_with_surfaces(boundary, fillet_surfaces, brep)?;
    let n = ring.n();
    let pole = pole_from_boundary(&ring);
    let sigma = sigma_center_jet_default(n)?;
    let beta = BetaReparam::bi4(n);
    let input = CapAssembleInput {
        boundary: ring,
        beta,
        sigma,
        pole,
        kappa: DEFAULT_KAPPA,
    };
    let cap = assemble_pns_bi4_cap(&input)?;

    let surfaces: Vec<SurfaceGeometry> = cap
        .patches
        .iter()
        .map(|patch| SurfaceGeometry::Nurbs(pns_bi4_to_nurbs(patch)))
        .collect();
    Ok(surfaces)
}

/// End-to-end PnS insertion: promote a blend boundary to N PnS
/// cap surfaces, poke the `blend_face` into N sub-faces meeting
/// at the central pole, and assign each sub-face its corresponding
/// PnS NURBS surface.
///
/// Returns the N sub-face handles + the central vertex.
///
/// The `sub_faces[k]` entry covers `blend_boundary.segments[2k]`
/// (FilletEnd arc) + `segments[2k+1]` (FaceEdge bridge) — i.e.,
/// the outer boundary of sub-face `k` is the concatenation of arc
/// `k` and bridge `k`.
pub fn insert_pns_cap_into_blend_face(
    brep: &mut BRep,
    blend_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<PnsInsertResult, PnsPromoteError> {
    // The poke op expects either 2N half-edges (even, N ≥ 3) OR N
    // half-edges (N ∈ [5, 16]). Cube corners at N=3, 4 where the
    // bridge edges were consumed by a prior kef end up with N HEs
    // outside that range. Each such HE carries a curve that spans
    // (arc_k + bridge_k); we split it at the arc/bridge junction so
    // the resulting 2N HEs correspond to BlendBoundary segments
    // 1-to-1 (poke then pairs them into sectors aligned with the
    // cap's bi-2 strip midpoint).
    let expected_n = boundary.segments.len() / 2;
    normalize_low_valence_face_for_poke(brep, blend_face, expected_n, boundary)?;

    // The committed loop can traverse the rim in the OPPOSITE
    // rotational direction to the boundary ring, and can start at any
    // segment. The poke pairs loop half-edges (2k, 2k+1) into sector
    // k, so a misaligned ring builds every cap sector for a different
    // span than the sub-face it lands on (measured 2026-07-26:
    // half-sector mispairing at n=3, cap_curve_gap 8.3e-2 once the
    // rim itself was exact). Align the boundary to the loop before
    // anything consumes it, then materialize the rim from it — at the
    // low valences the committed loop is a scaffold of stale Lines
    // whose curves and vertices only the serving route rebuilds.
    let aligned = if (3..=4).contains(&expected_n) {
        align_boundary_to_loop(brep, blend_face, boundary)
    } else {
        None
    };
    if let Some(aligned) = aligned.as_ref() {
        materialize_pns_rim_from_boundary(brep, blend_face, aligned);
    }
    let boundary = aligned.as_ref().unwrap_or(boundary);

    let surfaces = promote_blend_boundary_to_pns_surfaces(boundary, fillet_surfaces, brep)?;

    // Derive the pole the same way `promote_...` did, so the
    // topology's central vertex sits at the cap's geometric pole.
    let ring = blend_boundary_to_pns_ring_with_surfaces(boundary, fillet_surfaces, brep)?;
    let pole = pole_from_boundary(&ring);

    let poke = poke_n_gon_face_into_n_quads(brep, blend_face, pole)?;
    if poke.sub_faces.len() != surfaces.len() {
        return Err(PnsPromoteError::SubFaceCountMismatch {
            n_sub_faces: poke.sub_faces.len(),
            expected: surfaces.len(),
        });
    }

    // Assign each sub-face its patch. On the aligned low-valence path
    // the correspondence is verified by measurement (the sector pairing
    // has been observed rotated; the true patch reads ~1e-3 against the
    // face's rim, every other patch at sector scale ~5e-2); the n >= 5
    // route keeps its long-standing positional zip untouched.
    if aligned.is_some() {
        let assignment = measure_sub_face_patch_assignment(brep, &poke.sub_faces, &surfaces);
        for (face_idx, surface_idx) in assignment.into_iter().enumerate() {
            let sub_face = poke.sub_faces[face_idx];
            brep.faces[sub_face].surface = Some(SurfaceOnFace {
                geometry: surfaces[surface_idx].clone(),
                orientation: true,
            });
        }
    } else {
        for (sub_face, geometry) in poke.sub_faces.iter().zip(surfaces) {
            brep.faces[*sub_face].surface = Some(SurfaceOnFace {
                geometry,
                orientation: true,
            });
        }
    }

    // The sub-faces inherited their boundary pcurves from the pre-poke
    // blend face, whose surface had an entirely different
    // parameterization. On the aligned low-valence path the chart
    // correspondence is analytic by construction, and general 2-D
    // projection cannot be trusted to rebuild the pcurves — measured
    // 2026-07-26, it lands interior preimages (cap_uv (0.512, 0.639)
    // for a rim sample) and turns a 5e-3 surface into a 4.4e-2 pcurve
    // reading — so the cap edge geometry is authored from the chart.
    // The n >= 5 route keeps its long-standing projection recompute.
    if aligned.is_some() {
        author_pns_cap_edge_geometry(brep, &poke.sub_faces, poke.central_vertex);
    } else {
        for sub_face in &poke.sub_faces {
            let _ = brep.recompute_face_pcurves(*sub_face);
        }
    }

    Ok(PnsInsertResult {
        sub_faces: poke.sub_faces,
        central_vertex: poke.central_vertex,
    })
}

/// Re-author the cap-side edge pcurves after a later pass has
/// rewritten them (the corner-patch support completion "repairs" cap
/// boundary pcurves with 2-D projection, which on these bi-4 patches
/// lands interior preimages and replaces a rate-matched 6.5e-3 pcurve
/// with a 3.1e-2 one). The pole vertex is recovered as the unique
/// vertex shared by every cap face's outer loop.
pub(crate) fn reauthor_pns_cap_edge_geometry(brep: &mut BRep, cap_faces: &[FaceHandle]) {
    use std::collections::BTreeMap;
    let mut counts: BTreeMap<VertexHandle, usize> = BTreeMap::new();
    for &face in cap_faces {
        let Some(ld) = brep
            .faces
            .get(face)
            .and_then(|f| brep.loops.get(f.outer_loop))
        else {
            continue;
        };
        let start = ld.half_edge;
        let mut he = start;
        for _ in 0..8 {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            *counts.entry(hd.origin).or_insert(0) += 1;
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    let Some(pole) = counts
        .into_iter()
        .find(|&(_, c)| c == cap_faces.len())
        .map(|(v, _)| v)
    else {
        return;
    };
    author_pns_cap_edge_geometry(brep, cap_faces, pole);
}

/// Measure which patch belongs on which sub-face: for every
/// (face, patch) pair, the worst nearest-point distance from the
/// face's rim edge curves to the patch's `v = 1` rim row. Greedy
/// bijective assignment over ascending distance; falls back to the
/// positional identity for anything left unpaired (degenerate input).
fn measure_sub_face_patch_assignment(
    brep: &BRep,
    sub_faces: &[FaceHandle],
    surfaces: &[SurfaceGeometry],
) -> Vec<usize> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let n = sub_faces.len();
    let mut scores: Vec<(f64, usize, usize)> = Vec::with_capacity(n * n);
    for (fi, &face) in sub_faces.iter().enumerate() {
        // Collect this face's rim sample points (edges whose curve is
        // a materialized boundary segment — i.e., every outer-loop
        // edge that predates the poke and still carries a curve).
        let mut rim_samples: Vec<popkern_core::geometry::Point3> = Vec::new();
        if let Some(ld) = brep
            .faces
            .get(face)
            .and_then(|f| brep.loops.get(f.outer_loop))
        {
            let start = ld.half_edge;
            let mut he = start;
            for _ in 0..8 {
                let Some(hd) = brep.half_edges.get(he) else {
                    break;
                };
                let next = hd.next;
                // Poke-created spoke edges carry Line placeholders to
                // the central vertex; rim edges carry the materialized
                // segment curves. Distinguish by whether the twin's
                // face is another cap sub-face-to-be — cheaper and
                // orientation-free: skip edges shorter than nothing;
                // simply sample every edge, spokes included, but only
                // rim samples separate the candidates (spokes meet at
                // the shared pole for every patch).
                if let Some(curve) = brep.edges.get(hd.edge).and_then(|e| e.curve.as_ref()) {
                    let g = &curve.geometry;
                    for k in 1..8 {
                        let t = curve.t_start + (curve.t_end - curve.t_start) * f64::from(k) / 8.0;
                        rim_samples.push((g as &dyn CurveEval).point_at(t));
                    }
                }
                he = next;
                if he == start {
                    break;
                }
            }
        }
        for (si, surface) in surfaces.iter().enumerate() {
            let mut worst = 0.0_f64;
            for p in &rim_samples {
                let mut nearest = f64::MAX;
                for j in 0..=32 {
                    let u = f64::from(j) / 32.0;
                    for v in [1.0, 0.0] {
                        nearest = nearest.min(surface.point_at(u, v).distance_to(p));
                    }
                    // Spoke isos too, so spoke samples don't penalise
                    // the true patch.
                    nearest = nearest.min(surface.point_at(0.0, u).distance_to(p));
                }
                worst = worst.max(nearest);
            }
            scores.push((worst, fi, si));
        }
    }
    scores.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut assignment: Vec<Option<usize>> = vec![None; n];
    let mut surface_used = vec![false; n];
    for (_, fi, si) in scores {
        if assignment[fi].is_none() && !surface_used[si] {
            assignment[fi] = Some(si);
            surface_used[si] = true;
        }
    }
    assignment
        .into_iter()
        .enumerate()
        .map(|(i, a)| a.unwrap_or(i))
        .collect()
}

/// Author the cap sub-faces' edge curves and cap-side pcurves from
/// the assigned bi-4 patches' own charts.
///
/// Chart contract (fixed by the aligned boundary + poke pairing):
/// sub-face `k`'s rim is the patch's `v = 1` row — arc on
/// `u ∈ [0, 0.5]`, bridge on `[0.5, 1]`, junction at `(0.5, 1)` —
/// and its spokes are the `u = 1` (junction → pole) and `u = 0`
/// (pole → next junction) iso-curves, with the pole at `v = 0`.
///
/// Rim pcurves are rate-matched by a 1-D solve constrained to the
/// `v = 1` row: 2-D closest-point projection on these patches falls
/// into interior minima and cannot be used (measured; see call site).
/// Spoke edges get their true curves — the patch iso-curves — in
/// place of the poke's placeholder straight lines, which are exactly
/// the "cap surface authoring downstream rewrites these" obligation
/// the poke documents and PnS never previously honoured.
fn author_pns_cap_edge_geometry(
    brep: &mut BRep,
    sub_faces: &[FaceHandle],
    central_vertex: VertexHandle,
) {
    use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
    use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::{CurveGeometry, Point3, fitting};

    let mut spoke_curve_authored: Vec<popkern_core::topology::EdgeHandle> = Vec::new();
    // Two phases: rim edges and boundary-side spokes author curves
    // first, so interior-side spoke fits see the true seam curve
    // rather than the poke's placeholder Line.
    for phase in 0..2 {
        for &face in sub_faces {
            let Some(surface) = brep
                .faces
                .get(face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.clone())
            else {
                continue;
            };
            let Some(ld) = brep
                .faces
                .get(face)
                .and_then(|f| brep.loops.get(f.outer_loop))
            else {
                continue;
            };
            let start = ld.half_edge;
            let mut he = start;
            for _ in 0..8 {
                let Some((he_next, he_twin, he_origin, he_edge)) = brep
                    .half_edges
                    .get(he)
                    .map(|hd| (hd.next, hd.twin, hd.origin, hd.edge))
                else {
                    break;
                };
                let next = he_next;
                let dest = brep.half_edges.get(he_twin).map(|t| t.origin);
                let is_spoke = he_origin == central_vertex || dest == Some(central_vertex);
                if is_spoke {
                    // Tiling contract, read from `pns_cap_assemble`
                    // (`cur[0][v_c] = nxt[v_c][0]`): the seam curve
                    // pole -> rim-start(A) is the u = 0 BOUNDARY of the
                    // patch whose rim starts at A, and an INTERIOR curve
                    // of the other adjacent patch. So the spoke edge
                    // takes the boundary patch's u = 0 iso as its 3-D
                    // curve; the boundary-side face gets the exact
                    // Linear u = 0 pcurve, and the interior-side face
                    // gets a fitted preimage seeded along the chart
                    // diagonal (its endpoints are (0,0) and (1,1) by the
                    // sharing relations). Authoring plain isos on both
                    // sides was tried first and produced a 5.9e-2 shell
                    // gap — that guessed contract is falsified.
                    let origin_pos = brep.vertices[he_origin].position;
                    let dest_pos = dest
                        .map(|d| brep.vertices[d].position)
                        .unwrap_or(origin_pos);
                    let rim_junction = if he_origin == central_vertex {
                        dest_pos
                    } else {
                        origin_pos
                    };
                    // Boundary side iff THIS face's chart has its rim
                    // start (S(0,1)) at the junction.
                    let rim_start = surface.point_at(0.0, 1.0);
                    let is_boundary_side = rim_start.distance_to(&rim_junction) < 1.0e-6;
                    if (phase == 0) == !is_boundary_side {
                        he = next;
                        if he == start {
                            break;
                        }
                        continue;
                    }
                    let Some((t_start, t_end)) = brep
                        .edges
                        .get(he_edge)
                        .and_then(|e| e.curve.as_ref())
                        .map(|c| (c.t_start, c.t_end))
                    else {
                        he = next;
                        if he == start {
                            break;
                        }
                        continue;
                    };
                    let canonical_origin = brep
                        .edges
                        .get(he_edge)
                        .and_then(|e| brep.half_edges.get(e.half_edge))
                        .map(|c| c.origin);
                    let canonical_starts_at_pole = canonical_origin == Some(central_vertex);
                    if is_boundary_side {
                        // Author the shared edge curve from this patch's
                        // u = 0 iso (degree-4 LSQ reproduces the
                        // polynomial iso to machine precision), oriented
                        // from the canonical origin.
                        if !spoke_curve_authored.contains(&he_edge) {
                            let samples: Vec<Point3> = (0..=12)
                                .map(|j| {
                                    let s = f64::from(j) / 12.0;
                                    let v = if canonical_starts_at_pole { s } else { 1.0 - s };
                                    surface.point_at(0.0, v)
                                })
                                .collect();
                            let params: Vec<f64> = (0..=12).map(|j| f64::from(j) / 12.0).collect();
                            if let Ok(iso) = fitting::fit_curve(&samples, &params, 4, 6) {
                                let (t0, t1) = iso.parameter_range();
                                if let Some(edge_data) = brep.edges.get_mut(he_edge) {
                                    edge_data.curve = Some(CurveOnEdge {
                                        geometry: CurveGeometry::Nurbs(iso),
                                        t_start: t0,
                                        t_end: t1,
                                    });
                                    spoke_curve_authored.push(he_edge);
                                }
                            }
                        }
                        let (v_from, v_to) = if canonical_starts_at_pole {
                            (0.0, 1.0)
                        } else {
                            (1.0, 0.0)
                        };
                        // Re-read the (possibly re-authored) t range.
                        let (t_start, t_end) = brep
                            .edges
                            .get(he_edge)
                            .and_then(|e| e.curve.as_ref())
                            .map(|c| (c.t_start, c.t_end))
                            .unwrap_or((t_start, t_end));
                        brep.half_edges[he].pcurve =
                            Some(Pcurve::linear((0.0, v_from), (0.0, v_to), t_start, t_end));
                    } else {
                        // Interior side: fit the seam curve's preimage in
                        // THIS chart, seeded along the diagonal (pole at
                        // (0,0), far corner at (1,1)). Local coordinate
                        // descent per sample keeps the projection out of
                        // the interior minima that defeat the generic
                        // 2-D projector on these patches.
                        let Some(curve) = brep
                            .edges
                            .get(he_edge)
                            .and_then(|e| e.curve.as_ref())
                            .cloned()
                        else {
                            he = next;
                            if he == start {
                                break;
                            }
                            continue;
                        };
                        let refine_axis =
                            |fix: f64, lo0: f64, hi0: f64, p: &Point3, along_u: bool| {
                                let (mut lo, mut hi) = (lo0, hi0);
                                for _ in 0..20 {
                                    let m1 = lo + 0.382 * (hi - lo);
                                    let m2 = lo + 0.618 * (hi - lo);
                                    let d1 = if along_u {
                                        surface.point_at(m1, fix).distance_to(p)
                                    } else {
                                        surface.point_at(fix, m1).distance_to(p)
                                    };
                                    let d2 = if along_u {
                                        surface.point_at(m2, fix).distance_to(p)
                                    } else {
                                        surface.point_at(fix, m2).distance_to(p)
                                    };
                                    if d1 <= d2 {
                                        hi = m2;
                                    } else {
                                        lo = m1;
                                    }
                                }
                                0.5 * (lo + hi)
                            };
                        let n_samples = 12usize;
                        let (ct0, ct1) = (curve.t_start, curve.t_end);
                        // Sample in canonical curve direction; seed s from
                        // pole (0,0) to far corner (1,1) matching whether
                        // the canonical direction starts at the pole.
                        let mut uv_points = Vec::with_capacity(n_samples + 1);
                        let mut params = Vec::with_capacity(n_samples + 1);
                        let (lo_t, hi_t) = (ct0.min(ct1), ct0.max(ct1));
                        for j in 0..=n_samples {
                            let s = j as f64 / n_samples as f64;
                            let t = lo_t + s * (hi_t - lo_t);
                            let p = {
                                let g = &curve.geometry;
                                (g as &dyn CurveEval).point_at(t)
                            };
                            // Which end of the seed diagonal is this t?
                            let curve_start_is_pole = {
                                let g = &curve.geometry;
                                (g as &dyn CurveEval)
                                    .point_at(lo_t)
                                    .distance_to(&brep.vertices[central_vertex].position)
                                    < 1.0e-6
                            };
                            let seed = if curve_start_is_pole { s } else { 1.0 - s };
                            let (mut u, mut v) = (seed, seed);
                            for _ in 0..4 {
                                u = refine_axis(
                                    v,
                                    (u - 0.25).max(0.0),
                                    (u + 0.25).min(1.0),
                                    &p,
                                    true,
                                );
                                v = refine_axis(
                                    u,
                                    (v - 0.25).max(0.0),
                                    (v + 0.25).min(1.0),
                                    &p,
                                    false,
                                );
                            }
                            uv_points.push(Point3::new(u, v, 0.0));
                            params.push(t);
                        }
                        if std::env::var("POPKERN_PNS_WIRE_TRACE").is_ok_and(|v| v == "1") {
                            // Attribution: is the residual the fit, or
                            // a genuine inter-patch gap? Dense 2-D
                            // scan against the worst seam sample.
                            let mut worst_fit = 0.0_f64;
                            let mut worst_p = origin_pos;
                            for (q, t) in uv_points.iter().zip(&params) {
                                let img = surface.point_at(q.x, q.y);
                                let p = {
                                    let g = &curve.geometry;
                                    (g as &dyn CurveEval).point_at(*t)
                                };
                                let d = img.distance_to(&p);
                                if d > worst_fit {
                                    worst_fit = d;
                                    worst_p = p;
                                }
                            }
                            let mut dense_min = f64::MAX;
                            for iu in 0..=48 {
                                for iv in 0..=48 {
                                    let d = surface
                                        .point_at(f64::from(iu) / 48.0, f64::from(iv) / 48.0)
                                        .distance_to(&worst_p);
                                    dense_min = dense_min.min(d);
                                }
                            }
                            eprintln!(
                                "[pns-diag-seam] face={face:?} edge={he_edge:?} worst_fit={worst_fit:.4e} dense_min_at_worst={dense_min:.4e}"
                            );
                        }
                        if let Ok(fit) = fitting::fit_curve(&uv_points, &params, 3, 6) {
                            brep.half_edges[he].pcurve = Some(Pcurve {
                                geometry: PcurveGeometry::Nurbs(fit),
                                t_start: lo_t,
                                t_end: hi_t,
                            });
                        }
                    }
                } else {
                    if phase == 1 {
                        he = next;
                        if he == start {
                            break;
                        }
                        continue;
                    }
                    // Rim half-edge: 1-D solve along the v = 1 row, then a
                    // degree-3 fit through the recovered (u, 1) samples —
                    // rate-matched to the edge curve by sampling at edge
                    // parameters.
                    let Some(curve) = brep
                        .edges
                        .get(he_edge)
                        .and_then(|e| e.curve.as_ref())
                        .cloned()
                    else {
                        he = next;
                        if he == start {
                            break;
                        }
                        continue;
                    };
                    // Fit over an INCREASING parameter interval — the
                    // NURBS pcurve form evaluates at `t` directly, so the
                    // fit's parameters must be the declared t range.
                    let (t_from, t_to) = (
                        curve.t_start.min(curve.t_end),
                        curve.t_start.max(curve.t_end),
                    );
                    let solve_u = |p: &Point3| -> f64 {
                        let mut best = (0.0_f64, f64::MAX);
                        for j in 0..=64 {
                            let u = f64::from(j) / 64.0;
                            let d = surface.point_at(u, 1.0).distance_to(p);
                            if d < best.1 {
                                best = (u, d);
                            }
                        }
                        // Golden-section refinement around the dense winner.
                        let (mut lo, mut hi) = (
                            (best.0 - 1.0 / 64.0).max(0.0),
                            (best.0 + 1.0 / 64.0).min(1.0),
                        );
                        for _ in 0..24 {
                            let m1 = lo + 0.382 * (hi - lo);
                            let m2 = lo + 0.618 * (hi - lo);
                            if surface.point_at(m1, 1.0).distance_to(p)
                                <= surface.point_at(m2, 1.0).distance_to(p)
                            {
                                hi = m2;
                            } else {
                                lo = m1;
                            }
                        }
                        0.5 * (lo + hi)
                    };
                    let n_samples = 24usize;
                    let mut uv_points = Vec::with_capacity(n_samples + 1);
                    let mut params = Vec::with_capacity(n_samples + 1);
                    for j in 0..=n_samples {
                        let s = j as f64 / n_samples as f64;
                        let t = t_from + s * (t_to - t_from);
                        let p = {
                            let g = &curve.geometry;
                            (g as &dyn CurveEval).point_at(t)
                        };
                        uv_points.push(Point3::new(solve_u(&p), 1.0, 0.0));
                        params.push(t);
                    }
                    let fit = fitting::fit_curve(&uv_points, &params, 3, 5);
                    if let Ok(fit) = fit {
                        brep.half_edges[he].pcurve = Some(Pcurve {
                            geometry: PcurveGeometry::Nurbs(fit),
                            t_start: t_from,
                            t_end: t_to,
                        });
                    }
                }
                he = next;
                if he == start {
                    break;
                }
            }
        }
    }
}

/// If the face has exactly `n` half-edges and `n ∈ {3, 4}` (too
/// few for the poke's normal `N ∈ [5, 16]` split path), split each
/// HE at the arc/bridge junction (from `boundary`) to reach `2n`
/// half-edges. Splitting at the junction — rather than a naive
/// `t = 0.5` — keeps the sub-face loop aligned with the cap's
/// bi-2 strip midpoint (strip[1] = arc/bridge junction point),
/// so the cap's `uv = (0.5, 1)` boundary midpoint matches the
/// sub-face's mid-split vertex position.
fn normalize_low_valence_face_for_poke(
    brep: &mut BRep,
    face: FaceHandle,
    expected_n: usize,
    boundary: &BlendBoundary,
) -> Result<(), CoreError> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::topology::EdgeHandle;
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| CoreError::InvalidHandle("face not in BRep".into()))?
        .outer_loop;
    let loop_start = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("loop not in BRep".into()))?
        .half_edge;
    let mut count = 0usize;
    let mut he_edges: Vec<EdgeHandle> = Vec::new();
    let mut cursor = loop_start;
    loop {
        count += 1;
        he_edges.push(brep.half_edges[cursor].edge);
        cursor = brep.half_edges[cursor].next;
        if cursor == loop_start {
            break;
        }
        if count > 32 {
            return Err(CoreError::PreconditionFailed(
                "face loop exceeds 32 half-edges — likely unterminated".into(),
            ));
        }
    }
    if count != expected_n || !(3..=4).contains(&expected_n) {
        return Ok(());
    }
    let junctions: Vec<popkern_core::geometry::Point3> = (0..expected_n)
        .map(|k| {
            let arc = match &boundary.segments[2 * k] {
                super::boundary::BoundarySegment::FilletEnd { curve, .. } => curve,
                super::boundary::BoundarySegment::FaceEdge { curve, .. } => curve,
            };
            let (_, t_max) = arc.parameter_range();
            arc.point_at(t_max)
        })
        .collect();
    for edge in he_edges {
        let curve_on_edge = brep
            .edges
            .get(edge)
            .and_then(|e| e.curve.as_ref())
            .ok_or_else(|| CoreError::InvalidHandle("edge has no curve".into()))?
            .clone();
        let he_curve = &curve_on_edge.geometry;
        let (t_min, t_max) = he_curve.parameter_range();
        let mut best: Option<(f64, f64)> = None;
        for jp in &junctions {
            if let Some(t) = he_curve.closest_point(jp) {
                if t <= t_min + 1e-9 || t >= t_max - 1e-9 {
                    continue;
                }
                let d = he_curve.point_at(t).distance_to(jp);
                if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some((t, d));
                }
            }
        }
        let t_split = best.map(|(t, _)| t).unwrap_or(0.5 * (t_min + t_max));
        brep.split_edge(edge, t_split)?;
    }
    Ok(())
}

/// Walk `face`'s outer loop from its stored half-edge (the same
/// enumeration the poke uses) and pair each half-edge with the
/// boundary segment that LIES ON the surface across it.
///
/// Pairing is semantic, not proximity: the measured separation
/// between a segment and its own support surface is ~1e-5, against
/// ~1e-2 for every other surface at the corner — an earlier
/// endpoint-proximity matcher mis-paired arcs across strips and made
/// the rim worse. Returns `None` (nothing consumed) unless every
/// half-edge pairs with a distinct segment and consecutive loop
/// edges carry ring-adjacent segments.
fn pair_loop_half_edges_to_segments(
    brep: &BRep,
    face: FaceHandle,
    boundary: &BlendBoundary,
) -> Option<(Vec<popkern_core::topology::HalfEdgeHandle>, Vec<usize>)> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    /// True pairs measure ~1e-5 (station drift included); false
    /// pairs ~1e-2.
    const LIES_ON: f64 = 1.0e-3;
    let outer_loop = brep.faces.get(face)?.outer_loop;
    let start = brep.loops.get(outer_loop)?.half_edge;
    let mut loop_hes = Vec::new();
    let mut cursor = start;
    for _ in 0..32 {
        loop_hes.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            break;
        }
    }
    let n_seg = boundary.segments.len();
    if loop_hes.len() != n_seg {
        return None;
    }
    let mut assignment = Vec::with_capacity(n_seg);
    let mut used = vec![false; n_seg];
    for &he in &loop_hes {
        let across_surface = brep
            .half_edges
            .get(he)
            .and_then(|hd| brep.half_edges.get(hd.twin))
            .and_then(|t| brep.loops.get(t.loop_))
            .and_then(|l| brep.faces.get(l.face))
            .and_then(|f| f.surface.as_ref())?;
        let mut best: Option<(usize, f64)> = None;
        for (k, seg_used) in used.iter().enumerate() {
            if *seg_used {
                continue;
            }
            let curve = pns_segment_curve(boundary, k);
            let (t0, t1) = curve.parameter_range();
            let mut worst = 0.0_f64;
            for s in 0..=8 {
                let p = curve.point_at(t0 + (t1 - t0) * f64::from(s) / 8.0);
                let Some((u, v)) = across_surface.geometry.closest_point(&p) else {
                    worst = f64::MAX;
                    break;
                };
                worst = worst.max(across_surface.geometry.point_at(u, v).distance_to(&p));
            }
            if best.map(|(_, bd)| worst < bd).unwrap_or(true) {
                best = Some((k, worst));
            }
        }
        let (k, d) = best?;
        if d > LIES_ON {
            return None;
        }
        assignment.push(k);
        used[k] = true;
    }
    for i in 0..assignment.len() {
        let a = assignment[i] as isize;
        let b = assignment[(i + 1) % assignment.len()] as isize;
        let diff = (a - b).rem_euclid(n_seg as isize);
        if diff != 1 && diff != n_seg as isize - 1 {
            return None;
        }
    }
    Some((loop_hes, assignment))
}

/// The segment's curve regardless of variant.
fn pns_segment_curve(
    boundary: &BlendBoundary,
    k: usize,
) -> &popkern_core::geometry::nurbs_curve::NurbsCurve {
    match &boundary.segments[k] {
        super::boundary::BoundarySegment::FilletEnd { curve, .. } => curve,
        super::boundary::BoundarySegment::FaceEdge { curve, .. } => curve,
    }
}

/// Reverse a NURBS curve's direction (same image, opposite
/// parameterization). Mirror of `popkern_core`'s private
/// `reverse_nurbs_curve`.
fn reverse_pns_segment_curve(
    nurbs: &popkern_core::geometry::nurbs_curve::NurbsCurve,
) -> popkern_core::geometry::nurbs_curve::NurbsCurve {
    use popkern_core::geometry::curve::CurveEval;
    let (t_min, t_max) = nurbs.parameter_range();
    let mut knots = nurbs
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect::<Vec<_>>();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points = nurbs.control_points.clone();
    control_points.reverse();
    let mut weights = nurbs.weights.clone();
    weights.reverse();
    popkern_core::geometry::nurbs_curve::NurbsCurve {
        degree: nurbs.degree,
        knots,
        control_points,
        weights,
    }
}

/// Rebuild the boundary in the committed loop's own traversal order:
/// segment `i` of the result is the segment across loop half-edge
/// `i`, its curve oriented along the loop direction, starting from
/// the loop's stored half-edge (which the poke also starts from).
///
/// Sector parity is part of the contract: the poke pairs loop
/// half-edges `(2k, 2k+1)` into sector `k` and the cap adapter pairs
/// segments the same way, so segment 0 must be a `FilletEnd`. When
/// the loop's stored half-edge sits on a bridge instead, the loop
/// entry is rotated forward by one half-edge — an entry-pointer
/// update, not a topology change.
fn align_boundary_to_loop(
    brep: &mut BRep,
    face: FaceHandle,
    boundary: &BlendBoundary,
) -> Option<BlendBoundary> {
    let (mut loop_hes, mut assignment) = pair_loop_half_edges_to_segments(brep, face, boundary)?;
    // Rotate the loop entry so it starts on a FilletEnd (even index).
    if assignment.first().is_some_and(|&k| k % 2 == 1) {
        let outer_loop = brep.faces.get(face)?.outer_loop;
        let next = brep.half_edges.get(*loop_hes.first()?)?.next;
        brep.loops[outer_loop].half_edge = next;
        loop_hes.rotate_left(1);
        assignment.rotate_left(1);
    }
    if assignment.first().is_some_and(|&k| k % 2 == 1) {
        return None;
    }
    let n_seg = boundary.segments.len();
    let reversed = {
        let a = assignment[0] as isize;
        let b = assignment[1] as isize;
        (a - b).rem_euclid(n_seg as isize) == 1
    };
    let segments = assignment
        .iter()
        .map(|&k| {
            let curve = if reversed {
                reverse_pns_segment_curve(pns_segment_curve(boundary, k))
            } else {
                pns_segment_curve(boundary, k).clone()
            };
            match &boundary.segments[k] {
                super::boundary::BoundarySegment::FilletEnd { fillet_index, .. } => {
                    super::boundary::BoundarySegment::FilletEnd {
                        fillet_index: *fillet_index,
                        curve,
                    }
                }
                super::boundary::BoundarySegment::FaceEdge { face, .. } => {
                    super::boundary::BoundarySegment::FaceEdge { face: *face, curve }
                }
            }
        })
        .collect();
    let mut aligned = boundary.clone();
    aligned.segments = segments;
    Some(aligned)
}
/// Materialize the PnS rim: rewrite the blend face's outer-loop edge
/// curves and vertex positions from the loop-ALIGNED boundary ring.
///
/// At the low valences the committed corner loop is a mid-pipeline
/// scaffold: every edge carries a stale `Line` and consecutive edge
/// curves do not even meet (measured 2026-07-26 on the n=3 mixed
/// pyramid: 4.632e-2 curve-off-support, vertex positions that chain
/// nowhere). The incumbent direct-apex route rebuilds this rim inside
/// its own materialization step, which is why its finished bodies
/// measure clean; PnS previously authored only cap surfaces on top of
/// the scaffold, and the force-accept probe showed the resulting body
/// keeps the full defect. The boundary is exact — each segment lies
/// on its support surface at ~1e-5 or better — so it is the authority
/// for the rim.
///
/// Takes the output of [`align_boundary_to_loop`]: segment `i`
/// corresponds to loop half-edge `i` and is oriented along the loop,
/// so junction `i` (start of segment `i`) is the position of half-
/// edge `i`'s origin vertex. All-or-nothing: if the aligned segments
/// do not chain end-to-start the rewrite is abandoned before touching
/// anything, and the acceptance gate measures the scaffold as before.
fn materialize_pns_rim_from_boundary(brep: &mut BRep, face: FaceHandle, boundary: &BlendBoundary) {
    use popkern_core::geometry::CurveGeometry;
    use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};

    /// Aligned segments must chain within this slack; the ring closes
    /// to ~1e-12, so anything looser means the alignment is wrong.
    const CHAIN_SLACK: f64 = 1.0e-6;

    let Some(outer_loop) = brep.faces.get(face).map(|f| f.outer_loop) else {
        return;
    };
    let Some(start) = brep.loops.get(outer_loop).map(|l| l.half_edge) else {
        return;
    };
    let mut loop_hes = Vec::new();
    let mut cursor = start;
    for _ in 0..32 {
        loop_hes.push(cursor);
        let Some(hd) = brep.half_edges.get(cursor) else {
            return;
        };
        cursor = hd.next;
        if cursor == start {
            break;
        }
    }
    let n_seg = boundary.segments.len();
    if loop_hes.len() != n_seg {
        return;
    }
    // Verify the aligned ring chains before mutating anything.
    for i in 0..n_seg {
        let cur = pns_segment_curve(boundary, i);
        let next = pns_segment_curve(boundary, (i + 1) % n_seg);
        let (_, ct1) = cur.parameter_range();
        let (nt0, _) = next.parameter_range();
        if cur.point_at(ct1).distance_to(&next.point_at(nt0)) > CHAIN_SLACK {
            return;
        }
    }
    for (i, &he) in loop_hes.iter().enumerate() {
        let curve = pns_segment_curve(boundary, i).clone();
        let (t0, t1) = curve.parameter_range();
        let junction = curve.point_at(t0);
        let Some(hd) = brep.half_edges.get(he) else {
            continue;
        };
        let (edge, origin) = (hd.edge, hd.origin);
        // CurveOnEdge's t_start belongs to the edge's CANONICAL
        // half-edge origin, which may be the twin of `he`.
        let canonical_is_he = brep.edges.get(edge).map(|e| e.half_edge) == Some(he);
        let (t_start, t_end) = if canonical_is_he { (t0, t1) } else { (t1, t0) };
        if let Some(edge_data) = brep.edges.get_mut(edge) {
            edge_data.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve),
                t_start,
                t_end,
            });
        }
        brep.vertices[origin].position = junction;
    }
    // The moved vertices are shared with neighbouring faces whose
    // edge curves still end at the old positions — pull those curve
    // endpoints to the new vertices so every incident loop stays
    // closed. The rim curves just written end exactly at the moved
    // vertices, so this pass leaves them untouched. The threshold is
    // a gap-detection floor (same floor the pipeline's own coherence
    // pass uses), not a geometric tolerance.
    let _ = brep.repair_edge_curve_endpoint_coherence(1.0e-8);
}

/// PnS dispatch is on by default. The `POPKERN_PNS` env var is an
/// *off* switch — set to `0` to disable PnS and force every corner
/// onto the legacy rungs (rational-Gregory ladder, polar-NURBS
/// fallback, etc.). Useful for A/B comparisons or for dodging a
/// suspected PnS regression in the field without a recompile.
pub fn pns_env_enabled() -> bool {
    std::env::var("POPKERN_PNS")
        .map(|v| v != "0")
        .unwrap_or(true)
}

/// Dispatch entry point: try PnS promotion on a blend face, fall
/// through to the caller's next rung on failure. Returns
/// `Some(corner_patches)` on success; `None` if any step rejected
/// the input (invalid valence, segment count, assembly failure,
/// etc.).
///
/// This is gated by `pns_env_enabled()` at the caller's discretion.
/// When invoked and PnS succeeds, it bypasses the three-quad
/// decomposition / K-P promotion ladder entirely.
pub fn try_pns_dispatch(
    brep: &mut BRep,
    blend_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    vertex: VertexHandle,
    classification: crate::vertex_classifier::VertexBlendClassification,
    bridge_fallbacks: Vec<&'static str>,
) -> Result<Vec<crate::types::CornerPatchInfo>, PnsPromoteError> {
    use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo};

    // This step used to swallow its failure reason: it returned
    // `Option` and printed the cause behind a `POPKERN_PNS_TRACE`
    // env var, so in any normal run a PnS decline was invisible and
    // the corner silently fell through to the next rung. The reason
    // is now typed and returned, and the caller turns it into a
    // structured `pns.dispatch.declined` diagnostic — which is how
    // the n=3 / n=4 dispatch gaps were finally located.
    let result = insert_pns_cap_into_blend_face(brep, blend_face, boundary, fillet_surfaces)?;

    let n_sides = result.sub_faces.len();
    let patches = result
        .sub_faces
        .into_iter()
        .map(|face| CornerPatchInfo {
            vertex,
            n_sides,
            method: CornerMethod::PolyhedralNetSpline,
            achieved_continuity: ContinuityLevel::G1,
            classification,
            face,
            bridge_fallbacks: bridge_fallbacks.clone(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: None,
            cap_sphere: None,
            sub_patch_provenance: None,
            setback_cap: false,
            exact_boundary_cap: false,
        })
        .collect();
    Ok(patches)
}

/// Whether N=3 corners are released from direct-apex promotion so the
/// PnS cap can bid for them first.
///
/// Default **off**. This moves a whole valence class onto a different
/// construction route, which is the largest behavioural change in the
/// corner pipeline, so it stays opt-in until the corpus, the polar
/// gates and the direct-apex matrix have all been measured against it.
/// `POPKERN_PNS_FIRST_N3=1` enables it.
pub fn pns_first_at_low_valence() -> bool {
    std::env::var("POPKERN_PNS_FIRST_N3")
        .map(|v| v != "0")
        .unwrap_or(false)
}
