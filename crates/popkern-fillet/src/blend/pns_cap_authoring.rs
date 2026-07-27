//! Analytic cap edge/pcurve authoring for the low-valence PnS caps,
//! plus the measured sub-face/patch assignment. Split from
//! `pns_rim_materialize.rs`; see PK-192's ledger for the measured
//! provenance of every mechanism here.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

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
pub(crate) fn measure_sub_face_patch_assignment(
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
pub(crate) fn author_pns_cap_edge_geometry(
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
