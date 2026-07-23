//! PK-141 corner-cap layer, smooth-through class: absorb the
//! unconsumed remnants at valence-2 chain junctions into the
//! adjacent strips.
//!
//! On a region-path body (curved supports), a chain vertex where two
//! fillet sub-edges meet tangent-continuously (the wall-seam
//! crossing of the boss chain) leaves CORNER TERRITORY between the
//! two Phase A anchor cuts: a support sliver straddling the junction
//! (the boss-top plane quad) and a band remnant on the wall side
//! (fenced off by the anchor cut). The vertex blend has no corner to
//! build — the chain continues smoothly — so the territory belongs
//! to the strips themselves. This pass:
//!
//! 1. finds each valence-2 junction vertex of the operation whose
//!    two sub-edges are tangent-continuous there;
//! 2. collects the remnant faces at the vertex: faces carrying a
//!    SUPPORT surface (not any of the operation's fillet surfaces)
//!    whose loop stays within the junction ball (3 x radius);
//! 3. splits any remnant straddling the junction cross-section
//!    plane (mef between its two on-plane vertices);
//! 4. absorbs each piece into the strip face on its side with the
//!    clear-surface / kef / restore pattern (same-surface merges are
//!    the region path's own mechanism).
//!
//! Transactional per junction: candidate clone, validate()-gated.
//! Declines are silent per junction (the remnants simply remain, as
//! they do today); the pass returns how many pieces were absorbed.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::EdgeGeometry;

/// Absorb smooth-through junction remnants. Returns absorbed count.
pub(crate) fn absorb_junction_remnants(
    brep: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> usize {
    let trace = std::env::var_os("PK_JUNCTION_TRACE").is_some();
    let mut absorbed = 0usize;
    for (vertex, a, b) in smooth_junctions(brep, edge_geometries) {
        let Some(v_pos) = brep.vertices.get(vertex).map(|v| v.position) else {
            continue;
        };
        let radius = edge_geometries[a]
            .effective_radius
            .max(edge_geometries[b].effective_radius);
        let (Some(ta), Some(tb)) = (
            away_tangent_at(edge_geometries, a, v_pos),
            away_tangent_at(edge_geometries, b, v_pos),
        ) else {
            continue;
        };
        let Some(normal) = (tb - ta).normalized() else {
            continue;
        };
        let (ta, tb) = (ta, tb);
        let ball = radius * 3.0;
        let mut candidate = brep.clone();
        let mut moved = 0usize;
        let remnants = junction_remnants(&candidate, edge_geometries, vertex, v_pos, ball);
        if trace {
            eprintln!(
                "[junction] {vertex:?} egs ({a},{b}) remnants {remnants:?} normal \
                 ({:.2},{:.2},{:.2})",
                normal.x, normal.y, normal.z
            );
        }
        // True-corner junctions (the tangents turn): solve the
        // corner-cap geometry (ball center, contacts, extension
        // paths). The topology carve consumes this; until it lands
        // the solve is traced so every corner's geometry is
        // measurable on real bodies.
        if trace
            && ta.dot(&tb).abs() < 0.99
            && let Some(cap) = solve_corner_cap_for_junction(edge_geometries, a, b, v_pos, radius)
        {
            {
                eprintln!(
                    "[junction] {vertex:?} corner-cap geometry: c=({:.4},{:.4},{:.4}) shared                      contact ({:.4},{:.4},{:.4})",
                    cap.center.x,
                    cap.center.y,
                    cap.center.z,
                    cap.shared_contact.x,
                    cap.shared_contact.y,
                    cap.shared_contact.z
                );
            }
        }
        for face in remnants {
            let on_plane_tol = (tolerance.linear() * 100.0).max(radius * 1e-3);
            let pieces =
                split_straddling_remnant(&mut candidate, face, v_pos, normal, on_plane_tol);
            let Some(pieces) = pieces else {
                moved = 0;
                break;
            };
            let mut ok = true;
            for piece in pieces {
                match absorb_piece(&mut candidate, piece, edge_geometries, v_pos, normal, a, b) {
                    Some(()) => moved += 1,
                    None => {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok {
                moved = 0;
                break;
            }
        }
        if moved > 0 && candidate.validate().is_ok() {
            *brep = candidate;
            absorbed += moved;
            if trace {
                eprintln!("[junction] {vertex:?} absorbed {moved} piece(s)");
            }
        } else if trace {
            eprintln!("[junction] {vertex:?} declined (moved {moved})");
        }
    }
    absorbed
}

/// Valence-2 junction vertices of the operation: smooth-through
/// (tangent-continuous) junctions AND true corners both partition
/// their remnant territory with the same dividing plane — normal =
/// difference of the two AWAY-oriented chain tangents, which is the
/// chain tangent at a smooth junction and the diagonal bisector at
/// a corner.
fn smooth_junctions(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
) -> Vec<(VertexHandle, usize, usize)> {
    use std::collections::BTreeMap;
    let mut by_vertex: BTreeMap<VertexHandle, Vec<usize>> = BTreeMap::new();
    for (i, eg) in edge_geometries.iter().enumerate() {
        by_vertex
            .entry(eg.analysis.vertex_start)
            .or_default()
            .push(i);
        by_vertex.entry(eg.analysis.vertex_end).or_default().push(i);
    }
    let trace = std::env::var_os("PK_JUNCTION_TRACE").is_some();
    let mut out = Vec::new();
    for (v, egs) in by_vertex {
        let [a, b] = egs.as_slice() else {
            if trace {
                eprintln!("[junction] {v:?} skipped: valence {}", egs.len());
            }
            continue;
        };
        let Some(v_pos) = brep.vertices.get(v).map(|d| d.position) else {
            if trace {
                eprintln!("[junction] {v:?} skipped: vertex gone");
            }
            continue;
        };
        let (Some(ta), Some(tb)) = (
            away_tangent_at(edge_geometries, *a, v_pos),
            away_tangent_at(edge_geometries, *b, v_pos),
        ) else {
            continue;
        };
        let _ = (ta, tb);
        out.push((v, *a, *b));
    }
    out
}

/// Chain tangent of `eg` at the junction, sampled at the trim end
/// NEAREST the junction position (the rails run parallel to the
/// chain; sampling the far end folds the arc's curvature into the
/// tangent — measured cos 0.968 at the top seam junction, under the
/// old smoothness gate, purely from far-end sampling), oriented to
/// point AWAY from the junction.
fn away_tangent_at(edge_geometries: &[EdgeGeometry], idx: usize, at: Point3) -> Option<Vec3> {
    use popkern_core::geometry::curve::CurveEval;
    let eg = edge_geometries.get(idx)?;
    let (t0, t1) = eg.trim_3d_1.parameter_range();
    let p_start = eg.trim_3d_1.point_at(t0);
    let p_end = eg.trim_3d_1.point_at(t1);
    let step = (t1 - t0) * 0.05;
    let (a, b) = if at.distance_to(&p_start) <= at.distance_to(&p_end) {
        (p_start, eg.trim_3d_1.point_at(t0 + step))
    } else {
        (eg.trim_3d_1.point_at(t1 - step), p_end)
    };
    let t = (b - a).normalized()?;
    // Orient away: the far sample must lie on the +t side of `at`.
    let far = if at.distance_to(&p_start) <= at.distance_to(&p_end) {
        p_end
    } else {
        p_start
    };
    Some(if (far - at).dot(&t) >= 0.0 {
        t
    } else {
        t * -1.0
    })
}

/// Remnant faces at the junction: support-surface faces whose whole
/// outer loop lies inside the junction ball. Fillet-surface faces
/// (the strips) and big support faces are excluded.
fn junction_remnants(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
    vertex: VertexHandle,
    v_pos: Point3,
    ball: f64,
) -> Vec<FaceHandle> {
    let mut out = Vec::new();
    for (fh, fd) in brep.faces.iter() {
        let Some(sof) = fd.surface.as_ref() else {
            continue;
        };
        let this = format!("{:?}", sof.geometry);
        if edge_geometries
            .iter()
            .any(|eg| format!("{:?}", eg.surface) == this)
        {
            continue;
        }
        let positions = brep.face_outer_loop_positions(fh);
        if positions.is_empty() || positions.iter().any(|p| p.distance_to(&v_pos) > ball) {
            continue;
        }
        // Must actually touch the junction vertex.
        let touches = loop_half_edges(brep, fh)
            .into_iter()
            .any(|h| brep.half_edges.get(h).map(|d| d.origin) == Some(vertex));
        if touches {
            out.push(fh);
        }
    }
    out
}

fn loop_half_edges(brep: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let mut out = Vec::new();
    let Some(fd) = brep.faces.get(face) else {
        return out;
    };
    let mut loops = vec![fd.outer_loop];
    loops.extend(fd.inner_loops.iter().copied());
    for lh in loops {
        let Some(ld) = brep.loops.get(lh) else {
            continue;
        };
        let start = ld.half_edge;
        let mut cur = start;
        for _ in 0..brep.half_edges.len() + 1 {
            if !brep.half_edges.contains_key(cur) {
                break;
            }
            out.push(cur);
            cur = brep.half_edges[cur].next;
            if cur == start {
                break;
            }
        }
    }
    out
}

/// Split a remnant straddling the junction plane along its two
/// on-plane loop vertices. Returns the face pieces (1 when the
/// remnant is one-sided, 2 after a successful split), or None when
/// the remnant straddles but no clean on-plane vertex pair exists.
fn split_straddling_remnant(
    brep: &mut BRep,
    face: FaceHandle,
    v_pos: Point3,
    normal: Vec3,
    on_plane_tol: f64,
) -> Option<Vec<FaceHandle>> {
    let side = |p: &Point3| (*p - v_pos).dot(&normal);
    let positions = brep.face_outer_loop_positions(face);
    let has_pos = positions.iter().any(|p| side(p) > on_plane_tol);
    let has_neg = positions.iter().any(|p| side(p) < -on_plane_tol);
    if !(has_pos && has_neg) {
        return Some(vec![face]);
    }
    // On-plane loop vertices.
    let hes = loop_half_edges(brep, face);
    let mut on_plane: Vec<HalfEdgeHandle> = Vec::new();
    for &h in &hes {
        let v = brep.half_edges.get(h)?.origin;
        let p = brep.vertices.get(v)?.position;
        if side(&p).abs() <= on_plane_tol
            && !on_plane
                .iter()
                .any(|&o| brep.half_edges.get(o).map(|d| d.origin) == Some(v))
        {
            on_plane.push(h);
        }
    }
    let [ha, hb] = on_plane.as_slice() else {
        return None;
    };
    if brep.half_edges.get(*ha)?.loop_ != brep.half_edges.get(*hb)?.loop_ {
        return None;
    }
    let r = brep.mef(*ha, *hb).ok()?;
    // The mef edge runs across the junction: a straight chord is the
    // correct planar-support geometry (the remnants carry plane
    // surfaces; curved remnants keep pending markers via reconcile).
    let (va, vb) = (
        brep.half_edges
            .get(brep.edges.get(r.edge)?.half_edge)?
            .origin,
        {
            let he = brep.edges.get(r.edge)?.half_edge;
            let tw = brep.half_edges.get(he)?.twin;
            brep.half_edges.get(tw)?.origin
        },
    );
    let (pa, pb) = (
        brep.vertices.get(va)?.position,
        brep.vertices.get(vb)?.position,
    );
    if let Some(c) = popkern_core::geometry::curve::CurveOnEdge::line_segment(pa, pb, 1e-12) {
        brep.edges[r.edge].curve = Some(c);
    }
    Some(vec![face, r.face])
}

/// Absorb one remnant piece into the strip on its side of the
/// junction plane (clear-surface / kef / restore).
fn absorb_piece(
    brep: &mut BRep,
    piece: FaceHandle,
    edge_geometries: &[EdgeGeometry],
    v_pos: Point3,
    normal: Vec3,
    a: usize,
    b: usize,
) -> Option<()> {
    let positions = brep.face_outer_loop_positions(piece);
    if positions.is_empty() {
        return None;
    }
    let mut mean = 0.0;
    for p in &positions {
        mean += (*p - v_pos).dot(&normal);
    }
    mean /= positions.len() as f64;
    let side_eg = if mean >= 0.0 { a } else { b };
    // The strip face carrying this eg's surface, adjacent to piece.
    let target_signature = format!("{:?}", edge_geometries[side_eg].surface);
    let other_signature = format!(
        "{:?}",
        edge_geometries[if side_eg == a { b } else { a }].surface
    );
    let mut kef_he = None;
    for h in loop_half_edges(brep, piece) {
        let tw = brep.half_edges.get(h)?.twin;
        let Some(twd) = brep.half_edges.get(tw) else {
            continue;
        };
        let f = brep.loops.get(twd.loop_)?.face;
        let sig = brep
            .faces
            .get(f)
            .and_then(|d| d.surface.as_ref())
            .map(|s| format!("{:?}", s.geometry));
        if sig.as_deref() == Some(target_signature.as_str())
            || sig.as_deref() == Some(other_signature.as_str())
        {
            // Prefer the same-side strip; accept the other side only
            // if no same-side adjacency exists.
            if sig.as_deref() == Some(target_signature.as_str()) {
                kef_he = Some((h, f));
                break;
            }
            if kef_he.is_none() {
                kef_he = Some((h, f));
            }
        }
    }
    let (he, strip) = kef_he?;
    // kef kills the twin side's face: pass the half-edge ON the
    // strip (twin on the piece) so the piece merges into the strip.
    let strip_side = brep.half_edges.get(he)?.twin;
    let strip_surface = brep.faces.get(strip)?.surface.clone();
    brep.faces[strip].surface = None;
    brep.faces[piece].surface = None;
    let merged = brep.kef(strip_side).is_ok();
    brep.faces[strip].surface = strip_surface;
    if !merged {
        return None;
    }
    Some(())
}

/// Solve the corner-cap geometry for a true-corner junction: find
/// the SHARED support (present in both fillets' support pairs, by
/// structural identity), each fillet's own support and corner-side
/// trim endpoints, reconstruct the end stations, and run the
/// rolling-ball corner solve.
fn solve_corner_cap_for_junction(
    edge_geometries: &[EdgeGeometry],
    a: usize,
    b: usize,
    v_pos: Point3,
    radius: f64,
) -> Option<super::corner_cap_geometry::CornerCapGeometry> {
    use popkern_core::geometry::curve::CurveEval;
    let eg_a = edge_geometries.get(a)?;
    let eg_b = edge_geometries.get(b)?;
    let sig = |g: &popkern_core::geometry::surface::SurfaceGeometry| format!("{g:?}");
    let pairs_a = [&eg_a.analysis.surface_1, &eg_a.analysis.surface_2];
    let pairs_b = [&eg_b.analysis.surface_1, &eg_b.analysis.surface_2];
    let (shared, own_a, own_b) = pairs_a.iter().find_map(|&sa| {
        pairs_b.iter().find_map(|&sb| {
            (sig(sa) == sig(sb)).then(|| {
                let own_a = if sig(pairs_a[0]) == sig(sa) {
                    pairs_a[1]
                } else {
                    pairs_a[0]
                };
                let own_b = if sig(pairs_b[0]) == sig(sb) {
                    pairs_b[1]
                } else {
                    pairs_b[0]
                };
                (sa, own_a, own_b)
            })
        })
    })?;
    // Corner-side trim endpoints per fillet: the endpoint of each
    // trim nearest the junction vertex.
    let near_end = |trim: &popkern_core::geometry::NurbsCurve| {
        let (t0, t1) = trim.parameter_range();
        let p0 = trim.point_at(t0);
        let p1 = trim.point_at(t1);
        if v_pos.distance_to(&p0) <= v_pos.distance_to(&p1) {
            p0
        } else {
            p1
        }
    };
    let tol = 1e-9;
    let end_a = super::corner_cap_geometry::end_station_center(
        &eg_a.analysis.surface_1,
        &eg_a.analysis.surface_2,
        near_end(&eg_a.trim_3d_1),
        near_end(&eg_a.trim_3d_2),
        radius,
        tol,
    )?;
    let end_b = super::corner_cap_geometry::end_station_center(
        &eg_b.analysis.surface_1,
        &eg_b.analysis.surface_2,
        near_end(&eg_b.trim_3d_1),
        near_end(&eg_b.trim_3d_2),
        radius,
        tol,
    )?;
    super::corner_cap_geometry::solve_corner_cap(
        shared,
        [own_a, own_b],
        [end_a, end_b],
        radius,
        tol,
    )
}
