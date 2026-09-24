//! Rim refit for the setback cap fan: move rim vertices onto the cap
//! boundary rows and assign each rim edge the exact degree-9 Bezier
//! segment of its piece's boundary row (PK-096; row pick by
//! endpoint-corner identity since PK-225 — midpoint proximity was
//! ambiguous on wide N=3 pieces and rewrote spoke seams as WRONG grid
//! rows, fan seams g0 2.3e-2 / g1 1.0).
//!
//! Nested inside `setback_vertex_topology` (shares its private
//! bezier9 helpers and topology utilities).

use super::*;

/// Move rim vertices onto the cap boundary rows and assign each rim
/// edge (cap piece ↔ non-cap face) the exact degree-9 Bezier segment
/// of the piece's boundary row between its (re-projected) endpoints.
/// Non-rim curves that end at a moved vertex have their nearest
/// endpoint control point translated to match.
pub(super) fn refit_cap_rim(brep: &mut BRep, cap_faces: &[FaceHandle]) -> Result<(), String> {
    // Collect (rim half-edge on cap side, row, ..) jobs first.
    struct RimJob {
        edge: EdgeHandle,
        /// Cap-side half-edge (the loop side on the cap face).
        cap_he: popkern_core::topology::HalfEdgeHandle,
        /// Which grid side of the piece (0: u=0, 1: u=9-row, 2: v=0
        /// col, 3: v=9 col) — determines the pcurve's iso line.
        side: usize,
        row: [Point3; 10],
        va: VertexHandle,
        vb: VertexHandle,
    }
    let mut jobs: Vec<RimJob> = Vec::new();
    for &fh in cap_faces {
        // Boundary rows are plain cells on every bi-nonic-family
        // piece (the dual-jet blend affects interior cells only), so
        // the refit reads promoted corner pieces the same way.
        let Some(grid) = brep
            .faces
            .get(fh)
            .and_then(|f| f.surface.as_ref())
            .and_then(|s| match &s.geometry {
                popkern_core::geometry::SurfaceGeometry::Gregory(gp) => gp.binonic_grid(),
                _ => None,
            })
        else {
            continue;
        };
        let grid = *grid;
        let rows: [[Point3; 10]; 4] = [
            core::array::from_fn(|j| grid[0][j]),
            core::array::from_fn(|j| grid[9][j]),
            core::array::from_fn(|i| grid[i][0]),
            core::array::from_fn(|i| grid[i][9]),
        ];
        for he in face_loop_half_edges(brep, fh)? {
            let data = &brep.half_edges[he];
            let twin_face = brep.loops[brep.half_edges[data.twin].loop_].face;
            let twin_is_cap = cap_faces.contains(&twin_face);
            if twin_is_cap && fh > twin_face {
                // Interior cap↔cap edges are refit once (from the
                // lower-handle side); both grids agree on the shared
                // row after sweep reconciliation.
                continue;
            }
            let va = data.origin;
            let vb = brep.half_edges[data.next].origin;
            let pa = brep.vertices[va].position;
            let pb = brep.vertices[vb].position;
            // Pick the row by ENDPOINT-CORNER identity first: fan
            // edges terminate on grid corners (measured 1e-15 at
            // assignment), and two adjacent corners name the side
            // exactly. Midpoint proximity — the previous sole
            // criterion — is ambiguous on wide sectors: at N = 3 a
            // spoke's midpoint sits nearer ANOTHER row of the same
            // 120-degree piece, and the refit moved the spoke-mid
            // vertex 4.31e-3 onto the wrong row and rewrote the seam
            // curve as that row (fan seam g0 2.3e-2, g1 1.0 — the
            // whole N=3 setback decline class).
            let corners = [grid[0][0], grid[0][9], grid[9][0], grid[9][9]];
            let side_of_corner_pair = |a: usize, b: usize| -> Option<usize> {
                match (a.min(b), a.max(b)) {
                    (0, 1) => Some(0),
                    (2, 3) => Some(1),
                    (0, 2) => Some(2),
                    (1, 3) => Some(3),
                    _ => None,
                }
            };
            let nearest_corner = |p: Point3| -> (usize, f64) {
                corners
                    .iter()
                    .enumerate()
                    .map(|(i, c)| (i, p.distance_to(c)))
                    .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap_or((0, f64::INFINITY))
            };
            let (ca, da) = nearest_corner(pa);
            let (cb, db) = nearest_corner(pb);
            let corner_side = side_of_corner_pair(ca, cb).filter(|&ri| {
                let side_span = corners[ca].distance_to(&corners[cb]);
                let gate = 0.25 * side_span;
                ca != cb && da <= gate && db <= gate && {
                    let _ = ri;
                    true
                }
            });
            let ri = match corner_side {
                Some(ri) => ri,
                None => {
                    // Fallback: midpoint proximity (the historical
                    // criterion) for edges not corner-terminated.
                    let mid = match brep.edges[data.edge].curve.as_ref() {
                        Some(c) => {
                            let eval: &dyn CurveEval = &c.geometry;
                            eval.point_at(0.5 * (c.t_start + c.t_end))
                        }
                        None => Point3::new(
                            0.5 * (pa.x + pb.x),
                            0.5 * (pa.y + pb.y),
                            0.5 * (pa.z + pb.z),
                        ),
                    };
                    let mut best: Option<(f64, usize)> = None;
                    for (ri, row) in rows.iter().enumerate() {
                        let d = bezier9_nearest(row, mid).1;
                        if best.map(|(bd, _)| d < bd).unwrap_or(true) {
                            best = Some((d, ri));
                        }
                    }
                    best.ok_or("no rim row")?.1
                }
            };
            jobs.push(RimJob {
                edge: data.edge,
                cap_he: he,
                side: ri,
                row: rows[ri],
                va,
                vb,
            });
        }
    }

    // Move rim vertices onto their rows (average when a vertex is
    // claimed by multiple rows — e.g. c between two sectors; the rows
    // agree at corners to extraction accuracy).
    let mut targets: std::collections::HashMap<VertexHandle, Vec<Point3>> =
        std::collections::HashMap::new();
    for job in &jobs {
        for v in [job.va, job.vb] {
            let p = brep.vertices[v].position;
            let (t, _) = bezier9_nearest(&job.row, p);
            targets
                .entry(v)
                .or_default()
                .push(bezier9_eval(&job.row, t));
        }
    }
    // A vertex may move ONLY if every incident non-rim curve can
    // follow it: full-window NURBS whose end control point sits at
    // the vertex. Split trims carry WINDOWED curves (their endpoint
    // is curve(t_start), not a control point) — translating a CP
    // would deform the wrong span, and skipping it breaks loop
    // closure (measured 4.7e-4 gaps), so those vertices stay put and
    // the rim segments snap to them instead.
    let rim_edge_set: std::collections::HashSet<EdgeHandle> = jobs.iter().map(|j| j.edge).collect();
    let movable = |v: VertexHandle| -> bool {
        let old = brep.vertices[v].position;
        for (eh, edge) in brep.edges.iter() {
            let he = edge.half_edge;
            let va = brep.half_edges[he].origin;
            let vb = brep.half_edges[brep.half_edges[he].twin].origin;
            if va != v && vb != v {
                continue;
            }
            if rim_edge_set.contains(&eh) {
                continue;
            }
            let Some(curve) = edge.curve.as_ref() else {
                return false;
            };
            let popkern_core::geometry::CurveGeometry::Nurbs(n) = &curve.geometry else {
                return false;
            };
            let cp_end_at_v = n
                .control_points
                .first()
                .map(|p| p.distance_to(&old) < 1.0e-6)
                .unwrap_or(false)
                || n.control_points
                    .last()
                    .map(|p| p.distance_to(&old) < 1.0e-6)
                    .unwrap_or(false);
            if !cp_end_at_v {
                return false;
            }
        }
        true
    };
    let movable_set: std::collections::HashSet<VertexHandle> =
        targets.keys().copied().filter(|&v| movable(v)).collect();
    let mut moves: std::collections::HashMap<VertexHandle, (Point3, Point3)> =
        std::collections::HashMap::new();
    let mut moved_windowed_neighbour = false;
    for (v, list) in targets {
        let old = brep.vertices[v].position;
        let mut acc = [0.0_f64; 3];
        for p in &list {
            acc[0] += p.x;
            acc[1] += p.y;
            acc[2] += p.z;
        }
        let n = list.len() as f64;
        let new = Point3::new(acc[0] / n, acc[1] / n, acc[2] / n);
        if old.distance_to(&new) > 1.0e-12 {
            brep.vertices[v].position = new;
            moves.insert(v, (old, new));
            if !movable_set.contains(&v) {
                moved_windowed_neighbour = true;
            }
        }
    }
    // Vertices whose incident curves are WINDOWED splits used to stay
    // put, and the rim segments were ramped onto them — deforming the
    // exact row segment by the row-movement (~1.5e-3) right at the
    // fan corners. That deformation was the last author of the
    // G2-request fan's off-surface readings (PK-208: the tessellator
    // read a TRUE 1.48e-3 curve-vs-cap deviation on the rimsplit
    // stubs and fell back to the cracking vertex-blend path). Move those
    // vertices onto the rows too and let the dedicated endpoint-
    // coherence repair make every incident curve follow — windowed
    // curves included, which is exactly the case it exists for. The
    // surgery's commit postcondition + validate + rollback bound the
    // blast radius.
    if moved_windowed_neighbour {
        let repaired = brep.repair_edge_curve_endpoint_coherence(1.0e-8);
        if popkern_core::trace::trace_enabled("setback") {
            eprintln!(
                "  [rim-refit] windowed-neighbour vertices moved; {repaired} incident curves re-cohered"
            );
        }
    }

    // Assign exact row segments to the rim edges.
    let rim_edges: std::collections::HashSet<EdgeHandle> = jobs.iter().map(|j| j.edge).collect();
    for job in &jobs {
        let pa = brep.vertices[job.va].position;
        let pb = brep.vertices[job.vb].position;
        let (ta, _) = bezier9_nearest(&job.row, pa);
        let (tb, _) = bezier9_nearest(&job.row, pb);
        // The row's parameterization may run against the half-edge:
        // extract the ordered span and reverse the control points.
        let mut seg = if ta <= tb {
            bezier9_extract(&job.row, ta, tb)
        } else {
            let mut rev = bezier9_extract(&job.row, tb, ta);
            rev.reverse();
            rev
        };
        // Endpoints exactly at the vertices, with the deltas ramped
        // linearly across all control points (a degree-1 correction
        // in Bernstein form) so the segment stays smooth.
        let da = pa - bezier9_eval(&seg, 0.0);
        let db = pb - bezier9_eval(&seg, 1.0);
        for (i, cp) in seg.iter_mut().enumerate() {
            let f = i as f64 / 9.0;
            *cp = Point3::new(
                cp.x + da.x * (1.0 - f) + db.x * f,
                cp.y + da.y * (1.0 - f) + db.y * f,
                cp.z + da.z * (1.0 - f) + db.z * f,
            );
        }
        let mut knots = Vec::with_capacity(20);
        knots.extend(std::iter::repeat_n(0.0, 10));
        knots.extend(std::iter::repeat_n(1.0, 10));
        // Edge-curve direction convention (PK-303): the CANONICAL
        // half-edge's origin sits at t_start. `seg` runs va→vb along
        // the CAP-side half-edge; when that half-edge is the twin,
        // store the segment reversed so fraction-sampling consumers
        // (deviation metric, mesher boundary conformity, declared
        // accuracy, STEP trace sync) read the curve forward. Measured
        // before the fix: every rim edge whose cap side was the twin
        // carried a reversed curve — 300/528 half-edges persistently
        // "stale" at exactly the rim chord length (4.96e-2) on the
        // 12-edge G2 cube, unrepairable because the geometry itself
        // was exact.
        let canonical_forward = brep
            .edges
            .get(job.edge)
            .is_none_or(|e| e.half_edge == job.cap_he);
        let stored_cps: Vec<Point3> = if canonical_forward {
            seg.to_vec()
        } else {
            seg.iter().rev().copied().collect()
        };
        if let Some(edge) = brep.edges.get_mut(job.edge) {
            edge.curve = Some(popkern_core::geometry::CurveOnEdge {
                geometry: popkern_core::geometry::CurveGeometry::Nurbs(NurbsCurve {
                    degree: 9,
                    knots,
                    control_points: stored_cps,
                    weights: vec![1.0; 10],
                }),
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        // Cap-side pcurve: the row segment IS the iso line of the
        // piece domain between ta and tb, and de Casteljau extraction
        // parameterizes it linearly — so an exact, parameter-
        // synchronized pcurve is the straight domain segment. (The
        // STEP exporter demands trace synchronization; asynchronous
        // pcurves force face-local edges and an OPEN_SHELL export.)
        // grid[i][j] evaluates with B_i(u)·B_j(v): sides 0/1 are u
        // iso lines parameterized by v, sides 2/3 v iso lines by u.
        let (uv_a, uv_b) = match job.side {
            0 => ((0.0, ta), (0.0, tb)),
            1 => ((1.0, ta), (1.0, tb)),
            2 => ((ta, 0.0), (tb, 0.0)),
            _ => ((ta, 1.0), (tb, 1.0)),
        };
        if let Some(he_mut) = brep.half_edges.get_mut(job.cap_he) {
            he_mut.pcurve = Some(popkern_core::geometry::pcurve::Pcurve {
                geometry: popkern_core::geometry::pcurve::PcurveGeometry::Linear {
                    start: uv_a,
                    end: uv_b,
                },
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        // Strip-side pcurve: the twin keeps its pre-refit image (the
        // original arc path), stale against the rewritten curve by
        // the reconcile move size (measured 4.204e-4 on the unit
        // cube's worst rim seam — the whole-body G0 constant of
        // PK-392). The trim contract's authority is the served curve:
        // reproject it onto the strip surface.
        let twin = brep.half_edges[job.cap_he].twin;
        let strip_face = brep.loops[brep.half_edges[twin].loop_].face;
        if !cap_faces.contains(&strip_face)
            && let Some(surface) = brep
                .faces
                .get(strip_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.clone())
            && let Some(popkern_core::geometry::CurveGeometry::Nurbs(n)) = brep
                .edges
                .get(job.edge)
                .and_then(|e| e.curve.as_ref())
                .map(|c| c.geometry.clone())
            && let Some(fitted) =
                crate::construct::pcurve_attach_compute::compute_pcurve_on_surface(&n, &surface)
            && let Some(he_mut) = brep.half_edges.get_mut(twin)
        {
            he_mut.pcurve = Some(popkern_core::geometry::pcurve::Pcurve::from_nurbs(fitted));
        }
    }

    // Repair endpoint coherence for NON-rim curves at moved vertices:
    // translate whichever curve end is nearest the old position.
    let edge_keys: Vec<EdgeHandle> = brep.edges.keys().collect();
    for eh in edge_keys {
        if rim_edges.contains(&eh) {
            continue;
        }
        let he = brep.edges[eh].half_edge;
        let va = brep.half_edges[he].origin;
        let vb = brep.half_edges[brep.half_edges[he].twin].origin;
        for v in [va, vb] {
            let Some(&(old, new)) = moves.get(&v) else {
                continue;
            };
            let Some(curve) = brep.edges.get_mut(eh).and_then(|e| e.curve.as_mut()) else {
                continue;
            };
            if let popkern_core::geometry::CurveGeometry::Nurbs(n) = &mut curve.geometry {
                let front = n.control_points.first().map(|p| p.distance_to(&old));
                let back = n.control_points.last().map(|p| p.distance_to(&old));
                match (front, back) {
                    (Some(f), Some(b)) if f <= b && f < 1.0e-6 => {
                        *n.control_points.first_mut().ok_or("empty curve")? = new;
                    }
                    (Some(_), Some(b)) if b < 1.0e-6 => {
                        *n.control_points.last_mut().ok_or("empty curve")? = new;
                    }
                    _ => {}
                }
            }
        }
    }
    Ok(())
}
