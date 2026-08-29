//! PK-159 parallel-seam weld: collapse the duplicated mid-face seam
//! left by interference-completed parallel fillets.
//!
//! When `interference_seam` upgrades a mid-face collision, each
//! fillet's apply splits the shared support face along the seam's
//! face-plane PROJECTION, so the face survives as a zero-width
//! sliver between two coincident split edges — and the split
//! vertices sit on the face plane, not on the true seam (the strips'
//! intersection ridge lies one sagitta inside the material; measured
//! 8.4e-3 on the r=0.6 unit-cube pair, with the four support corners
//! carrying the (√2−1)·r stale stubs the terminal settle then
//! declines on).
//!
//! This pass recognizes the sliver — a created plane face whose
//! every boundary edge's far side is one of two parallel-axis
//! cylinder strips or a zero-length end pair — and:
//!
//! 1. `kef`-merges the sliver into one strip (the other long edge
//!    becomes the strip×strip seam edge);
//! 2. `kve`-contracts the zero-length end edges (welding the
//!    duplicated seam endpoints);
//! 3. re-authors the seam edge as the exact parallel-cylinder
//!    intersection line (closed-form circle×circle in the
//!    cross-section), moves its endpoint vertices onto it, rebuilds
//!    the adjacent chords, and writes iso pcurves in both strip
//!    charts.
//!
//! Runs BEFORE the terminal settle, whose corner machinery then
//! closes the four support-corner stubs against the now-exact seam.
//! Transactional per sliver: candidate clone, `validate()` gate.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Weld every recognized duplicated parallel seam. Returns the
/// number of seams welded.
pub(crate) fn weld_parallel_seams(
    source: &BRep,
    result: &mut BRep,
    tolerance: &Tolerance,
) -> usize {
    let trace = popkern_core::trace::trace_enabled("interference");
    let pre: std::collections::HashSet<FaceHandle> = source.faces.keys().collect();
    let mut welded = 0usize;
    // Slivers are independent; re-scan after each weld.
    for _ in 0..8 {
        let Some(c) = recognize_sliver(result, &pre, tolerance) else {
            break;
        };
        let mut cand = result.clone();
        match weld_one(&mut cand, &c, tolerance) {
            Ok(()) => {
                if cand.validate().is_ok() {
                    *result = cand;
                    welded += 1;
                    if trace {
                        eprintln!(
                            "[iseam] welded sliver {:?} seam {:?}",
                            c.sliver, c.seam_edge
                        );
                    }
                } else {
                    if trace {
                        eprintln!("[iseam] sliver {:?} weld failed validation", c.sliver);
                    }
                    break;
                }
            }
            Err(reason) => {
                if trace {
                    eprintln!(
                        "[iseam] sliver {sliver:?} weld declined: {reason}",
                        sliver = c.sliver
                    );
                }
                break;
            }
        }
    }
    welded
}

struct SliverCandidate {
    sliver: FaceHandle,
    /// Long edge to kef (half-edge on the STRIP side, which survives).
    kef_he: HalfEdgeHandle,
    /// The surviving long edge (becomes the seam).
    seam_edge: EdgeHandle,
    /// Zero-length end edges to contract after the merge.
    end_edges: Vec<EdgeHandle>,
    cyl_a: (Point3, Vec3, f64),
    cyl_b: (Point3, Vec3, f64),
}

fn cylinder_of(brep: &BRep, face: FaceHandle) -> Option<(Point3, Vec3, f64)> {
    match brep.faces.get(face)?.surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) => Some((*origin, *axis, *radius)),
        _ => None,
    }
}

fn recognize_sliver(
    brep: &BRep,
    pre: &std::collections::HashSet<FaceHandle>,
    tolerance: &Tolerance,
) -> Option<SliverCandidate> {
    let tol = tolerance.linear().max(1e-9) * 100.0;
    'faces: for (fh, f) in brep.faces.iter() {
        // The sliver can KEEP the pre-operation support-face handle
        // (the apply splits the support and leaves the remainder on
        // the original handle) — the created-face requirement sits
        // on the two strips instead.
        if !f.inner_loops.is_empty() {
            continue;
        }
        if !matches!(
            f.surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Plane { .. })
        ) {
            continue;
        }
        // Walk the boundary: collect long edges (far side = cylinder
        // strip) and near-zero-length edges.
        let mut long: Vec<(HalfEdgeHandle, EdgeHandle, FaceHandle)> = Vec::new();
        let mut ends: Vec<EdgeHandle> = Vec::new();
        let l = brep.loops.get(f.outer_loop)?;
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..16 {
            let he = brep.half_edges.get(cur)?;
            let tw = brep.half_edges.get(he.twin)?;
            let p0 = brep.vertices.get(he.origin)?.position;
            let p1 = brep.vertices.get(tw.origin)?.position;
            let len = p0.distance_to(&p1);
            let far = brep.loops.get(tw.loop_)?.face;
            if len <= tol {
                ends.push(he.edge);
            } else if cylinder_of(brep, far).is_some() && !pre.contains(&far) {
                long.push((he.twin, he.edge, far));
            } else {
                continue 'faces;
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
        if long.len() != 2 || long[0].2 == long[1].2 {
            continue;
        }
        let (cyl_a, cyl_b) = (cylinder_of(brep, long[0].2)?, cylinder_of(brep, long[1].2)?);
        // Parallel axes only (the analytic seam class).
        if cyl_a.1.cross(&cyl_b.1).length() > 1e-6 * cyl_b.1.length() {
            continue;
        }
        // Zero width: the two long edges coincide within tolerance.
        let he0 = brep.half_edges.get(long[0].0.min(long[1].0))?;
        let _ = he0;
        let mid = |eh: EdgeHandle| -> Option<Point3> {
            let e = brep.edges.get(eh)?;
            let he = brep.half_edges.get(e.half_edge)?;
            let tw = brep.half_edges.get(he.twin)?;
            let p0 = brep.vertices.get(he.origin)?.position;
            let p1 = brep.vertices.get(tw.origin)?.position;
            Some(Point3::new(
                0.5 * (p0.x + p1.x),
                0.5 * (p0.y + p1.y),
                0.5 * (p0.z + p1.z),
            ))
        };
        let (m0, m1) = (mid(long[0].1)?, mid(long[1].1)?);
        if m0.distance_to(&m1) > tol {
            continue;
        }
        return Some(SliverCandidate {
            sliver: fh,
            kef_he: long[0].0,
            seam_edge: long[1].1,
            end_edges: ends,
            cyl_a,
            cyl_b,
        });
    }
    None
}

fn weld_one(
    brep: &mut BRep,
    c: &SliverCandidate,
    tolerance: &Tolerance,
) -> Result<(), &'static str> {
    // 1. Merge the sliver into the strip across the kef edge.
    brep.kef(c.kef_he).map_err(|_| "kef")?;
    // 2. Contract the zero-length end pairs.
    for eh in &c.end_edges {
        let Some(e) = brep.edges.get(*eh) else {
            continue;
        };
        let he = e.half_edge;
        if brep.kve(he).is_err() {
            return Err("kve");
        }
    }
    // 3. True seam: closed-form parallel-cylinder intersection.
    let e = brep.edges.get(c.seam_edge).ok_or("seam gone")?;
    let he = brep.half_edges.get(e.half_edge).ok_or("seam he")?;
    let tw = brep.half_edges.get(he.twin).ok_or("seam twin")?;
    let (v0, v1) = (he.origin, tw.origin);
    let p0 = brep.vertices.get(v0).ok_or("v0")?.position;
    let p1 = brep.vertices.get(v1).ok_or("v1")?.position;
    let axis = c.cyl_a.1.normalized().ok_or("axis")?;
    let sep = c.cyl_b.0 - c.cyl_a.0;
    let w = sep - axis * sep.dot(&axis);
    let d = w.length();
    let (r_a, r_b) = (c.cyl_a.2, c.cyl_b.2);
    if d <= 1e-12 || d >= r_a + r_b || d <= (r_a - r_b).abs() {
        return Err("degenerate circles");
    }
    let w_hat = w.normalized().ok_or("w")?;
    let u_hat = axis.cross(&w_hat);
    let alpha = (d * d + r_a * r_a - r_b * r_b) / (2.0 * d);
    let h = (r_a * r_a - alpha * alpha).sqrt();
    let base = c.cyl_a.0 + w_hat * alpha;
    let mid_old = Point3::new(
        0.5 * (p0.x + p1.x),
        0.5 * (p0.y + p1.y),
        0.5 * (p0.z + p1.z),
    );
    let cand_p = base + u_hat * h;
    let cand_m = base - u_hat * h;
    let seam_pt = if cand_p.distance_to(&mid_old) <= cand_m.distance_to(&mid_old) {
        cand_p
    } else {
        cand_m
    };
    let project = |p: Point3| -> Point3 { seam_pt + axis * (p - seam_pt).dot(&axis) };
    let (q0, q1) = (project(p0), project(p1));
    // Move the endpoints onto the seam and rebuild adjacent chords.
    super::edge_reconcile::snap_vertex_and_adjacent_chords(brep, v0, q0, c.seam_edge);
    super::edge_reconcile::snap_vertex_and_adjacent_chords(brep, v1, q1, c.seam_edge);
    // Author the seam line (primary half-edge direction).
    let e = brep.edges.get(c.seam_edge).ok_or("seam gone")?;
    let primary_from_v0 = {
        let phe = brep.half_edges.get(e.half_edge).ok_or("seam he")?;
        phe.origin == v0
    };
    let (a, b) = if primary_from_v0 { (q0, q1) } else { (q1, q0) };
    let line = popkern_core::geometry::CurveOnEdge::line_segment(a, b, tolerance.linear())
        .ok_or("degenerate seam span")?;
    if let Some(em) = brep.edges.get_mut(c.seam_edge) {
        em.curve = Some(line);
    }
    // Iso pcurves in both strip charts via projection.
    let faces: Vec<FaceHandle> = {
        let e = brep.edges.get(c.seam_edge).ok_or("seam gone")?;
        let he = brep.half_edges.get(e.half_edge).ok_or("seam he")?;
        let tw = brep.half_edges.get(he.twin).ok_or("seam twin")?;
        [he.loop_, tw.loop_]
            .into_iter()
            .filter_map(|l| brep.loops.get(l).map(|l| l.face))
            .collect()
    };
    for f in faces {
        let _ = brep.recompute_face_pcurves(f);
    }
    Ok(())
}
