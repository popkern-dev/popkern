//! Chart handedness from a face's stored trim ring (PK-480), split
//! from `chart_branch_clustering.rs` under the file-size rule.

use popkern_math::PortableMath;

/// Chart handedness from the face's OWN stored trim ring (PK-480):
/// the walk's material-left rotation sense must follow the STORED
/// outer loop's winding about the trim interior — the same ring
/// every containment parity test reads. The graph-registered
/// boundary embeddings are an unreliable proxy for that winding:
/// split boundary pieces register in graph half-edge direction,
/// which follows the EDGE, not the face's loop traversal (measured
/// on the angled tap's crest strip: the registered sum said +1 while
/// the stored ring winds −1; every wedge fan at the entry-cone
/// circle's junctions then paired the arc arrival with its own twin,
/// the circle closed as a zero-area slit, and the merge stranded one
/// half-claimed side). Segments are re-oriented to chain head-to-tail
/// so a pcurve stored against its half-edge's direction cannot poison
/// the sign. `None` on missing pcurves or a degenerate area.
pub(crate) fn stored_ring_sigma(
    brep: &popkern_core::topology::BRep,
    face: popkern_core::topology::FaceHandle,
) -> Option<f64> {
    let ring = crate::clip::face_uv_boundary_chained(brep, face)?;
    // Scope: rings of OPEN edges only. A closed edge (rim circle,
    // both ends one vertex) carries a direction-untrustworthy pcurve
    // (older producers store edge-aligned pcurves — the PK-478
    // closed-side lesson), so its segment can silently mirror the
    // shoelace (measured: the angled tap's bore wall read +1 from its
    // rims while the registered runs correctly said −1, and flipping
    // σ broke the thread stage's bottom rim assembly). Faces with
    // closed boundary pieces keep their orientation flag.
    for (heh, _) in &ring {
        let he = brep.half_edges.get(*heh)?;
        let twin = he.mate().and_then(|m| brep.half_edges.get(m))?;
        if he.origin == twin.origin {
            return None;
        }
    }
    // A seam's two copies on one chart line: the ring runs up the seam
    // and back down it with no period between, so a pole it passes
    // round has no row run and the shoelace counts only the other
    // pieces (a ball's face round its pole, fused to a rod, the seam up
    // to the north pole and back at u = 0.4 both ways, read 2.44 for a
    // face of the whole chart, flipping σ against an honest flag:
    // PK-1187). Such a ring does not bound the face in the chart.
    let mid = |pc: &popkern_core::geometry::pcurve::Pcurve| {
        let (t0, t1) = pc.parameter_range();
        pc.point_at(0.5 * (t0 + t1))
    };
    for (i, (ha, pa)) in ring.iter().enumerate() {
        let ea = brep.half_edges.get(*ha)?.edge;
        for (hb, pb) in &ring[i + 1..] {
            if brep.half_edges.get(*hb)?.edge == ea && mid(pa).distance_to(&mid(pb)) < 1e-9 {
                return None;
            }
        }
    }
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for (_, pc) in &ring {
        let (t0, t1) = pc.parameter_range();
        let seg: Vec<(f64, f64)> = (0..=8)
            .map(|k| {
                let q = pc.point_at(t0 + (t1 - t0) * k as f64 / 8.0);
                (q.x, q.y)
            })
            .collect();
        // The ring must CHAIN head-to-tail in loop order: a segment
        // that would need re-orientation means the pcurve directions
        // are not loop-faithful and the shoelace sign is untrustworthy.
        if let (Some(&prev), Some(&s0), Some(&s1)) = (pts.last(), seg.first(), seg.last()) {
            let d_fwd = (s0.0 - prev.0).phypot(s0.1 - prev.1);
            let d_rev = (s1.0 - prev.0).phypot(s1.1 - prev.1);
            if d_rev < d_fwd {
                return None;
            }
        }
        pts.extend(seg);
    }
    if pts.len() < 3 {
        return None;
    }
    let mut area2 = 0.0;
    let mut scale = 0.0_f64;
    for i in 0..pts.len() {
        let p = pts[i];
        let q = pts[(i + 1) % pts.len()];
        area2 += p.0 * q.1 - q.0 * p.1;
        scale = scale.max(p.0.abs() + p.1.abs());
    }
    if area2.abs() < 1e-12 * (1.0 + scale * scale) {
        None
    } else {
        Some(area2.signum())
    }
}
