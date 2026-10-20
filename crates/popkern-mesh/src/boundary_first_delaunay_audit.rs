//! The boundary-first tessellator's Delaunay validity audit (split
//! out of `boundary_first_edge_samples.rs` for the 500-line cap).

/// Env-gated Delaunay validity audit (the delaunay-audit trace topic): count
/// circumcircle violations in a Bowyer-Watson output. Proved the
/// PK-146 seam-band fans were LEGAL (zero violations) — the defect
/// was the cocircular lattice, not the triangulation.
pub(crate) fn delaunay_audit(
    face_handle: popkern_core::topology::FaceHandle,
    bw_pts: &[(f64, f64)],
    tris: &[[usize; 3]],
) {
    if !popkern_core::trace::trace_enabled("delaunay-audit") || bw_pts.len() <= 500 {
        return;
    }
    {
        let mut violations = 0usize;
        let mut worst = 0.0_f64;
        for t in tris {
            let (a, b, c) = (bw_pts[t[0]], bw_pts[t[1]], bw_pts[t[2]]);
            let d2 = 2.0 * (a.0 * (b.1 - c.1) + b.0 * (c.1 - a.1) + c.0 * (a.1 - b.1));
            if d2.abs() < 1e-30 {
                continue;
            }
            let ux = ((a.0 * a.0 + a.1 * a.1) * (b.1 - c.1)
                + (b.0 * b.0 + b.1 * b.1) * (c.1 - a.1)
                + (c.0 * c.0 + c.1 * c.1) * (a.1 - b.1))
                / d2;
            let uy = ((a.0 * a.0 + a.1 * a.1) * (c.0 - b.0)
                + (b.0 * b.0 + b.1 * b.1) * (a.0 - c.0)
                + (c.0 * c.0 + c.1 * c.1) * (b.0 - a.0))
                / d2;
            let r2 = (a.0 - ux).powi(2) + (a.1 - uy).powi(2);
            for (k, p) in bw_pts.iter().enumerate() {
                if k == t[0] || k == t[1] || k == t[2] {
                    continue;
                }
                let dd = (p.0 - ux).powi(2) + (p.1 - uy).powi(2);
                if dd < r2 * (1.0 - 1e-9) {
                    violations += 1;
                    let depth = (r2.sqrt() - dd.sqrt()) / r2.sqrt();
                    if depth > worst {
                        worst = depth;
                    }
                    break;
                }
            }
        }
        eprintln!(
            "  [bf-delaunay-audit] face {face_handle:?}: {} tris, {violations} circumcircle violations (worst depth {worst:.3})",
            tris.len()
        );
    }
}
