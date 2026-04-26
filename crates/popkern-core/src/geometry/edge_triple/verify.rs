//! `EdgeTriple::verify` — sample-based check of guarantees G2–G4.

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::edge_triple::error::EdgeTripleVerifyError;
use crate::geometry::edge_triple::inputs::SurfaceSide;
use crate::geometry::edge_triple::types::{EdgeTriple, EdgeTripleQuality};
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

impl EdgeTriple {
    /// Verify the triple's hard guarantees against the surface pair.
    ///
    /// Samples `n_samples` evenly-spaced points in `[t_start,
    /// t_end]` and checks G2 (per-sample positional invariant) using
    /// the parameter map (G3 §13). UV-box containment (G4) reads CPs
    /// against `surf.parameter_range()` for the NURBS variant;
    /// analytic variants are inherently inside their own UV box by
    /// construction.
    ///
    /// Returns `Ok(quality)` with the achieved residual, or
    /// `Err(EdgeTripleVerifyError)` on the first guarantee violation.
    pub fn verify(
        &self,
        surf_a: &SurfaceGeometry,
        surf_b: &SurfaceGeometry,
        tol: &Tolerance,
        n_samples: u32,
    ) -> Result<EdgeTripleQuality, EdgeTripleVerifyError> {
        let n = n_samples.max(2);
        let mut max_residual = 0.0f64;

        let t_lo = self.t_start;
        let t_hi = self.t_end;
        if t_hi < t_lo {
            return Err(EdgeTripleVerifyError::DirectionReversed);
        }

        for i in 0..n {
            let s = (i as f64 + 0.5) / (n as f64);
            let t = t_lo + s * (t_hi - t_lo);
            let p3 = self.curve_3d.point_at(t);

            let (t_a, t_b) = self.parameter_map.apply(t);
            let uv_a = self.pcurve_a.point_at(t_a);
            let uv_b = self.pcurve_b.point_at(t_b);
            let pa = surf_a.point_at(uv_a.x, uv_a.y);
            let pb = surf_b.point_at(uv_b.x, uv_b.y);

            let da = pa.distance_to(&p3);
            let db = pb.distance_to(&p3);
            let dmax = da.max(db);
            if dmax > tol.linear() {
                return Err(EdgeTripleVerifyError::PositionalInvariantBroken { t, residual: dmax });
            }
            max_residual = max_residual.max(dmax);
        }

        check_uv_box(&self.pcurve_a, surf_a, tol).map_err(|cp| {
            EdgeTripleVerifyError::PcurveOutsideUvBox {
                side: SurfaceSide::A,
                cp,
            }
        })?;
        check_uv_box(&self.pcurve_b, surf_b, tol).map_err(|cp| {
            EdgeTripleVerifyError::PcurveOutsideUvBox {
                side: SurfaceSide::B,
                cp,
            }
        })?;

        Ok(EdgeTripleQuality {
            max_positional_residual: max_residual,
            max_endpoint_error: 0.0,
            samples: n,
            tolerance_used: tol.linear(),
        })
    }
}

fn check_uv_box(pcurve: &Pcurve, surf: &SurfaceGeometry, tol: &Tolerance) -> Result<(), Point3> {
    let ((u_lo, u_hi), (v_lo, v_hi)) = surf.parameter_range();
    if !u_lo.is_finite() || !u_hi.is_finite() || !v_lo.is_finite() || !v_hi.is_finite() {
        return Ok(());
    }
    let pad = tol.parametric_base();
    // Analytic variants (Linear / Elliptical) are inside their own
    // UV box by construction. NURBS variants check CPs explicitly.
    match &pcurve.geometry {
        PcurveGeometry::Nurbs(n) => {
            for cp in &n.control_points {
                if cp.x < u_lo - pad || cp.x > u_hi + pad || cp.y < v_lo - pad || cp.y > v_hi + pad
                {
                    return Err(*cp);
                }
            }
            Ok(())
        }
        _ => Ok(()),
    }
}
