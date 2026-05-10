use crate::diagnostics::PcurveEscape;
use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, TopologyRole};

/// Detect half-edge pcurves that, when sampled, produce (u, v) values
/// outside the surface's valid parameter range. Slicers evaluating
/// such pcurves get garbage from the surface and treat the face as
/// non-manifold. Limited to surfaces with bounded UV (NURBS, sphere,
/// torus); planes and cylinders have unbounded ranges where escape
/// isn't well-defined.
pub(super) fn detect_pcurve_escapes(brep: &BRep) -> Vec<PcurveEscape> {
    let mut out = Vec::new();
    let n_samples = 11;
    for (heh, he) in &brep.half_edges {
        if brep
            .edges
            .get(he.edge)
            .is_some_and(|edge| edge.role == TopologyRole::ConstructionResidual)
        {
            continue;
        }
        let Some(pcurve) = he.pcurve.as_ref() else {
            continue;
        };
        let Some(loop_) = brep.loops.get(he.loop_) else {
            continue;
        };
        let face = loop_.face;
        let Some(face_data) = brep.faces.get(face) else {
            continue;
        };
        if face_data.role == TopologyRole::ConstructionResidual {
            continue;
        }
        let Some(surf) = face_data.surface.as_ref() else {
            continue;
        };
        let bounded = matches!(
            &surf.geometry,
            SurfaceGeometry::Nurbs(_)
                | SurfaceGeometry::Sphere { .. }
                | SurfaceGeometry::Torus { .. }
                | SurfaceGeometry::Gregory(_)
        );
        if !bounded {
            continue;
        }
        let ((u_lo, u_hi), (v_lo, v_hi)) = surf.geometry.parameter_range();
        if !u_lo.is_finite() || !u_hi.is_finite() || !v_lo.is_finite() || !v_hi.is_finite() {
            continue;
        }
        let (t_min, t_max) = pcurve.parameter_range();
        let mut max_overshoot = 0.0_f64;
        let mut worst_uv = Point3::new(0.0, 0.0, 0.0);
        let mut worst_t = t_min;
        for i in 0..n_samples {
            let s = i as f64 / (n_samples - 1) as f64;
            let t = t_min + s * (t_max - t_min);
            let uv = pcurve.point_at(t);
            let mut over = 0.0_f64;
            if uv.x < u_lo {
                over = over.max(u_lo - uv.x);
            }
            if uv.x > u_hi {
                over = over.max(uv.x - u_hi);
            }
            if uv.y < v_lo {
                over = over.max(v_lo - uv.y);
            }
            if uv.y > v_hi {
                over = over.max(uv.y - v_hi);
            }
            if over > max_overshoot {
                max_overshoot = over;
                worst_uv = uv;
                worst_t = t;
            }
        }
        if max_overshoot > 1e-9 {
            out.push(PcurveEscape {
                half_edge: heh,
                edge: he.edge,
                face,
                parameter: worst_t,
                uv: worst_uv,
                uv_range: ((u_lo, u_hi), (v_lo, v_hi)),
                max_overshoot,
            });
        }
    }
    out
}
