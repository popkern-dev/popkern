//! The variational-fit cap's own pcurves (PK-574): exact Linear iso
//! pcurves on every cap half-edge. Split from `setback_fit_fan` for the
//! file budget; it shares that module's helpers.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::variational_cap::VariationalCap;
use popkern_core::topology::{BRep, FaceHandle};

use super::param_on;
use crate::construct::setback_vertex_topology::face_loop_half_edges;

/// Exact Linear iso pcurves on every cap half-edge: each edge lies on
/// one boundary iso of its sector (o at (0,0), spokes on u = 0 / v = 0,
/// rim halves on u = 1 / v = 1). The iso is chosen per edge, the one
/// that carries both ends, and the fixed coordinate is written exactly:
/// a corner vertex projected on its own fell onto the other rim at
/// u = 1 − 1e-16, and the certifier's boundary surrogate, which needs
/// exact membership (PK-541), then refused a dual-jet rim (PK-595).
pub(super) fn author_cap_pcurves(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
    cap: &VariationalCap,
) -> Result<(), String> {
    for (s, &face) in cap_faces.iter().enumerate() {
        let sec = &cap.sectors[s];
        let isos = [
            sec.iso_curve_at_u(1.0),
            sec.iso_curve_at_v(1.0),
            sec.iso_curve_at_v(0.0),
            sec.iso_curve_at_u(0.0),
        ];
        // The parameter along iso `k`, snapped to its ends.
        let along = |k: usize, p: Point3| -> f64 {
            let t = param_on(&isos[k], p);
            if t < 1e-9 {
                0.0
            } else if t > 1.0 - 1e-9 {
                1.0
            } else {
                t
            }
        };
        let uv = |k: usize, t: f64| -> (f64, f64) {
            match k {
                0 => (1.0, t),
                1 => (t, 1.0),
                2 => (t, 0.0),
                _ => (0.0, t),
            }
        };
        let miss = |k: usize, a: Point3, b: Point3| {
            isos[k]
                .point_at(along(k, a))
                .distance_to(&a)
                .max(isos[k].point_at(along(k, b)).distance_to(&b))
        };
        for he in face_loop_half_edges(brep, face)? {
            let a = brep.vertices[brep.half_edges[he].origin].position;
            let b = brep.vertices[brep.half_edges[brep.half_edges[he].twin].origin].position;
            let k = (0..4)
                .min_by(|&x, &y| miss(x, a, b).total_cmp(&miss(y, a, b)))
                .unwrap_or(0);
            let (uv_a, uv_b) = (uv(k, along(k, a)), uv(k, along(k, b)));
            brep.half_edges[he].pcurve = Some(Pcurve::linear(uv_a, uv_b, 0.0, 1.0));
        }
    }
    Ok(())
}
