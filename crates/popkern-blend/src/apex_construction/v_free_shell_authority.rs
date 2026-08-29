//! Face-use trim-authority hygiene for the V-free consolidation.

use popkern_core::topology::BRep;

/// Face-use pcurve authorities can survive a weld with a stale
/// image. STEP export prefers them over the healed half-edge
/// pcurves (measured: one plane-face cube edge exported 2.9e-2
/// off after the aggregate consolidation while the body-side
/// pcurve was exact), so every authority whose pairing died or
/// whose image no longer tracks its edge curve is dropped here —
/// the healed half-edge pcurve is the fallback authority. The
/// gate is the body's own linear tolerance ×10, matching the
/// stale-pcurve repair convention.
pub(super) fn retire_stale_face_use_pcurve_authorities(brep: &mut BRep) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    let gate = (brep.tolerance.linear() * 10.0).max(1e-9);
    let keys: Vec<popkern_core::topology::FaceUseTrimKey> =
        brep.face_use_pcurves.keys().copied().collect();
    for key in keys {
        let pairing_alive = brep
            .half_edges
            .get(key.half_edge)
            .and_then(|he| brep.loops.get(he.loop_))
            .is_some_and(|loop_data| loop_data.face == key.face);
        if !pairing_alive {
            brep.face_use_pcurves.remove(&key);
            continue;
        }
        let Some(surface) = brep
            .faces
            .get(key.face)
            .and_then(|face| face.surface.as_ref())
            .map(|s| s.geometry.clone())
        else {
            continue;
        };
        let Some(curve) = brep
            .half_edges
            .get(key.half_edge)
            .and_then(|he| brep.edges.get(he.edge))
            .and_then(|edge| edge.curve.clone())
        else {
            continue;
        };
        let Some(pcurve) = brep.face_use_pcurves.get(&key).cloned() else {
            continue;
        };
        let (p_lo, p_hi) = pcurve.parameter_range();
        let (t_lo, t_hi) = (curve.t_start, curve.t_end);
        // Same-fraction comparison in both orientations: the
        // authority may be wound opposite to the edge curve.
        let mut worst_forward = 0.0_f64;
        let mut worst_reverse = 0.0_f64;
        for k in 0..=4 {
            let f = k as f64 / 4.0;
            let uv = pcurve.point_at(p_lo + f * (p_hi - p_lo));
            let on_surface = surface.point_at(uv.x, uv.y);
            let fwd = curve.geometry.point_at(t_lo + f * (t_hi - t_lo));
            let rev = curve.geometry.point_at(t_hi - f * (t_hi - t_lo));
            worst_forward = worst_forward.max(on_surface.distance_to(&fwd));
            worst_reverse = worst_reverse.max(on_surface.distance_to(&rev));
        }
        if worst_forward.min(worst_reverse) > gate {
            brep.face_use_pcurves.remove(&key);
        }
    }
}
