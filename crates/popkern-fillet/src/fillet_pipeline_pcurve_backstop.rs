//! Displaced-pcurve backstop (PK-211), run at pipeline end.
//!
//! A pcurve authored against TRANSITIONAL mid-carve geometry (the kef
//! merge hands the fillet face half-edges whose 3-D edges are not yet
//! on its surface) can survive every repair pass — the attach fit was
//! made against a surface the edge was 1.8 off of, and by the time
//! the geometry settles, the guards see a present pcurve and keep it
//! (measured: offset_bracket's strip-end rulings shipped images 1.84
//! off their own lines, read by the blend evidence as a phantom 0.43
//! rail error). At pipeline end the geometry is final: any pcurve
//! whose image leaves its edge curve MACROSCOPICALLY (beyond 1e-3 and
//! the imprecision budget) flags its face for a recompute, whose
//! survival rule only ever replaces with a strictly better-fitting
//! candidate.

use popkern_core::topology::BRep;

/// Flag faces carrying a macroscopically displaced pcurve and
/// recompute their pcurves against the final geometry. See module
/// docs for the defect class this closes.
pub(crate) fn recompute_displaced_pcurves(brep: &mut BRep) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let gate = 1.0e-3_f64.max(2.0 * brep.max_imprecision);
    let mut flagged: Vec<popkern_core::topology::FaceHandle> = Vec::new();
    for (_eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let ce: &dyn CurveEval = &curve.geometry;
        let refs: Vec<_> = (0..=16)
            .map(|i| ce.point_at(curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 16.0))
            .collect();
        let Some(hd) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        for he in [edge.half_edge, hd.twin] {
            let Some(h) = brep.half_edges.get(he) else {
                continue;
            };
            let Some(pc) = h.pcurve.as_ref() else {
                continue;
            };
            let Some(face) = brep.loops.get(h.loop_).map(|l| l.face) else {
                continue;
            };
            if flagged.contains(&face) {
                continue;
            }
            let Some(sf) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let displaced = [0.25, 0.5, 0.75].iter().any(|f| {
                let uv = pc.point_at(pc.t_start + f * (pc.t_end - pc.t_start));
                let img = sf.geometry.point_at(uv.x, uv.y);
                refs.iter()
                    .map(|r| img.distance_to(r))
                    .fold(f64::MAX, f64::min)
                    > gate
            });
            if displaced {
                flagged.push(face);
            }
        }
    }
    for f in flagged {
        let _ = brep.recompute_face_pcurves(f);
    }
}

/// Final pcurve backstops, in order: the PK-211 displaced-pcurve
/// recompute, then the PK-294 periodic-loop winding repair (late
/// surgery can leave a loop's stored pcurve chain winding ≥2 periods
/// — the flagship boss plate's mitred seam sliver welded 98 open
/// edges without it). Returns the number of winding-repaired faces.
pub(crate) fn final_pcurve_backstops(brep: &mut BRep) -> usize {
    recompute_displaced_pcurves(brep);
    brep.repair_periodic_loop_winding()
}
