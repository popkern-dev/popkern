//! Per-pair clip-and-classify tail of the parallel face-pair scan —
//! split from `ff_parallel.rs` under the 500-line budget (PK-511).

use popkern_core::geometry::intersection::SurfaceSurfaceIntersection;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::ff_parallel_residue;
use popkern_core::geometry::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;

use super::ff_parallel_trace::{TraceFlags, trace_clip_input, trace_clip_output};
use super::{FacePairTask, PairInputs};
use crate::clip::clip_curve_to_faces_cached;
use crate::intersect::IntersectionRecord;

/// Clip one pair's shared SSI curves to the two faces and classify
/// tangency — the per-pair tail of the retired `process_one_pair`.
pub(super) fn clip_and_classify(
    task: &FacePairTask,
    inp: &PairInputs<'_>,
    result: &popkern_core::geometry::intersection::SurfaceSurfaceIntersectionResult,
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
    trace: TraceFlags,
) -> Vec<IntersectionRecord> {
    let FacePairTask { face_a, face_b, .. } = *task;
    let mut records = Vec::new();
    // Pair-level pierce junction inventory, shared by every record of
    // this pair (PK-511): the inventory is record-independent and a
    // multi-wrap pair carries one record per wrap.
    let mut pierce_cache = crate::clip::PierceInventoryCache::default();
    for isect in &result.intersections {
        match isect {
            SurfaceSurfaceIntersection::Curve(curve) => {
                let mut residues: Vec<crate::clip::TrimResidue> = Vec::new();
                let clipped_list = clip_curve_to_faces_cached(
                    a,
                    face_a,
                    b,
                    face_b,
                    curve,
                    Some(&mut residues),
                    &mut pierce_cache,
                );
                ff_parallel_residue::push_residue_records(
                    &mut records,
                    residues,
                    face_a,
                    face_b,
                    trace.clip,
                );
                if trace.clip {
                    trace_clip_input(
                        face_a,
                        face_b,
                        inp.surf_a,
                        inp.surf_b,
                        curve,
                        clipped_list.len(),
                    );
                }
                for mut clipped in clipped_list {
                    if trace.clip {
                        use popkern_core::geometry::curve::CurveEval as _;
                        let (t0, t1) = clipped.curve_3d.parameter_range();
                        let p0 = clipped.curve_3d.point_at(t0);
                        let p1 = clipped.curve_3d.point_at(t1);
                        eprintln!(
                            "  [clip-piece] pre-refine t[{t0:.3},{t1:.3}] ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                            p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                        );
                    }
                    // Boundary-crossing endpoints refine onto the
                    // exact edge-curve x surface junction so the
                    // graph mints ONE node per junction (see
                    // clip_endpoint_refine).
                    crate::clip_endpoint_refine::refine_endpoints_onto_boundaries(
                        a,
                        face_a,
                        inp.surf_a,
                        (inp.r1.0, inp.r1.1),
                        b,
                        face_b,
                        inp.surf_b,
                        (inp.r2.0, inp.r2.1),
                        &mut clipped,
                        curve,
                        tolerance,
                    );
                    if trace.clip {
                        trace_clip_output(face_a, face_b, &clipped);
                    }
                    // Zero-extent clip residue: a record whose whole
                    // clipped span is one point (seam-touch junction
                    // windows on periodic pairs) mints a self-loop
                    // shared edge in the merge graph that the fragment
                    // walk chews into mixed loops (measured: torus ×
                    // cylinder seam-tangent junction; the walk closed
                    // one 10-he loop over both fragments of the face).
                    {
                        use popkern_core::geometry::curve::CurveEval as _;
                        let (c0, c1) = clipped.curve_3d.parameter_range();
                        let chord_span = clipped
                            .curve_3d
                            .point_at(c0)
                            .distance_to(&clipped.curve_3d.point_at(c1));
                        let mid_span = clipped
                            .curve_3d
                            .point_at(c0)
                            .distance_to(&clipped.curve_3d.point_at(0.5 * (c0 + c1)));
                        if chord_span.max(mid_span) < tolerance.linear() {
                            continue;
                        }
                    }
                    // Surface honesty: a kept window's interior must
                    // lie ON both carrier surfaces. A near-tangent
                    // pair's marched trace can bridge two hyperbola
                    // branches through the tangency region where no
                    // real intersection exists; the tiny bridge
                    // windows survive the UV clip (their pcurve
                    // crossings look legitimate) and strand in the
                    // merge (measured, PK-464: the hex-corner
                    // chamfer ring cone grazing the across-flats
                    // plane at its base — two 0.007-span windows per
                    // wall, 0.65 OFF the cone, 4 half-claimed record
                    // sides). The gate sits far above the record-fit
                    // slop (~1e-4) and far below any bridge.
                    let both_analytic = !matches!(
                        inp.surf_a,
                        popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
                    ) && !matches!(
                        inp.surf_b,
                        popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
                    );
                    {
                        use popkern_core::geometry::curve::CurveEval as _;
                        use popkern_core::geometry::surface::SurfaceEval as _;
                        let (c0, c1) = clipped.curve_3d.parameter_range();
                        // Analytic-analytic pairs use the tight band
                        // (their projected distances are exact). A
                        // NURBS carrier's multisheet distance can
                        // read a phantom ~0.08 for an on-surface
                        // point (the deep-grind pin), so mixed pairs
                        // keep only a COARSE structural-lie gate: a
                        // record piece sweeping the wrong way round a
                        // chart sits ~the feature radius off its
                        // carrier (measured 1.4 on the M5 tip-band
                        // arc that overwrapped the chamfer cone).
                        let band = if both_analytic {
                            tolerance.linear() * 1000.0
                        } else {
                            0.1_f64.max(tolerance.linear() * 1000.0)
                        };
                        let mut worst = 0.0_f64;
                        for k in 1..4 {
                            let t = c0 + (c1 - c0) * k as f64 / 4.0;
                            let p = clipped.curve_3d.point_at(t);
                            for (s, uv_curve) in [
                                (inp.surf_a, &clipped.curve_uv1),
                                (inp.surf_b, &clipped.curve_uv2),
                            ] {
                                let projected =
                                    popkern_core::geometry::uv_projection::project_point_to_uv_multisheet(
                                        s, &p,
                                    )
                                    .map(|(u, v)| s.point_at(u, v).distance_to(&p))
                                    .unwrap_or(0.0);
                                // The piece's own PCURVE is the
                                // chart-true witness (PK-474): on a
                                // multi-turn helical NURBS the cold
                                // multisheet projection converges
                                // into the WRONG TURN and reads a
                                // phantom ~pitch (measured 0.87 for
                                // the quarter-pitch rod's on-surface
                                // seam helix — the honest record was
                                // dropped and the seam ladder lost
                                // one link). A tangency BRIDGE stays
                                // caught either way: its 3D curve
                                // sits off the surface its own
                                // pcurve evaluates onto (0.65
                                // measured on the hex-corner cone).
                                let q = uv_curve.point_at(t);
                                let via_pcurve = s.point_at(q.x, q.y).distance_to(&p);
                                worst = worst.max(projected.min(via_pcurve));
                            }
                        }
                        if worst > band {
                            if trace.clip {
                                eprintln!(
                                    "  [clip-drop] A{face_a:?} x B{face_b:?}: window interior {worst:.4} off a carrier surface (tangency bridge)"
                                );
                            }
                            continue;
                        }
                        if trace.clip && worst > tolerance.linear() * 1000.0 {
                            eprintln!(
                                "  [clip-masked-phantom] A{face_a:?} x B{face_b:?}: interior {worst:.4} off a carrier (under the mixed-pair gate {band:.3})"
                            );
                        }
                    }
                    // NURBS carriers keep tuned behavior: the
                    // multisheet distance on a coiled wall can read a
                    // phantom deviation for an on-surface point, and
                    // the gate dropped a load-bearing deep-grind
                    // window (measured: the z=0.45 band-grind pin
                    // tore at 1021 open edges).
                    // Classify tangency at curve midpoint.
                    // Uses popkern_core::adaptive_precision::
                    // adaptive_cross_ratio_sq_lt so near-parallel
                    // normals produce deterministic decisions —
                    // the f64 fast path handles clear cases,
                    // and double-precision recomputation kicks
                    // in only when the comparison lands inside
                    // the ambiguity margin.
                    let (t0, t1) = clipped.curve_3d.parameter_range();
                    let t_mid = (t0 + t1) * 0.5;
                    let mid_uv1 = clipped.curve_uv1.point_at(t_mid);
                    let mid_uv2 = clipped.curve_uv2.point_at(t_mid);
                    let n1 = inp.surf_a.normal_at(mid_uv1.x, mid_uv1.y);
                    let n2 = inp.surf_b.normal_at(mid_uv2.x, mid_uv2.y);
                    let is_tangent = popkern_core::adaptive_precision::adaptive_cross_ratio_sq_lt(
                        (n1.x, n1.y, n1.z),
                        (n2.x, n2.y, n2.z),
                        tolerance.angular(),
                    );

                    records.push(IntersectionRecord {
                        residue_owner_a: None,
                        face_a,
                        face_b,
                        curve: clipped,
                        is_tangent,
                        is_coincident: false,
                        start_boundary: Default::default(),
                        end_boundary: Default::default(),
                    });
                }
            }
            SurfaceSurfaceIntersection::CoincidentRegion { .. } => {
                // Coincident face pair — skip. Coincidence engine
                // records the overlap separately.
            }
        }
    }
    // PK-313: land clip-crossing endpoints exactly on their face
    // iso-boundaries (mechanism and honesty bounds in
    // `ff_junction_refine`).
    crate::parallel::ff_junction_refine::refine_record_junction_endpoints(
        &mut records,
        a,
        face_a,
        inp.surf_a,
        b,
        face_b,
        inp.surf_b,
        tolerance,
    );
    records
}
