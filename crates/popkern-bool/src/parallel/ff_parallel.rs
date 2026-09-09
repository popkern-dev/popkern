//! Parallel face-face intersection dispatch (PK-277).
//!
//! Holds the per-pair body of Stage 1 intersection: each candidate face
//! pair is resolved to intersection records independently, so pairs can
//! be processed on separate threads. Determinism is by construction:
//! results are collected in the caller's candidate order (rayon's
//! ordered `collect`), so the record sequence is positionally identical
//! to the sequential loop regardless of thread count or scheduling.
//!
//! The per-pair computation reads `&BRep` only (plain `SlotMap` storage,
//! no interior mutability), so sharing both solids across threads is
//! sound without locks.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::intersection::{SurfaceSurfaceIntersection, intersect_surfaces};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::clip_curve_to_faces;
use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::intersect::{IntersectionRecord, bounded_parameter_range};

#[path = "ff_parallel_trace.rs"]
mod ff_parallel_trace;
use ff_parallel_trace::{trace_clip_input, trace_clip_output, trace_pair_surfaces};

/// A face pair selected by BVH overlap filtering, awaiting intersection.
#[derive(Debug, Clone, Copy)]
pub(crate) struct FacePairTask {
    /// Face from solid A.
    pub face_a: FaceHandle,
    /// Face from solid B.
    pub face_b: FaceHandle,
}

/// Debug-trace gates for the pair loop, read from the environment once
/// per Stage 1 call instead of once per pair (`env::var_os` takes the
/// process environment lock).
#[derive(Debug, Clone, Copy)]
pub(crate) struct TraceFlags {
    /// the `pair-scan` trace topic: candidate/surface dumps.
    pub pair: bool,
    /// the `clip` trace topic: clip input/output dumps.
    pub clip: bool,
}

impl TraceFlags {
    /// Read both trace gates from the environment.
    pub(crate) fn from_env() -> Self {
        Self {
            pair: popkern_core::trace::trace_enabled("pair-scan"),
            clip: popkern_core::trace::trace_enabled("clip"),
        }
    }

    /// True if any trace stream is active. Tracing forces sequential
    /// processing so interleaved thread output cannot garble the dump.
    /// (Only consulted by the rayon dispatch arm — the sequential-only
    /// build has no scheduling decision to make.)
    #[cfg_attr(not(feature = "parallel"), allow(dead_code))]
    pub(crate) fn any(&self) -> bool {
        self.pair || self.clip
    }
}

/// A surface-surface intersection failure on one face pair.
///
/// Carries everything the Stage 1 caller needs to build the same
/// [`BoolDiagnostics`] the sequential loop produced: the failing pair,
/// the solver's message, and a 3D location proxy (§13.1 — the midpoint
/// of face A's bounded UV range).
#[derive(Debug)]
pub(crate) struct PairFailure {
    /// Face from solid A.
    pub face_a: FaceHandle,
    /// Face from solid B.
    pub face_b: FaceHandle,
    /// The intersection solver's error message.
    pub message: String,
    /// 3D proxy for where the solver gave up.
    pub location: popkern_core::geometry::Point3,
}

impl PairFailure {
    /// Build the Stage 1 failure diagnostics, preserving the records
    /// accumulated from pairs processed before this failure (in
    /// candidate order), exactly as the sequential loop reported them.
    pub(crate) fn into_diagnostics(self, records: Vec<IntersectionRecord>) -> BoolDiagnostics {
        BoolDiagnostics {
            error: Some(BoolError::IntersectionFailed(format!(
                "face_a={:?} face_b={:?}: {}",
                self.face_a, self.face_b, self.message
            ))),
            intersection_curves: records.into_iter().map(|r| r.curve).collect(),
            face_classifications: Vec::new(),
            failed_stage: Some(Stage::Intersection),
            entries: vec![crate::error::DiagnosticEntry {
                severity: crate::error::Severity::Error,
                phase: Stage::Intersection,
                message: format!("Surface intersection failed: {}", self.message),
                entities: vec![crate::error::EntityRef::FacePair(self.face_a, self.face_b)],
                location: Some(self.location),
                measured_value: None,
                remediation: Some(crate::error::Remediation {
                    action: crate::error::RemediationAction::RetryWithPerturbation,
                    parameter: None,
                    value: None,
                    rationale: "Surface intersection solver failed. \
                        Try boolean_with_perturbation() to slightly \
                        offset the geometry and retry."
                        .into(),
                }),
            }],
            phase_metrics: Vec::new(),
        }
    }
}

/// Process candidate face pairs into intersection records.
///
/// With the `parallel` feature (default) pairs run on rayon's pool;
/// results are collected in input order, so the output is positionally
/// identical to sequential processing at any thread count. Tracing
/// forces the sequential path to keep dump ordering readable. Pairs
/// after a failing pair are still computed in parallel mode (their
/// results are discarded by the caller's first-failure fold), which
/// wastes bounded work on the error path but changes no output.
pub(crate) fn process_pairs(
    tasks: &[FacePairTask],
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
    trace: TraceFlags,
) -> Vec<Result<Vec<IntersectionRecord>, PairFailure>> {
    // Phase 1 (sequential, cheap): resolve each pair's surfaces and
    // bounded UV windows.
    let inputs: Vec<Option<PairInputs<'_>>> = tasks
        .iter()
        .map(|task| resolve_pair_inputs(task, a, b, trace))
        .collect();

    // SSI memo over EXACTLY-equal inputs (PK-304). Analytic
    // primitives arrive as N trimmed patches of ONE surface with
    // identical full-period windows — a 16-patch torus × cylinder
    // boolean re-marched the same complete intersection per pair
    // (measured 48 s of the 54 s repro). Exact structural equality is
    // the consumption gate (the PK-278 pattern): the memo only fires
    // when the solver inputs are bit-identical, so the records are
    // byte-for-byte what the per-pair solves produced.
    let mut unique_of_task: Vec<Option<usize>> = vec![None; tasks.len()];
    let mut unique: Vec<&PairInputs<'_>> = Vec::new();
    for (ti, inp) in inputs.iter().enumerate() {
        if let Some(inp) = inp {
            let idx = match unique.iter().position(|u| u.same_solver_inputs(inp)) {
                Some(existing) => existing,
                None => {
                    unique.push(inp);
                    unique.len() - 1
                }
            };
            unique_of_task[ti] = Some(idx);
        }
    }
    let solve = |inp: &PairInputs<'_>| {
        intersect_surfaces(
            inp.surf_a, inp.r1.0, inp.r1.1, inp.surf_b, inp.r2.0, inp.r2.1, tolerance,
        )
        .map_err(|e| e.to_string())
    };
    let ssi: Vec<
        Result<popkern_core::geometry::intersection::SurfaceSurfaceIntersectionResult, String>,
    > = {
        #[cfg(feature = "parallel")]
        {
            if !trace.any() {
                use rayon::prelude::*;
                unique.par_iter().map(|u| solve(u)).collect()
            } else {
                unique.iter().map(|u| solve(u)).collect()
            }
        }
        #[cfg(not(feature = "parallel"))]
        {
            unique.iter().map(|u| solve(u)).collect()
        }
    };

    // Per-task clip + classify from the shared solve, in task
    // order (rayon's ordered collect keeps the output positionally
    // identical to sequential processing at any thread count). Pairs
    // after a failing pair still compute in parallel mode (discarded
    // by the caller's first-failure fold).
    let one_task = |ti: usize| -> Result<Vec<IntersectionRecord>, PairFailure> {
        let (Some(ui), Some(inp)) = (unique_of_task[ti], inputs[ti].as_ref()) else {
            return Ok(Vec::new());
        };
        let task = &tasks[ti];
        match &ssi[ui] {
            Err(message) => {
                // §13.1: report a 3D location for the failure. The
                // midpoint of face_a's bounded UV range is the closest
                // cheap proxy for "where" the solver gave up.
                let u_mid = 0.5 * (inp.r1.0.0 + inp.r1.0.1);
                let v_mid = 0.5 * (inp.r1.1.0 + inp.r1.1.1);
                Err(PairFailure {
                    face_a: task.face_a,
                    face_b: task.face_b,
                    message: message.clone(),
                    location: inp.surf_a.point_at(u_mid, v_mid),
                })
            }
            Ok(result) => Ok(clip_and_classify(task, inp, result, a, b, tolerance, trace)),
        }
    };
    #[cfg(feature = "parallel")]
    {
        if !trace.any() {
            use rayon::prelude::*;
            return (0..tasks.len()).into_par_iter().map(one_task).collect();
        }
    }
    (0..tasks.len()).map(one_task).collect()
}

/// Resolved SSI solver inputs for one pair: surface references plus
/// bounded UV windows. [`Self::same_solver_inputs`] is the exact
/// structural equality that keys the SSI memo.
struct PairInputs<'x> {
    surf_a: &'x popkern_core::geometry::surface::SurfaceGeometry,
    r1: ((f64, f64), (f64, f64)),
    surf_b: &'x popkern_core::geometry::surface::SurfaceGeometry,
    r2: ((f64, f64), (f64, f64)),
}

impl PairInputs<'_> {
    /// Bit-exact equality of everything `intersect_surfaces` reads.
    fn same_solver_inputs(&self, other: &Self) -> bool {
        self.r1 == other.r1
            && self.r2 == other.r2
            && self.surf_a == other.surf_a
            && self.surf_b == other.surf_b
    }
}

/// Surface + window resolution for one pair. Faces without surface
/// geometry yield `None` (the sequential loop skipped them).
fn resolve_pair_inputs<'x>(
    task: &FacePairTask,
    a: &'x BRep,
    b: &'x BRep,
    trace: TraceFlags,
) -> Option<PairInputs<'x>> {
    let FacePairTask { face_a, face_b } = *task;
    if trace.pair {
        eprintln!("  [pair] {face_a:?} x {face_b:?}");
    }
    let surf_a = &a.faces[face_a].surface.as_ref()?.geometry;
    let surf_b = &b.faces[face_b].surface.as_ref()?.geometry;

    let r1 = bounded_parameter_range(a, face_a, surf_a);
    let r2 = bounded_parameter_range(b, face_b, surf_b);

    if trace.pair {
        trace_pair_surfaces(face_a, face_b, surf_a, surf_b, r1, r2);
    }
    Some(PairInputs {
        surf_a,
        r1,
        surf_b,
        r2,
    })
}

/// Clip one pair's shared SSI curves to the two faces and classify
/// tangency — the per-pair tail of the retired `process_one_pair`.
fn clip_and_classify(
    task: &FacePairTask,
    inp: &PairInputs<'_>,
    result: &popkern_core::geometry::intersection::SurfaceSurfaceIntersectionResult,
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
    trace: TraceFlags,
) -> Vec<IntersectionRecord> {
    let FacePairTask { face_a, face_b } = *task;
    let mut records = Vec::new();
    for isect in &result.intersections {
        match isect {
            SurfaceSurfaceIntersection::Curve(curve) => {
                let clipped_list = clip_curve_to_faces(a, face_a, b, face_b, curve);
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
                            let p = clipped.curve_3d.point_at(c0 + (c1 - c0) * k as f64 / 4.0);
                            for s in [inp.surf_a, inp.surf_b] {
                                let d =
                                    popkern_core::geometry::uv_projection::project_point_to_uv_multisheet(
                                        s, &p,
                                    )
                                    .map(|(u, v)| s.point_at(u, v).distance_to(&p))
                                    .unwrap_or(0.0);
                                worst = worst.max(d);
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
    super::ff_junction_refine::refine_record_junction_endpoints(
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
