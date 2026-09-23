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

use popkern_core::geometry::intersection::intersect_surfaces_with_seeds;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::intersect::{IntersectionRecord, bounded_parameter_range};

#[path = "ff_parallel_clip.rs"]
mod ff_parallel_clip;
#[path = "ff_parallel_residue.rs"]
mod ff_parallel_residue;
#[path = "ff_parallel_trace.rs"]
mod ff_parallel_trace;
pub(crate) use ff_parallel_trace::TraceFlags;
use ff_parallel_trace::trace_pair_surfaces;

/// A face pair selected by BVH overlap filtering, awaiting intersection.
#[derive(Debug, Clone)]
pub(crate) struct FacePairTask {
    /// Face from solid A.
    pub face_a: FaceHandle,
    /// Face from solid B.
    pub face_b: FaceHandle,
    /// Edge-face pierce points adjacent to this pair: known-exact
    /// junction points ON the pair's intersection curves, fed to the
    /// SSI as march seeds (rescues narrow near-tangent lenses the
    /// interior seeding misses — PK-465).
    pub seeds: Vec<popkern_core::geometry::Point3>,
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
        intersect_surfaces_with_seeds(
            inp.surf_a, inp.r1.0, inp.r1.1, inp.surf_b, inp.r2.0, inp.r2.1, tolerance, inp.seeds,
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
            Ok(result) => Ok(ff_parallel_clip::clip_and_classify(
                task, inp, result, a, b, tolerance, trace,
            )),
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
    seeds: &'x [popkern_core::geometry::Point3],
}

impl PairInputs<'_> {
    /// Bit-exact equality of everything the seeded SSI reads. Seeds
    /// participate: two pairs sharing surfaces and windows but carrying
    /// different pierce seeds may march different curve sets.
    fn same_solver_inputs(&self, other: &Self) -> bool {
        self.r1 == other.r1
            && self.r2 == other.r2
            && self.surf_a == other.surf_a
            && self.surf_b == other.surf_b
            && self.seeds.len() == other.seeds.len()
            && self
                .seeds
                .iter()
                .zip(other.seeds)
                .all(|(p, q)| p.x == q.x && p.y == q.y && p.z == q.z)
    }
}

/// Surface + window resolution for one pair. Faces without surface
/// geometry yield `None` (the sequential loop skipped them).
fn resolve_pair_inputs<'x>(
    task: &'x FacePairTask,
    a: &'x BRep,
    b: &'x BRep,
    trace: TraceFlags,
) -> Option<PairInputs<'x>> {
    let FacePairTask { face_a, face_b, .. } = *task;
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
        seeds: &task.seeds,
    })
}
