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
    /// `PK_PAIR_TRACE`: candidate/surface dumps.
    pub pair: bool,
    /// `PK_CLIP_TRACE`: clip input/output dumps.
    pub clip: bool,
}

impl TraceFlags {
    /// Read both trace gates from the environment.
    pub(crate) fn from_env() -> Self {
        Self {
            pair: std::env::var_os("PK_PAIR_TRACE").is_some(),
            clip: std::env::var_os("PK_CLIP_TRACE").is_some(),
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
    #[cfg(feature = "parallel")]
    {
        if !trace.any() {
            use rayon::prelude::*;
            return tasks
                .par_iter()
                .map(|task| process_one_pair(task, a, b, tolerance, trace))
                .collect();
        }
    }
    tasks
        .iter()
        .map(|task| process_one_pair(task, a, b, tolerance, trace))
        .collect()
}

/// Intersect one face pair and classify each clipped curve.
///
/// Faces without surface geometry yield no records (the sequential loop
/// skipped them). Coincident regions are skipped here as well — the
/// coincidence engine records overlaps separately.
fn process_one_pair(
    task: &FacePairTask,
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
    trace: TraceFlags,
) -> Result<Vec<IntersectionRecord>, PairFailure> {
    let FacePairTask { face_a, face_b } = *task;
    if trace.pair {
        eprintln!("  [pair] {face_a:?} x {face_b:?}");
    }
    let surf_a = match a.faces[face_a].surface.as_ref() {
        Some(s) => &s.geometry,
        None => return Ok(Vec::new()),
    };
    let surf_b = match b.faces[face_b].surface.as_ref() {
        Some(s) => &s.geometry,
        None => return Ok(Vec::new()),
    };

    let ((u1_min, u1_max), (v1_min, v1_max)) = bounded_parameter_range(a, face_a, surf_a);
    let ((u2_min, u2_max), (v2_min, v2_max)) = bounded_parameter_range(b, face_b, surf_b);

    if trace.pair {
        trace_pair_surfaces(
            face_a,
            face_b,
            surf_a,
            surf_b,
            ((u1_min, u1_max), (v1_min, v1_max)),
            ((u2_min, u2_max), (v2_min, v2_max)),
        );
    }

    let result = match intersect_surfaces(
        surf_a,
        (u1_min, u1_max),
        (v1_min, v1_max),
        surf_b,
        (u2_min, u2_max),
        (v2_min, v2_max),
        tolerance,
    ) {
        Ok(result) => result,
        Err(e) => {
            // §13.1: report a 3D location for the failure. The midpoint
            // of face_a's bounded UV range is the closest cheap proxy
            // for "where" the intersection solver gave up.
            let u_mid = 0.5 * (u1_min + u1_max);
            let v_mid = 0.5 * (v1_min + v1_max);
            return Err(PairFailure {
                face_a,
                face_b,
                message: e.to_string(),
                location: surf_a.point_at(u_mid, v_mid),
            });
        }
    };

    let mut records = Vec::new();
    for isect in result.intersections {
        match isect {
            SurfaceSurfaceIntersection::Curve(curve) => {
                let clipped_list = clip_curve_to_faces(a, face_a, b, face_b, &curve);
                if trace.clip {
                    trace_clip_input(face_a, face_b, surf_a, surf_b, &curve, clipped_list.len());
                }
                for clipped in clipped_list {
                    if trace.clip {
                        trace_clip_output(face_a, face_b, &clipped);
                    }
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
                    let n1 = surf_a.normal_at(mid_uv1.x, mid_uv1.y);
                    let n2 = surf_b.normal_at(mid_uv2.x, mid_uv2.y);
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
    Ok(records)
}

fn trace_pair_surfaces(
    face_a: FaceHandle,
    face_b: FaceHandle,
    surf_a: &popkern_core::geometry::surface::SurfaceGeometry,
    surf_b: &popkern_core::geometry::surface::SurfaceGeometry,
    ((u1_min, u1_max), (v1_min, v1_max)): ((f64, f64), (f64, f64)),
    ((u2_min, u2_max), (v2_min, v2_max)): ((f64, f64), (f64, f64)),
) {
    let kind = |s: &popkern_core::geometry::surface::SurfaceGeometry| match s {
        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
        popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. } => "Cylinder",
        popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. } => "Sphere",
        popkern_core::geometry::surface::SurfaceGeometry::Cone { .. } => "Cone",
        popkern_core::geometry::surface::SurfaceGeometry::Torus { .. } => "Torus",
        _ => "Nurbs/other",
    };
    eprintln!(
        "  [pair-surf] {face_a:?}({}) u[{u1_min:.1},{u1_max:.1}] v[{v1_min:.1},{v1_max:.1}] x {face_b:?}({}) u[{u2_min:.1},{u2_max:.1}] v[{v2_min:.1},{v2_max:.1}]",
        kind(surf_a),
        kind(surf_b)
    );
    if let popkern_core::geometry::surface::SurfaceGeometry::Plane { u_axis, v_axis, .. } = surf_a {
        eprintln!(
            "    [pair-axes] A ua=({:.2},{:.2},{:.2}) va=({:.2},{:.2},{:.2})",
            u_axis.x, u_axis.y, u_axis.z, v_axis.x, v_axis.y, v_axis.z
        );
    }
}

fn trace_clip_input(
    face_a: FaceHandle,
    face_b: FaceHandle,
    surf_a: &popkern_core::geometry::surface::SurfaceGeometry,
    surf_b: &popkern_core::geometry::surface::SurfaceGeometry,
    curve: &popkern_core::geometry::curve_on_surface::CurveOnSurface,
    kept: usize,
) {
    let (t0, t1) = curve.curve_3d.parameter_range();
    let p0 = curve.curve_3d.point_at(t0);
    let p1 = curve.curve_3d.point_at(t1);
    let plane_dist = |s: &popkern_core::geometry::surface::SurfaceGeometry,
                      p: &popkern_core::geometry::Point3|
     -> f64 {
        match s {
            popkern_core::geometry::surface::SurfaceGeometry::Plane { origin, normal, .. } => {
                ((p.x - origin.x) * normal.x
                    + (p.y - origin.y) * normal.y
                    + (p.z - origin.z) * normal.z)
                    .abs()
            }
            _ => f64::NAN,
        }
    };
    eprintln!(
        "  [clip] A{face_a:?} x B{face_b:?}: t[{t0:.3},{t1:.3}] ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) dA={:.1e} dB={:.1e} kept={}",
        p0.x,
        p0.y,
        p0.z,
        p1.x,
        p1.y,
        p1.z,
        plane_dist(surf_a, &p0),
        plane_dist(surf_b, &p0),
        kept
    );
}

fn trace_clip_output(
    face_a: FaceHandle,
    face_b: FaceHandle,
    clipped: &popkern_core::geometry::curve_on_surface::CurveOnSurface,
) {
    let (t0, t1) = clipped.curve_3d.parameter_range();
    let p0 = clipped.curve_3d.point_at(t0);
    let p1 = clipped.curve_3d.point_at(t1);
    eprintln!(
        "  [clip-out] A{face_a:?} x B{face_b:?}: t[{t0:.3},{t1:.3}] ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) closed={}",
        p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, clipped.is_closed
    );
}
