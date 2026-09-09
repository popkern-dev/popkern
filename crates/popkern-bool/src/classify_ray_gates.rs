//! Boundary-derived geometric hit gates for the classify ray
//! (split from `classify_ray.rs` under the 500-line file budget).
//!
//! Each gate judges one candidate ray hit against the face's
//! boundary-sample geometry in the surface's own frame — never the
//! solver window or stored pcurves, so parameter conventions and
//! seam branches cannot bias the verdict.

use popkern_core::{
    geometry::{Point3, Vec3, surface::SurfaceGeometry},
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

/// Verdict of a per-hit geometric gate.
pub(super) enum HitGate {
    /// The hit survives this gate.
    Accept,
    /// The hit is a phantom for this face — skip it, the ray stays valid.
    SkipHit,
    /// The hit is ambiguous — the whole ray is inconclusive.
    PoisonRay,
}

/// CYLINDER hits: accept only within the face's boundary-derived
/// angular span (+tol). The padded parameter window is solver
/// headroom, and rays passing BEHIND a partial arc wall "hit" the
/// phantom half of the cylinder (measured: an arm-socket tool's
/// discarded fragments sampled near a 180° arc front classified
/// Inside and the assembly went non-manifold). Full closed rims have
/// degenerate spans and keep the padded behavior; FINAL bodies only
/// (the caller guards on `skip_hole_hits`).
///
/// Parameter conventions differ between the solver window and the
/// stored chart, so the gate is GEOMETRIC: the hit's 3D point (from
/// the ray, not the u/v) must fall within the face's
/// boundary-derived angular span in the surface's own basis. The
/// span comes from CURVE SAMPLES, not vertices alone: a face
/// spanning exactly half the cylinder has all its vertices ON the
/// atan2 ±π cut, where the covered half is decided by the sign of a
/// zero (measured: a boolean emitted −0.0 in a vertex coordinate and
/// the span flipped to the phantom half, rejecting every interior
/// hit — a prow-rear block read Outside at an interior point and a
/// seat-cutter wall was dropped). The span is the minimal circular
/// arc covering the samples (largest-gap), so cut placement cannot
/// bias it; a full rim's gap is one sample spacing, which yields
/// span ≈ τ and keeps the padded behavior.
pub(super) fn cylinder_angular_gate(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    origin: &Point3,
    direction: &Vec3,
    t: f64,
) -> HitGate {
    let SurfaceGeometry::Cylinder {
        origin: cyl_origin,
        axis: cyl_axis,
        ..
    } = surface
    else {
        return HitGate::Accept;
    };
    let samples = brep.face_boundary_sample_points(face_handle, 8);
    if samples.len() < 3 {
        return HitGate::Accept;
    }
    let (perp1, perp2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(*cyl_axis);
    let ang = |p: &Point3| -> f64 {
        let d = *p - *cyl_origin;
        d.dot(&perp2).atan2(d.dot(&perp1))
    };
    let mut angs: Vec<f64> = samples.iter().map(&ang).collect();
    angs.sort_by(|a, b| a.total_cmp(b));
    let tau = std::f64::consts::TAU;
    let mut max_gap = 0.0;
    let mut span_start = angs[0];
    for i in 0..angs.len() {
        let a = angs[i];
        let b = if i + 1 < angs.len() {
            angs[i + 1]
        } else {
            angs[0] + tau
        };
        if b - a > max_gap {
            max_gap = b - a;
            span_start = b;
        }
    }
    let span = tau - max_gap;
    if span > 1.0e-6 && span <= std::f64::consts::PI + 1.0e-6 {
        let hit3d = *origin + *direction * t;
        let rel = (ang(&hit3d) - span_start).rem_euclid(tau);
        let margin = 1.0e-6;
        if rel > span + margin && rel < tau - margin {
            if popkern_core::trace::trace_enabled("ray-classify") {
                eprintln!(
                    "  [ray-rej-cyl] {face_handle:?} rel={rel:.3} span={span:.3} from {span_start:.3}"
                );
            }
            return HitGate::SkipHit;
        }
    }
    HitGate::Accept
}

/// SPHERE-CAP gate (PK-345): a partial sphere face (a filleted-cube
/// corner octant) accepts only hits whose direction-from-center
/// falls inside the boundary-derived cap — the padded parameter box
/// covers the WHOLE sphere, and rays from a point just outside a
/// corner hit the phantom 7/8 (measured: 4 of 8 identical corner
/// probes classified Inside). Geometric like the angular gate: the
/// cap is min over boundary-sample dots with the mean direction; a
/// full sphere's seam boundary yields a non-positive min and keeps
/// the padded behavior. FINAL bodies only (caller-guarded).
pub(super) fn sphere_cap_gate(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    origin: &Point3,
    direction: &Vec3,
    t: f64,
    tolerance: &Tolerance,
) -> HitGate {
    let SurfaceGeometry::Sphere {
        center: sph_center, ..
    } = surface
    else {
        return HitGate::Accept;
    };
    let samples = brep.face_boundary_sample_points(face_handle, 8);
    if samples.len() < 3 {
        return HitGate::Accept;
    }
    let dirs: Vec<Vec3> = samples
        .iter()
        .filter_map(|p| {
            let d = *p - *sph_center;
            let len = d.length();
            (len > tolerance.linear()).then(|| d * (1.0 / len))
        })
        .collect();
    if dirs.len() < 3 {
        return HitGate::Accept;
    }
    let mut mean = Vec3::new(0.0, 0.0, 0.0);
    for d in &dirs {
        mean = mean + *d;
    }
    let mean_len = mean.length();
    if mean_len <= 1.0e-9 {
        return HitGate::Accept;
    }
    let mean = mean * (1.0 / mean_len);
    let min_dot = dirs
        .iter()
        .map(|d| d.dot(&mean))
        .fold(f64::INFINITY, f64::min);
    if min_dot > 0.05 {
        let hit3d = *origin + *direction * t;
        let hd = hit3d - *sph_center;
        let hd_len = hd.length();
        if hd_len > tolerance.linear() && (hd * (1.0 / hd_len)).dot(&mean) < min_dot - 0.05 {
            if popkern_core::trace::trace_enabled("ray-classify") {
                eprintln!(
                    "  [ray-rej-sph] {face_handle:?} dot={:.3} min={min_dot:.3}",
                    (hd * (1.0 / hd_len)).dot(&mean)
                );
            }
            return HitGate::SkipHit;
        }
    }
    HitGate::Accept
}

/// AXIAL gate for cylinder and cone hits: the hit's coordinate along
/// the surface axis must fall within the face's boundary-derived
/// axial span (+one relative margin). The padded parameter box is
/// solver headroom, and near a chamfered rim a near-horizontal ray
/// legitimately crosses the cap and cone, then ALSO "hits" the
/// wall's infinite cylinder BELOW the weld (measured: rays from
/// under a rim-chamfered column counted cap + cone + phantom wall at
/// v=0.29 against a wall trimmed to [1,59] — every below-body tool
/// fragment classified Inside and the subtract kept a skirt; the
/// mirror case accepted an extended-cone hit at z=2.3 on a cone
/// banded [0,1]). Geometric, like the angular gate: boundary curve
/// samples, not pcurves, so seam branches cannot bias it. FINAL
/// bodies only (caller-guarded).
pub(super) fn axial_span_gate(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    origin: &Point3,
    direction: &Vec3,
    t: f64,
    tolerance: &Tolerance,
) -> HitGate {
    let axis_frame = match surface {
        SurfaceGeometry::Cylinder {
            origin: o, axis, ..
        } => Some((*o, *axis)),
        SurfaceGeometry::Cone { apex, axis, .. } => Some((*apex, *axis)),
        _ => None,
    };
    let Some((axis_origin, axis_dir)) = axis_frame else {
        return HitGate::Accept;
    };
    let samples = brep.face_boundary_sample_points(face_handle, 8);
    if samples.len() < 3 {
        return HitGate::Accept;
    }
    let axial = |p: &Point3| -> f64 { (*p - axis_origin).dot(&axis_dir) };
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for sp in &samples {
        let a = axial(sp);
        lo = lo.min(a);
        hi = hi.max(a);
    }
    let hit3d = *origin + *direction * t;
    let a = axial(&hit3d);
    // Margin: 5% of the span (+tolerance). The coincident-overlap
    // machinery probes at ~1e-3 of flush interfaces that sit exactly
    // ON a face's axial end — a tight margin rejected those
    // knife-edge hits and re-tuned flush verdicts (measured: the
    // flush gusset union lost a cooriented OnBoundary and kept a
    // duplicate skin). The rim phantoms this gate exists for sit
    // 20–130% outside their spans.
    let margin = tolerance.linear().max(0.05 * (hi - lo));
    if a < lo - margin || a > hi + margin {
        if popkern_core::trace::trace_enabled("ray-classify") {
            eprintln!("  [ray-rej-axial] {face_handle:?} a={a:.3} span=[{lo:.3},{hi:.3}]");
        }
        return HitGate::SkipHit;
    }
    // A hit OUTSIDE the span but inside the margin is on the
    // surface's extension right at the face's end — it can neither
    // be counted (an extended-cone phantom 0.045 past a chamfer
    // ring's top edge broke the parity and kept a hex wall fragment
    // inside the tool, PK-460) nor silently dropped (the flush-probe
    // knife hits the margin exists for). The ray is inconclusive;
    // another direction decides.
    let slop = tolerance.linear() * 10.0;
    if a < lo - slop || a > hi + slop {
        if popkern_core::trace::trace_enabled("ray-classify") {
            eprintln!(
                "  [ray-band-axial] {face_handle:?} a={a:.3} span=[{lo:.3},{hi:.3}] — ray inconclusive"
            );
        }
        return HitGate::PoisonRay;
    }
    HitGate::Accept
}

/// PLANE hits: accept only within the face's boundary-derived
/// in-plane bounding box (+10% margin). The padded parameter window
/// is solver headroom; a ray passing BEYOND a bounded plane face
/// "hits" its infinite plane (measured, PK-341: rays from inside a
/// tube counted the annular cap at r = 1.5 against an outer rim of
/// 0.7 — parity flipped, every interior sample read Outside, both
/// groove-wall fragments were dropped, and the subtraction served a
/// non-closed body). The box is a strict SUPERSET of the face trim
/// (boundary samples), so rejection can never lose a real hit;
/// geometric in the plane's own frame, so pcurve/window conventions
/// cannot bias it. Planes are non-periodic — no seam-branch hazard —
/// so this gate runs in BOTH modes, unlike the cylinder and axial
/// gates.
pub(super) fn plane_box_gate(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    origin: &Point3,
    direction: &Vec3,
    t: f64,
    tolerance: &Tolerance,
) -> HitGate {
    let SurfaceGeometry::Plane {
        origin: pl_origin,
        normal: pl_normal,
        ..
    } = surface
    else {
        return HitGate::Accept;
    };
    let samples = brep.face_boundary_sample_points(face_handle, 8);
    if samples.len() < 3 {
        return HitGate::Accept;
    }
    let (e1, e2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(*pl_normal);
    let coord = |p: &Point3| -> (f64, f64) {
        let d = *p - *pl_origin;
        (d.dot(&e1), d.dot(&e2))
    };
    let (mut lo1, mut hi1) = (f64::MAX, f64::MIN);
    let (mut lo2, mut hi2) = (f64::MAX, f64::MIN);
    for sp in &samples {
        let (a, b) = coord(sp);
        lo1 = lo1.min(a);
        hi1 = hi1.max(a);
        lo2 = lo2.min(b);
        hi2 = hi2.max(b);
    }
    let hit3d = *origin + *direction * t;
    let (a, b) = coord(&hit3d);
    // 10% covers the chordal underestimate of circular rims sampled
    // 8-per-edge (sagitta ≈ 7.6% of the radius); the phantom class
    // sits ≥ 100% outside.
    let m1 = tolerance.linear().max(0.1 * (hi1 - lo1));
    let m2 = tolerance.linear().max(0.1 * (hi2 - lo2));
    if a < lo1 - m1 || a > hi1 + m1 || b < lo2 - m2 || b > hi2 + m2 {
        if popkern_core::trace::trace_enabled("ray-classify") {
            eprintln!(
                "  [ray-rej-plane] {face_handle:?} ({a:.3},{b:.3}) box=[{lo1:.3},{hi1:.3}]x[{lo2:.3},{hi2:.3}]"
            );
        }
        return HitGate::SkipHit;
    }
    HitGate::Accept
}
