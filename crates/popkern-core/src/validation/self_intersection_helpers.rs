use crate::geometry::Point3;
use crate::geometry::analytic_surfaces::analytic_signed_distance;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, FaceHandle};

#[allow(clippy::too_many_arguments)]
pub(super) fn cylinder_cylinder_near_samples(
    surface_a: &SurfaceGeometry,
    u_range_a: (f64, f64),
    v_range_a: (f64, f64),
    surface_b: &SurfaceGeometry,
    u_range_b: (f64, f64),
    v_range_b: (f64, f64),
    linear_tol: f64,
) -> Option<Vec<Point3>> {
    let SurfaceGeometry::Cylinder {
        origin: origin_a,
        axis: axis_a,
        radius: radius_a,
    } = surface_a
    else {
        return None;
    };
    let SurfaceGeometry::Cylinder {
        origin: origin_b,
        axis: axis_b,
        radius: radius_b,
    } = surface_b
    else {
        return None;
    };

    const GRID: usize = 17;
    let threshold = (linear_tol * 10.0).max(1.0e-8);
    let mut samples = Vec::new();
    // Validation containment for unsupported skew cylinder-cylinder pairs:
    // use bounded implicit-distance sampling plus trim classification instead
    // of the generic marcher, whose 4-D subdivision can be disproportionate
    // for small fillet cylinders that only meet at trimmed boundaries.
    sample_cylinder_against_cylinder(
        surface_a,
        u_range_a,
        v_range_a,
        CylinderDistanceTarget {
            origin: *origin_b,
            axis: *axis_b,
            radius: *radius_b,
        },
        threshold,
        &mut samples,
    );
    sample_cylinder_against_cylinder(
        surface_b,
        u_range_b,
        v_range_b,
        CylinderDistanceTarget {
            origin: *origin_a,
            axis: *axis_a,
            radius: *radius_a,
        },
        threshold,
        &mut samples,
    );

    if samples.is_empty() {
        for i in 0..GRID {
            let alpha = i as f64 / (GRID - 1) as f64;
            for j in 0..GRID {
                let beta = j as f64 / (GRID - 1) as f64;
                let u = u_range_a.0 + (u_range_a.1 - u_range_a.0) * alpha;
                let v = v_range_a.0 + (v_range_a.1 - v_range_a.0) * beta;
                let p = surface_a.point_at(u, v);
                let distance = distance_to_cylinder(p, *origin_b, *axis_b, *radius_b);
                if distance <= threshold {
                    samples.push(p);
                }
            }
        }
    }

    Some(samples)
}

#[derive(Clone, Copy)]
struct CylinderDistanceTarget {
    origin: Point3,
    axis: crate::geometry::Vec3,
    radius: f64,
}

fn sample_cylinder_against_cylinder(
    sampled_surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    other: CylinderDistanceTarget,
    threshold: f64,
    samples: &mut Vec<Point3>,
) {
    const GRID: usize = 17;
    for i in 0..GRID {
        let alpha = i as f64 / (GRID - 1) as f64;
        for j in 0..GRID {
            let beta = j as f64 / (GRID - 1) as f64;
            let u = u_range.0 + (u_range.1 - u_range.0) * alpha;
            let v = v_range.0 + (v_range.1 - v_range.0) * beta;
            let p = sampled_surface.point_at(u, v);
            let distance = distance_to_cylinder(p, other.origin, other.axis, other.radius);
            if distance <= threshold {
                samples.push(p);
            }
        }
    }
}

fn distance_to_cylinder(
    point: Point3,
    origin: Point3,
    axis: crate::geometry::Vec3,
    radius: f64,
) -> f64 {
    let offset = point - origin;
    let axial = offset.dot(&axis);
    let radial = offset - axis * axial;
    (radial.length() - radius).abs()
}

pub(super) fn trace_surface_kind(brep: &BRep, face: FaceHandle) -> &'static str {
    match brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Plane { .. }) => "Plane",
        Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder",
        Some(SurfaceGeometry::Cone { .. }) => "Cone",
        Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
        Some(SurfaceGeometry::Torus { .. }) => "Torus",
        Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
        Some(SurfaceGeometry::Gregory(_)) => "Gregory",
        Some(SurfaceGeometry::PolyhedralSpline(_)) => "PolyhedralSpline",
        None => "None",
    }
}

#[derive(Clone, Copy)]
pub(super) struct TraceSurfacePatch {
    pub(super) face: FaceHandle,
    pub(super) u_range: (f64, f64),
    pub(super) v_range: (f64, f64),
}

pub(super) fn trace_self_intersection_candidate(
    brep: &BRep,
    index: usize,
    patch_a: TraceSurfacePatch,
    patch_b: TraceSurfacePatch,
) {
    if std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_none() {
        return;
    }
    let face_a = patch_a.face;
    let face_b = patch_b.face;
    eprintln!(
        "[self-intersection candidate] #{index} {face_a:?}/{} u=({:.6},{:.6}) v=({:.6},{:.6}) \
         {face_b:?}/{} u=({:.6},{:.6}) v=({:.6},{:.6})",
        trace_surface_kind(brep, patch_a.face),
        patch_a.u_range.0,
        patch_a.u_range.1,
        patch_a.v_range.0,
        patch_a.v_range.1,
        trace_surface_kind(brep, patch_b.face),
        patch_b.u_range.0,
        patch_b.u_range.1,
        patch_b.v_range.0,
        patch_b.v_range.1,
    );
    if let Some(detail) = trace_cylinder_detail(brep, patch_a.face) {
        eprintln!(
            "[self-intersection candidate]   {:?} {detail}",
            patch_a.face
        );
    }
    if let Some(detail) = trace_cylinder_detail(brep, patch_b.face) {
        eprintln!(
            "[self-intersection candidate]   {:?} {detail}",
            patch_b.face
        );
    }
}

fn trace_cylinder_detail(brep: &BRep, face: FaceHandle) -> Option<String> {
    match &brep.faces.get(face)?.surface.as_ref()?.geometry {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some(format!(
            "Cylinder origin=({:.6},{:.6},{:.6}) axis=({:.6},{:.6},{:.6}) r={:.6}",
            origin.x, origin.y, origin.z, axis.x, axis.y, axis.z, radius
        )),
        _ => None,
    }
}

/// One near-crossing sample between an analytic surface and a
/// free-form face: the 3D point plus its UV on the FREE-FORM surface
/// (known exactly from the sampling grid — no iterative projection
/// involved).
pub(in crate::validation) struct AnalyticFreeFormSample {
    pub(in crate::validation) point: Point3,
    pub(in crate::validation) free_form_uv: (f64, f64),
}

/// PK-251 validation fast path for analytic × free-form pairs: locate
/// the crossing curve by grid-sampling the free-form surface's signed
/// implicit distance to the analytic surface and bisecting grid edges
/// with a sign change — the marching solver's adaptive step collapses
/// on grazing fillet-seam contact (measured: a single `validate_all`
/// spent 25 minutes marching pyramid3 cap/strip and cap/plane pairs),
/// while the implicit distance is O(1) per sample. Returns `None`
/// when the pair is not analytic × free-form (caller falls back to
/// the generic solver); the samples (possibly empty) otherwise, in
/// (analytic-operand, free-form-operand) classification order.
#[allow(clippy::too_many_arguments)]
pub(in crate::validation) fn analytic_free_form_near_samples(
    analytic: &SurfaceGeometry,
    free_form: &SurfaceGeometry,
    free_form_u_range: (f64, f64),
    free_form_v_range: (f64, f64),
    linear_tol: f64,
) -> Option<Vec<AnalyticFreeFormSample>> {
    analytic_signed_distance(analytic, &Point3::new(0.0, 0.0, 0.0))?;
    if !matches!(
        free_form,
        SurfaceGeometry::Nurbs(_) | SurfaceGeometry::Gregory(_)
    ) {
        return None;
    }
    // Out-resolve the free-form control nets this path serves
    // (fillet cap sectors run 38..67 control points per direction).
    const GRID: usize = 96;
    let threshold = (linear_tol * 10.0).max(1.0e-8);
    let signed_distance = |u: f64, v: f64| -> f64 {
        let p = free_form.point_at(u, v);
        analytic_signed_distance(analytic, &p).unwrap_or(f64::INFINITY)
    };
    let parameter = |i: usize, j: usize| -> (f64, f64) {
        (
            free_form_u_range.0
                + (free_form_u_range.1 - free_form_u_range.0) * i as f64 / (GRID - 1) as f64,
            free_form_v_range.0
                + (free_form_v_range.1 - free_form_v_range.0) * j as f64 / (GRID - 1) as f64,
        )
    };
    let mut grid = vec![0.0_f64; GRID * GRID];
    for i in 0..GRID {
        for j in 0..GRID {
            let (u, v) = parameter(i, j);
            grid[i * GRID + j] = signed_distance(u, v);
        }
    }
    let mut samples = Vec::new();
    let mut push_sample = |u: f64, v: f64| {
        samples.push(AnalyticFreeFormSample {
            point: free_form.point_at(u, v),
            free_form_uv: (u, v),
        });
    };
    let bisect = |mut a: (f64, f64), mut sa: f64, mut b: (f64, f64)| {
        for _ in 0..32 {
            let mid = ((a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5);
            let sm = signed_distance(mid.0, mid.1);
            if (sm <= 0.0) == (sa <= 0.0) {
                a = mid;
                sa = sm;
            } else {
                b = mid;
            }
        }
        ((a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5)
    };
    for i in 0..GRID {
        for j in 0..GRID {
            let sd = grid[i * GRID + j];
            let (u, v) = parameter(i, j);
            // Grazing contact without a sign change still counts.
            if sd.abs() <= threshold {
                push_sample(u, v);
                continue;
            }
            // Sign changes along grid edges localize the crossing
            // curve; bisection on the implicit distance lands on it.
            if i + 1 < GRID {
                let sd_right = grid[(i + 1) * GRID + j];
                if sd_right.abs() > threshold && (sd <= 0.0) != (sd_right <= 0.0) {
                    let (bu, bv) = bisect((u, v), sd, parameter(i + 1, j));
                    push_sample(bu, bv);
                }
            }
            if j + 1 < GRID {
                let sd_up = grid[i * GRID + j + 1];
                if sd_up.abs() > threshold && (sd <= 0.0) != (sd_up <= 0.0) {
                    let (bu, bv) = bisect((u, v), sd, parameter(i, j + 1));
                    push_sample(bu, bv);
                }
            }
        }
    }
    Some(samples)
}
