//! §5.3 item #11 — zebra-stripe regression diagnostics for blend
//! faces.
//!
//! Zebra analysis shades the surface with alternating black/white
//! bands based on a parallel-stripe reflection function; G2
//! discontinuities manifest as **zig-zag** kinks in the stripe
//! boundaries (where they cross patch edges). Companion Class-A
//! tool to reflection lines.
//!
//! * [`zebra_stripes_on_face`] — compute the zebra stripe field
//!   on a single face (thin wrapper over
//!   [`popkern_core::healing::compute_zebra_stripes`]).
//! * [`zebra_boundary_kink`] — sample across the shared edge
//!   between two faces and measure the maximum offset between
//!   their zebra-stripe transitions. A G2-continuous boundary
//!   produces `≈ 0`; breaks produce non-zero offsets
//!   proportional to the surface kink angle.
//! * [`blend_zebra_kink_max`] — the one-call Class-A metric for
//!   a blend face against all its support faces.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::healing::{ZebraStripeField, compute_zebra_stripes};
use popkern_core::topology::{BRep, FaceHandle};

use crate::geometry_compute::face_param_range;

/// Result of running zebra-stripe analysis on a single face.
#[derive(Debug, Clone)]
pub struct ZebraStripeReport {
    /// The face that was analysed.
    pub face: FaceHandle,
    /// Grid density used for the sampling.
    pub grid_density: usize,
    /// Raw zebra stripe field (grid of 0.0/1.0 values + UV per cell).
    pub field: ZebraStripeField,
}

/// Compute the zebra-stripe field on a single face of the BRep.
pub fn zebra_stripes_on_face(
    brep: &BRep,
    face: FaceHandle,
    view_dir: Vec3,
    stripe_dir: Vec3,
    stripe_spacing: f64,
    grid_density: usize,
) -> Option<ZebraStripeReport> {
    let face_ref = brep.faces.get(face)?;
    let surface_on_face = face_ref.surface.as_ref()?;
    let (u_range, v_range) = face_param_range(brep, face, &surface_on_face.geometry);
    let field = compute_zebra_stripes(
        &surface_on_face.geometry,
        u_range,
        v_range,
        view_dir,
        stripe_dir,
        stripe_spacing,
        grid_density,
    );
    Some(ZebraStripeReport {
        face,
        grid_density,
        field,
    })
}

/// §5.3 item #11 — zebra-boundary kink between two adjacent faces.
///
/// Samples each face on a `grid_density × grid_density` grid of
/// surface points, computes the zebra-stripe value at each, and
/// for pairs whose 3D positions are within `boundary_tol` (the
/// "boundary samples") compares the stripe values on either
/// side. For a G2-continuous boundary the stripes transition
/// smoothly across the edge — paired samples should have
/// matching values. For a G0/G1 break, some paired samples
/// disagree.
///
/// Returns the **mismatch ratio** — the fraction of boundary
/// sample pairs whose zebra values disagree (scaled into
/// `[0, 1]`). A UI can surface this as `0 % = Class A G2`,
/// `50 % = random half-stripes = full G0 break`.
///
/// Returns `None` if no boundary pairs are identified (faces
/// aren't actually adjacent or their 3D sample grids don't
/// overlap within tolerance).
#[allow(clippy::too_many_arguments)]
pub fn zebra_boundary_kink(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    view_dir: Vec3,
    stripe_dir: Vec3,
    stripe_spacing: f64,
    grid_density: usize,
    boundary_tol: f64,
) -> Option<f64> {
    let report_a = zebra_stripes_on_face(
        brep,
        face_a,
        view_dir,
        stripe_dir,
        stripe_spacing,
        grid_density,
    )?;
    let report_b = zebra_stripes_on_face(
        brep,
        face_b,
        view_dir,
        stripe_dir,
        stripe_spacing,
        grid_density,
    )?;

    let samples_a = sample_field_in_3d(brep, face_a, &report_a);
    let samples_b = sample_field_in_3d(brep, face_b, &report_b);
    if samples_a.is_empty() || samples_b.is_empty() {
        return None;
    }

    let tol_sq = boundary_tol * boundary_tol;
    let mut matched = 0usize;
    let mut mismatched = 0usize;
    for (pt_a, val_a) in &samples_a {
        for (pt_b, val_b) in &samples_b {
            let dx = pt_a.x - pt_b.x;
            let dy = pt_a.y - pt_b.y;
            let dz = pt_a.z - pt_b.z;
            if dx * dx + dy * dy + dz * dz > tol_sq {
                continue;
            }
            matched += 1;
            if (val_a - val_b).abs() > 0.5 {
                mismatched += 1;
            }
        }
    }
    if matched == 0 {
        return None;
    }
    Some(mismatched as f64 / matched as f64)
}

/// One-call Class-A zebra metric for a blend face: max
/// mismatch-ratio across any blend-vs-**non-Gregory**-support
/// boundary. Gregory-vs-Gregory spoke boundaries are excluded
/// (item #3's scope).
#[allow(clippy::too_many_arguments)]
pub fn blend_zebra_kink_max(
    brep: &BRep,
    blend_face: FaceHandle,
    view_dir: Vec3,
    stripe_dir: Vec3,
    stripe_spacing: f64,
    grid_density: usize,
    boundary_tol: f64,
) -> Option<f64> {
    let supports = crate::reflection_line::adjacent_support_faces(brep, blend_face)?;
    let mut worst: Option<f64> = None;
    for support in supports {
        if brep
            .faces
            .get(support)
            .and_then(|f| f.surface.as_ref())
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Gregory(_)))
        {
            continue;
        }
        if let Some(ratio) = zebra_boundary_kink(
            brep,
            blend_face,
            support,
            view_dir,
            stripe_dir,
            stripe_spacing,
            grid_density,
            boundary_tol,
        ) {
            worst = Some(match worst {
                Some(m) => m.max(ratio),
                None => ratio,
            });
        }
    }
    worst
}

/// Evaluate the face's surface at each of `report.field.params`
/// to get a 3D position, paired with the zebra stripe value.
fn sample_field_in_3d(
    brep: &BRep,
    face: FaceHandle,
    report: &ZebraStripeReport,
) -> Vec<(Point3, f64)> {
    use popkern_core::geometry::surface::SurfaceEval;
    let Some(face_ref) = brep.faces.get(face) else {
        return Vec::new();
    };
    let Some(surface_on_face) = face_ref.surface.as_ref() else {
        return Vec::new();
    };
    let geom: &SurfaceGeometry = &surface_on_face.geometry;
    let mut out = Vec::with_capacity(
        report.field.params.len() * report.field.params.first().map_or(0, |r| r.len()),
    );
    for (param_row, val_row) in report.field.params.iter().zip(report.field.grid.iter()) {
        for (&(u, v), &val) in param_row.iter().zip(val_row.iter()) {
            out.push((geom.point_at(u, v), val));
        }
    }
    out
}
