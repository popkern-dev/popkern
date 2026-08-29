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

    // SEAM-SAMPLED continuous stripe phase (PK-133). Two prior
    // constructions were falsified with measurement:
    // - the binary field paired at 0.04 had a random-agreement FLOOR
    //   (~0.5) on curved caps — interior fairing and strip-threshold
    //   tightening moved nothing;
    // - continuous phase paired at 0.01 still decorrelates: the
    //   stripe coordinate's surface gradient is ~2κ/spacing (≈ 270
    //   per unit at these fillet curvatures), so even 0.01 of
    //   positional pairing drift spans multiple stripes (measured
    //   0.88–0.97 while the seam normal residual is ≤ 1e-2 rad).
    // The artifact-free measure evaluates BOTH surfaces' normals at
    // the SAME boundary point: sample the shared edge curves and
    // compare the continuous stripe coordinates there. Any mismatch
    // is then a genuine reflected-phase jump across the seam.
    let _ = (&report_a, &report_b, boundary_tol);
    let shared_points = shared_edge_sample_points(brep, face_a, face_b, grid_density.max(8));
    if shared_points.is_empty() {
        return None;
    }
    let spacing = if stripe_spacing.abs() < 1e-30 {
        1.0
    } else {
        stripe_spacing
    };
    let phase_on = |face: FaceHandle, point: &Point3, seed: Option<(f64, f64)>| -> Option<f64> {
        use popkern_core::geometry::surface::SurfaceEval;
        use popkern_core::geometry::uv_projection::{
            project_point_to_uv, project_point_to_uv_near,
        };
        let geom = &brep.faces.get(face)?.surface.as_ref()?.geometry;
        // Global inversion on the cap NURBS lands up to 4e-2 off the
        // seam station (measured) — a pcurve-derived seed +
        // Gauss-Newton refinement is the honest station.
        let (u, v) = match seed {
            Some((su, sv)) => project_point_to_uv_near(geom, point, su, sv)?,
            None => project_point_to_uv(geom, point)?,
        };
        if popkern_core::trace::trace_enabled("zebra") {
            let back = geom.point_at(u, v);
            eprintln!(
                "  [zebra] face {face:?} uv ({u:.4},{v:.4}) inversion residual {:.3e}",
                back.distance_to(point)
            );
        }
        let n = geom.normal_at(u, v).normalized()?;
        let reflected = n * (2.0 * n.dot(&view_dir)) - view_dir;
        Some(reflected.dot(&stripe_dir) / spacing)
    };
    let mut matched = 0usize;
    let mut mismatched = 0usize;
    for sample in &shared_points {
        let (Some(pa), Some(pb)) = (
            phase_on(face_a, &sample.point, sample.seed_a),
            phase_on(face_b, &sample.point, sample.seed_b),
        ) else {
            continue;
        };
        matched += 1;
        if (pa - pb).abs() > 0.5 {
            mismatched += 1;
        }
    }
    if matched == 0 {
        return None;
    }
    Some(mismatched as f64 / matched as f64)
}

/// One seam station: the 3D point plus per-side UV SEEDS from the
/// half-edges' pcurves (refined by the caller — the pcurve station
/// is parameter-speed-approximate but an excellent Newton seed).
struct SharedEdgeSample {
    point: Point3,
    seed_a: Option<(f64, f64)>,
    seed_b: Option<(f64, f64)>,
}

/// Sample stations along every edge shared by the two faces.
fn shared_edge_sample_points(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    samples_per_edge: usize,
) -> Vec<SharedEdgeSample> {
    use popkern_core::geometry::curve::CurveEval;
    let mut out = Vec::new();
    for (_, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let fa = brep.loops.get(he.loop_).map(|l| l.face);
        let fb = brep.loops.get(twin.loop_).map(|l| l.face);
        let (side_a, side_b) = if fa == Some(face_a) && fb == Some(face_b) {
            (edge.half_edge, he.twin)
        } else if fa == Some(face_b) && fb == Some(face_a) {
            (he.twin, edge.half_edge)
        } else {
            continue;
        };
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        // Pcurves are loop-aligned (start at the half-edge origin);
        // map the shared 3D fraction onto each side's own traversal
        // direction of the edge curve.
        let start_3d = curve.geometry.point_at(curve.t_start);
        let seed_at = |heh: popkern_core::topology::HalfEdgeHandle, f: f64| -> Option<(f64, f64)> {
            let data = brep.half_edges.get(heh)?;
            let pc = data.pcurve.as_ref()?;
            let origin = brep.vertices.get(data.origin)?.position;
            let forward = origin.distance_to(&start_3d) <= 1.0e-9
                || origin.distance_to(&start_3d)
                    <= origin.distance_to(&curve.geometry.point_at(curve.t_end));
            let (p0, p1) = pc.parameter_range();
            let fp = if forward { f } else { 1.0 - f };
            let q = pc.point_at(p0 + (p1 - p0) * fp);
            Some((q.x, q.y))
        };
        for i in 0..samples_per_edge {
            // Interior stations only: endpoint normals belong to the
            // corner-junction class, measured by the corner gates.
            let f = (i as f64 + 0.5) / samples_per_edge as f64;
            let t = curve.t_start + (curve.t_end - curve.t_start) * f;
            out.push(SharedEdgeSample {
                point: curve.geometry.point_at(t),
                seed_a: seed_at(side_a, f),
                seed_b: seed_at(side_b, f),
            });
        }
    }
    out
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
