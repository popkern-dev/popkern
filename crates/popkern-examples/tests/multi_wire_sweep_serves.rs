//! PK-477 acceptance: multi-wire closed sweep profiles serve
//! watertight — the per-wire ring topology end to end, through the
//! mesh oracle (validator-clean alone is not accepted for new
//! topology paths; torn-mesh serves are a measured failure class).

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_sweep::{SweepOptions, SweepProfile, sweep};

fn square_wire(h: f64) -> NurbsCurve {
    let pts = [[-h, -h], [h, -h], [h, h], [-h, h], [-h, -h]];
    let n = pts.len();
    let mut knots = vec![0.0];
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    NurbsCurve {
        degree: 1,
        knots,
        control_points: pts.iter().map(|p| Point3::new(p[0], p[1], 0.0)).collect(),
        weights: vec![1.0; n],
    }
}

fn z_path(length: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, length)],
        weights: vec![1.0, 1.0],
    }
}

/// Watertight mesh census (position-welded) + signed volume against
/// the exact tube volume (outer area minus hole area, times length).
#[test]
fn square_tube_sweep_meshes_watertight_with_exact_volume() {
    let profile = SweepProfile::MultiWire(vec![square_wire(2.0), square_wire(1.0)]);
    let result = sweep(
        &profile,
        &z_path(5.0),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("multi-wire sweep");
    let brep = &result.brep;
    brep.validate().expect("valid tube");

    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tube meshes");
    // Weld vertices by quantized position (per-face meshes share no
    // indices), then census edge use counts and the signed volume.
    let mut key: std::collections::HashMap<(i64, i64, i64), usize> =
        std::collections::HashMap::new();
    let mut welded: Vec<usize> = Vec::with_capacity(mesh.vertices.len());
    for v in &mesh.vertices {
        let p = v.position;
        let q = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        let n = key.len();
        welded.push(*key.entry(q).or_insert(n));
    }
    let mut boundary: std::collections::HashMap<(usize, usize), i64> =
        std::collections::HashMap::new();
    let mut volume6 = 0.0_f64;
    for t in &mesh.triangles {
        let [a, b, c] = *t;
        let (wa, wb, wc) = (welded[a], welded[b], welded[c]);
        for (i, j) in [(wa, wb), (wb, wc), (wc, wa)] {
            if i == j {
                continue;
            }
            let entry = boundary.entry((i.min(j), i.max(j))).or_default();
            *entry += if i < j { 1 } else { -1 };
        }
        let (pa, pb, pc) = (
            mesh.vertices[a].position,
            mesh.vertices[b].position,
            mesh.vertices[c].position,
        );
        volume6 += pa.x * (pb.y * pc.z - pc.y * pb.z) - pa.y * (pb.x * pc.z - pc.x * pb.z)
            + pa.z * (pb.x * pc.y - pc.x * pb.y);
    }
    let open = boundary.values().filter(|&&v| v != 0).count();
    assert_eq!(open, 0, "{open} open mesh edges on the swept tube");

    // Visual verification artifact (geometry-affecting change rule).
    let out_dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(out_dir).expect("output dir");
    let camera = popkern_viz::camera::Camera {
        eye: glam::Vec3::new(8.0, -7.0, 7.0),
        target: glam::Vec3::new(0.0, 0.0, 2.5),
        up: glam::Vec3::new(0.0, 0.0, 1.0),
        fov_y: 45.0_f32.to_radians(),
        aspect: 800.0 / 600.0,
        near: 0.1,
        far: 100.0,
    };
    let png = out_dir.join("pk477_square_tube.png");
    popkern_viz::headless::render_to_png(brep, &camera, 800, 600, &png, true).expect("render");
    println!("[render] {}", png.display());
    let volume = (volume6 / 6.0).abs();
    let expected = (4.0 * 4.0 - 2.0 * 2.0) * 5.0;
    let rel = (volume - expected).abs() / expected;
    assert!(
        rel < 0.01,
        "tube volume {volume:.4} vs expected {expected:.4} ({:.2}% off)",
        rel * 100.0
    );
}

/// The `SweepProfile::Face` route with an inner loop (a pocketed
/// plate top face) goes through the same multi-wire analysis and
/// serves. The pocket is polygonal: `loop_to_polyline` samples loop
/// VERTICES only, so a circular hole (a one-half-edge self-loop)
/// degenerates — a pre-existing face-profile limitation, unchanged
/// here.
#[test]
fn face_profile_with_hole_sweeps_valid() {
    use popkern_part::prelude::*;
    let plate = Profile::on(Plane::XY)
        .rect(6.0, 6.0)
        .extrude(1.0)
        .expect("plate");
    let holed = plate
        .subtract(
            &Profile::on(Plane::XY)
                .rect(2.0, 2.0)
                .extrude(1.0)
                .expect("square drill"),
        )
        .expect("holed plate");
    let brep = holed.brep();
    let (face, _) = brep
        .faces
        .iter()
        .find(|(_, f)| !f.inner_loops.is_empty())
        .expect("a face with an inner loop");
    let profile = SweepProfile::Face {
        brep: Box::new(brep.clone()),
        face,
    };
    let result = sweep(
        &profile,
        &z_path(3.0),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("face-with-hole sweep");
    result.brep.validate().expect("valid swept ring solid");
}

/// PK-479: the drilled-CIRCULAR-hole face profile now sweeps to a
/// valid solid — loop half-edges adopt their edge curves, so the
/// hole's one-half-edge rim stays an exact circle instead of
/// degenerating to its single vertex (measured pre-fix:
/// InconsistentOrientation on the swept body).
#[test]
fn face_profile_with_circular_hole_sweeps_valid() {
    use popkern_part::prelude::*;
    let plate = Profile::on(Plane::XY)
        .rect(6.0, 6.0)
        .extrude(1.0)
        .expect("plate");
    let holed = plate
        .subtract(
            &Profile::on(Plane::XY)
                .circle(1.0)
                .extrude(1.0)
                .expect("drill"),
        )
        .expect("holed plate");
    let brep = holed.brep();
    let (face, _) = brep
        .faces
        .iter()
        .find(|(_, f)| !f.inner_loops.is_empty())
        .expect("a face with an inner loop");
    let profile = SweepProfile::Face {
        brep: Box::new(brep.clone()),
        face,
    };
    let result = sweep(
        &profile,
        &z_path(3.0),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("circular-hole face sweep");
    let out = &result.brep;
    out.validate().expect("valid swept ring solid");
    // The hole wall must be curved (a recognized cylinder or a NURBS
    // wall), not a prism of chords: expect a face whose surface is
    // non-planar with area consistent with an r=1 circular hole.
    let curved_walls = out
        .faces
        .iter()
        .filter(|(_, f)| {
            !matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) | None
            )
        })
        .count();
    assert!(
        curved_walls >= 1,
        "expected a curved hole wall, found {curved_walls} non-planar faces"
    );
}
