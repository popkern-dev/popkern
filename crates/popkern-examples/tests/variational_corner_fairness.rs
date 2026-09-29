//! PK-574 row 879: fairness of the variational corner fit against the
//! KP20 cap on the same corners.
//!
//! Each corner is built twice: as dispatched (KP20 serves, except where
//! it declines and the fit is the served construction) and with the KP20
//! cap forced to decline (the fit serves). Per build, over the cap faces:
//! - certified curvature enclosures (PK-510 interval bounds; a face the
//!   enclosure cannot bound is counted, never guessed);
//! - the sampled largest principal curvature times the smallest fillet
//!   radius (the overshoot relative to the strips' 1/r);
//! - zebra seam readings on the surface normals (cap↔neighbour and cap↔cap
//!   classes: 95th-percentile normal and slope jumps), with a PNG per
//!   build for the report.
//!
//! Gate (design doc E3): the fit's curvature overshoot is not worse than
//! KP20's by more than 25 % on any corner both serve; past that the
//! Moreton–Séquin polish becomes planned work.

use std::collections::HashSet;

use popkern_blend::dispatch_overrides::with_kp20_cap_forced_to_decline;
use popkern_blend::{
    BlendResult, ContinuityLevel, CornerMethod, FilletOptions, FilletRadius, FilletSpec,
    fillet_edges,
};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::healing::interval_curvature::guaranteed_curvature_bounds;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};
use popkern_examples::zebra_attribution::{SeamClass, SeamSummary, classify_crossings, summarize};
use popkern_measure::{CurvatureKind, curvature_at};
use popkern_viz::camera::Camera;
use popkern_viz::zebra::{ZebraNormals, rasterize_zebra_frame};

struct Corner {
    name: &'static str,
    body: BRep,
    vertex: VertexHandle,
    radii: Vec<f64>,
}

fn vertex_at(body: &BRep, at: Point3) -> VertexHandle {
    body.vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&at) < 1e-6)
        .map(|(h, _)| h)
        .expect("vertex")
}

fn apex(body: &BRep) -> VertexHandle {
    body.vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex")
}

fn corners() -> Vec<Corner> {
    let tol = Tolerance::default();
    let pyramid = |n, h| BRep::n_gonal_pyramid(n, 1.0, h, tol).expect("pyramid");
    let cube = BRep::unit_cube(tol).expect("cube");
    let notch = Transform3::translation(0.6, 0.6, 0.6).compose(&Transform3::scale_uniform(0.6));
    let tool = BRep::unit_cube(tol)
        .expect("cube")
        .transform(&notch)
        .expect("notch");
    let notched = popkern_bool::subtraction(BRep::unit_cube(tol).expect("cube"), tool, tol, None)
        .expect("notched cube");
    let mut out = Vec::new();
    for (name, n, h, radii) in [
        ("pyramid3 h1.2", 3, 1.2, vec![0.16]),
        ("pentagon flagship", 5, 1.2, vec![0.16]),
        ("pentagon5 blunt", 5, 3.0, vec![0.05]),
        (
            "pentagon5 variable",
            5,
            0.8,
            vec![0.02, 0.04, 0.06, 0.04, 0.02],
        ),
    ] {
        let body = pyramid(n, h);
        let vertex = apex(&body);
        out.push(Corner {
            name,
            body,
            vertex,
            radii,
        });
    }
    for (name, radii) in [
        ("cube corner 0.2", vec![0.2]),
        ("cube corner 0.05/0.2/0.3", vec![0.05, 0.2, 0.3]),
    ] {
        out.push(Corner {
            name,
            body: cube.clone(),
            vertex: vertex_at(&cube, Point3::new(1.0, 1.0, 1.0)),
            radii,
        });
    }
    out.push(Corner {
        name: "notched saddle",
        vertex: vertex_at(&notched, Point3::new(0.6, 0.6, 0.6)),
        body: notched,
        radii: vec![0.08],
    });
    out
}

fn fillet(corner: &Corner) -> BlendResult {
    let mut edges: Vec<_> = corner
        .body
        .half_edges
        .values()
        .filter(|he| he.origin == corner.vertex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .into_iter()
        .zip(corner.radii.iter().cycle())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    fillet_edges(&corner.body, &specs, &options, Tolerance::default(), None).expect("fillet")
}

/// Fairness readings of one build's cap.
#[derive(Debug, Clone)]
struct Reading {
    construction: String,
    level: ContinuityLevel,
    faces: usize,
    /// Certified `max |H|` over the cap, `None` when a face could not
    /// be enclosed.
    mean_bound: Option<f64>,
    /// Faces whose enclosure has degenerate cells.
    unbounded_faces: usize,
    /// Sampled `max |κ|` × smallest radius.
    overshoot: f64,
    cap_strip: SeamSummary,
    cap_cap: SeamSummary,
}

fn finite_range(r: (f64, f64)) -> (f64, f64) {
    if r.0.is_finite() && r.1.is_finite() {
        r
    } else {
        (0.0, 1.0)
    }
}

/// A camera on the cap: target its centroid, look along its mean
/// outward normal (tilted a little so no seam runs along the view),
/// framed by its bounding diagonal.
fn cap_camera(brep: &BRep, caps: &[FaceHandle]) -> Camera {
    let (mut sum, mut normal, mut count) =
        (Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0), 0.0);
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for &face in caps {
        let Some(s) = brep.faces[face].surface.as_ref() else {
            continue;
        };
        let ((u0, u1), (v0, v1)) = {
            let (u, v) = s.geometry.parameter_range();
            (finite_range(u), finite_range(v))
        };
        for i in 0..=4 {
            for j in 0..=4 {
                let (u, v) = (
                    u0 + (u1 - u0) * i as f64 / 4.0,
                    v0 + (v1 - v0) * j as f64 / 4.0,
                );
                let p = s.geometry.point_at(u, v);
                let n = s.geometry.normal_at(u, v);
                sum = sum + Vec3::new(p.x, p.y, p.z);
                normal = normal + if s.orientation { n } else { n * -1.0 };
                count += 1.0;
                for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                    lo[k] = lo[k].min(c);
                    hi[k] = hi[k].max(c);
                }
            }
        }
    }
    let centre = sum * (1.0 / count);
    // Outward: away from the body's vertex centroid (face orientation
    // conventions differ between cap constructions).
    let body_centre = {
        let (mut acc, mut k) = (Vec3::new(0.0, 0.0, 0.0), 0.0);
        for (_, v) in brep.vertices.iter() {
            acc = acc + Vec3::new(v.position.x, v.position.y, v.position.z);
            k += 1.0;
        }
        acc * (1.0 / k)
    };
    let mut n = normal.normalized().expect("cap normal");
    if n.dot(&(centre - body_centre)) < 0.0 {
        n = n * -1.0;
    }
    let up = if n.z.abs() > 0.9 {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let tilt = up.cross(&n).normalized().expect("tilt");
    let diag = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
    let eye = centre + (n + tilt * 0.3) * (3.5 * diag);
    let f = |v: Vec3| glam::Vec3::new(v.x as f32, v.y as f32, v.z as f32);
    Camera {
        eye: f(eye),
        target: f(centre),
        up: f(up),
        fov_y: 0.5,
        aspect: 1.0,
        near: 0.001,
        far: 100.0,
    }
}

fn read(result: &BlendResult, radii: &[f64], camera: &Camera, png: &std::path::Path) -> Reading {
    let brep = &result.brep;
    let caps: Vec<FaceHandle> = result.corner_patches.iter().map(|p| p.face).collect();
    let construction = match result.corner_patches.first().map(|p| p.method) {
        Some(CornerMethod::PlateMin) => "variational fit".to_string(),
        Some(m) => format!("{m:?}"),
        None => "none".to_string(),
    };
    let level = result
        .corner_patches
        .iter()
        .map(|p| p.achieved_continuity)
        .min_by_key(|l| *l as u8)
        .unwrap_or(ContinuityLevel::G0);
    let r_min = radii.iter().copied().fold(f64::INFINITY, f64::min);
    let (mut mean_bound, mut unbounded_faces, mut kmax) = (Some(0.0f64), 0usize, 0.0f64);
    for &face in &caps {
        let Some(s) = brep.faces[face].surface.as_ref() else {
            continue;
        };
        let (u, v) = s.geometry.parameter_range();
        let (u, v) = (finite_range(u), finite_range(v));
        // PK-510 cannot enclose a Gregory patch's corner blend (0/0), so
        // a Gregory face is recorded unbounded without refining.
        let enclosable = matches!(
            s.geometry,
            popkern_core::geometry::SurfaceGeometry::Nurbs(_)
        );
        let bounds = enclosable.then(|| guaranteed_curvature_bounds(&s.geometry, u, v, 1.0, 6));
        if bounds.is_none_or(|b| b.degenerate_cells > 0) {
            unbounded_faces += 1;
            mean_bound = None;
        } else if let (Some(m), Some(b)) = (mean_bound.as_mut(), bounds) {
            *m = m.max(b.mean.lo.abs().max(b.mean.hi.abs()));
        }
        // Interior samples: a Gregory patch's rational blend is 0/0 at
        // its parameter corners, which read spurious curvature.
        for i in 1..32 {
            for j in 1..32 {
                let (uu, vv) = (
                    u.0 + (u.1 - u.0) * i as f64 / 32.0,
                    v.0 + (v.1 - v.0) * j as f64 / 32.0,
                );
                if let Some(k) = curvature_at(&s.geometry, uu, vv, CurvatureKind::MaxPrincipalAbs) {
                    kmax = kmax.max(k);
                }
            }
        }
    }
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::display())
        .expect("mesh");
    let frame = rasterize_zebra_frame(brep, &mesh, camera, 480, 480);
    popkern_viz::headless::write_png(png, 480, 480, &frame.shade(ZebraNormals::Surface, 12))
        .expect("png");
    let cap_set: HashSet<FaceHandle> = caps.iter().copied().collect();
    let is_plane = |f: FaceHandle| {
        matches!(
            brep.faces
                .get(f)
                .and_then(|x| x.surface.as_ref())
                .map(|s| &s.geometry),
            Some(popkern_core::geometry::SurfaceGeometry::Plane { .. })
        )
    };
    let role = |f: FaceHandle| {
        if cap_set.contains(&f) {
            SeamClass::CAP
        } else if is_plane(f) {
            SeamClass::SUPPORT
        } else {
            SeamClass::STRIP
        }
    };
    let classes = classify_crossings(&frame.seam_crossings(ZebraNormals::Surface), role);
    let summary = |label: &str| {
        classes
            .iter()
            .find(|(c, _)| c.label() == label)
            .map(|(_, rows)| summarize(rows))
            .unwrap_or_default()
    };
    Reading {
        construction,
        level,
        faces: caps.len(),
        mean_bound,
        unbounded_faces,
        overshoot: kmax * r_min,
        cap_strip: {
            let (rim, support) = (summary("cap-rim"), summary("cap-support"));
            SeamSummary {
                rows: rim.rows + support.rows,
                seams: rim.seams + support.seams,
                value_p50: rim.value_p50.max(support.value_p50),
                value_p95: rim.value_p95.max(support.value_p95),
                normal_p50: rim.normal_p50.max(support.normal_p50),
                normal_p95: rim.normal_p95.max(support.normal_p95),
                slope_p50: rim.slope_p50.max(support.slope_p50),
                slope_p95: rim.slope_p95.max(support.slope_p95),
            }
        },
        cap_cap: summary("cap-cap"),
    }
}

fn row(name: &str, r: &Reading) -> String {
    format!(
        "| {name} | {} | {:?} | {} | {} | {:.2} | {:.1e} / {:.1e} | {:.1e} / {:.1e} |",
        r.construction,
        r.level,
        r.faces,
        match r.mean_bound {
            Some(m) => format!("{m:.1}"),
            None => format!("unbounded ({} faces)", r.unbounded_faces),
        },
        r.overshoot,
        r.cap_strip.normal_p95,
        r.cap_strip.slope_p95,
        r.cap_cap.normal_p95,
        r.cap_cap.slope_p95,
    )
}

#[test]
fn row_879_fit_fairness_against_kp20() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("variational_corner_fairness");
    std::fs::create_dir_all(&dir).expect("dir");
    let mut table = vec![
        "| corner | construction | level | faces | certified max abs H | max abs k * r_min | cap-rim normal / slope p95 | cap-cap normal / slope p95 |".to_string(),
        "|---|---|---|---|---|---|---|---|".to_string(),
    ];
    let mut violations = Vec::new();
    let mut kp20: Vec<(&'static str, f64)> = Vec::new();
    for corner in corners() {
        let dispatched = fillet(&corner);
        let fitted = with_kp20_cap_forced_to_decline(|| fillet(&corner));
        let caps: Vec<FaceHandle> = dispatched.corner_patches.iter().map(|p| p.face).collect();
        let camera = cap_camera(&dispatched.brep, &caps);
        let stem = corner.name.replace([' ', '/'], "_");
        let a = read(
            &dispatched,
            &corner.radii,
            &camera,
            &dir.join(format!("{stem}_dispatched.png")),
        );
        let b = read(
            &fitted,
            &corner.radii,
            &camera,
            &dir.join(format!("{stem}_fit.png")),
        );
        table.push(row(corner.name, &a));
        table.push(row(corner.name, &b));
        if a.construction != "variational fit" {
            kp20.push((corner.name, a.overshoot));
        }
        // PK-605: the tangent (G1) fit re-freezes its bending reference
        // three times; its max k*r_min is pinned where it serves at G1
        // (before: blunt 5.10, variable 5.73, flagship 1.48,
        // pyramid3 1.17).
        let fit_pin = [
            ("pentagon5 blunt", 1.76),
            ("pentagon5 variable", 2.20),
            ("pentagon flagship", 1.40),
            ("pyramid3 h1.2", 1.22),
        ];
        if let Some((_, pin)) = fit_pin.iter().find(|(n, _)| *n == corner.name) {
            assert!(
                b.overshoot <= *pin,
                "{}: fit k*r {:.2} over its pin {pin:.2}",
                corner.name,
                b.overshoot
            );
        }
        if a.construction != "variational fit" && b.overshoot > 1.25 * a.overshoot {
            violations.push(format!(
                "{}: fit {:.2} vs {} {:.2}",
                corner.name, b.overshoot, a.construction, a.overshoot
            ));
        }
    }
    eprintln!("{}\nPNGs: {}", table.join("\n"), dir.display());
    // Max k*r_min of the dispatched KP20 cap, pinned per corner. PK-599
    // (bi-9 resampling), PK-600 (guide third-order fair interior) and
    // PK-601 (curvature-variation fair interior, the fairer of the two
    // served): cube and saddle 6.77 -> 2.39 -> 1.76, pyramid3 2.78 ->
    // 1.74 -> 1.23, flagship 2.07 -> 1.84 -> 1.84, blunt 6.43 -> 5.19 ->
    // 1.70, variable 4.72 -> 4.71 -> 3.75.
    for (name, bound) in [
        ("cube corner 0.2", 1.82),
        ("notched saddle", 1.82),
        ("pyramid3 h1.2", 1.28),
        ("pentagon flagship", 1.90),
        ("pentagon5 blunt", 1.76),
        ("pentagon5 variable", 3.85),
    ] {
        let reading = kp20.iter().find(|(n, _)| *n == name).map(|(_, r)| *r);
        let reading = reading.unwrap_or_else(|| panic!("{name}: KP20 did not serve"));
        assert!(
            reading <= bound,
            "{name}: KP20 cap k*r {reading:.2} over {bound:.2}"
        );
    }
    assert!(
        violations.is_empty(),
        "E3 gate: fit overshoot > 1.25x KP20: {violations:?}"
    );
}
