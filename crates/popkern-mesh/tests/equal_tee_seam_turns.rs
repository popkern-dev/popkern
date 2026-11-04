//! PK-765: the equal-radius tees from the kernel API, the branch turned
//! about its own axis so its seam runs through each quarter, including
//! the singular point S+ where the two intersection ellipses cross.
//! With the seam through S+ the surface intersection traced a phantom
//! half-ellipse beyond the branch's end on the chart branch one period
//! below the face's, the clip's hull test skipped it as "incomplete",
//! and the union declined on its curveless edge. Every turn of the axis
//! tee (x0 = 0) and the x0 = -20 tee, union, subtract and intersect,
//! serves on its charts (no graph-merge fallback, no face walked
//! off-chart) with the analytic volume (main x branch overlap by
//! quadrature; the mesh volume is within its chord error). The x0 = -20
//! difference is pinched at S± and declines (PK-992; split serves it).

use popkern_bool::chart_walk_fallback::take_chart_walk_fallback_events;
use popkern_bool::merge_fallback::take_merge_fallback_events;
use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, ContactKind, NonManifoldPolicy, boolean,
    boolean_with_options,
};
use popkern_core::geometry::Vec3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const R: f64 = 40.0;
const H: f64 = 120.0;
const LEN: f64 = 100.0;

/// The main cylinder: radius R along +z, z in [0, H].
fn main_body() -> BRep {
    BRep::cylinder(R, H, Tolerance::default()).expect("main")
}

/// The branch: radius R along +x from x = x0, axis at z = 60, turned
/// `k` quarters about its own axis (its seam through S+ for k = 1).
fn branch(x0: f64, k: usize) -> BRep {
    let rod = BRep::cylinder(R, LEN, Tolerance::default()).expect("branch");
    let t = Transform3::translation(x0, 0.0, 60.0)
        .compose(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2))
        .compose(&Transform3::rotation_z(
            k as f64 * std::f64::consts::FRAC_PI_2,
        ));
    rod.transform(&t).expect("placed")
}

/// main ∩ branch for a branch over `x >= x0` (quadrature over the
/// branch's cross-section).
fn overlap(x0: f64) -> f64 {
    let n = 200_000;
    let mut sum = 0.0;
    for i in 0..n {
        let y = -R + 2.0 * R * (f64::from(i) + 0.5) / f64::from(n);
        let hc = (R * R - y * y).max(0.0).sqrt();
        let depth = (hc - x0.max(-hc)).max(0.0);
        sum += depth * 2.0 * (R * R - y * y).max(0.0).sqrt() * 2.0 * R / f64::from(n);
    }
    sum
}

fn mesh_volume(body: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(body, &popkern_mesh::TessellationOptions::default())
        .expect("tessellation");
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.vertices[i].position);
            let (a, b, c) = (
                Vec3::new(a.x, a.y, a.z),
                Vec3::new(b.x, b.y, b.z),
                Vec3::new(c.x, c.y, c.z),
            );
            a.dot(&b.cross(&c)) / 6.0
        })
        .sum()
}

/// The pinched difference: declined by default with point contacts at
/// S±, served by the split policy with both recorded and the analytic
/// volume.
fn pinched_difference(x0: f64, k: usize, want: f64) -> Result<(), String> {
    let tol = Tolerance::default();
    match boolean(
        BoolOperation::Subtraction,
        main_body(),
        branch(x0, k),
        tol,
        None,
    ) {
        Err(BoolError::NonManifoldResult { contacts }) => {
            let at_singular = |c: &popkern_bool::NonManifoldContact| {
                c.kind == ContactKind::Vertex
                    && c.location.x.abs() < 1e-3
                    && (c.location.y.abs() - R).abs() < 1e-3
                    && (c.location.z - 60.0).abs() < 1e-3
            };
            if contacts.len() != 2 || !contacts.iter().all(at_singular) {
                return Err(format!("declined with {contacts:?}, want S+ and S-"));
            }
        }
        other => return Err(format!("want NonManifoldResult, got {other:?}")),
    }
    let split = boolean_with_options(
        BoolOperation::Subtraction,
        main_body(),
        branch(x0, k),
        tol,
        None,
        &BoolOptions {
            non_manifold: NonManifoldPolicy::SplitToManifold,
        },
    )
    .map_err(|e| format!("split: {e}"))?;
    if split.topology_contacts.len() != 2 {
        return Err(format!(
            "split recorded {} contacts",
            split.topology_contacts.len()
        ));
    }
    let got = mesh_volume(&split);
    if (got - want).abs() > 2e-3 * want {
        return Err(format!("split volume {got} vs {want}"));
    }
    Ok(())
}

#[test]
fn equal_tees_serve_every_seam_turn_and_operation() {
    let mut failures = Vec::new();
    for x0 in [0.0, -20.0] {
        let shared = overlap(x0);
        let main = std::f64::consts::PI * R * R * H;
        for (op, want) in [
            (
                BoolOperation::Union,
                main + std::f64::consts::PI * R * R * LEN - shared,
            ),
            (BoolOperation::Subtraction, main - shared),
            (BoolOperation::Intersection, shared),
        ] {
            for k in 0..4 {
                let name = format!("x0 {x0} {op:?} seam turn {k}");
                let _ = take_merge_fallback_events();
                let _ = take_chart_walk_fallback_events();
                // A branch through the main cylinder (x0 < 0) leaves the
                // main body's material in two lobes, |z − 60| > |x|,
                // touching only at the singular points S± = (0, ±40,
                // 60): the difference is non-manifold there (PK-992).
                // It declines with both points, and the split policy
                // serves it with them recorded.
                let pinched = x0 < 0.0 && op == BoolOperation::Subtraction;
                if pinched {
                    if let Err(e) = pinched_difference(x0, k, want) {
                        failures.push(format!("{name}: {e}"));
                    }
                    continue;
                }
                let body = match boolean(op, main_body(), branch(x0, k), Tolerance::default(), None)
                {
                    Ok(b) => b,
                    Err(e) => {
                        failures.push(format!("{name}: {e}"));
                        continue;
                    }
                };
                if let Err(e) = body.validate() {
                    failures.push(format!("{name}: invalid {e}"));
                    continue;
                }
                let merges = take_merge_fallback_events();
                let walks = take_chart_walk_fallback_events();
                if !merges.is_empty() || !walks.is_empty() {
                    failures.push(format!("{name}: fallbacks {merges:?} {walks:?}"));
                }
                let got = mesh_volume(&body);
                if (got - want).abs() > 2e-3 * want {
                    failures.push(format!("{name}: volume {got} vs {want}"));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
