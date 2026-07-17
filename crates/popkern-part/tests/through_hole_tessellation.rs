//! PK-130: through-hole tessellation class gates.
//!
//! The historical gap (PK-100 era: interior hole walls never
//! tessellated — periodic wall ring u-winding ±2 periods,
//! boundary-first declined) closed en route through the chart-walk
//! and mesher campaigns; the old `through_hole_tessellation_gap_is_
//! pinned` fixture passes at HEAD. These gates widen the claim from
//! one fixture to the drilled-part class: axis-aligned, off-center,
//! tilted, coaxial tube, and cross-drilled cylinder — ALL watertight.
//!
//! Residual defect (measured 2026-07-17, open): cube-bore volume
//! plateaus ~1.3–2.0e-3 relative ABOVE analytic regardless of chord
//! tolerance. Localized: a seam-adjacent band of bore-wall triangles
//! (fine mesh: all 220 deep triangles at θ ∈ [−0.79, 0.79] around
//! the cylinder seam, chords up to 0.63 rad, r down to 0.191) while
//! boundary rings are exact (257/257 on-circle at 1e-5). Falsified:
//! fixed-N ring sampling, chord-lying T-points, cap membranes over
//! the bore mouth, and the seam-spine synthetic loop (bore walls are
//! not seam-spine faces). Suspect: boundary-first covering near the
//! seam columns (diametral-disk evacuation × interior grid). The
//! three bore gates below pin the measured plateau (3e-3) — tighten
//! to 1e-3 when the seam band refines. Tube and cross-drill are
//! strict.

use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_part::Part;

fn drilled(base: BRep, drill: BRep) -> Part {
    let out = popkern_bool::boolean(
        popkern_bool::BoolOperation::Subtraction,
        base,
        drill,
        Tolerance::default(),
        None,
    )
    .expect("subtract");
    Part::from_brep(out, Tolerance::default(), "drilled body").expect("validates (genus-aware)")
}

fn signed_mesh_volume(part: &Part) -> f64 {
    // Fine chord tolerance so cylinder faceting error stays well
    // under the volume gates (default 0.005 chord height costs
    // ~4e-3 relative on an r=0.2 bore — measured, not missing walls).
    let options = popkern_mesh::TessellationOptions {
        chord_height: 2.0e-4,
        ..popkern_mesh::TessellationOptions::default()
    };
    let mesh = part.mesh(&options).expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    v6 / 6.0
}

fn assert_watertight_and_volume(part: &Part, expected: f64, rel_tol: f64, label: &str) {
    let opens = part.open_edge_count().expect("mesh");
    assert_eq!(opens, 0, "[{label}] drilled body must mesh watertight");
    let got = signed_mesh_volume(part);
    assert!(
        ((got - expected) / expected).abs() < rel_tol,
        "[{label}] volume oracle: got {got:.6}, expected {expected:.6} (rel {:.2e})",
        (got - expected) / expected
    );
}

/// Axis-aligned centered through-hole (the original pin's fixture)
/// with the volume oracle the pin never had.
#[test]
fn centered_through_hole_is_watertight_volume_pinned() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let drill = BRep::cylinder(0.2, 2.0, Tolerance::default())
        .expect("cyl")
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .expect("move");
    let part = drilled(cube, drill);
    let expected = 1.0 - std::f64::consts::PI * 0.2 * 0.2;
    // 3e-3 bound = the measured seam-band plateau class (see module
    // header); a missing wall shows as tens of percent.
    assert_watertight_and_volume(&part, expected, 3.0e-3, "centered");
}

/// Off-center hole 0.12 from a side wall (thin remaining web) — the
/// near-wall class where trim charts get asymmetric.
#[test]
fn off_center_through_hole_is_watertight_volume_pinned() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let drill = BRep::cylinder(0.18, 2.0, Tolerance::default())
        .expect("cyl")
        .transform(&Transform3::translation(0.3, 0.7, -0.5))
        .expect("move");
    let part = drilled(cube, drill);
    let expected = 1.0 - std::f64::consts::PI * 0.18 * 0.18;
    assert_watertight_and_volume(&part, expected, 3.0e-3, "off-center");
}

/// Tilted through-hole: drill yawed 15° from the cube normal, still
/// entering the top face and exiting the bottom face. Removed volume
/// is exactly π r² / cos θ while the slanted core stays inside the
/// side walls.
#[test]
fn tilted_through_hole_is_watertight_volume_pinned() {
    let theta = 15.0_f64.to_radians();
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let drill = BRep::cylinder(0.15, 3.0, Tolerance::default())
        .expect("cyl")
        .transform(&Transform3::translation(0.0, 0.0, -1.5))
        .expect("center")
        .transform(&Transform3::rotation_x(theta))
        .expect("tilt")
        .transform(&Transform3::translation(0.5, 0.5, 0.5))
        .expect("place");
    let part = drilled(cube, drill);
    let expected = 1.0 - std::f64::consts::PI * 0.15 * 0.15 / theta.cos();
    assert_watertight_and_volume(&part, expected, 3.0e-3, "tilted");
}

/// Coaxial tube: cylinder minus centered cylinder — BOTH walls are
/// full-period periodic faces (the annular caps carry two closed
/// rings each).
#[test]
fn coaxial_tube_is_watertight_and_exact() {
    let outer = BRep::cylinder(0.5, 1.0, Tolerance::default()).expect("outer");
    let drill = BRep::cylinder(0.2, 2.0, Tolerance::default())
        .expect("cyl")
        .transform(&Transform3::translation(0.0, 0.0, -0.5))
        .expect("move");
    let part = drilled(outer, drill);
    let expected = std::f64::consts::PI * (0.5 * 0.5 - 0.2 * 0.2);
    assert_watertight_and_volume(&part, expected, 1.5e-3, "tube");
}

/// Cross-drilled cylinder: radial hole through a curved wall — the
/// cylinder×cylinder rim class. No closed-form volume at default
/// quality budget; gates are watertightness and genus-1 topology.
#[test]
fn cross_drilled_cylinder_is_watertight() {
    let body = BRep::cylinder(0.4, 1.0, Tolerance::default()).expect("body");
    let drill = BRep::cylinder(0.15, 2.0, Tolerance::default())
        .expect("drill")
        .transform(&Transform3::translation(0.0, 0.0, -1.0))
        .expect("center")
        .transform(&Transform3::rotation_y(90.0_f64.to_radians()))
        .expect("radial")
        .transform(&Transform3::translation(0.0, 0.0, 0.5))
        .expect("place");
    let part = drilled(body, drill);
    let opens = part.open_edge_count().expect("mesh");
    assert_eq!(opens, 0, "cross-drilled cylinder must mesh watertight");
}

/// Diagnostic probe: count bore-ring vertices at both chord options.
#[test]
#[ignore = "PK-130 probe — run with --ignored --nocapture"]
fn bore_ring_sampling_probe() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let drill = BRep::cylinder(0.2, 2.0, Tolerance::default())
        .expect("cyl")
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .expect("move");
    let part = drilled(cube, drill);
    for ch in [5.0e-3, 2.0e-4] {
        let options = popkern_mesh::TessellationOptions {
            chord_height: ch,
            ..popkern_mesh::TessellationOptions::default()
        };
        let mesh = part.mesh(&options).expect("mesh");
        let mut top_ring = 0usize;
        let mut _wall = 0usize;
        let mut on_circle = 0usize;
        for v in &mesh.vertices {
            let r = ((v.position.x - 0.5).powi(2) + (v.position.y - 0.5).powi(2)).sqrt();
            if (r - 0.2).abs() < 5.0e-3 {
                _wall += 1;
                if (v.position.z - 1.0).abs() < 1.0e-6 {
                    top_ring += 1;
                    if (r - 0.2).abs() < 1.0e-5 {
                        on_circle += 1;
                    }
                }
            }
        }
        let mut v6 = 0.0;
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]].position,
                mesh.vertices[t[1]].position,
                mesh.vertices[t[2]].position,
            );
            v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x);
        }
        let mut cap_over_hole = 0usize;
        let mut cap_over_hole_area = 0.0;
        let mut wall_low = 0usize;
        let mut inside_stats: Vec<(f64, f64, f64, f64)> = Vec::new();
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]].position,
                mesh.vertices[t[1]].position,
                mesh.vertices[t[2]].position,
            );
            let cx = (a.x + b.x + c.x) / 3.0 - 0.5;
            let cy = (a.y + b.y + c.y) / 3.0 - 0.5;
            let cz = (a.z + b.z + c.z) / 3.0;
            let cr = (cx * cx + cy * cy).sqrt();
            let flat = (a.z - c.z).abs() < 1.0e-9 && (b.z - c.z).abs() < 1.0e-9;
            if flat && (cz < 1.0e-9 || (cz - 1.0).abs() < 1.0e-9) && cr < 0.195 {
                cap_over_hole += 1;
                let ux = b.x - a.x;
                let uy = b.y - a.y;
                let vx = c.x - a.x;
                let vy = c.y - a.y;
                cap_over_hole_area += 0.5 * (ux * vy - uy * vx).abs();
            }
            if !flat && cr < 0.199 && cz > 0.01 && cz < 0.99 {
                wall_low += 1;
                let th = cy.atan2(cx);
                let dth_ab = ((b.y - 0.5).atan2(b.x - 0.5) - (a.y - 0.5).atan2(a.x - 0.5)).abs();
                inside_stats.push((th, cz, cr, dth_ab.min(std::f64::consts::TAU - dth_ab)));
            }
        }
        eprintln!(
            "chord_height={ch:.0e}: verts={} top-ring={top_ring} ON-CIRCLE={on_circle} vol={:.6} exp={:.6} cap-tris-over-hole={cap_over_hole} (area {cap_over_hole_area:.6}) wall-tris-inside={wall_low}",
            mesh.vertices.len(),
            v6 / 6.0,
            1.0 - std::f64::consts::PI * 0.04
        );
        let mut theta_bins = [0usize; 8];
        for (th, _, _, _) in &inside_stats {
            let b = (((th + std::f64::consts::PI) / std::f64::consts::TAU) * 8.0) as usize;
            theta_bins[b.min(7)] += 1;
        }
        let mut worst = inside_stats.clone();
        worst.sort_by(|x, y| x.2.partial_cmp(&y.2).expect("finite"));
        eprintln!("  theta bins: {theta_bins:?}");
        for w in worst.iter().take(4) {
            eprintln!(
                "  deep tri: theta={:.2} z={:.3} r={:.4} dtheta={:.3}",
                w.0, w.1, w.2, w.3
            );
        }
    }
}
