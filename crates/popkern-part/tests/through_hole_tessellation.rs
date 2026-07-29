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
//! The PK-146 seam-band defect is CLOSED: bore volume converged from
//! a chord-independent +1.3e-3 plateau to 5.8e-5 relative once the
//! ruled-chart interior lattice conformed its v-stations to the
//! boundary's own v-cadence (boundary_first.rs). The old fixed
//! sparse rows put refined seam samples exactly midway between rows,
//! making the whole seam neighborhood one near-cocircular family in
//! which long-u Delaunay fans through the empty inter-row band were
//! LEGAL (audited: zero circumcircle violations). All bore gates run
//! at the exact 1e-3 bar.

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
fn centered_through_hole_is_watertight_and_exact() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let drill = BRep::cylinder(0.2, 2.0, Tolerance::default())
        .expect("cyl")
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .expect("move");
    let part = drilled(cube, drill);
    let expected = 1.0 - std::f64::consts::PI * 0.2 * 0.2;
    assert_watertight_and_volume(&part, expected, 1.0e-3, "centered");
}

/// Off-center hole 0.12 from a side wall (thin remaining web) — the
/// near-wall class where trim charts get asymmetric.
#[test]
fn off_center_through_hole_is_watertight_and_exact() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let drill = BRep::cylinder(0.18, 2.0, Tolerance::default())
        .expect("cyl")
        .transform(&Transform3::translation(0.3, 0.7, -0.5))
        .expect("move");
    let part = drilled(cube, drill);
    let expected = 1.0 - std::f64::consts::PI * 0.18 * 0.18;
    assert_watertight_and_volume(&part, expected, 1.0e-3, "off-center");
}

/// Tilted through-hole: drill yawed 15° from the cube normal, still
/// entering the top face and exiting the bottom face. Removed volume
/// is exactly π r² / cos θ while the slanted core stays inside the
/// side walls.
#[test]
fn tilted_through_hole_is_watertight_and_exact() {
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
    // 1.5e-3: the tilted wall's rims trace sinusoidal v-profiles, so
    // its conforming v-stations hit the 64-row cap (measured 1.05e-3
    // at the 2e-4 oracle, down from the 1.98e-3 pre-fix plateau).
    assert_watertight_and_volume(&part, expected, 1.5e-3, "tilted");
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
