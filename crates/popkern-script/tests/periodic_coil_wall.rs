//! A coil spring's tube wall is written in periodic form (PK-892): the
//! wire circle closes smoothly, so the wall is a ring in u like a
//! cylinder's lateral face, not a clamped chart holding its seam twice.
//! Its volume is the tube's exact one, pi r^2 times the helix length
//! (the wire radius is far under the helix's curvature radius), and it
//! meshes watertight.

use popkern_core::geometry::SurfaceGeometry;
use popkern_script::Plane;
use popkern_script::body_sweep::HelixParams;
use popkern_script::prelude::*;

#[test]
fn a_coil_springs_wall_is_periodic_and_its_volume_exact() {
    let (big_r, pitch, wire_r, turns) = (0.3_f64, 0.15_f64, 0.04_f64, 6.5_f64);
    let mut m = Model::new();
    let z_axis = AxisSpec::Literal {
        origin: [0.0, 0.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    };
    let spring = m
        .profile(Plane::XY)
        .circle(wire_r)
        .sketch("wire_section")
        .expect("wire sketch")
        .helical_sweep(&mut m, z_axis, HelixParams::new(big_r, pitch, turns))
        .expect("coil");
    let brep = m.brep_of(spring).expect("brep").clone();
    let walls: Vec<_> = brep
        .faces
        .values()
        .filter_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Nurbs(n)) => Some(n),
            _ => None,
        })
        .collect();
    assert_eq!(walls.len(), 1, "one tube wall");
    assert_eq!(walls[0].u_period(), Some(1.0), "the wall is a ring in u");
    let length = turns * ((std::f64::consts::TAU * big_r).powi(2) + pitch * pitch).sqrt();
    let exact = std::f64::consts::PI * wire_r * wire_r * length;
    let measured = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
    assert!(
        ((measured - exact) / exact).abs() < 1e-4,
        "volume {measured} vs exact {exact}"
    );
    let mesh = popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(popkern_mesh::welded_open_edges(&mesh), 0, "watertight mesh");
}
