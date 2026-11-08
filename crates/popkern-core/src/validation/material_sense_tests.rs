//! Material-sense resolution pins.

use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn unit_cube_material_sense_is_determinate_and_flips_with_body_orientation() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let seed = brep.faces.keys().next().unwrap();
    let sense = brep
        .component_material_sense(seed)
        .expect("cube sense determinate");

    let mut reversed = brep.clone();
    reversed.reverse_body_orientation();
    let seed_r = reversed.faces.keys().next().unwrap();
    let sense_r = reversed
        .component_material_sense(seed_r)
        .expect("reversed cube sense determinate");
    assert_ne!(
        sense, sense_r,
        "reversing the body orientation must flip the material sense"
    );
}

/// The sense does not depend on where the body sits: a thin bore far
/// from the origin reads as it does at the origin (the coarse grid's
/// residual used to be multiplied by the placement offset).
#[test]
fn material_sense_is_placement_invariant_for_a_thin_bore() {
    use crate::geometry::transform::Transform3;
    let tol = Tolerance::default();
    let bore = BRep::cylinder(0.04, 0.28, tol).unwrap();
    let seed = bore.faces.keys().next().unwrap();
    let home = bore.component_material_sense(seed).expect("determinate");
    let mut flipped = Vec::new();
    for i in -8..=8 {
        for j in -8..=8 {
            let (x, y) = (0.61 * f64::from(i) / 4.0, 0.31 * f64::from(j) / 4.0);
            let placed = bore
                .transform(&Transform3::translation(x, y, -0.1))
                .unwrap();
            let seed = placed.faces.keys().next().unwrap();
            if placed.component_material_sense(seed) != Some(home) {
                flipped.push((x, y));
            }
        }
    }
    assert!(flipped.is_empty(), "sense changed at {flipped:?}");
}

/// The sense is determinate for a body far larger than the tolerance
/// cube at a coarse tolerance: the floor is relative to the integral,
/// not (100 · linear)³ (which read this cone, volume 0.34, as
/// indeterminate at linear = 1e-2).
#[test]
fn material_sense_is_determinate_at_a_coarse_tolerance() {
    let tol = Tolerance::default().with_linear(1e-2);
    for body in [
        BRep::cone(0.4, 2.0, tol).unwrap(),
        BRep::cylinder(0.15, 0.5, tol).unwrap(),
        BRep::unit_cube(tol).unwrap(),
    ] {
        let seed = body.faces.keys().next().unwrap();
        let fine = {
            let mut b = body.clone();
            b.tolerance = Tolerance::default();
            b.component_material_sense(seed)
        };
        assert!(fine.is_some());
        assert_eq!(body.component_material_sense(seed), fine);
    }
}

/// A seam pcurve that stops √ε short of a pole, as a pole vertex
/// projected through acos does (v = π − 2.6e-8 on a 0.072 sphere, a
/// row 3e-9 across), still closes the seam pair across the pole: the
/// face is the whole sphere, not a zero-area slit (which classified
/// every point Outside and left the sense unknown).
#[test]
fn a_seam_short_of_the_pole_still_bounds_the_whole_sphere() {
    use crate::geometry::pcurve::Pcurve;
    let short = std::f64::consts::PI - 2.6e-8;
    let mut sphere = BRep::sphere(0.072, Tolerance::default()).unwrap();
    let hes: Vec<_> = sphere.half_edges.keys().collect();
    for h in hes {
        let he = &sphere.half_edges[h];
        let Some(p) = he.pcurve.as_ref() else {
            continue;
        };
        let (a, b) = (p.point_at(p.t_start), p.point_at(p.t_end));
        let pull = |v: f64| if v > short { short } else { v };
        let (t0, t1) = (p.t_start, p.t_end);
        sphere.half_edges[h].pcurve =
            Some(Pcurve::linear((a.x, pull(a.y)), (b.x, pull(b.y)), t0, t1));
    }
    let face = sphere.faces.keys().next().unwrap();
    let surface = sphere.faces[face]
        .surface
        .as_ref()
        .unwrap()
        .geometry
        .clone();
    assert_eq!(
        sphere.face_uv_trim_containment(face, &surface, 1.0, 1.0),
        Some(crate::validation::FaceUvContainment::Interior)
    );
    assert!(sphere.component_material_sense(face).is_some());
}

/// A seam-only face whose seam projections scatter by ulps (a window
/// 2.2e-15 wide, past an ε·period cutoff) is still the full wrap: the
/// sense is determinate wherever the sphere sits (measured: two of a
/// pipped die's 21 pip spheres read unknown at their placements).
#[test]
fn a_seam_only_sphere_is_determinate_at_every_placement() {
    use crate::geometry::transform::Transform3;
    let sphere = BRep::sphere(0.0021, Tolerance::default().with_linear(1e-9)).unwrap();
    let mut unknown = Vec::new();
    for i in -6..=6 {
        for j in -6..=6 {
            let (x, y) = (0.0044 * f64::from(i), 0.0095 * f64::from(j) / 3.0);
            let placed = sphere
                .transform(&Transform3::translation(x, y, -0.0023))
                .unwrap();
            let seed = placed.faces.keys().next().unwrap();
            if placed.component_material_sense(seed).is_none() {
                unknown.push((x, y));
            }
        }
    }
    assert!(unknown.is_empty(), "sense unknown at {unknown:?}");
}

/// The enclosed volume is the volume: a third of the wedge integral
/// (it read three times it). Planar faces integrate exactly on the
/// midpoint grid (a unit cube, one scaled by two).
#[test]
fn the_enclosed_volume_of_boxes_is_their_volume() {
    use crate::geometry::transform::Transform3;
    let tol = crate::tolerance::Tolerance::default();
    let cube = crate::topology::BRep::unit_cube(tol).expect("cube");
    let long = crate::topology::BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::scale_uniform(2.0)))
        .expect("box");
    for (body, want) in [(cube, 1.0), (long, 8.0)] {
        let faces: Vec<_> = body.faces.keys().collect();
        let v = body.enclosed_folded_volume(&faces).expect("volume");
        assert!((v - want).abs() < 1e-9, "{v} != {want}");
    }
}
