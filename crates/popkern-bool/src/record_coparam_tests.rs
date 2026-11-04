use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Every record's uv images against its 3D curve at matched fractions.
fn worst_misfit(a: &BRep, b: &BRep) -> f64 {
    let tol = Tolerance::default();
    let ix = crate::intersect::find_intersections(a, b, &tol).expect("intersections");
    assert!(!ix.records.is_empty());
    let mut worst = 0.0_f64;
    for r in &ix.records {
        for (brep, face, uv) in [
            (a, r.face_a, &r.curve.curve_uv1),
            (b, r.face_b, &r.curve.curve_uv2),
        ] {
            let surface = &brep.faces[face].surface.as_ref().expect("surface").geometry;
            worst = worst.max(super::misfit(&r.curve.curve_3d, uv, surface));
        }
    }
    worst
}

fn turned(b: BRep, r: &Transform3) -> BRep {
    b.transform(r).expect("turned")
}

#[test]
fn record_uv_curves_follow_their_3d_curve_at_every_fraction() {
    let tol = Tolerance::default();
    let gate = super::GATE_FRACTION * tol.linear();
    let generic = [
        Transform3::rotation_x(0.37).compose(&Transform3::rotation_z(1.1)),
        Transform3::rotation_y(2.3).compose(&Transform3::rotation_x(-0.8)),
        Transform3::rotation_z(0.2).compose(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2)),
    ];
    for r in &generic {
        // A sphere drilled through, and a cone standing on a cylinder.
        let ball = turned(BRep::sphere(0.5, tol).expect("sphere"), r);
        let bore = BRep::cylinder(0.2, 2.0, tol)
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -1.0)))
            .expect("bore");
        let bore = turned(bore, r);
        let w = worst_misfit(&ball, &bore);
        assert!(w <= gate, "drilled sphere {r:?}: {w:.3e}");
        let cone = BRep::cone(0.6, 1.0, tol)
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 0.3)))
            .expect("cone");
        let rod = BRep::cylinder(0.4, 2.0, tol)
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -1.0)))
            .expect("rod");
        let w = worst_misfit(&turned(cone, r), &turned(rod, r));
        assert!(w <= gate, "cone on rod {r:?}: {w:.3e}");
    }
}
