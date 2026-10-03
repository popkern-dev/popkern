use super::*;

#[test]
fn rigid_motions_are_accepted_bit_for_bit() {
    let turn = Placement::rotation(Point3::new(1.0, 2.0, 3.0), Vec3::new(0.0, 0.0, 1.0), 90.0)
        .expect("rotation");
    let back = Placement::from_matrix(*turn.transform(), "x").expect("rigid");
    assert!(orthonormality_defect(turn.transform()) <= 4e-16);
    assert_eq!(back.origin(), turn.origin());
    let moved = turn.then(&Placement::translation(1.0, 0.0, 0.0));
    let p = moved.apply_point(&Point3::new(2.0, 2.0, 3.0));
    assert!(
        (p.x - 2.0).abs() < 1e-12 && (p.y - 3.0).abs() < 1e-12,
        "{p:?}"
    );
}

#[test]
fn scale_shear_and_mirror_are_refused() {
    for (matrix, mirror) in [
        (Transform3::scale(2.0, 2.0, 2.0), ""),
        (Transform3::scale(1.0, 1.0, -1.0), ", mirrored"),
    ] {
        match Placement::from_matrix(matrix, "bolt_1") {
            Err(AssemblyError::NonRigidPlacement {
                instance,
                mirror: m,
                ..
            }) => {
                assert_eq!(instance, "bolt_1");
                assert_eq!(m, mirror);
            }
            other => panic!("{other:?}"),
        }
    }
    let mut shear = Transform3::identity();
    shear.matrix[0][1] = 0.1;
    assert!(Placement::from_matrix(shear, "s").is_err());
}

#[test]
fn round_off_is_reorthonormalised_and_larger_defects_are_refused() {
    let mut nearly =
        *Placement::rotation(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 0.0), 30.0)
            .expect("rotation")
            .transform();
    nearly.matrix[0][0] += 1e-11;
    let fixed = Placement::from_matrix(nearly, "n").expect("repaired");
    assert!(orthonormality_defect(fixed.transform()) < 1e-15);
    nearly.matrix[0][0] += 1e-6;
    assert!(Placement::from_matrix(nearly, "n").is_err());
}

#[test]
fn inverse_undoes_the_motion() {
    let p = Placement::rotation(Point3::new(1.0, -2.0, 0.5), Vec3::new(0.3, 1.0, 0.2), 47.0)
        .expect("rotation")
        .then(&Placement::translation(3.0, 4.0, 5.0));
    let q = Point3::new(0.7, 0.1, -0.4);
    let r = p.inverse().apply_point(&p.apply_point(&q));
    assert!((r.x - q.x).abs() < 1e-14 && (r.y - q.y).abs() < 1e-14 && (r.z - q.z).abs() < 1e-14);
}

#[test]
fn serde_form_is_the_matrix_and_validates_on_read() {
    let p = Placement::translation(1.0, 2.0, 3.0);
    let json = serde_json::to_string(&p).expect("json");
    assert_eq!(
        json,
        "[[1.0,0.0,0.0,1.0],[0.0,1.0,0.0,2.0],[0.0,0.0,1.0,3.0],[0.0,0.0,0.0,1.0]]"
    );
    let back: Placement = serde_json::from_str(&json).expect("read");
    assert_eq!(back, p);
    assert!(
        serde_json::from_str::<Placement>("[[2,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]").is_err()
    );
}
