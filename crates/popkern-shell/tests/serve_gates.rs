//! Serve-gate and decline-path tests for popkern-shell (PK-321/PK-330):
//! collapse declines with binding entities, curved-body validation
//! declines (un-narrowed by PK-324), open-sheet containment, and
//! approximation-honesty warnings.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn shell_past_feature_size_declines_offset_too_large() {
    // A 0.6 wall on a unit cube crosses the opposite walls — a pure
    // plane-plane collision the retired midpoint-normal detector was
    // structurally blind to. The collapse gate declines it with the
    // binding edge and a usable max-safe estimate.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = shell(
        &brep,
        0.6,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    match result {
        Err(ShellError::OffsetTooLarge {
            requested,
            max_safe,
            binding,
        }) => {
            assert!((requested - 0.6).abs() < 1e-12);
            assert!(max_safe > 0.0 && max_safe < 0.6, "max_safe {max_safe}");
            assert!(binding.contains("edge"), "{binding}");
        }
        other => panic!("expected OffsetTooLarge decline, got {other:?}"),
    }
}

#[test]
fn offset_shrink_past_feature_size_declines() {
    // Shrinking a unit cube by 0.6 per side inverts it — the collapse
    // gate must decline rather than serve an inside-out body.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = popkern_shell::offset::offset_solid(&brep, -0.6, &ShellOptions::default(), tol());
    assert!(
        matches!(result, Err(ShellError::OffsetTooLarge { .. })),
        "expected OffsetTooLarge decline, got {result:?}"
    );
}

#[test]
fn shell_cylinder_serves_with_exact_curved_recompute() {
    // PK-322: analytic recompute arms (coaxial circles, straight-ruling
    // seams, curved vertex projection) make the cylinder shell exact —
    // the stale-clone decline pin flips to a measured serving test.
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let cap = brep
        .faces
        .iter()
        .find(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            )
        })
        .map(|(h, _)| h)
        .unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[cap],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();

    assert!(
        result.brep.validate_all().is_empty(),
        "shelled cylinder must be validate_all-clean: {:?}",
        result.brep.validate_all()
    );
    // Inner void skin: lateral radius shrunk to 0.9 (material-inward),
    // remaining cap moved 0.1 into the body.
    let mut radii: Vec<f64> = result
        .brep
        .faces
        .values()
        .filter_map(|f| f.surface.as_ref())
        .filter_map(|s| match &s.geometry {
            SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    radii.sort_by(f64::total_cmp);
    assert_eq!(radii, vec![0.9, 1.0], "outer lateral kept, inner at 0.9");
}

#[test]
fn shell_thin_prism_declines_offset_too_large() {
    // A 0.2 wall in a prism only 0.05 thick in Z cannot exist: the skin
    // inverts. The collapse gate declines with the binding edge (the old
    // contract accepted "any valid result or error"; the honest outcome
    // is exactly one thing).
    use popkern_core::geometry::{Point3, Vec3};

    let profile = &[
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let brep = BRep::prism(profile, Vec3::new(0.0, 0.0, 1.0), 0.05, tol()).unwrap();
    let result = shell(
        &brep,
        0.2,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    assert!(
        matches!(result, Err(ShellError::OffsetTooLarge { .. })),
        "expected OffsetTooLarge decline, got {result:?}"
    );
}

#[test]
fn shell_midplane_past_feature_size_declines() {
    // 1.2 total wall on a unit cube: the inward half (0.6) inverts.
    // The decline reports the caller's full-thickness units.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = shell(
        &brep,
        1.2,
        ShellDirection::Midplane,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    match result {
        Err(ShellError::OffsetTooLarge {
            requested,
            max_safe,
            ..
        }) => {
            assert!((requested - 1.2).abs() < 1e-9, "requested {requested}");
            assert!(max_safe > 0.9 && max_safe < 1.2, "max_safe {max_safe}");
        }
        other => panic!("expected OffsetTooLarge decline, got {other:?}"),
    }
}

#[test]
fn thicken_open_sheet_serves_closed_solid_with_rim_walls() {
    // An open sheet (cube minus one face) thickens into a CLOSED solid:
    // both skins plus one wall face per free rim edge (PK-324; the old
    // pipeline first served two disconnected skins while claiming to
    // build side walls, then declined with OpenBoundary as containment).
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    let outer_loop = brep.faces[fh].outer_loop;
    brep.loops.remove(outer_loop);
    brep.faces.remove(fh);

    let result = popkern_shell::thicken::thicken(&brep, 0.05, &ShellOptions::default(), tol())
        .expect("open-sheet thicken serves");
    // PK-326 rim corners: the outer skin serves rounded (8 strips on
    // the convex edges, 4 sphere triangles at its closed corners, rim
    // arcs at the 4 boundary corners), so the count is 17 outer +
    // 5 inner + 4 rim wall quads + 4 rim corner sectors.
    assert_eq!(result.brep.faces.len(), 30);
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| matches!(w, popkern_shell::ShellWarning::JoinStyleDowngraded { .. })),
        "rim corners must serve without downgrade: {:?}",
        result.warnings
    );
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "served invalid: {errors:?}");
    // Closed: no loopless half-edges survive.
    let free = result
        .brep
        .half_edges
        .iter()
        .filter(|(_, he)| !result.brep.loops.contains_key(he.loop_))
        .count();
    assert_eq!(free, 0, "free half-edges remain after rim closing");
}

#[test]
fn cone_offset_is_exact_and_silent() {
    // PK-322: the cone gained a closed-form offset arm (apex translated
    // d/sin(half_angle) along the axis) — its offset surface is a Cone,
    // never a NURBS approximation, and carries no approximation warning.
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_shell::classify::classify_faces;
    use popkern_shell::surface_offset::compute_offset_surfaces;

    let brep = BRep::cone(1.0, 2.0, tol()).unwrap();
    let classification = classify_faces(&brep, 0.05, &[], ShellDirection::Inward, &[]).unwrap();
    let result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
    assert!(result.failed.is_empty(), "failed: {:?}", result.failed);
    assert!(
        result.warnings.is_empty(),
        "exact offsets carry no warnings: {:?}",
        result.warnings
    );
    assert!(
        result
            .surfaces
            .values()
            .any(|s| matches!(s, SurfaceGeometry::Cone { .. })),
        "cone lateral face must offset to a Cone"
    );
}

#[test]
fn nurbs_offset_path_emits_approximation_warning() {
    // A genuinely free-form face (NURBS surface) has no closed form, so
    // its offset goes through the degree-3 least-squares fit. The
    // dispatch must SAY so: every fitted face carries
    // FaceOffsetApproximated with the measured deviation (the warning
    // existed for years and was never emitted before PK-321).
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    use popkern_core::geometry::{NurbsSurface, Point3};
    use popkern_shell::classify::FaceClassification;
    use popkern_shell::surface_offset::compute_offset_surfaces;

    // A gently curved quadratic patch grafted onto a cube face — only
    // the offset dispatch runs here, so the body needn't remain
    // watertight around the graft.
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    brep.faces[fh].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(NurbsSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 0.0, 0.1),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 0.5, 0.1),
                Point3::new(0.5, 0.5, 0.2),
                Point3::new(1.0, 0.5, 0.1),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.5, 1.0, 0.1),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 3,
            num_v: 3,
            weights: vec![1.0; 9],
        }),
        orientation: true,
    });
    let classification = FaceClassification {
        offset_faces: vec![(fh, 0.02, true)],
        removed_faces: vec![],
    };
    let result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
    let deviations: Vec<f64> = result
        .warnings
        .iter()
        .filter_map(|w| match w {
            ShellWarning::FaceOffsetApproximated {
                achieved_deviation, ..
            } => Some(*achieved_deviation),
            _ => None,
        })
        .collect();
    assert_eq!(
        deviations.len(),
        1,
        "the fitted face must warn with its measured deviation: {:?}",
        result.warnings
    );
    assert!(deviations[0].is_finite() && deviations[0] >= 0.0);
}

#[test]
fn shell_cone_through_apex_declines_never_serves() {
    // An inward wall so thick the offset cone's apex passes through the
    // base plane (d/sin(half_angle) > height): the exact offset SURFACE
    // still exists (translated cone — through-apex is a trim concern,
    // never an approximation), but the void does not — the offset base
    // plane no longer meets the offset cone's material nappe, the
    // vertex solve declines, and the serve gate contains the stale skin
    // as a transactional ValidationFailed. Upgrading this decline to
    // OffsetTooLarge with the exact bound (height·sin(half_angle)) is
    // max_safe_offset scope (PK-328).
    let brep = BRep::cone(1.0, 2.0, tol()).unwrap();
    let result = shell(
        &brep,
        1.0,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    match result {
        Err(ShellError::ValidationFailed { errors }) => {
            assert!(!errors.is_empty());
        }
        Err(ShellError::OffsetTooLarge { .. }) => {
            // The PK-328 upgrade landed: also acceptable, strictly better.
        }
        other => panic!("expected a decline, got {other:?}"),
    }
}
