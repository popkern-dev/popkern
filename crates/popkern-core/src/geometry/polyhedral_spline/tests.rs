use super::*;
use crate::geometry::Point3;
use crate::geometry::SurfaceEval;

fn point(i: usize) -> Point3 {
    Point3::new(i as f64, (i % 5) as f64 * 0.25, (i % 7) as f64 * 0.125)
}

fn points(n: usize) -> Vec<Point3> {
    (0..n).map(point).collect()
}

#[test]
fn regular_quad_tiles_the_dual_cell_by_exact_quarters() {
    let stencil: [Point3; 9] = points(9).try_into().unwrap();
    let spline = construct_regular_quad(stencil).unwrap();
    // The biquadratic mask patch spans the vertex's dual cell; the
    // construction ships it split at the parametric midlines so that
    // regular components share the extraordinary sectors' quarter-cell
    // boundary granularity.
    assert_eq!(spline.patches.len(), 4);
    assert!(
        spline
            .patches
            .iter()
            .all(|p| (p.degree_u, p.degree_v) == (2, 2))
    );
    // De Casteljau subdivision is exact: the quarters reproduce the
    // undivided mask patch pointwise.
    let whole = super::mask::apply_inline_mask_to_patch(
        &super::tables::REGULAR_MASK,
        PolyhedralSplineFamily::RegularQuad,
        &stencil,
    )
    .unwrap();
    let quarters = [(0, 0.0, 0.0), (1, 0.0, 0.5), (2, 0.5, 0.0), (3, 0.5, 0.5)];
    for (patch_index, u0, v0) in quarters {
        let quarter = &spline.patches[patch_index];
        for iu in 0..=4 {
            for iv in 0..=4 {
                let (s, t) = (iu as f64 / 4.0, iv as f64 / 4.0);
                let got = quarter.evaluate(s, t);
                let want = whole.evaluate(u0 + 0.5 * s, v0 + 0.5 * t);
                assert!(
                    got.distance_to(&want) < 1e-12,
                    "quarter {patch_index} deviates at ({s},{t})"
                );
            }
        }
    }
}

#[test]
fn all_table_driven_families_construct_expected_patch_counts() {
    for (n, count, stencil_len) in [(3, 3, 7), (5, 5, 11), (6, 24, 13), (7, 28, 15), (8, 32, 17)] {
        let spline = construct_extraordinary_star(n, &points(stencil_len)).unwrap();
        assert_eq!(spline.patches.len(), count, "extraordinary n={n}");
        assert!(
            spline
                .patches
                .iter()
                .all(|p| (p.degree_u, p.degree_v) == (3, 3))
        );
    }

    for (n, count, stencil_len) in [
        (3, 3, 12),
        (5, 5, 20),
        (6, 24, 24),
        (7, 28, 28),
        (8, 32, 32),
    ] {
        let spline = construct_ngon(n, &points(stencil_len)).unwrap();
        assert_eq!(spline.patches.len(), count, "ngon n={n}");
        assert!(
            spline
                .patches
                .iter()
                .all(|p| (p.degree_u, p.degree_v) == (3, 3))
        );
    }

    for n in 3..=8 {
        let spline = construct_polar(n, &points(n + 1)).unwrap();
        assert_eq!(spline.patches.len(), n, "polar n={n}");
        assert!(
            spline
                .patches
                .iter()
                .all(|p| (p.degree_u, p.degree_v) == (3, 2))
        );
    }

    assert_eq!(
        construct_t0(points(14).try_into().unwrap())
            .unwrap()
            .patches
            .len(),
        4
    );
    assert_eq!(
        construct_t1(points(18).try_into().unwrap())
            .unwrap()
            .patches
            .len(),
        8
    );
    assert_eq!(
        construct_t2(points(20).try_into().unwrap())
            .unwrap()
            .patches
            .len(),
        16
    );
}

#[test]
fn unsupported_valence_declines_with_structured_error() {
    let err = construct_ngon(4, &points(16)).unwrap_err();
    assert!(matches!(
        err,
        PolyhedralSplineError::UnsupportedValence {
            family: "n-gon",
            value: 4,
            ..
        }
    ));
}

#[test]
fn wrong_stencil_count_declines_before_mask_application() {
    let err = construct_polar(5, &points(5)).unwrap_err();
    assert!(matches!(
        err,
        PolyhedralSplineError::StencilPointCount {
            family: PolyhedralSplineFamily::Polar { sectors: 5 },
            expected: 6,
            got: 5,
        }
    ));
}

#[test]
fn degree_raise_preserves_surface_evaluation() {
    let patch = &construct_polar(5, &points(6)).unwrap().patches[0];
    let raised = patch.degree_raised_to_bicubic();
    assert_eq!((raised.degree_u, raised.degree_v), (3, 3));
    for (u, v) in [(0.0, 0.0), (0.25, 0.75), (0.5, 0.5), (1.0, 1.0)] {
        let a = patch.evaluate(u, v);
        let b = raised.evaluate(u, v);
        assert!(a.approx_eq(&b, 1e-10), "u={u}, v={v}: {a:?} != {b:?}");
    }
}

#[test]
fn bezier_nurbs_conversion_preserves_evaluation() {
    let patch = &construct_t0(points(14).try_into().unwrap())
        .unwrap()
        .patches[0];
    let nurbs = patch.to_bezier_nurbs();
    for (u, v) in [(0.0, 0.0), (0.2, 0.8), (0.5, 0.5), (1.0, 1.0)] {
        let a = patch.evaluate(u, v);
        let b = nurbs.point_at(u, v);
        assert!(a.approx_eq(&b, 1e-10), "u={u}, v={v}: {a:?} != {b:?}");
    }
}

#[test]
fn c0_edge_validation_accepts_identical_edges() {
    let patch = construct_regular_quad(points(9).try_into().unwrap())
        .unwrap()
        .patches
        .remove(0);
    let report = validate_c0_edge_match(&patch, PatchEdge::UMin, &patch, PatchEdge::UMin, 9, 1e-12);
    assert!(report.accepted);
    assert_eq!(report.max_position_gap, 0.0);
}

#[test]
fn control_net_adapter_constructs_regular_quad_vertex() {
    let verts = points(9);
    let faces = vec![
        vec![0, 1, 4, 3],
        vec![1, 2, 5, 4],
        vec![3, 4, 7, 6],
        vec![4, 5, 8, 7],
    ];
    let net = ControlNet::new(verts, faces).unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    assert_eq!(splines.len(), 1);
    assert_eq!(splines[0].family, PolyhedralSplineFamily::RegularQuad);
}

#[test]
fn control_net_adapter_constructs_polar_fan() {
    let verts = points(4);
    let faces = vec![vec![0, 1, 2], vec![0, 2, 3], vec![0, 3, 1]];
    let net = ControlNet::new(verts, faces).unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    assert_eq!(splines.len(), 1);
    assert_eq!(
        splines[0].family,
        PolyhedralSplineFamily::Polar { sectors: 3 }
    );
    assert_eq!(splines[0].patches.len(), 3);
}

#[test]
fn reference_obj_polar_fixture_constructs_expected_family() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar3sct.obj"
    ))
    .unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    assert_eq!(splines.len(), 1);
    assert_eq!(
        splines[0].family,
        PolyhedralSplineFamily::Polar { sectors: 3 }
    );
}

#[test]
fn reference_obj_extraordinary_fixture_constructs_expected_family() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop3sct.obj"
    ))
    .unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    assert_eq!(splines.len(), 1);
    assert_eq!(
        splines[0].family,
        PolyhedralSplineFamily::ExtraordinaryStar { valence: 3 }
    );
}

#[test]
fn reference_obj_ngon_fixture_constructs_expected_family() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon3.obj"
    ))
    .unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    assert_eq!(splines.len(), 1);
    assert_eq!(splines[0].family, PolyhedralSplineFamily::Ngon { sides: 3 });
}

#[test]
fn reference_obj_t_fixtures_construct_expected_families() {
    for (fixture, family, patch_count) in [
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"),
            PolyhedralSplineFamily::T0,
            4,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T1.obj"),
            PolyhedralSplineFamily::T1,
            8,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T2.obj"),
            PolyhedralSplineFamily::T2,
            16,
        ),
    ] {
        let net = ControlNet::from_obj_str(fixture).unwrap();
        let splines = construct_from_control_net(&net).unwrap();
        assert_eq!(splines.len(), 1, "{family:?}");
        assert_eq!(splines[0].family, family);
        assert_eq!(splines[0].patches.len(), patch_count);
    }
}

#[test]
fn atlas_detects_reference_fixture_internal_seams() {
    for (fixture, expected_seams) in [
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar3sct.obj"
            ),
            3,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop3sct.obj"
            ),
            3,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon3.obj"),
            3,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"),
            4,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T1.obj"),
            8,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T2.obj"),
            24,
        ),
    ] {
        let net = ControlNet::from_obj_str(fixture).unwrap();
        let splines = construct_from_control_net(&net).unwrap();
        let atlas = build_patch_atlas(&splines[0], 9, 1e-9);
        assert_eq!(atlas.patch_count, splines[0].patches.len());
        assert!(
            atlas.seams.len() >= expected_seams,
            "{:?}: seams={}",
            splines[0].family,
            atlas.seams.len()
        );
    }
}

#[test]
fn seam_validation_accepts_reference_fixtures() {
    for fixture in [
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar3sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop3sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon3.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T1.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T2.obj"),
    ] {
        let net = ControlNet::from_obj_str(fixture).unwrap();
        let splines = construct_from_control_net(&net).unwrap();
        let reports = validate_spline_atlas(&splines[0], 9, 1e-8, 1e-6)
            .unwrap_or_else(|err| panic!("{:?}: {err}", splines[0].family));
        assert!(!reports.is_empty(), "{:?}", splines[0].family);
    }
}

#[test]
fn seam_validation_accepts_reference_valence_matrix() {
    for (fixture, family, patch_count) in [
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar3sct.obj"
            ),
            PolyhedralSplineFamily::Polar { sectors: 3 },
            3,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar4sct.obj"
            ),
            PolyhedralSplineFamily::Polar { sectors: 4 },
            4,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar5sct.obj"
            ),
            PolyhedralSplineFamily::Polar { sectors: 5 },
            5,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar6sct.obj"
            ),
            PolyhedralSplineFamily::Polar { sectors: 6 },
            6,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar7sct.obj"
            ),
            PolyhedralSplineFamily::Polar { sectors: 7 },
            7,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar8sct.obj"
            ),
            PolyhedralSplineFamily::Polar { sectors: 8 },
            8,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop3sct.obj"
            ),
            PolyhedralSplineFamily::ExtraordinaryStar { valence: 3 },
            3,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop5sct.obj"
            ),
            PolyhedralSplineFamily::ExtraordinaryStar { valence: 5 },
            5,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop6sct.obj"
            ),
            PolyhedralSplineFamily::ExtraordinaryStar { valence: 6 },
            24,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop7sct.obj"
            ),
            PolyhedralSplineFamily::ExtraordinaryStar { valence: 7 },
            28,
        ),
        (
            include_str!(
                "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop8sct.obj"
            ),
            PolyhedralSplineFamily::ExtraordinaryStar { valence: 8 },
            32,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon3.obj"),
            PolyhedralSplineFamily::Ngon { sides: 3 },
            3,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon5.obj"),
            PolyhedralSplineFamily::Ngon { sides: 5 },
            5,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon6.obj"),
            PolyhedralSplineFamily::Ngon { sides: 6 },
            24,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon7.obj"),
            PolyhedralSplineFamily::Ngon { sides: 7 },
            28,
        ),
        (
            include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon8.obj"),
            PolyhedralSplineFamily::Ngon { sides: 8 },
            32,
        ),
    ] {
        let net = ControlNet::from_obj_str(fixture).unwrap();
        let splines = construct_from_control_net(&net).unwrap();
        assert_eq!(splines.len(), 1, "{family:?}");
        assert_eq!(splines[0].family, family);
        assert_eq!(splines[0].patches.len(), patch_count, "{family:?}");
        let reports = validate_spline_atlas(&splines[0], 11, 1e-8, 1e-6)
            .unwrap_or_else(|err| panic!("{family:?}: {err}"));
        assert!(!reports.is_empty(), "{family:?}");
    }
}

#[test]
fn commercial_validation_accepts_reference_valence_matrix() {
    for fixture in [
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar3sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar4sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar5sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar6sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar7sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar8sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop3sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop5sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop6sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop7sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop8sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon3.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon5.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon6.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon7.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon8.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T1.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T2.obj"),
    ] {
        let net = ControlNet::from_obj_str(fixture).unwrap();
        let splines = construct_from_control_net(&net).unwrap();
        let report = validate_polyhedral_spline(&splines[0], 11, 1e-8, 1e-6, 1e-12)
            .unwrap_or_else(|err| panic!("{:?}: {err}", splines[0].family));
        assert_eq!(report.patch_count, splines[0].patches.len());
        assert!(!report.seam_reports.is_empty(), "{:?}", splines[0].family);
        assert_eq!(report.regularity_reports.len(), splines[0].patches.len());
        assert_eq!(report.brep_readiness.patch_count, splines[0].patches.len());
    }
}
