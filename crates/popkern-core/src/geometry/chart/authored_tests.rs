use super::*;
use crate::geometry::{Point3, Vec3};

fn sphere() -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    }
}

fn cone() -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 1.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: 0.4,
    }
}

#[test]
fn a_sphere_pole_stated_as_a_chart_row_is_the_measured_one() {
    let d = ChartDescriptor::of(&sphere(), &Tolerance::default());
    let row = d.singular_rows[0];
    let got = authored_row(
        &sphere(),
        &Tolerance::default(),
        (row.axis, row.at),
        Some(SingularKind::Chart),
    )
    .expect("the pole");
    assert_eq!(got, row);
}

#[test]
fn a_row_stated_where_none_collapses_is_missing() {
    let d = ChartDescriptor::of(&sphere(), &Tolerance::default());
    let row = d.singular_rows[0];
    let middle = match row.axis {
        Axis::V => 0.5 * (d.domain.1.0 + d.domain.1.1),
        Axis::U => 0.5 * (d.domain.0.0 + d.domain.0.1),
    };
    let err = authored_row(
        &sphere(),
        &Tolerance::default(),
        (row.axis, middle),
        Some(SingularKind::Chart),
    )
    .expect_err("no row there");
    assert!(matches!(err, AuthoredRowMismatch::Missing { .. }), "{err}");
}

#[test]
fn an_apex_stated_as_a_smooth_hub_differs_in_kind() {
    let d = ChartDescriptor::of(&cone(), &Tolerance::default());
    let Some(row) = d.singular_rows.first().copied() else {
        panic!("a cone's apex row: {d:?}");
    };
    assert_eq!(row.kind, SingularKind::Apex);
    let err = authored_row(
        &cone(),
        &Tolerance::default(),
        (row.axis, row.at),
        Some(SingularKind::Chart),
    )
    .expect_err("an apex");
    assert!(
        matches!(
            err,
            AuthoredRowMismatch::KindDiffers {
                measured: SingularKind::Apex,
                ..
            }
        ),
        "{err}"
    );
}

#[test]
fn a_row_stated_of_either_kind_matches_an_apex() {
    let d = ChartDescriptor::of(&cone(), &Tolerance::default());
    let row = d.singular_rows[0];
    let got = authored_row(&cone(), &Tolerance::default(), (row.axis, row.at), None)
        .expect("the apex, of whichever kind");
    assert_eq!(got.kind, SingularKind::Apex);
}
