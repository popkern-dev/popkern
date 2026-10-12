//! PK-761: a record lying on a face's own boundary edge with the
//! touching-only shortcut out of reach. A pin stands with its whole
//! bottom rim flush on a base's top face (the rim is a closed record on
//! the pin wall's own boundary edge, PK-756's class) while a bar of the
//! base runs through the pin at mid-height, so the bodies overlap and
//! the full pipeline runs. A split along the rim would cut a
//! zero-height sliver off the pin wall whose sample lies on the base's
//! plane; the union must serve on its charts with no such fragment
//! (the pin wall keeps its main region and the bar's two windows, the
//! flush bottom disc classifies OnBoundary). The same body built in
//! popkern live reads in OpenCASCADE as valid with volume 569.15391936,
//! the analytic 569.15391936. A sliver would show as an extra wall face
//! or a decline.

use popkern_bool::chart_walk_fallback::take_chart_walk_fallback_events;
use popkern_bool::merge_fallback::take_merge_fallback_events;
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{AnalyticProfile, ProfileShape, SketchPlane};

fn extrude(shape: ProfileShape<f64>, center: [f64; 2], z: f64, height: f64) -> BRep {
    let sheet = AnalyticProfile { shape, center }
        .to_brep(
            SketchPlane {
                origin: Point3::new(0.0, 0.0, z),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
            },
            Tolerance::default(),
            None,
        )
        .expect("sheet");
    let face = sheet.regions[0].face;
    let mut b = sheet.brep;
    b.extrude(face, Vec3::new(0.0, 0.0, 1.0), height, None)
        .expect("extrude");
    b
}

fn rect(w: f64, h: f64, center: [f64; 2], z: f64, height: f64) -> BRep {
    extrude(
        ProfileShape::Rect {
            width: w,
            height: h,
        },
        center,
        z,
        height,
    )
}

fn union_on_charts(a: BRep, b: BRep, name: &str) -> BRep {
    let _ = take_merge_fallback_events();
    let _ = take_chart_walk_fallback_events();
    let out = boolean(BoolOperation::Union, a, b, Tolerance::default(), None)
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    out.validate().unwrap_or_else(|e| panic!("{name}: {e}"));
    out.validate_curve_surface()
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    let merges = take_merge_fallback_events();
    assert!(merges.is_empty(), "{name}: graph-merge fallback {merges:?}");
    let walks = take_chart_walk_fallback_events();
    assert!(walks.is_empty(), "{name}: faces walked off-chart {walks:?}");
    out
}

#[test]
fn a_pin_flush_on_a_base_and_crossed_by_its_bar_unions_without_a_rim_sliver() {
    // Posts rooted inside the base and a bar overhanging them: no face
    // of the supports is flush with another, so only the pin's rim is.
    let base = rect(10.0, 10.0, [0.0, 0.0], 0.0, 5.0);
    let base = union_on_charts(base, rect(0.4, 0.6, [-4.4, 0.0], 4.0, 2.5), "post 1");
    let base = union_on_charts(base, rect(0.4, 0.6, [4.4, 0.0], 4.0, 2.5), "post 2");
    let base = union_on_charts(base, rect(9.8, 1.0, [0.0, 0.0], 6.0, 1.0), "bar");
    let pin = extrude(ProfileShape::Circle { radius: 2.0 }, [0.0, 0.0], 5.0, 5.0);
    let out = union_on_charts(base, pin, "pin");
    let walls = out
        .faces
        .values()
        .filter(|f| {
            f.surface.as_ref().is_some_and(|s| {
                matches!(
                    s.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }
                )
            })
        })
        .count();
    assert_eq!(walls, 2, "the pin wall's regions, no rim sliver");
    assert_eq!(out.faces.len(), 29, "faces");
}
