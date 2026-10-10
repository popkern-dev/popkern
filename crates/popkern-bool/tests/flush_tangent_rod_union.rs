//! PK-750 at the boolean boundary: a box whose caps are flush with a
//! rod's and whose long faces are tangent to it along the box's own
//! vertical edges. The rod is built as popkern live builds it (an
//! analytic circle profile whose seam sits at the top of the circle,
//! exactly on the tangent line), so the box's x-face cut line takes
//! over the rod's seam and the rod's near rim arc is a coincident
//! record. The union must serve on the cellular path with every face
//! walked on its chart: no graph-merge fallback and no 3D-Newell walk.

use popkern_bool::chart_walk_fallback::take_chart_walk_fallback_events;
use popkern_bool::merge_fallback::take_merge_fallback_events;
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{AnalyticProfile, ProfileShape, SketchPlane};

const R: f64 = 4.5;
const H: f64 = 16.0;

fn extrude(shape: ProfileShape<f64>, plane: SketchPlane) -> BRep {
    let sheet = AnalyticProfile {
        shape,
        center: [0.0, 0.0],
    }
    .to_brep(plane, Tolerance::default(), None)
    .expect("sheet");
    let face = sheet.regions[0].face;
    let mut b = sheet.brep;
    b.extrude(face, Vec3::new(0.0, 0.0, 1.0), H, None)
        .expect("extrude");
    b
}

/// A rod of radius R on the z axis through (x, 0), its circle starting
/// at +y (the seam line at (x, R)).
fn rod(x: f64) -> BRep {
    extrude(
        ProfileShape::Circle { radius: R },
        SketchPlane {
            origin: Point3::new(x, 0.0, 0.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(-1.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        },
    )
}

/// The 12 × 2R box between x = ±6, flush with the rods' caps.
fn slab() -> BRep {
    extrude(
        ProfileShape::Rect {
            width: 12.0,
            height: 2.0 * R,
        },
        SketchPlane {
            origin: Point3::new(0.0, 0.0, 0.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        },
    )
}

fn serve_on_charts(a: BRep, b: BRep, name: &str) -> BRep {
    let _ = take_merge_fallback_events();
    let _ = take_chart_walk_fallback_events();
    let out = boolean(BoolOperation::Union, a, b, Tolerance::default(), None)
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    out.validate().unwrap_or_else(|e| panic!("{name}: {e}"));
    let merges = take_merge_fallback_events();
    assert!(merges.is_empty(), "{name}: graph-merge fallback {merges:?}");
    let walks = take_chart_walk_fallback_events();
    assert!(walks.is_empty(), "{name}: faces walked off-chart {walks:?}");
    out
}

#[test]
fn a_flush_tangent_box_unions_with_a_rod_on_its_charts() {
    let half = serve_on_charts(rod(-6.0), slab(), "half stadium");
    assert_eq!(half.faces.len(), 6, "half stadium faces");
}

#[test]
fn a_flush_tangent_box_unions_with_two_rods_on_their_charts() {
    let rods = serve_on_charts(rod(-6.0), rod(6.0), "two rods");
    let stadium = serve_on_charts(rods, slab(), "stadium");
    assert_eq!(stadium.faces.len(), 6, "stadium faces");
}
