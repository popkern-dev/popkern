//! The serve tail's locality (PK-815, docs/architecture/local-boolean-
//! scope.md "the touched region"): a boolean normalizes the structure
//! it created or changed, and structure it carried unchanged from `a`
//! keeps its form.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn axis_box(min: Point3, max: Point3) -> BRep {
    let profile = [
        Point3::new(min.x, min.y, min.z),
        Point3::new(max.x, min.y, min.z),
        Point3::new(max.x, max.y, min.z),
        Point3::new(min.x, max.y, min.z),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), max.z - min.z, tol()).expect("box")
}

/// Faces on the plane y = 0 with outward normal -y.
fn faces_on_y0(b: &BRep) -> usize {
    b.faces
        .values()
        .filter(|f| {
            f.surface.as_ref().is_some_and(|s| match &s.geometry {
                SurfaceGeometry::Plane { origin, normal, .. } => {
                    origin.y.abs() < 1e-9 && normal.x.abs() < 1e-9 && normal.z.abs() < 1e-9
                }
                _ => false,
            })
        })
        .count()
}

/// A 4 x 2 x 1 block whose y = 0 side is two coplanar faces (the
/// profile carries a collinear vertex at x = 2): a co-surfacial pair
/// the block arrives with.
fn split_sided_block() -> BRep {
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(4.0, 2.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), 1.0, tol()).expect("block")
}

#[test]
fn a_pair_the_body_arrived_with_keeps_its_form_away_from_the_cut() {
    let block = split_sided_block();
    assert_eq!(faces_on_y0(&block), 2, "the fixture's y = 0 side is split");
    // A pocket in the far corner, clear of the y = 0 side.
    let pocket = axis_box(Point3::new(3.0, 1.5, 0.5), Point3::new(3.5, 2.5, 2.0));
    let out = boolean(BoolOperation::Subtraction, block, pocket, tol(), None).expect("pocket");
    out.validate().expect("valid");
    assert_eq!(
        faces_on_y0(&out),
        2,
        "the untouched pair is not merged by the cut"
    );
}

#[test]
fn a_pair_the_boolean_creates_is_merged() {
    // Two flush boxes side by side: their y = 0 faces meet edge to edge
    // on one plane only in the union, which merges them.
    let left = axis_box(Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 1.0));
    let right = axis_box(Point3::new(2.0, 0.0, 0.0), Point3::new(4.0, 2.0, 1.0));
    let out = boolean(BoolOperation::Union, left, right, tol(), None).expect("union");
    out.validate().expect("valid");
    assert_eq!(
        faces_on_y0(&out),
        1,
        "the boolean's own co-surfacial pair is merged"
    );
    assert_eq!(out.faces.len(), 6, "the union is one box");
}
