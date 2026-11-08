use crate::error::ValidationError;
use crate::geometry::transform::Transform3;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// The outer shell's folded-normal volume of `body`'s only solid region.
fn outer_volume(body: &BRep) -> f64 {
    let region = body.solid_regions()[0];
    let shell = body.region_outer_shell(region).expect("an outer shell");
    body.enclosed_folded_volume(&body.shell_faces(shell))
        .expect("determinate")
}

/// Every primitive faces out of its material: the planar Euler builds
/// (cube, prisms of either profile winding, pyramids, frustum) and the
/// curved ones (sphere, cylinder, cone, torus).
#[test]
fn every_primitive_faces_out_of_its_material() {
    let square = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let mut clockwise = square;
    clockwise.reverse();
    // An L whose loops start at its reflex corner, and a block whose
    // side starts across a collinear corner: three corners read the
    // first's cap the wrong way and the second's not at all.
    let ell = [
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 2.0, 0.0),
        Point3::new(1.0, 2.0, 0.0),
    ];
    let split_side = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(4.0, 2.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
    ];
    let up = Vec3::new(0.0, 0.0, 1.0);
    let bodies = [
        ("cube", BRep::unit_cube(tol())),
        ("ccw prism", BRep::prism(&square, up, 2.0, tol())),
        ("cw prism", BRep::prism(&clockwise, up, 2.0, tol())),
        ("reflex-start prism", BRep::prism(&ell, up, 1.0, tol())),
        ("split-side prism", BRep::prism(&split_side, up, 1.0, tol())),
        ("square pyramid", BRep::square_pyramid(2.0, 1.0, tol())),
        ("hex pyramid", BRep::n_gonal_pyramid(6, 1.0, 2.0, tol())),
        ("sphere", BRep::sphere(1.5, tol())),
        ("cylinder", BRep::cylinder(0.5, 2.0, tol())),
        ("cone", BRep::cone(0.5, 2.0, tol())),
        ("torus", BRep::torus(2.0, 0.5, tol())),
    ];
    for (name, body) in bodies {
        let body = body.unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(outer_volume(&body) > 0.0, "{name} faces into its material");
        let errors = body.validate_all();
        assert!(errors.is_empty(), "{name}: {errors:?}");
    }
}

/// The negative control: a cube turned inside out is reported, naming
/// its region and outer shell.
#[test]
fn a_solid_facing_into_its_material_is_reported() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    cube.reverse_body_orientation();
    let region = cube.solid_regions()[0];
    let outer = cube.region_outer_shell(region);
    let errors = cube.validate_region_sense_all();
    assert!(
        matches!(
            errors.as_slice(),
            [ValidationError::RegionSense { region: r, shell, outer: true, volume }]
                if *r == region && Some(*shell) == outer && *volume < 0.0
        ),
        "{errors:?}"
    );
    assert!(
        cube.validate_all()
            .iter()
            .any(|e| matches!(e, ValidationError::RegionSense { .. }))
    );
}

/// A small cube inside a big one, turned out of its own material and
/// inverted into the host, is a cavity whose faces point into the
/// hollow (negative volume); without the turn its faces point into the
/// host's material and the cavity shell is reported.
#[test]
fn a_cavity_faces_into_the_hollow() {
    let hollow = |turn: bool| {
        let mut host = BRep::unit_cube(tol()).expect("cube");
        let mut tool = BRep::unit_cube(tol())
            .and_then(|t| t.transform(&Transform3::scale_uniform(0.5)))
            .and_then(|t| t.transform(&Transform3::translation(0.25, 0.25, 0.25)))
            .expect("small cube");
        if turn {
            tool.reverse_body_orientation();
        }
        let solid = host.solid_regions()[0];
        let map = host.append_body(&tool);
        host.invert_regions_into(&[map.regions[&tool.solid_regions()[0]]], solid);
        host
    };
    let good = hollow(true);
    let errors = good.validate_region_sense_all();
    assert!(errors.is_empty(), "{errors:?}");

    let bad = hollow(false);
    let errors = bad.validate_region_sense_all();
    assert!(
        matches!(
            errors.as_slice(),
            [ValidationError::RegionSense { outer: false, volume, .. }] if *volume > 0.0
        ),
        "{errors:?}"
    );
}

/// A sheet bounds no solid region, so it is not judged either way.
#[test]
fn a_sheet_is_not_judged() {
    let mut open = BRep::open_box(tol()).expect("open box");
    assert!(open.validate_region_sense_all().is_empty());
    open.reverse_body_orientation();
    assert!(open.validate_region_sense_all().is_empty());
}

/// A hollow cube turned inside out after its regions were laid down
/// reads both shells the wrong way; exchanging the outer shell with
/// the cavity's makes the inner cube's skin the outer boundary and the
/// region reads true again, partners and all.
#[test]
fn exchanging_the_outer_shell_fixes_a_body_turned_inside_out() {
    let mut host = BRep::unit_cube(tol()).expect("cube");
    let mut tool = BRep::unit_cube(tol())
        .and_then(|t| t.transform(&Transform3::scale_uniform(0.5)))
        .and_then(|t| t.transform(&Transform3::translation(0.25, 0.25, 0.25)))
        .expect("small cube");
    tool.reverse_body_orientation();
    let solid = host.solid_regions()[0];
    let map = host.append_body(&tool);
    host.invert_regions_into(&[map.regions[&tool.solid_regions()[0]]], solid);
    host.reverse_body_orientation();
    assert_eq!(host.validate_region_sense_all().len(), 2);

    let outer = host.region_outer_shell(solid).expect("outer");
    let cavity = host.region_shells(solid)[1];
    assert!(!host.exchange_outer_shell(solid, outer), "not a cavity");
    assert!(host.exchange_outer_shell(solid, cavity));
    assert_eq!(host.region_outer_shell(solid), Some(cavity));
    let errors = host.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    let hollow = host.cavity_regions()[0];
    let hollow_shell = host.region_outer_shell(hollow).expect("cavity shell");
    assert!(
        host.shell_faces(hollow_shell)
            .iter()
            .all(|&f| host.face_back_shell(f) == Some(outer))
    );
}
