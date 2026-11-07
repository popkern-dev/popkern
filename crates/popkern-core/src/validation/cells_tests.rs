use crate::constructors::BRepBuilder;
use crate::error::{RegionDefect, ValidationError};
use crate::geometry::Point3;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceSide, RegionKind, Shell};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Every cell-structure finding on `body`.
fn cells(body: &BRep) -> Vec<ValidationError> {
    let mut errors = body.validate_ring_membership_all();
    errors.extend(body.validate_wedges_all());
    errors.extend(body.validate_shell_cells_all());
    errors.extend(body.validate_regions_all());
    errors
}

fn small_cube_at(offset: f64) -> BRep {
    BRep::unit_cube(tol())
        .and_then(|t| t.transform(&Transform3::scale_uniform(0.5)))
        .and_then(|t| t.transform(&Transform3::translation(offset, offset, offset)))
        .expect("small cube")
}

/// A cube with a cube-shaped cavity in its middle.
fn hollow_cube() -> BRep {
    let mut host = BRep::unit_cube(tol()).expect("cube");
    let mut tool = small_cube_at(0.25);
    tool.reverse_body_orientation();
    let solid = host.solid_regions()[0];
    let map = host.append_body(&tool);
    host.invert_regions_into(&[map.regions[&tool.solid_regions()[0]]], solid);
    host
}

/// Solids (planar, curved, of genus one, hollow, two in one body),
/// sheets and wires all have a sound cell structure.
#[test]
fn every_kind_of_body_has_a_sound_cell_structure() {
    let mut wire = BRepBuilder::new(tol());
    let a = wire.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let b = wire.add_vertex(Point3::new(1.0, 0.0, 0.0));
    wire.add_edge(a, b);
    let mut pair = BRep::unit_cube(tol()).expect("cube");
    pair.append_body(&small_cube_at(3.0));
    let bodies = [
        ("cube", BRep::unit_cube(tol()).expect("cube")),
        ("sphere", BRep::sphere(1.0, tol()).expect("sphere")),
        (
            "cylinder",
            BRep::cylinder(0.5, 2.0, tol()).expect("cylinder"),
        ),
        ("cone", BRep::cone(0.5, 2.0, tol()).expect("cone")),
        ("torus", BRep::torus(2.0, 0.5, tol()).expect("torus")),
        ("hollow cube", hollow_cube()),
        ("two cubes", pair),
        ("open box", BRep::open_box(tol()).expect("open box")),
        ("wire", wire.build_unchecked()),
    ];
    for (name, body) in bodies {
        let errors = cells(&body);
        assert!(errors.is_empty(), "{name}: {errors:?}");
    }
}

/// A cube's mates run against each other; a cylinder rim's mate (a
/// closed edge) runs against its first fin.
#[test]
fn mates_run_against_each_other() {
    for body in [
        BRep::unit_cube(tol()).expect("cube"),
        BRep::cylinder(0.5, 2.0, tol()).expect("cylinder"),
    ] {
        for (edge, data) in &body.edges {
            let first = data.half_edge;
            assert_eq!(body.fin_is_forward(first), Some(true));
            let mate = body.mate(first).expect("a manifold edge");
            assert_eq!(body.fin_is_forward(mate), Some(false), "{edge:?}");
            let wedges = body.edge_wedges(edge).expect("decidable");
            assert_eq!(wedges.len(), 2);
            assert_eq!(wedges[0].from_side, FaceSide::Front);
            assert_eq!(wedges[0].to_side, FaceSide::Front);
            assert_eq!(wedges[1].from_side, FaceSide::Back);
            assert_eq!(wedges[1].to_side, FaceSide::Back);
        }
    }
}

/// Negative control: a fin unlinked from its ring is reported, and
/// the opening it makes in the solid puts each free edge's one wedge
/// across the solid's shell and the outside's.
#[test]
fn a_fin_off_its_ring_and_an_opened_solid_are_reported() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let (edge, first) = cube
        .edges
        .iter()
        .map(|(e, d)| (e, d.half_edge))
        .next()
        .expect("an edge");
    let mate = cube.mate(first).expect("mate");
    cube.half_edges[first].clear_mate();
    cube.half_edges[mate].clear_mate();
    let errors = cube.validate_ring_membership_all();
    assert!(
        matches!(errors.as_slice(), [ValidationError::RingMembership { edge: e, half_edge }]
            if *e == edge && *half_edge == mate),
        "{errors:?}"
    );
    let errors = cube.validate_wedges_all();
    assert!(
        matches!(errors.as_slice(), [ValidationError::WedgeShellMismatch { from, to, .. }]
            if from == to && *from == first),
        "{errors:?}"
    );
}

/// Negative control: a face with its sides exchanged puts both its
/// sides in the wrong shells, so all four of its wedges span two
/// shells.
#[test]
fn a_face_with_its_sides_exchanged_is_reported() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let face = cube.faces.keys().next().expect("a face");
    let data = &mut cube.faces[face];
    std::mem::swap(&mut data.back, &mut data.front);
    let errors = cube.validate_wedges_all();
    assert_eq!(errors.len(), 8, "{errors:?}");
    assert!(
        errors
            .iter()
            .all(|e| matches!(e, ValidationError::WedgeShellMismatch { .. }))
    );
}

/// Negative control: two cubes' solid sides put in one shell (and
/// their outside sides in another) make shells of two pieces, each of
/// two spheres' Euler count.
#[test]
fn a_shell_of_two_pieces_is_reported() {
    let mut pair = BRep::unit_cube(tol()).expect("cube");
    let map = pair.append_body(&small_cube_at(3.0));
    let first = pair.solid_regions()[0];
    let (solid_shell, partner) = {
        let f = pair.shell_faces(pair.region_outer_shell(first).expect("outer"))[0];
        (pair.faces[f].back, pair.faces[f].front)
    };
    let moved: Vec<_> = map.faces.values().copied().collect();
    for f in moved {
        pair.faces[f].back = solid_shell;
        pair.faces[f].front = partner;
    }
    let errors = pair.validate_shell_cells_all();
    for shell in [solid_shell, partner] {
        assert!(
            errors.iter().any(|e| matches!(e,
                ValidationError::ShellDisconnected { shell: s, pieces: 2 } if *s == shell)),
            "{errors:?}"
        );
        assert!(
            errors.iter().any(|e| matches!(e,
                ValidationError::ShellEuler { shell: s, euler: 4 } if *s == shell)),
            "{errors:?}"
        );
    }
}

/// Negative control: a loop listed twice on its face counts twice, and
/// the shell's count is odd.
#[test]
fn a_wrong_euler_count_is_reported() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let face = cube.faces.keys().next().expect("a face");
    let outer = cube.faces[face].outer_loop;
    cube.faces[face].inner_loops.push(outer);
    let errors = cube.validate_shell_cells_all();
    assert_eq!(errors.len(), 2, "{errors:?}");
    assert!(
        errors
            .iter()
            .all(|e| matches!(e, ValidationError::ShellEuler { euler: 1, .. }))
    );
}

/// Negative controls for the region structure: an outside that is a
/// solid, a region with no shell, a shell of a region that is gone.
#[test]
fn a_broken_region_structure_is_reported() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let solid = cube.solid_regions()[0];

    let mut solid_outside = cube.clone();
    let outside = solid_outside.outside;
    solid_outside.regions[outside].kind = RegionKind::Solid;
    assert!(
        solid_outside
            .validate_regions_all()
            .iter()
            .any(|e| matches!(
                e,
                ValidationError::RegionStructure {
                    defect: RegionDefect::OutsideNotVoid,
                    ..
                }
            ))
    );

    let mut shell_less = cube.clone();
    shell_less.regions[solid].shell = Default::default();
    assert!(matches!(
        shell_less.validate_regions_all().as_slice(),
        [ValidationError::RegionStructure { region, defect: RegionDefect::NoOuterShell }]
            if *region == solid
    ));

    let mut orphan = cube.clone();
    let gone = orphan.regions.insert(crate::topology::Region {
        kind: RegionKind::Void,
        shell: Default::default(),
    });
    orphan.regions.remove(gone);
    orphan.shells.insert(Shell::new(Default::default(), gone));
    assert!(orphan.validate_regions_all().iter().any(|e| matches!(
        e,
        ValidationError::RegionStructure {
            defect: RegionDefect::MissingRegion,
            ..
        }
    )));
}
