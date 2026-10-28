//! Cell complex model and validator (PK-969).
//!
//! The fixture lifts a whole closed BRep into a complex (one complex
//! face per BRep face, one edge per BRep edge). Each invariant then
//! has a negative control: a fault that breaks only it must be
//! reported as exactly its finding.

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::artifact::complex_artifact_json;
use super::ids::{EdgeId, FaceId, LoopId, ShellId};
use super::model::{CellComplex, FaceOwner, SideIndex};
use super::operand_set::OperandSet;
use super::validate::{ComplexFinding, validate_complex};
use super::*;

/// Lift every face of `brep` (operand 0) into a fresh complex.
pub(crate) fn complex_of_brep(brep: &BRep) -> CellComplex {
    let mut c = CellComplex::new(1);
    let mut vmap = HashMap::new();
    for (vh, v) in &brep.vertices {
        vmap.insert(vh, c.add_vertex(v.position, 1e-6, None));
    }
    let mut emap = HashMap::new();
    for (eh, e) in &brep.edges {
        let he = &brep.half_edges[e.half_edge];
        let tw = &brep.half_edges[he.twin];
        let id = c.add_edge(
            vmap[&he.origin],
            vmap[&tw.origin],
            e.curve.clone(),
            1e-6,
            None,
        );
        emap.insert(eh, (id, e.half_edge));
    }
    for (fh, f) in &brep.faces {
        let face = c.add_face(
            vec![FaceOwner {
                operand: 0,
                face: fh,
                same_sense: true,
            }],
            true,
            Vec::new(),
        );
        for (k, &lh) in std::iter::once(&f.outer_loop)
            .chain(f.inner_loops.iter())
            .enumerate()
        {
            let start = brep.loops[lh].half_edge;
            let mut uses = Vec::new();
            let mut he = start;
            loop {
                let h = &brep.half_edges[he];
                let (edge, rep) = emap[&h.edge];
                uses.push((edge, he == rep, None));
                he = h.next;
                if he == start {
                    break;
                }
            }
            c.add_loop(face, &uses, k == 0);
        }
    }
    c
}

fn cube_complex() -> CellComplex {
    complex_of_brep(&BRep::unit_cube(Tolerance::default()).expect("cube builds"))
}

#[test]
fn lifted_cube_is_well_formed() {
    let c = cube_complex();
    assert_eq!(
        (c.vertices.len(), c.edges.len(), c.faces.len(), c.fins.len()),
        (8, 12, 6, 24)
    );
    assert!(c.edges.iter().all(|e| e.radial.len() == 2));
    assert_eq!(validate_complex(&c), Vec::new());
}

#[test]
fn lifted_cylinder_seam_carries_two_fins_of_one_face() {
    let c = complex_of_brep(&BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder"));
    assert_eq!(validate_complex(&c), Vec::new());
    let seam_like = c.edges.iter().any(|e| {
        e.radial.len() == 2 && c.fins[e.radial[0].index()].face == c.fins[e.radial[1].index()].face
    });
    assert!(seam_like, "the cylinder seam is used twice by its wall");
}

#[test]
fn broken_loop_link_is_an_open_loop() {
    let mut c = cube_complex();
    let first = c.loops[0].first;
    let second = c.fins[first.index()].next;
    c.fins[first.index()].next = first; // short-circuit the ring
    let findings = validate_complex(&c);
    assert!(
        findings
            .iter()
            .any(|f| matches!(f, ComplexFinding::OpenLoop { .. })),
        "{findings:?}"
    );
    let _ = second;
}

#[test]
fn reversed_fin_is_a_loop_gap() {
    let mut c = cube_complex();
    let fin = c.loops[0].first;
    c.fins[fin.index()].forward = !c.fins[fin.index()].forward;
    let findings = validate_complex(&c);
    assert!(
        findings
            .iter()
            .any(|f| matches!(f, ComplexFinding::LoopGap { .. })),
        "{findings:?}"
    );
}

#[test]
fn dropped_radial_entry_is_reported_with_the_free_edge() {
    let mut c = cube_complex();
    let fin = c.edges[0].radial.pop().expect("a fin");
    let findings = validate_complex(&c);
    assert!(findings.contains(&ComplexFinding::FreeEdge {
        edge: EdgeId::from_index(0)
    }));
    assert!(
        findings
            .iter()
            .any(|f| matches!(f, ComplexFinding::RadialMembership { fin: x, .. } if *x == fin)),
        "{findings:?}"
    );
}

#[test]
fn face_without_outer_loop_first_is_reported() {
    let mut c = cube_complex();
    let l = c.faces[0].loops[0];
    c.loops[l.index()].outer = false;
    assert!(validate_complex(&c).contains(&ComplexFinding::FaceLoops {
        face: FaceId::from_index(0)
    }));
}

#[test]
fn dangling_fin_edge_is_reported_before_anything_else() {
    let mut c = cube_complex();
    c.fins[0].edge = EdgeId::from_index(999);
    let findings = validate_complex(&c);
    assert!(
        findings
            .iter()
            .all(|f| matches!(f, ComplexFinding::DanglingId { .. })),
        "{findings:?}"
    );
    assert!(!findings.is_empty());
}

#[test]
fn unassigned_side_and_missing_membership_are_reported() {
    let mut c = cube_complex();
    c.shells.push(model::ComplexShell {
        sides: vec![(FaceId::from_index(0), SideIndex::Ahead)],
        membership: Some(OperandSet::new()),
    });
    c.shells.push(model::ComplexShell::default());
    c.faces[0].sides[0].shell = Some(ShellId::from_index(0));
    let findings = validate_complex(&c);
    assert!(findings.contains(&ComplexFinding::UnassignedSide {
        face: FaceId::from_index(0),
        side: SideIndex::Behind
    }));
    assert!(findings.contains(&ComplexFinding::MissingMembership {
        shell: ShellId::from_index(1)
    }));
}

#[test]
fn artifact_is_deterministic_and_complete() {
    let a = complex_artifact_json(&cube_complex());
    let b = complex_artifact_json(&cube_complex());
    assert_eq!(a, b);
    for key in [
        "\"vertices\"",
        "\"edges\"",
        "\"fins\"",
        "\"faces\"",
        "\"shells\"",
    ] {
        assert!(a.contains(key), "{key} missing");
    }
    assert_eq!(a.matches("\"radial\"").count(), 12);
}

#[test]
fn operand_set_grows_past_64_operands() {
    let mut s = OperandSet::new();
    assert!(s.is_empty());
    s.insert(3);
    s.insert(70);
    assert!(s.contains(3) && s.contains(70) && !s.contains(4));
    assert_eq!(s.len(), 2);
    assert_eq!(s.iter().collect::<Vec<_>>(), vec![3, 70]);
}

/// Shells and membership of a lifted closed body: exactly two shells,
/// the one behind every face's outward normal inside the operand.
fn assert_inside_outside(brep: &BRep) {
    let mut c = complex_of_brep(brep);
    let seed = brep.faces.keys().next().expect("a face");
    let outward = matches!(
        brep.component_material_sense(seed),
        Some(popkern_core::validation::MaterialSense::FoldedOutward)
    );
    for f in &mut c.faces {
        f.folded_is_outward = outward;
    }
    let geo = radial::FaceGeometry {
        operands: [brep, brep],
    };
    let sides = radial::order_radially(&mut c, &geo, 1e-4, 1e-6).expect("radial order");
    regions::assemble_shells(&mut c, &sides).expect("shells");
    assert_eq!(c.shells.len(), 2, "inside and outside");
    let classified = regions::assign_membership(&mut c, |_, _, _| None).expect("membership");
    assert_eq!(classified, 0, "every shell touches an owned face");
    assert_eq!(validate_complex(&c), Vec::new());
    for face in &c.faces {
        let behind = face.sides[SideIndex::Behind as usize]
            .shell
            .expect("assigned");
        let ahead = face.sides[SideIndex::Ahead as usize]
            .shell
            .expect("assigned");
        let m = |s: ShellId| c.shells[s.index()].membership.clone().expect("membership");
        assert!(m(behind).contains(0) && !m(ahead).contains(0));
    }
}

#[test]
fn lifted_cube_has_an_inside_and_an_outside() {
    assert_inside_outside(&BRep::unit_cube(Tolerance::default()).expect("cube"));
}

#[test]
fn lifted_cylinder_has_an_inside_and_an_outside() {
    assert_inside_outside(&BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder"));
}

/// Negative control: one face declared with the wrong material sense
/// puts an "inside" side in the outside shell, and membership must
/// report the conflict rather than pick a side.
#[test]
fn one_wrong_face_sense_is_a_membership_conflict() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let mut c = complex_of_brep(&brep);
    let outward = matches!(
        brep.component_material_sense(brep.faces.keys().next().expect("face")),
        Some(popkern_core::validation::MaterialSense::FoldedOutward)
    );
    for f in &mut c.faces {
        f.folded_is_outward = outward;
    }
    c.faces[0].folded_is_outward = !outward; // the fault
    let geo = radial::FaceGeometry {
        operands: [&brep, &brep],
    };
    let sides = radial::order_radially(&mut c, &geo, 1e-4, 1e-6).expect("radial order");
    regions::assemble_shells(&mut c, &sides).expect("shells");
    let result = regions::assign_membership(&mut c, |_, _, _| None);
    assert!(
        matches!(
            result,
            Err(regions::RegionDecline::MembershipConflict { operand: 0, .. })
        ),
        "{result:?}"
    );
}

#[test]
fn loop_fins_terminates_on_a_corrupt_cycle() {
    let mut c = cube_complex();
    let first = c.loops[0].first;
    let second = c.fins[first.index()].next;
    c.fins[second.index()].next = second; // a self-cycle not through `first`
    let fins = c.loop_fins(LoopId::from_index(0));
    assert!(fins.len() <= c.fins.len() + 1);
    assert!(
        validate_complex(&c)
            .iter()
            .any(|f| matches!(f, ComplexFinding::OpenLoop { .. }))
    );
}
